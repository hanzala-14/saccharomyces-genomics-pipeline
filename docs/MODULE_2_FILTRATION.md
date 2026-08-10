
---

# 📖 Module 2: Read Filtration

* **Repository:** `hanzala-14/saccharomyces-genomics-pipeline`
* **Module Source:** `modules/filter_reads.nf`
* **Workflow Source:** `main.nf`
* **Pipeline Configuration:** `nextflow.config`
* **Pipeline Version:** `2.1.2`
* **Last Updated:** `02 August 2026`


---

## 1) Purpose & Scope

Module 2 performs quality-based read filtration on the merged strain-level FASTQs produced by Module 1. It uses **fastp** to:

1. Trim low-quality bases from read ends (sliding window).
2. Remove reads shorter than a configurable minimum length.
3. Disable polyG tail trimming (`-G`) — appropriate for non-NovaSeq/NextSeq platforms.
4. Detect and trim adapters automatically (PE adapter detection enabled).
5. Perform overrepresented sequence analysis.
6. Generate per-strain JSON and HTML quality reports.
7. Produce a cohort-wide `filtration_summary.tsv` aggregating key metrics across all strains.

The module handles **paired-end** and **single-end** streams independently through two dedicated processes, preserving the channel separation established in Module 1.

---

## 2) Design Principles

### 2.1 Separate processes for PE and SE

Rather than using conditional logic inside a single process, Module 2 uses two distinct processes (`FILTER_PE` and `FILTER_SE`). This keeps each process simple, testable, and correctly parameterized for its read layout.

### 2.2 Channel-driven routing (no filename parsing)

The PE/SE decision is made entirely by which Nextflow channel carries the data. `MERGE_PE.out.merged_pe` feeds `FILTER_PE`; `MERGE_SE.out.merged_se` feeds `FILTER_SE`. Filenames are not inspected to determine read layout.

### 2.3 Empty-channel safety

If all strains are PE-only (no SE data exists), `MERGE_SE.out.merged_se` is an empty channel. Nextflow natively handles this: `FILTER_SE` simply never executes. The pipeline does not hang or error. The same applies in reverse (all SE, no PE).

### 2.4 Flat output directory

All filtered files from all strains are published into a single flat directory (`filtration/filtered/`). Strain identity is encoded in the filename prefix, making it easy to glob or iterate without traversing subdirectories.

### 2.5 Dynamic Parameter Injection

Filtration parameters (read length, quality thresholds, and base trimming) are injected dynamically from `nextflow.config`. If a parameter (like `trim_front1`) is set to `0` or `null`, Nextflow safely ignores it and omits the flag from the command line. This allows rapid adaptation to noisy datasets without altering the core `.nf` script.

### 2.6 Post-filter validation (fail-safe)

Every filtered output undergoes:

1. **Gzip integrity check** (`gzip -t`) — catches truncated or corrupted output.
2. **Non-empty validation** — verifies at least 1 read survives filtering. If all reads are filtered out, the task fails with a clear error message rather than silently passing an empty file downstream.

### 2.7 Traceability

Read counts after filtering are logged to stdout, captured by Nextflow's trace system, and visible in the execution log.

### 2.8 Reports alongside data

JSON reports (machine-readable, MultiQC-compatible) and HTML reports (human-readable) are published to a dedicated `filtration/reports/` directory for centralized QC review.

---

## 3) Input Contract

### 3.1 From Module 1 → FILTER_PE

| Channel | Tuple shape | Source |
| --- | --- | --- |
| `MERGE_PE.out.merged_pe` | `(strain_id, PE_R1.fastq.gz, PE_R2.fastq.gz)` | `MERGE_PE` process |

### 3.2 From Module 1 → FILTER_SE

| Channel | Tuple shape | Source |
| --- | --- | --- |
| `MERGE_SE.out.merged_se` | `(strain_id, SE.fastq.gz)` | `MERGE_SE` process |

---

## 4) Process Specification

### 4.1 `FILTER_PE`

**Role:** Quality-filter paired-end merged reads for a strain.

* **Label:** `base`
* **Container:** `biocontainers/fastp:1.3.6--h5f740d0_0`
* **Input tuple:** `(strain_id, path(r1), path(r2))`
* **Output tuples:**

| Emit | Shape | Description |
| --- | --- | --- |
| `filtered_pe` | `(strain_id, <strain>.filt.R1.fastq.gz, <strain>.filt.R2.fastq.gz)` | Filtered paired reads |
| `json_report` | `(strain_id, <strain>.fastp.json)` | Machine-readable QC metrics |
| `html_report` | `(strain_id, <strain>.fastp.html)` | Human-readable QC report |

#### fastp flags (Static & Dynamic)

| Flag | Value | Rationale |
| --- | --- | --- |
| `-i` / `-I` | Input R1 / R2 | Paired-end input |
| `-o` / `-O` | Output R1 / R2 | Filtered output |
| `-G` | — | Disable polyG trimming (not NovaSeq/NextSeq) |
| `-l` | `params.min_read_length` | Minimum read length after trimming |
| `--detect_adapter_for_pe` | — | Auto-detect and trim PE adapters |
| `-q` | `params.qualified_quality` | Enforce Q-score threshold *(Injected if >0)* |
| `--cut_right ...` | `params.cut_mean_quality` | Sliding window quality *(Injected if >0)* |
| `--trim_front1` | `params.trim_front1` | Trim N bases from R1 start *(Injected if >0)* |
| `--adapter_fasta` | `params.adapter_fasta` | Custom adapters *(Injected if not null)* |

#### Post-filter checks

1. `gzip -t` on both R1 and R2 output files.
2. Verify at least 4 lines (1 FASTQ record) exist in each output.
3. Log final read counts to stdout.

---

### 4.2 `FILTER_SE`

**Role:** Quality-filter single-end merged reads for a strain.

* **Label:** `base`
* **Container:** `biocontainers/fastp:1.3.6--h5f740d0_0`
* **Input tuple:** `(strain_id, path(se))`
* **Output tuples:**

| Emit | Shape | Description |
| --- | --- | --- |
| `filtered_se` | `(strain_id, <strain>.filt.SE.fastq.gz)` | Filtered single-end reads |
| `json_report` | `(strain_id, <strain>.SE.fastp.json)` | Machine-readable QC metrics |
| `html_report` | `(strain_id, <strain>.SE.fastp.html)` | Human-readable QC report |

#### fastp flags

Same as `FILTER_PE`, except:

* Only `-i` / `-o` are used (no mate).
* `--detect_adapter_for_pe` is omitted (not applicable to SE).
* R2-specific trimming parameters are omitted.

#### Post-filter checks

1. `gzip -t` on output file.
2. Verify at least 4 lines (1 FASTQ record) exist.
3. Log final read count to stdout.

---

### 4.3 `WRITE_FILTER_SUMMARY`

**Role:** Parse all fastp JSON reports and produce a cohort-wide TSV summarizing filtration metrics.

* **Label:** `tiny`
* **Input:** Collected PE JSON paths + collected SE JSON paths
* **Output:** `filtration_summary.tsv`

#### TSV columns

| Column | Description |
| --- | --- |
| `strain_id` | Strain identifier |
| `mode` | `PE` or `SE` |
| `reads_before` | Total reads before filtering |
| `reads_after` | Total reads after filtering |
| `reads_passed_pct` | Percentage of reads surviving filtration |
| `q20_rate_before` | % bases ≥ Q20 before filtering |
| `q30_rate_before` | % bases ≥ Q30 before filtering |
| `adapter_trimmed_pct` | % reads with adapter detected and trimmed |
| `duplicate_rate` | Estimated duplication rate |

---

## 5) Wiring in `main.nf`

```groovy
include {
    FILTER_PE
    FILTER_SE
    WRITE_FILTER_SUMMARY
} from './modules/filter_reads.nf'

// =========================================================================
// MODULE 2: Read Filtration (fastp)
// =========================================================================

// Both channels are safe to pass even if empty
FILTER_PE(MERGE_PE.out.merged_pe)
FILTER_SE(MERGE_SE.out.merged_se)

// --- Filtration Summary Report ---
pe_jsons = FILTER_PE.out.json_report.map { id, json -> json }.collect().ifEmpty([])
se_jsons = FILTER_SE.out.json_report.map { id, json -> json }.collect().ifEmpty([])

WRITE_FILTER_SUMMARY(pe_jsons, se_jsons)

```

---

## 6) Configuration

In `nextflow.config`:

```groovy
params {
    // --- Module 2: Read Filtration (fastp) ---
    adapter_fasta     = null                     // Custom adapter FASTA (optional, fastp auto-detects)
    trim_front1       = 0                        // Off
    trim_tail1        = 0                        // Off
    trim_front2       = 0                        // Off
    trim_tail2        = 0                        // Off
    min_read_length   = 30                       // Discard reads shorter than 30bp
    qualified_quality = 20                       // Enforce Q20 quality threshold
    cut_mean_quality  = 20                       // Enforce Q20 sliding window mean
}

```

*Note: Any additional custom fastp flags can still be passed via `ext.args` in the process block if required.*

---

## 7) Output Layout

```text
results/
└── filtration/
    ├── filtered/
    │   ├── strainA.filt.R1.fastq.gz
    │   ├── strainA.filt.R2.fastq.gz
    │   ├── strainB.filt.SE.fastq.gz
    │   └── ...
    └── reports/
        ├── strainA.fastp.json
        ├── strainA.fastp.html
        ├── strainB.SE.fastp.json
        └── ...

```

---

## 8) Downstream Channels (for Module 3+)

| Emit | Channel shape | Typical next step |
| --- | --- | --- |
| `FILTER_PE.out.filtered_pe` | `(strain_id, filt.R1.fastq.gz, filt.R2.fastq.gz)` | Alignment (BWA, Bowtie2), assembly |
| `FILTER_SE.out.filtered_se` | `(strain_id, filt.SE.fastq.gz)` | SE alignment, assembly |
| `FILTER_PE.out.json_report` | `(strain_id, fastp.json)` | MultiQC aggregation |
| `FILTER_SE.out.json_report` | `(strain_id, SE.fastp.json)` | MultiQC aggregation |

---

## 9) Tool Reference

| Tool | Version | Source | Purpose |
| --- | --- | --- | --- |
| fastp | 1.3.6 | bioconda | All-in-one FASTQ preprocessing |

---

## 10) Change Management Rule

Any logic change affecting filtration parameters, process signatures, output filenames, or publishDir paths must update:

1. `docs/MODULE_2_FILTRATION.md` (this file)
2. `modules/filter_reads.nf`
3. The Module 2 section of `main.nf`

This guarantees implementation and documentation remain synchronized.