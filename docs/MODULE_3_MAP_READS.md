
---

# 📖 Module 3: BWA Mapping, Duplicate Marking & Alignment QC

* **Repository:** `hanzala-14/saccharomyces-genomics-pipeline`
* **Module Source:** `modules/map_reads.nf`
* **Workflow Source:** `main.nf`
* **Pipeline Configuration:** `nextflow.config`
* **Pipeline Version:** `2.1.3`
* **Last Updated:** `05 August 2026`

---

## 1) Purpose & Scope

Module 3 takes the quality-filtered FASTQs from Module 2, aligns them to a reference genome, processes alignments to remove PCR/optical duplicates, and collects rigorous alignment statistics.

It handles:

1. Reference genome acquisition (local FASTA or automated fallback via NCBI Assembly accession).
2. Reference indexing (`bwa index`, `samtools faidx`, Picard sequence dictionary creation).
3. Paired-end (`MAP_PE`) and single-end (`MAP_SE`) read mapping using `bwa mem` piped directly into `samtools view` and `samtools sort` to eliminate intermediate disk I/O bottlenecks.
4. Read group (`@RG`) injection per run to ensure compatibility with downstream variant callers (e.g., GATK).
5. High-performance duplicate marking using a pure-samtools pipeline (`collate` $\rightarrow$ `fixmate` $\rightarrow$ `sort` $\rightarrow$ `markdup`), bypassing heavy Java/Picard memory overhead.
6. Generation of comprehensive alignment statistics (`flagstat`, `idxstats`, `samtools stats`).
7. Cohort-wide mapping summary reporting (`mapping_summary.tsv`) with automated low-mapping rate warnings.

---

## 2) Design Principles

### 2.1 Zero Intermediate Disk Bloat (Piped Streaming)

To prevent massive SAM or raw coordinate-sorted BAM files from clogging disk storage and slowing down I/O, alignment streams are piped directly:
`bwa mem ... | samtools view -Sb -F 4 - | samtools sort ...`
Only the final sorted BAM is written to disk.

### 2.2 Pure-Samtools Duplicate Marking

Rather than relying on memory-intensive Java/Picard tools (which frequently cause heap-space crashes on standard hardware), duplicate marking is executed via a highly optimized, multi-step `samtools` workflow (`collate` $\rightarrow$ `fixmate` $\rightarrow$ `sort` $\rightarrow$ `markdup`). This provides massive speedups and keeps RAM usage exceptionally low.

### 2.3 Flexible Reference Resolution

The pipeline prioritizes a local reference genome specified in `params.reference`. If a local file is absent and an assembly accession is provided (`params.genome_id`), the pipeline automatically queries the NCBI FTP server, downloads the latest genomic fasta, and prepares it for indexing.

### 2.4 Fail-Safe BAM Validation

Every mapping process verifies output integrity. Tasks that produce empty BAM files (0 mapped reads) fail immediately rather than passing corrupted or blank alignments downstream.

### 2.5 Automated Audit & Compliance

Post-alignment QC metrics are harvested into a single cohort-wide `mapping_summary.tsv`, making it straightforward to audit total reads, mapping percentages, duplicate rates, and insert sizes at a glance.

---

## 3) Input Contract

### 3.1 From Module 2

| Channel | Tuple Shape | Source |
| --- | --- | --- |
| `FILTER_PE.out.filtered_pe` | `(strain_id, filt.R1.fastq.gz, filt.R2.fastq.gz)` | `FILTER_PE` |
| `FILTER_SE.out.filtered_se` | `(strain_id, filt.SE.fastq.gz)` | `FILTER_SE` |

### 3.2 Reference Configuration

* **Local Mode:** `params.reference` (Path to local `.fasta` file)
* **Remote Mode:** `params.genome_id` (NCBI Assembly Accession, e.g., `GCF_000146045.2`)

---

## 4) Process Specification

### 4.1 `FETCH_REFERENCE`

**Role:** Automatically fetch reference genome from NCBI FTP if no local reference path is provided.

* **Label:** `base`
* **Input:** `val(genome_id)`
* **Output:** `path("*.fasta")`

### 4.2 `BWA_INDEX`

**Role:** Build indexes required for alignment and downstream analyses.

* **Label:** `medium`
* **Input:** `path(fasta)`
* **Output:** `tuple path(fasta), path("${fasta}.*")`
* **Actions:** Executes `bwa index`, `samtools faidx`, and Picard `CreateSequenceDictionary`.

### 4.3 `MAP_PE` & `MAP_SE`

**Role:** Align paired-end or single-end reads to the indexed reference genome.

* **Label:** `low`
* **Input:** Filtered FASTQs + indexed reference tuple
* **Output:** `tuple val(strain_id), path("${strain_id}.sorted.bam")`
* **Key Features:** Dynamically constructs Illumina read groups (`@RG`) containing sample IDs (`SM`), library tags (`LB`), and platform units (`PU`). Filters unmapped reads (`-F 4`) on the fly.

### 4.4 `MARKDUP`

**Role:** Mark PCR and optical duplicates using pure `samtools`.

* **Label:** `medium`
* **Input:** `tuple val(strain_id), path(sorted_bam)`
* **Output:**
* `tuple val(strain_id), path("${strain_id}.mdup.bam"), path("${strain_id}.mdup.bam.bai")`
* `tuple val(strain_id), path("${strain_id}.mdup.metrics.txt")`


* **Workflow:** Name collate $\rightarrow$ Fixmate $\rightarrow$ Coordinate sort $\rightarrow$ Markdup $\rightarrow$ Index.

### 4.5 `MAPPING_STATS`

**Role:** Collect standard alignment metrics on the deduplicated BAM.

* **Label:** `tiny`
* **Input:** Deduplicated BAM and BAI
* **Output:** `flagstat.txt`, `idxstats.txt`, `stats.txt`

### 4.6 `WRITE_MAPPING_SUMMARY`

**Role:** Aggregate alignment stats, duplicate metrics, and insert sizes into a single cohort summary table.

* **Label:** `tiny`
* **Input:** Collected flagstat, metrics, and stats files
* **Output:** `mapping_summary.tsv`

#### TSV Columns

| Column | Description |
| --- | --- |
| `strain_id` | Strain identifier |
| `total_reads` | Total sequence reads processed |
| `mapped_reads` | Number of successfully mapped reads |
| `mapped_pct` | Percentage of successfully mapped reads |
| `properly_paired_pct` | Percentage of properly paired reads (PE only) |
| `singletons_pct` | Percentage of singleton reads |
| `duplicate_pct` | Estimated PCR/optical duplicate rate |
| `insert_size_mean` | Average insert size |
| `status` | `PASS` or `WARN_LOW_MAPPING` (based on `params.min_mapping_pct`) |

---

## 5) Configuration

In `nextflow.config`:

```groovy
params {
    // --- Module 3: Mapping ---
    reference         = 'data/references/cer_eub/cer_eub_reference.fasta' // Local FASTA path
    genome_id         = null                     // NCBI Assembly accession fallback
    save_reference    = false                    // Publish indexed reference to outdir?
    min_mapping_pct   = 0                        // Mapping rate warning threshold (0 = disabled)
}

```

---

## 6) Output Layout

```text
results/
├── reference/
│   ├── *.fasta
│   ├── *.fasta.bwt
│   ├── *.fasta.fai
│   └── *.fasta.dict
└── mapped/
    ├── mapping_summary.tsv
    ├── strains/
    │   ├── strainA.mdup.bam
    │   ├── strainA.mdup.bam.bai
    │   ├── strainB.mdup.bam
    │   └── ...
    └── stats/
        ├── strainA.mdup.metrics.txt
        ├── strainA.flagstat.txt
        ├── strainA.idxstats.txt
        ├── strainA.stats.txt
        └── ...

```

---

## 7) Tool Reference

| Tool | Version | Purpose |
| --- | --- | --- |
| BWA | 0.7.19 | Maximal Exact Match (MEM) read alignment |
| Samtools | 1.24 | BAM sorting, viewing, indexing, and duplicate marking |
| Picard | 3.5.0 | Sequence dictionary generation (`CreateSequenceDictionary`) |

---

## 8) Change Management Rule

Any logic change affecting mapping parameters, process signatures, duplicate marking steps, or output paths must update:

1. `docs/MODULE_3_MAPPING.md` (this file)
2. `modules/map_reads.nf`
3. The Module 3 section of `main.nf`