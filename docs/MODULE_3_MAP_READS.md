
---

# 📖 Module 3: BWA Mapping, Duplicate Marking & Alignment QC

* **Repository:** `hanzala-14/saccharomyces-genomics-pipeline`
* **Module Source:** `modules/map_reads.nf`
* **Workflow Source:** `main.nf`
* **Pipeline Configuration:** `nextflow.config`
* **Pipeline Version:** `2.1.3`
* **Last Updated:** `06 September 2026`

---

## 1) Purpose & Scope

Module 3 takes the quality-filtered FASTQs from Module 2, aligns them to a reference genome, consolidates all mappings belonging to the same biological strain, processes alignments to remove PCR/optical duplicates, and collects rigorous alignment statistics.

It handles:

1. Reference genome acquisition (local FASTA or automated fallback via NCBI Assembly accession).

2. Reference indexing (`bwa index`).

3. Paired-end (`MAP_PE`) and single-end (`MAP_SE`) read mapping using `bwa mem` piped directly into `samtools view` and `samtools sort` to eliminate unnecessary intermediate SAM I/O.

4. Read group (`@RG`) injection per sequencing run to ensure compatibility with downstream variant callers (e.g., GATK).

5. Strain-level BAM consolidation through `MERGE_STRAIN_BAMS`, ensuring that all PE and/or SE mappings belonging to the same biological strain are represented by a single BAM before duplicate marking.

6. High-performance duplicate marking using a pure-samtools pipeline (`collate` → `fixmate` → `sort` → `markdup`), bypassing heavy Java/Picard memory overhead.

7. Generation of comprehensive alignment statistics (`flagstat`, `idxstats`, `samtools stats`).

8. Cohort-wide mapping summary reporting (`mapping_summary.tsv`) with automated low-mapping rate warnings.

---

## 2) Design Principles

### 2.1 Zero Intermediate Disk Bloat (Piped Streaming)

To prevent massive SAM or unnecessary intermediate coordinate-sorted BAM files from clogging disk storage and slowing down I/O, alignment streams are piped directly:

```text
bwa mem ... | samtools view -Sb -F 4 - | samtools sort ...
```

Only the mapping output BAM is written by each mapping task.

### 2.2 Run-Level Mapping, Strain-Level Consolidation

PE and SE sequencing runs are mapped independently so that each sequencing run retains its own read-group information and traceability.

The resulting BAMs are then grouped by `strain_id` and consolidated by `MERGE_STRAIN_BAMS` before duplicate marking.

This guarantees:

```text
PE run(s) ─┐
           ├─→ MERGE_STRAIN_BAMS → one BAM per strain → MARKDUP
SE run(s) ─┘
```

A strain with multiple sequencing runs therefore produces **one downstream marked BAM**, regardless of whether its input consists of PE runs, SE runs, or a combination of both.

### 2.3 Pure-Samtools Duplicate Marking

Rather than relying on memory-intensive Java/Picard tools, duplicate marking is executed using a samtools workflow:

```text
collate → fixmate → coordinate sort → markdup → index
```

This reduces Java heap requirements and keeps the duplicate-marking stage lightweight.

### 2.4 Flexible Reference Resolution

The pipeline prioritizes a local reference genome specified in `params.reference`. If a local file is absent and an assembly accession is provided (`params.genome_id`), the pipeline automatically queries the NCBI FTP server, downloads the genomic FASTA, and prepares it for indexing.

### 2.5 Fail-Safe BAM Validation

Mapping and strain-level BAM consolidation validate their outputs before downstream processing. BAMs that are empty or otherwise invalid are prevented from silently propagating into duplicate marking and variant calling.

### 2.6 Automated Audit & Compliance

Post-alignment QC metrics are harvested into a single cohort-wide `mapping_summary.tsv`, making it straightforward to audit total reads, mapping percentages, duplicate rates, and insert sizes at a glance.

### 2.7 Reference Handling & The `save_reference` Toggle

The pipeline supports two modes for acquiring a reference genome:

1. **Local Mode:** Pointing `params.reference` to an existing local `.fasta` file.

2. **NCBI Fetch Mode:** Providing an NCBI assembly accession via `params.genome_id` to automatically download the genome.

Regardless of the mode, the pipeline automatically detects missing BWA indices and generates them.

**The `save_reference` Parameter:**

* **If `false` (Default):** The pipeline builds the indices in the temporary Nextflow work directory and does not publish the indexed reference bundle.

* **If `true`:** The pipeline publishes the FASTA and generated BWA index bundle into the configured reference output directory.

---

## 3) Input Contract

### 3.1 From Module 2

| Channel                     | Tuple Shape                                       | Source      |
| --------------------------- | ------------------------------------------------- | ----------- |
| `FILTER_PE.out.filtered_pe` | `(strain_id, filt.R1.fastq.gz, filt.R2.fastq.gz)` | `FILTER_PE` |
| `FILTER_SE.out.filtered_se` | `(strain_id, filt.SE.fastq.gz)`                   | `FILTER_SE` |

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

**Role:** Build BWA indices required for alignment.

* **Label:** `medium`
* **Input:** `path(fasta)`
* **Output:** `tuple path(fasta), path("${fasta}.*")`
* **Actions:** Executes `bwa index`.

### 4.3 `MAP_PE` & `MAP_SE`

**Role:** Align paired-end or single-end reads to the indexed reference genome while preserving strain and run-level provenance.

* **Label:** `low`
* **Input:** Filtered FASTQs + indexed reference tuple
* **Output:**

  * `MAP_PE`: `tuple val(strain_id), path("${strain_id}.PE.sorted.bam")`
  * `MAP_SE`: `tuple val(strain_id), path("${strain_id}.SE.sorted.bam")`

**Key Features:**

* Uses `bwa mem` with dynamically constructed Illumina read groups (`@RG`).
* Includes sample ID (`SM`), library (`LB`), and platform unit (`PU`) information.
* Filters unmapped reads (`-F 4`) during the mapping stream.
* Writes separate PE and SE BAM filenames so that different mapping branches cannot collide for the same biological strain.

### 4.4 `MERGE_STRAIN_BAMS`

**Role:** Consolidate all mapped BAMs belonging to the same biological strain into a single BAM before duplicate marking.

* **Label:** `base`
* **Input:** Grouped mapped BAMs by `strain_id`
* **Output:** `tuple val(strain_id), path("${strain_id}.sorted.bam")`

#### Execution Logic

1. Receives all mapped BAMs associated with a single `strain_id`.
2. If only one BAM is present, it is passed through as the strain-level BAM.
3. If multiple BAMs are present, they are merged using `samtools merge`.
4. The resulting strain-level BAM is validated before being passed downstream.

This process is essential for strains with multiple sequencing runs and especially for **mixed PE/SE strains**.

Without strain-level consolidation, the PE and SE mapping branches could independently produce BAMs for the same biological strain, resulting in duplicate same-strain inputs downstream.

### 4.5 `MARKDUP`

**Role:** Mark PCR and optical duplicates on the consolidated strain-level BAM using pure `samtools`.

* **Label:** `medium`
* **Input:** `tuple val(strain_id), path(sorted_bam)`
* **Output:**

```text
tuple val(strain_id),
      path("${strain_id}.mdup.bam"),
      path("${strain_id}.mdup.bam.bai")
```

and:

```text
tuple val(strain_id),
      path("${strain_id}.mdup.metrics.txt")
```

**Workflow:**

```text
samtools collate
      ↓
samtools fixmate
      ↓
samtools sort
      ↓
samtools markdup
      ↓
samtools index
```

### 4.6 `MAPPING_STATS`

**Role:** Collect standard alignment metrics on the deduplicated BAM.

* **Label:** `tiny`
* **Input:** Deduplicated BAM and BAI
* **Output:** `flagstat.txt`, `idxstats.txt`, `stats.txt`

### 4.7 `WRITE_MAPPING_SUMMARY`

**Role:** Aggregate alignment stats, duplicate metrics, and insert sizes into a single cohort summary table.

* **Label:** `tiny`
* **Input:** Collected flagstat, metrics, and stats files
* **Output:** `mapping_summary.tsv`

#### TSV Columns

| Column                | Description                                                    |
| --------------------- | -------------------------------------------------------------- |
| `strain_id`           | Strain identifier                                              |
| `total_reads`         | Total sequence reads processed                                 |
| `mapped_reads`        | Number of successfully mapped reads                            |
| `mapped_pct`          | Percentage of successfully mapped reads                        |
| `properly_paired_pct` | Percentage of properly paired reads (PE only)                  |
| `singletons_pct`      | Percentage of singleton reads                                  |
| `duplicate_pct`       | Estimated PCR/optical duplicate rate                           |
| `insert_size_mean`    | Average insert size                                            |
| `status`              | `PASS` or `WARN_LOW_MAPPING` based on `params.min_mapping_pct` |

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

├── reference/                  # Populated ONLY if save_reference = true
│   ├── *.fasta
│   ├── *.fasta.bwt
│   └── *.fasta.pac             # (and other BWA index files)
│
└── mapped/
    ├── mapping_summary.tsv
    │
    ├── strains/
    │   ├── strainA.mdup.bam
    │   ├── strainA.mdup.bam.bai
    │   ├── strainB.mdup.bam
    │   └── ...
    │
    └── stats/
        ├── strainA.mdup.metrics.txt
        ├── strainA.flagstat.txt
        ├── strainA.idxstats.txt
        ├── strainA.stats.txt
        └── ...
```

The intermediate mapping BAMs produced by `MAP_PE` and `MAP_SE` are temporary Nextflow work-directory artifacts and are consolidated by `MERGE_STRAIN_BAMS` before duplicate marking.

---

## 7) Tool Reference

| Tool     | Version | Purpose                                                        |
| -------- | ------- | -------------------------------------------------------------- |
| BWA      | 0.7.17  | Maximal Exact Match (MEM) read alignment                       |
| Samtools | 1.16.1  | BAM sorting, viewing, merging, indexing, and duplicate marking |

---

## 8) Change Management Rule

Any logic change affecting mapping parameters, process signatures, BAM consolidation, duplicate marking steps, or output paths must update:

1. `docs/MODULE_3_MAPPING.md` (this file)
2. `modules/map_reads.nf`
3. The Module 3 section of `main.nf`

---
