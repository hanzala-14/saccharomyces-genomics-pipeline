
---

# 📖 Module 4: Genome Coverage Profiling & Sliding Window Analysis

* **Repository:** `hanzala-14/saccharomyces-genomics-pipeline`
* **Module Source:** `modules/coverage.nf`
* **Workflow Source:** `main.nf`
* **Pipeline Configuration:** `nextflow.config`
* **Pipeline Version:** `2.1.4`
* **Last Updated:** `06 August 2026`

---

## 1) Purpose & Scope

Module 4 provides deep genomic coverage profiling and smoothing for mapped cohorts. Because depth profiling is computationally intensive and not always required for immediate variant-calling passes, this module is built with an **optional ON/OFF toggle**.

When enabled, it computes exact per-base depth from deduplicated BAMs, converts the metrics into formatted bedGraph files, sorts them against reference layouts, and calculates **sliding window median coverage** to eliminate local genomic noise (such as repetitive regions or GC bias).

---

## 2) Design Principles

### 2.1 Conditional Execution (ON/OFF Toggle)

To save compute time and storage when running large cohorts (like your 105 lager strains or future 5,000-genome compendium), Module 4 can be turned on or off globally via `params.run_coverage`.

### 2.2 Leveraging Nextflow Caching (`-resume`)

When `run_coverage` is toggled off, the pipeline bypasses coverage calculation entirely and generates downstream outputs immediately. If coverage data is needed later, flipping the parameter to `true` and running with `-resume` allows Nextflow to instantly ingest the cached BAMs from Module 3 and compute coverage without re-running upstream steps.

### 2.3 Concurrent Strain Processing

Unlike naive Bash scripts that loop sequentially through files, Module 4 processes each strain concurrently across available CPU cores, dramatically accelerating execution times for large datasets.

### 2.4 Noise Smoothing via Sliding Windows

Using `bedtools map` with a median calculation over pre-defined genomic sliding windows produces clean, normalized coverage profiles suitable for structural variant analysis, aneuploidy screening, and publication-ready visualization.

---

## 3) Input Contract

### 3.1 From Module 3

| Channel | Tuple Shape | Source |
| --- | --- | --- |
| `MARKDUP.out.markdup_bam` | `(strain_id, mdup.bam, mdup.bam.bai)` | `MARKDUP` |

### 3.2 Reference Files

* **Chromosome Size File (`params.genome_file`):** Tab-delimited reference layout (`chrom \t size`).
* **Sliding Windows BED File (`params.sliding_windows`):** Pre-defined window intervals across the reference genome.

---

## 4) Process Specification

### 4.1 `COMPUTE_COVERAGE`

**Role:** Calculate exact per-base depth, convert to sorted bedGraph, and compute sliding window median coverage.

* **Label:** `low`
* **Container:** `biocontainers/bedtools:2.31.1--hcb2786e_0`
* **Input tuple:** `tuple val(strain_id), path(bam), path(bai)` + reference path variables.
* **Output tuple:** `tuple val(strain_id), path("${strain_id}.slidingwindow.tab")`

#### Execution Steps

1. **Per-base depth:** `bedtools genomecov -ibam ${bam} -g ${genome_file} -d`
2. **BedGraph conversion:** `awk` formats depth output into standard coordinate intervals (`chrom start end depth`).
3. **Coordinate sorting:** `bedtools sort` organizes intervals according to the reference genome layout.
4. **Sliding window aggregation:** `bedtools map` computes the median coverage value (`-c 4 -o median`) across the sliding window BED intervals.
5. **Cleanup & Validation:** Temporary tab/bedGraph files are scrubbed to preserve disk space, and output size is validated to ensure non-empty results.

---

## 5) Wiring in `main.nf`

```groovy
include { COMPUTE_COVERAGE } from './modules/coverage.nf'

// =========================================================================
// MODULE 4: Coverage Profiling (Conditional Execution)
// =========================================================================
if (params.run_coverage) {
    COMPUTE_COVERAGE(
        MARKDUP.out.markdup_bam,
        file(params.genome_file),
        file(params.sliding_windows)
    )
}

```

---

## 6) Configuration

In `nextflow.config`:

```groovy
params {
    // --- Module 4: Coverage Profiling (Optional) ---
    run_coverage      = false                    // Toggle Module 4 ON (true) or OFF (false)
    genome_file       = 'data/references/sensustricto_genomefile.tab'        // Genome chromosome size file
    sliding_windows   = 'data/references/sensustrictoslidingwindows.bed'     // Sliding windows BED file
}

```

---

## 7) Output Layout

```text
results/
└── coverage/
    ├── strainA.slidingwindow.tab
    ├── strainB.slidingwindow.tab
    └── ...

```

---

## 8) Tool Reference

| Tool | Version | Purpose |
| --- | --- | --- |
| BEDTools | 2.31.1 | Genome coverage, interval sorting, and feature mapping |

---

## 9) Change Management Rule

Any logic change affecting coverage parameters, bedtools flags, reference mapping files, or publishDir paths must update:

1. `docs/MODULE_4_COVERAGE.md` (this file)
2. `modules/coverage.nf`
3. The Module 4 section of `main.nf`