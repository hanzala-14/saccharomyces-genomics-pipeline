
---

# 📖 Module 8: Phylogenomics (IBSx & Tree Construction)

* **Repository:** `hanzala-14/saccharomyces-genomics-pipeline`
* **Module Source:** `modules/phylogenomy.nf`
* **Custom Scripts:** `bin/run_ibs.R`, `bin/build_tree.R`, `bin/ibs.cpp`
* **Workflow Source:** `main.nf`
* **Pipeline Configuration:** `nextflow.config`
* **Pipeline Version:** `2.1.8`
* **Last Updated:** `14 August 2026`

---

## 1) Purpose & Scope

Module 8 is the analytical terminus of the pipeline, transforming pristine, hard-filtered VCFs into evolutionary insights. It is engineered specifically to unravel the complex reticulate evolutionary history of the *Saccharomyces Sensu Stricto* complex.

Because yeast species frequently undergo hybridization, introgression, and incomplete lineage sorting (ILS), genome-wide consensus trees are often biologically misleading. To counteract this, Module 8 executes a three-stage chromosomal phylogenomics workflow: (1) Vectorized haploid-equivalent allele splitting, (2) High-performance Identity-by-State (IBS) matrix calculation via an OpenMP C++ backend, and (3) Automated tree inference, splits-graph generation, and Robinson-Foulds topological cross-comparisons.

---

## 2) Algorithmic & Mathematical Architecture

### 2.1 Genotype Normalization & Allele Splitting (`PREP_ALLELES` / `awk`)

Network algorithms (like NeighborNet) cannot natively process unphased diploid representations (e.g., `0/1`). The module utilizes a highly optimized `awk` vectorization engine to normalize all VCF genotypes and mathematically split them into distinct `Allele 1` and `Allele 2` haploid-equivalent arrays.

* **Invariant Site Purging:** To reduce matrix bloat and computational overhead, monomorphic sites (where all samples share the reference allele) are actively evaluated and purged from the data stream.
* **Missing Data Handling:** Genotypes lacking sufficient depth or quality (e.g., `./.`) are explicitly converted to a null state (`..`) to ensure downstream denominator penalties are strictly avoided during IBS calculations.

### 2.2 High-Performance IBS Computation (`ibs.cpp` & `run_ibs.R`)

Calculating pairwise genetic distance across millions of variable sites requires significant computational scaling. `run_ibs.R` acts as an orchestrator, passing the haploid matrices into a custom, multithreaded C++ engine.

* **The Distance Metric:** The engine calculates the pairwise IBS distance $D_{ij}$ between any two yeast strains $i$ and $j$ based on the proportion of shared alleles:

$$D_{ij} = 1 - \frac{\text{Shared Alleles}_{ij}}{\text{Valid Alleles}_{ij}}$$


* **Dynamic Missingness Correction:** If a locus is missing in either strain $i$ or $j$, the C++ engine dynamically decrements the $\text{Valid Alleles}$ denominator. This prevents regions of low sequencing coverage from artificially inflating evolutionary distance.
* **Chunked OpenMP Threading:** To prevent out-of-memory (OOM) cluster crashes, the C++ engine chunks the genomic arrays (defined by `ibs_block_rows = 5000`) and processes the $N \times N$ matrix populating loops via OpenMP parallelization, strictly adhering to Nextflow's CPU allocations.

---

## 3) Process Specification & Data Flow

### Stage 1: `PREP_ALLELES`

* **Input:** Hard-filtered `${cohort_name}.final.vcf.gz`
* **Execution:** Dynamically queries the VCF header/data to extract all contigs, slicing the genome into independent chromosomal VCFs. Applies the AWK genotype splitter. Physically destroys resulting slices that possess 0 informative variants to shield the downstream C++ engine from null-pointer crashes.
* **Output:** Clean `.allele1.tsv` and `.allele2.tsv` matrices; `.meta.json` checksums tracking exact sample and variant counts for data provenance.

### Stage 2: `RUN_IBS`

* **Input:** Haploid matrices from Stage 1.
* **Execution:** Compiles and triggers `ibs.cpp`. Matrix block size and OpenMP thread count are explicitly passed via CLI arguments to guarantee stable cluster scaling.
* **Output:** `FINAL_combined_distance_matrix.tsv` + 16 chromosomal distance matrices.

---

## 4) Configuration Contract (`nextflow.config`)

The mathematical depth and structural behavior of the module are entirely controlled via the config file:

```groovy
    // IBS Compute Backend Limits
    ibs_block_rows       = 8192          // C++ chunk size for memory bounding
    ibs_report_blocks    = 5             // OpenMP logging frequency

```

---

## 5) Dependency Reference

| Tool | Version / Origin | Scientific Purpose |
| --- | --- | --- |
| `bcftools` | 1.19 | Contig extraction and variant isolation |
| `awk` | Native Linux | Vectorized genotype string manipulation |
| `ibs.cpp` | Custom Binary | High-speed, OpenMP-accelerated IBS distance math |
| `ape` | R Package | NJ algorithms, bootstrapping, and midpoint rooting |
| `phangorn` | R Package | NeighborNet algorithms and Robinson-Foulds metrics |


## 8) Change Management Rule

Any logic change affecting the `awk` string splitting, the OpenMP thread allocation, or the `saveAs` interceptor routing must update:

1. `docs/MODULE_8_PHYLOGENOMICS.md` (this file)
2. `modules/phylogenomy.nf`
3. The R/C++ binaries located in `bin/`

---
