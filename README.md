
---

# Saccharomyces Genomics Pipeline

**Current Pipeline Version:** `2.1.8`

A production-grade, fault-tolerant, and highly optimized Whole Genome Sequencing (WGS) variant calling and comparative phylogenomics pipeline engineered for the *Saccharomyces* genus, with a primary focus on historic lager yeast lineages and large-scale strain compendiums.

Designed to scale from local development to high-performance computing (HPC) environments for large *Saccharomyces* strain collections.

---

## Pipeline Architecture & Workflow

The pipeline is organized into modular, independently testable Nextflow sub-workflows, scaling from raw data acquisition to deep genomic profiling, population-scale joint genotyping, strict variant filtration, and topological phylogenomics.

```mermaid
graph TD
    classDef active fill:#e1f5fe,stroke:#0288d1,stroke-width:2px,color:#000;
    classDef optional fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px,color:#000;
    classDef future fill:#f5f5f5,stroke:#616161,stroke-width:2px,stroke-dasharray: 5 5,color:#555;

    A[Samplesheet / ENA Accessions] -->|Module 1| B(Data Acquisition & QC)
    B -->|Filtered Fastqs| C(Module 2: Read Filtration fastp)
    C -->|Clean Reads| D(Module 3: Piped BWA Mapping, Strain BAM Consolidation & Samtools Markdup)
    
    D -->|run_coverage = true| E(Module 4: Coverage Profiling)
    D -->|BAM Streams| F(Module 5: GATK Variant Calling)
    
    F -->|GVCFs| G(Module 6: Joint Genotyping & Merging)
    G -->|Master VCFs| H(Module 7: Variant Filtration & QC)
    
    H -->|Analysis-Ready VCFs| I(Module 8: Phylogenomics & Trees)

    class B,C,D,F,G,H,I active;
    class E optional;
````

---

## Implemented Core Modules (Modules 1–8)

### [Module 1: Data Acquisition & QC Classification]

* **Metadata-Driven Acquisition:** Determines whether each ENA/SRA run is paired-end (`PE`) or single-end (`SE`) using the ENA filereport API before downloading reads.
* **High-Speed Downloads:** Primary ENA FASTQ retrieval via multi-threaded **`aria2c`** (`-x 16 -s 16`), backed by an NCBI `prefetch` / `fasterq-dump` fallback.
* **Strict Integrity:** Validates downloaded FASTQs using non-empty file checks and `gzip -t`. Paired-end runs require valid R1 **and** R2 before acquisition can succeed.
* **Fault-Tolerant Acquisition:** Incomplete or failed ENA downloads trigger the configured fallback/retry mechanism instead of silently converting a paired-end run into single-end data.
* **Smart Classification:** Automatically classifies runs into Paired-End (`PE_PASS`), Single-End (`SE_INPUT` / `SE_FALLBACK`), or dropped (`DROP`) while generating cohort-wide TSV summaries.

### [Module 2: Read Filtration]

* **Dual Streams:** Handles Paired-End (`FILTER_PE`) and Single-End (`FILTER_SE`) layouts independently through dedicated channel-driven processes.
* **Post-Filter Safety:** Enforces strict non-empty read counts and gzip structural validation post-filtering.

### [Module 3: Mapping, Strain BAM Consolidation, Duplicate Marking & Alignment QC]

* **Zero Disk Bloat:** Implements piped streaming (`bwa mem ... | samtools view -Sb -F 4 - | samtools sort ...`) so intermediate uncompressed alignments never touch the disk.
* **Run-Level Mapping:** Paired-end and single-end sequencing runs are mapped independently while retaining run-specific read-group information.
* **Strain-Level BAM Consolidation:** All mapped BAMs belonging to the same biological strain are grouped and merged by `MERGE_STRAIN_BAMS` before duplicate marking. This prevents duplicate downstream BAM/GVCF inputs when a strain contains multiple sequencing runs or mixed PE/SE data.
* **Pure-Samtools Duplicate Marking:** Uses a `samtools collate → fixmate → sort → markdup` workflow instead of Java/Picard-based duplicate marking, reducing memory overhead.

### [Module 4: Genome Coverage Profiling (Optional Toggle)]

* **Toggleable Execution:** Controlled via `params.run_coverage = true / false`. Skip coverage calculation during fast iterations, then flip it on and use Nextflow's `-resume` to backfill depth data.
* **Sliding Window Normalization:** Converts to bedGraph and calculates sliding window median coverage via `bedtools map` to reduce local genomic noise.

### [Module 5: Variant Calling (GATK HaplotypeCaller)]

* **Dynamic Reference Agnosticism:** Completely eliminates hardcoded chromosome lists by dynamically parsing the `.fai` index on the fly, allowing the workflow to adapt to custom reference genomes.
* **Biological Ploidy Calibration:** Features a config-driven `ploidy_map` that identifies organelle contigs (Mitochondria/Plasmids) and calls them as haploid (`-ploidy 1`), while nuclear DNA remains diploid by default.

### [Module 6: Joint Genotyping & Subgenome Merging]

* **Crash-Proof Database Updates:** Utilizes an atomic Bash wrapper for `GenomicsDBImport` that safely isolates existing databases during updates, protecting database integrity in the event of compute or workflow failures.
* **Universal *Sensu Stricto* Routing:** Employs a config-based taxonomy dictionary (`species_map`) to automatically detect, isolate, and route individual subgenomes (e.g., *cerevisiae*, *eubayanus*) into clean species-specific master VCFs.

### [Module 7: Variant Filtration & Quality Control]

* **Parallel Execution Architecture:** Splits SNPs and INDELs into independent computational streams, reducing unnecessary sequential processing.
* **Config-Driven Mathematics:** Injects GATK hard-filtering thresholds directly from `nextflow.config`, avoiding hardcoded filtering parameters and improving portability.
* **Biological Accuracy:** Dynamically extracts species prefixes to mask unmappable repetitive regions (e.g., Ty elements, telomeres) and automatically left-aligns INDEL coordinates to standardize downstream analysis.

### [Module 8: Phylogenomics & Evolutionary Trees]

* **C++ Identity-by-State (IBS) Engine:** Uses vectorized, haploid-equivalent allele matrices with an OpenMP-accelerated C++ backend for high-speed genetic distance calculations.
* **Introgression Detection:** Automatically generates chromosomal-level NeighborNet networks (SplitsTree compatible) and Robinson-Foulds topological heatmaps to help identify hybridization and incomplete lineage sorting.

---

## Upcoming Modules & Roadmap

As the project scales toward comprehensive population genomics, the following modules are currently in development for subsequent releases:

* **Module 9: Functional Annotation** (Integration of `SnpEff` / `VEP` to predict the phenotypic impact of surviving variants).
* **Module 10: Admixture & Population Structure** (Automated PCA and ancestral sub-population modeling).

---

## 💻 Quick Start

### Prerequisites

* [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (`>=23.10.0`)
* Container engine ([Docker](https://www.docker.com/) or [Singularity/Apptainer](https://apptainer.org/)) OR [Conda](https://docs.conda.io/)

### 1. Clone the Repository

```bash
git clone https://github.com/hanzala-14/saccharomyces-genomics-pipeline.git
cd saccharomyces-genomics-pipeline
```

### 2. Prepare your Samplesheet (`samplesheet.csv`)

Create a CSV file defining your strains. You can mix ENA/SRA accessions and local FASTQ paths:

```csv
strain_id,accession,r1,r2
Strain_A,SRR1234567,,
Strain_B,,/path/to/local_R1.fastq.gz,/path/to/local_R2.fastq.gz
```

### 3. Run the Pipeline

```bash
# Local Development
nextflow run main.nf -profile conda

# Docker
nextflow run main.nf -profile docker

# Production / HPC Cluster
nextflow run main.nf -profile singularity
```

---

## Advanced Usage: Incremental Updates & Compendiums

This pipeline supports **incremental database ingestion**. If your lab maintains a massive compendium (e.g., a 5,000-strain GenomicsDB on an external office drive), you do **not** need to re-run historical data to joint-call new strains.

To append new samples to an existing database:

1. Create a samplesheet containing *only* your new strains.
2. Run the pipeline, passing the parent directory of your existing databases to the `genomicsdb_update_path` flag:

```bash
nextflow run main.nf \
  --samplesheet new_strains.csv \
  --genomicsdb_update_path "/media/Office_Drive/Yeast_Compendium/Joint_Genotyping" \
  -profile singularity \
  -resume
```

*Nextflow will fast-track the new strains to GVCFs, safely inject them into the existing compendium, and generate updated Master VCFs.*

---

## ⚙️ Configuration & Parameters

All operational settings can be tuned inside `nextflow.config`.

```groovy
params {

    samplesheet       = 'samplesheet.csv'
    outdir            = 'results'

    // Module 1: Data Acquisition
    max_Retries              = 5
    min_r2_len               = 30
    sra_max_forks            = 4
    aria2c_connections       = 16
    aria2c_min_split_size    = '1M'
    aria2c_connect_timeout   = 30
    aria2c_timeout           = 60
    aria2c_max_tries         = 3
    aria2c_retry_wait        = 10
    aria2c_summary_interval  = 10

    ena_api_max_time         = 60
    ena_api_retries          = 3
    ena_api_retry_wait       = 5
    fetch_retry_sleep        = 10

    prefetch_timeout         = 3600
    prefetch_max_size        = '50G'
    fasterq_timeout          = 1800

    // Module 3: Mapping
    reference         = 'data/references/cer_eub/cer_eub_reference.fasta'
    genome_id         = null
    save_reference    = false
    min_mapping_pct   = 0

    // Module 5 & 6
    ploidy                    = 2
    ploidy_map                = [ 'M': 1, 'P': 1 ]
    include_extrachromosomal = false
    merge_strategy            = 'separate'
    extrachromosomal_strategy = 'bundled'

    // Taxonomy routing dictionary
    species_map = [
        'c': 'cerevisiae', 'e': 'eubayanus', 'p': 'paradoxus',
        'm': 'mikatae', 'k': 'kudriavzevii', 'u': 'uvarum',
        'a': 'arboricola', 'j': 'jurei'
    ]

    // Module 7 Storage Toggles
    keep_snp_vcf         = true
    keep_indel_vcf       = true
    keep_merged_vcf      = true

    // Module 8 Phylogenomics Engine
    run_phylogeny        = true
    ibs_strategy         = 'chromosomal'
}
```

---

## Output Directory Layout

```text
results/
├── Pipeline_Info/          # Execution reports, timeline, and DAG
├── Strains/                # Strain lists, QC reports, and merged FASTQs
├── Filtration/             # Filtered FASTQs and QC outputs
├── Mapped/                 # Sorted BAMs, index files, and mapping summaries
├── Haplotype_Calling/      # Per-strain GVCFs and manifests
├── Joint_Genotyping/
│   ├── c_I/                # Stateful GenomicsDB workspaces
│   └── Final_Merged/       # Pre-filtration cohort VCFs
├── Variant_Filtration/     # Hard-filtered SNPs, INDELs, and QC reports
└── Phylogeny/
    └── Trees/              # Phylogenetic outputs and network files
```

---

## 📖 Module Documentation

For deep technical specifications, input/output contracts, and command flag rationales, check the module guides:

* [`docs/MODULE_1_DATA_ACQUISITION.md`](docs/MODULE_1_DATA_ACQUISITION.md)
* [`docs/MODULE_2_READ_FILTRATION.md`](docs/MODULE_2_READ_FILTRATION.md)
* [`docs/MODULE_3_MAPPING.md`](docs/MODULE_3_MAPPING.md)
* [`docs/MODULE_4_COVERAGE.md`](docs/MODULE_4_COVERAGE.md)
* [`docs/MODULE_5_VARIANT_CALLING.md`](docs/MODULE_5_VARIANT_CALLING.md)
* [`docs/MODULE_6_JOINT_GENOTYPING.md`](docs/MODULE_6_JOINT_GENOTYPING.md)
* [`docs/MODULE_7_VARIANT_FILTRATION.md`](docs/MODULE_7_VARIANT_FILTRATION.md)
* [`docs/MODULE_8_PHYLOGENOMICS.md`](docs/MODULE_8_PHYLOGENOMICS.md)

---

## 📜 License

Distributed under the MIT License. See `LICENSE` for more information.

