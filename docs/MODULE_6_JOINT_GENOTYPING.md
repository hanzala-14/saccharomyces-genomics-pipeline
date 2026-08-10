
---

# 📖 Module 6: Joint Genotyping & Merging

* **Repository:** `hanzala-14/saccharomyces-genomics-pipeline`
* **Module Source:** `modules/joint_genotyping.nf`
* **Workflow Source:** `main.nf`
* **Pipeline Configuration:** `nextflow.config`
* **Pipeline Version:** `2.1.6`
* **Last Updated:** `09 August 2026`

---

## 1) Purpose & Scope

Module 6 is the capstone of the pipeline. It aggregates the individual Genomic VCFs (GVCFs) produced in Module 5 and performs **Joint Genotyping** across the entire population cohort. This multi-sample joint calling dramatically increases statistical power, allowing GATK to accurately call rare variants, filter out noise, and definitively distinguish homozygous reference sites from missing data.

The final outputs of this module are fully merged, analysis-ready master VCFs tailored for downstream admixture, PCA, and maximum-likelihood phylogenomic tree construction.

---

## 2) Design Principles

### 2.1 Atomic Updates & Compendium Data Protection

Standard Nextflow processes are stateless, but GATK’s `GenomicsDB` is a highly complex, stateful data directory. To allow the pipeline to append new strains to a massive existing database (e.g., a 5,000-strain compendium on an external drive) without risking data loss, this module utilizes a strict **Source-Isolation Architecture**:
* **Read-Only Source:** The script safely copies the existing database from the external source into an isolated Nextflow work directory. The source database is never modified directly.
* **Isolated Update:** GATK updates this local, isolated copy with the new GVCFs.
* **Safe Output Overwrite:** To prevent Nextflow `publishDir` conflicts, a pure-Bash `rm -rf` command clears out the *local output directory* (`results/Joint_Genotyping/`) just before publishing. It never executes against the source path, guaranteeing zero data corruption to the original compendium.

### 2.2 Universal *Sensu Stricto* Taxonomy Routing

The pipeline completely eliminates hardcoded organism logic (e.g., `if (chrom == 'c_I')`). It relies on a `species_map` dictionary in `nextflow.config` to dynamically tokenize chromosome names at the underscore (`c_`, `e_`, `p_`, etc.) and automatically route them to the correct species cohort. It scales infinitely from simple 2-way *S. pastorianus* hybrids up to massive 8-species pan-genomic *Sensu Stricto* compendiums.

### 2.3 Extrachromosomal Modularity

Population geneticists often require strictly nuclear genomes to calculate accurate phylogenetic distances. By utilizing the `extrachromosomal_strategy` toggle, the pipeline can dynamically alter the output prefix before merging:

* **Bundled:** Mitochondria and plasmids are packaged directly into the main species VCF.
* **Isolated:** The pipeline splits the biological components, outputting cleanly separated Nuclear, Mitochondrial, and Plasmid VCFs for each species.

---

## 3) Input Contract

| Channel | Tuple Shape | Source |
| --- | --- | --- |
| `ch_for_genomicsdb` | `(chrom, list(gvcfs), list(tbis))` | `HAPLOTYPE_CALLER.out.gvcf` |
| `gatk_ref` | `(fasta, fai, dict)` | `PREPARE_GATK_REF.out` |

---

## 4) Process Specification

### 4.1 `GENOMICSDB_IMPORT`

**Role:** Builds or incrementally updates a multi-dimensional datastore for variant calling.

* **Label:** `high`
* **Concurrency:** Throttled to `maxForks = 2` to prevent memory bandwidth starvation.
* **Dynamic Feature:** Reads `params.genomicsdb_update_path`. If null, builds a fresh DB. If provided, copies and incrementally updates the database with new GVCFs, saving days of compute time for large compendiums.

### 4.2 `GENOTYPE_GVCFS`

**Role:** Executes the joint calling algorithms across the cohort database to produce a raw VCF per chromosome.

* **Label:** `high`
* **Output:** `cohort.${chrom}.vcf.gz` (e.g., `cohort.c_I.vcf.gz`)

### 4.3 `MERGE_VCFS`

**Role:** A Picard-driven process that stitches the scattered per-chromosome VCFs into the final unified subgenome/species files.

* **Label:** `medium`
* **Dynamic Routing:** This process is biologically naive; it simply takes whatever `prefix` channel it is passed by `main.nf` (e.g., `cerevisiae`, `eubayanus_nuclear`) and outputs `<prefix>_cohort.vcf.gz`.

---

## 5) Wiring in `main.nf`

```groovy
// =========================================================================
// MODULE 6: JOINT GENOTYPING
// =========================================================================

// Group Module 5 output by chromosome
ch_for_genomicsdb = HAPLOTYPE_CALLER.out.gvcf
    .map { _strain_id, chrom, gvcf, tbi -> tuple(chrom, gvcf, tbi) }    
    .groupTuple(by: 0)

// 6A: Build or Update the Database
GENOMICSDB_IMPORT(
    ch_for_genomicsdb,
    params.genomicsdb_update_path ?: ""
)

// 6B: Call the final variants using the Database from 6A
GENOTYPE_GVCFS(
    GENOMICSDB_IMPORT.out.db,
    PREPARE_GATK_REF.out.gatk_ref.collect()
)

// 6C: Intelligent Subgenome Merging (Sensu Stricto Universal)
if (params.merge_strategy == 'separate') {
    ch_for_merge = GENOTYPE_GVCFS.out.cohort_vcf
        .map { chrom, vcf, tbi -> 
            // Split the chromosome name by the underscore
            def parts = chrom.tokenize('_')
            def prefix = parts[0]
            def suffix = parts.size() > 1 ? parts[1] : ""
            
            // 1. Map the prefix to the species name
            def species_name = params.species_map.get(prefix, "${prefix}_other")
            def final_prefix = species_name // Default to bundled
            
            // 2. Apply the extrachromosomal isolation logic if requested
            if (params.extrachromosomal_strategy == 'isolated') {
                if (suffix.startsWith('M')) {
                    final_prefix = "${species_name}_mitochondria"
                } else if (suffix.startsWith('P')) {
                    final_prefix = "${species_name}_plasmid"
                } else {
                    final_prefix = "${species_name}_nuclear"
                }
            }
            
            tuple(final_prefix, vcf, tbi)
        }
        .groupTuple(by: 0)
} else {
    // Group all chromosomes together under one name
    ch_for_merge = GENOTYPE_GVCFS.out.cohort_vcf
        .map { _chrom, vcf, tbi -> tuple('sensu_stricto_combined', vcf, tbi) }
        .groupTuple(by: 0)
}

// Pass the intelligently grouped channel to Picard
MERGE_VCFS(ch_for_merge)

```

---

## 6) Configuration

In `nextflow.config`:

```groovy
params {
    // --- Module 6: Joint Genotyping ---
    batch_size             = 50   // Prevents GenomicsDBImport from consuming too much RAM at once
    
    // Master Path for Incremental Compendium Updates
    genomicsdb_update_path = null // Leave null for fresh build. Point to an existing DB parent path to update.
    
    // Merging strategy for joint-called VCFs
    merge_strategy = 'separate'  // Options: 'separate' (subgenome) or 'together' (all chromosomes combined)
    
    // Extrachromosomal strategy (only applies if merge_strategy is 'separate')
    extrachromosomal_strategy = 'bundled' // Options: 'bundled' (all in one) or 'isolated' (splits nuclear/mito/plasmid)

    // Universal Sensu Stricto Dictionary
    species_map = [
        'c': 'cerevisiae',
        'e': 'eubayanus',
        'p': 'paradoxus',
        'm': 'mikatae',
        'k': 'kudriavzevii',
        'u': 'uvarum',
        'a': 'arboricola',
        'j': 'jurei'
    ]
}

```

---

## 7) Output Layout

The database directories are safely persisted to allow future incremental updates, while the final `Final_Merged` folder contains the deliverable datasets.

```text
results/
└── Joint_Genotyping/
    ├── c_I/
    │   ├── genomicsdb_c_I/          <-- Stateful GATK Database (can be incrementally updated)
    │   ├── cohort.c_I.vcf.gz        <-- Raw per-chromosome joint VCF
    │   └── sample_map_c_I.txt
    ├── e_XVI/
    │   └── ...
    └── Final_Merged/
        ├── cerevisiae_cohort.vcf.gz         <-- Final, analysis-ready Master VCF
        └── eubayanus_cohort.vcf.gz          <-- Final, analysis-ready Master VCF

```

---

## 8) Tool Reference

| Tool | Version | Purpose |
| --- | --- | --- |
| GATK4 | 4.4.0.0 | GenomicsDBImport & GenotypeGVCFs |
| Picard | 3.5.0 | MergeVcfs |

## 9) Change Management Rule

Any logic change affecting bash collision-handling (`rm -rf` inside Nextflow blocks), the taxonomy dictionary, merging grouping logic, or GATK DB paths must update:

1. `docs/MODULE_6_JOINT_GENOTYPING.md` (this file)
2. `modules/joint_genotyping.nf`
3. The Module 6 routing block of `main.nf`