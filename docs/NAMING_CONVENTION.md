
---

# 📖 Pipeline Input Contracts & Naming Conventions

## Purpose

This document defines the naming conventions, input schema rules, taxonomy routing dictionaries, and database architecture required to successfully run the *Saccharomyces* Genomics Pipeline (Modules 1 through 6).

---

## 1. Samplesheet Schema (Module 1)

Required header:

```csv
strain_id,accession,r1,r2

```

### Row modes

* **ENA mode:** `accession` set, `r1/r2` empty
* **Local mode:** `r1` set, `r2` optional, `accession` empty

### Invalid rows

* Both accession and r1 missing
* Accession and local paths provided together in the same row

### Path resolution

Local paths are resolved at channel creation time:

* Absolute paths (e.g., `/home/user/data/sample_R1.fastq.gz`) are used as-is.
* Relative paths are resolved from the Nextflow launch directory.
* Missing or empty files cause an immediate pipeline abort with a clear error message.

---

## 2. Strain & Run Naming Conventions (Modules 1-4)

### `strain_id`

* **Required.**
* Sanitized to safe characters: `[A-Za-z0-9_.-]`
* Underscores replace any invalid characters to prevent downstream bash errors.

### `run_id`

* **ENA:** Accession ID (e.g., `SRR10047173`).
* **Local:** `LOCAL_<strain_id>_<r1_basename>` (sanitized).
* *Example:* File `WLP830.masodikkk.R1.fastq.gz` → `run_id = LOCAL_WLP830_WLP830.masodikkk.R1.fastq`
* The basename (without `.gz`) guarantees uniqueness per file.



### Merge Behavior

When a strain has multiple runs (from any source combination):

* All `PE_PASS` runs are concatenated into a single `<strain_id>_PE_R1.fastq.gz` / `<strain_id>_PE_R2.fastq.gz`.
* All `SE_INPUT` / `SE_FALLBACK` runs are concatenated into a single `<strain_id>_SE.fastq.gz`.
* Concatenation order is deterministic: sorted by `run_id` alphabetically.

---

## 3. Chromosome & Taxonomy Routing (Modules 5-6)

To support universal pan-genomic complexes (e.g., 8-species *Sensu Stricto* hybrid references) without hardcoding logic, the reference genome FASTA headers **must** follow a strict `prefix_suffix` naming convention.

### 3.1 The Prefix: Species Taxonomy

The string before the first underscore (`_`) dictates the biological species. The pipeline checks this prefix against `params.species_map` in `nextflow.config` to dynamically route merged VCFs to the correct subgenome cohort.

| Prefix | Species Mapping | Output VCF Prefix |
| --- | --- | --- |
| `c` | *cerevisiae* | `cerevisiae_cohort.vcf.gz` |
| `e` | *eubayanus* | `eubayanus_cohort.vcf.gz` |
| `p` | *paradoxus* | `paradoxus_cohort.vcf.gz` |
| `m` | *mikatae* | `mikatae_cohort.vcf.gz` |
| `k` | *kudriavzevii* | `kudriavzevii_cohort.vcf.gz` |
| `u` | *uvarum* | `uvarum_cohort.vcf.gz` |
| `a` | *arboricola* | `arboricola_cohort.vcf.gz` |
| `j` | *jurei* | `jurei_cohort.vcf.gz` |

*(Note: Unrecognized prefixes will default to `<prefix>_other_cohort.vcf.gz`)*

### 3.2 The Suffix: Organelle & Ploidy Detection

The string after the underscore dictates ploidy calling and optional extrachromosomal isolation.

* **Nuclear Chromosomes (e.g., `c_I`, `e_XVI`):** Default to `-ploidy 2`.
* **Mitochondria (`_M`, e.g., `c_M`):** Caught by `params.ploidy_map`. Forced to `-ploidy 1`.
* **Plasmids (`_P`, e.g., `c_P`):** Caught by `params.ploidy_map`. Forced to `-ploidy 1`.

If `params.extrachromosomal_strategy = 'isolated'` is set, the pipeline will output `cerevisiae_nuclear_cohort.vcf.gz` and `cerevisiae_mitochondria_cohort.vcf.gz` separately.

---

## 4. Compendium Database Architecture (Incremental Updates)

When pointing `params.genomicsdb_update_path` to an existing database (e.g., a 5,000-strain compendium on an external hard drive), the pipeline expects a highly specific folder architecture.

If the provided directory does not match this internal layout, the pipeline will assume the database does not exist and will crash or start a fresh build.

### 4.1 Required Directory Structure

The path you provide in the configuration flag must be the **parent directory** containing the individual chromosome folders.

If your execution command is:
`--genomicsdb_update_path "/media/Office_Drive/Yeast_Compendium/Joint_Genotyping"`

The pipeline will look for the database by automatically appending `/${chrom}/genomicsdb_${chrom}`. Therefore, the external drive **must** be structured exactly like this:

```text
/media/Office_Drive/Yeast_Compendium/Joint_Genotyping/
├── c_I/
│   └── genomicsdb_c_I/          <-- The actual GATK stateful database directory
├── c_II/
│   └── genomicsdb_c_II/
├── ...
└── e_XVI/
    └── genomicsdb_e_XVI/

```

### 4.2 Critical Rules for Compendium Updating

1. **Never rename the `genomicsdb_${chrom}` folders.** The naming convention is hardcoded into GATK's workspace schema.
2. **Do not point the path directly at a `.vcf` file.** You cannot incrementally update a flat text VCF file; you must point to the parent directory containing the `genomicsdb_` workspaces.
3. **Atomic Safety:** The pipeline copies these external folders to a local temporary working directory before updating them. Your external drive will not be corrupted if the pipeline crashes midway.