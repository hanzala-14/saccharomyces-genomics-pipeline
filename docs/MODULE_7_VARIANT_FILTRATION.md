
---

# 📖 Module 7: Variant Filtration & Quality Control

* **Repository:** `hanzala-14/saccharomyces-genomics-pipeline`
* **Module Source:** `modules/local/variant_filtration.nf`
* **Workflow Source:** `main.nf`
* **Pipeline Configuration:** `nextflow.config`
* **Pipeline Version:** `2.1.7`
* **Last Updated:** `12 August 2026`

---

## 1) Purpose & Scope

Module 7 acts as the critical bridge between upstream variant calling and downstream evolutionary analysis (phylogenomics, admixture, and functional annotation). Joint-genotyped VCFs are inherently noisy; this module applies strict mathematical thresholds (GATK Hard Filtering) to eliminate false positives, sequencing artifacts, and strand biases.

To maximize biological accuracy, the module strictly isolates SNPs from INDELs (as they possess distinct error profiles), aggressively filters out unmappable repetitive regions (such as telomeres and Ty elements), and standardizes deletion coordinates via left-alignment. The final output is a pristine, analysis-ready master VCF and a survivor QC report.

---

## 2) Design Principles

### 2.1 Parallel Execution Architecture

Traditional bash scripts process SNPs and INDELs sequentially. This module routes the raw cohort VCF into two independent, parallel channels (`FILTER_SNPS` and `FILTER_INDELS`). Nextflow executes these simultaneously, effectively cutting the total compute time for the filtration phase in half.

### 2.2 Config-Driven Thresholds (Zero Hardcoding)

Yeast genomics requires completely different filtration mathematics than human genomics. To ensure pipeline portability, all GATK mathematical thresholds (QD, FS, MQ, etc.) and biological annotation paths are extracted from the `.nf` scripts and injected dynamically via `nextflow.config`.

### 2.3 Dynamic Interval Slicing & Extrachromosomal Safety

Because the pipeline supports pan-genome joint calling across the *Sensu Stricto* complex, the module features a dynamic Groovy routing block. It automatically reads the species from the VCF cohort name, matches it against a config dictionary, and dynamically slices the master intervals file (e.g., extracting only `c_` chromosomes for *cerevisiae*). Furthermore, it includes a safety override (`SKIP_EXTRACHROMOSOMAL_`) to ensure that mitochondrial and plasmid VCFs are processed without crashing GATK.

### 2.4 Granular Storage Optimization

A standard GATK filtration workflow generates massive hard drive bloat. Module 7 utilizes highly granular, conditional Nextflow `publishDir` toggles (`keep_snp_vcf`, `keep_indel_vcf`, `keep_merged_vcf`). Intermediate files are relegated to temporary Nextflow work directories by default, giving the user absolute control over what is permanently saved to the results directory.

### 2.5 Downstream Readiness (Index Tracking)

Downstream R packages (like `SNPRelate`) and annotation tools (`SnpEff`) will critically fail if a VCF lacks a Tabix index. Module 7 natively tracks and emits `.tbi` files alongside every `.vcf.gz` output, guaranteeing perfect compatibility with downstream tools.

---

## 3) Input Contract

| Channel | Tuple Shape | Source |
| --- | --- | --- |
| `MERGE_VCFS.out` | `(cohort_name, vcf, tbi)` | `MERGE_VCFS` (Module 6) |
| `ch_fasta_clean` | `(fasta)` | `PREPARE_GATK_REF.out.gatk_ref` (first) |
| `ch_fai_clean` | `(fai)` | `PREPARE_GATK_REF.out.gatk_ref` (first) |
| `ch_dict_clean` | `(dict)` | `PREPARE_GATK_REF.out.gatk_ref` (first) |
| `ch_intervals` | `(intervals_file)` | `params.repetitive_intervals` |

---

## 4) Process Specification

### 4.1 `FILTER_SNPS`

**Role:** Isolates single nucleotide polymorphisms, dynamically masks repetitive genome intervals, and applies strict SNP-specific threshold math.

* **Label:** `high` (Dynamic Java heap allocation via `-Xmx`)
* **Output:** `${cohort_name}.snpvarfiltered.vcf.gz`

### 4.2 `FILTER_INDELS`

**Role:** Isolates insertions/deletions, dynamically masks repetitive intervals, applies INDEL-specific threshold math (excluding MQ filters), and left-aligns coordinates.

* **Label:** `high`
* **Output:** `${cohort_name}.indelvarfilteredleftal.vcf.gz`

### 4.3 `MERGE_AND_CLEAN`

**Role:** Recombines the filtered SNP and INDEL streams, physically drops any variant that failed the hard filters, and generates a survivor statistical report.

* **Label:** `high`
* **Outputs:**
* `${cohort_name}.final.vcf.gz`
* `${cohort_name}.variant_eval.txt` (Parsed by MultiQC)



---

## 5) Wiring in `main.nf`

```groovy
    // =========================================================================
    // MODULE 7: VARIANT FILTRATION & QC
    // =========================================================================
    
    // Load the intervals file (with a safety check so it crashes instantly if missing)
    ch_intervals = file(params.repetitive_intervals, checkIfExists: true)

    // Extract the individual reference files from Module 5's prepared GATK bundle
    // We use .first() so Nextflow knows it can reuse these files endlessly
    // Note: Linter warnings suppressed using '_' prefix for unused variables
    ch_fasta_clean = PREPARE_GATK_REF.out.gatk_ref.map { fasta, _fai, _dict -> fasta }.first()
    ch_fai_clean   = PREPARE_GATK_REF.out.gatk_ref.map { _fasta, fai, _dict -> fai }.first()
    ch_dict_clean  = PREPARE_GATK_REF.out.gatk_ref.map { _fasta, _fai, dict -> dict }.first()

    // 1. Launch SNP and INDEL filtering simultaneously on the merged cohort VCFs
    FILTER_SNPS(
        MERGE_VCFS.out,
        ch_fasta_clean,
        ch_fai_clean,
        ch_dict_clean,
        ch_intervals
    )

    FILTER_INDELS(
        MERGE_VCFS.out,
        ch_fasta_clean,
        ch_fai_clean,
        ch_dict_clean,
        ch_intervals
    )

    // 2. The Join Operator: Wait for both to finish, then match them by cohort name (by: 0)
    ch_filtered_joined = FILTER_SNPS.out.filtered_snps
        .join(FILTER_INDELS.out.filtered_indels, by: 0)

    // 3. Merge and clean the surviving variants
    MERGE_AND_CLEAN(
        ch_filtered_joined,
        ch_fasta_clean,
        ch_fai_clean,
        ch_dict_clean
    )

```

---

## 6) Configuration

In `nextflow.config`:

```groovy
    // ========================================================================
    // MODULE 7: VARIANT FILTRATION PARAMETERS
    // ========================================================================
    
    // 1. Output Storage Toggles (What to save to the hard drive)
    keep_snp_vcf         = true    // Save the SNPs-only filtered VCF?
    keep_indel_vcf       = true    // Save the INDELs-only left-aligned VCF?
    keep_merged_vcf      = true    // Save the final recombined (SNP+INDEL) VCF?

    // 2. Downstream Analysis Target (What to feed to Module 8)
    downstream_vcf       = 'merged' // Options: 'snps', 'indels', or 'merged'
    
    // 3. Biological Annotations
    repetitive_intervals = "${projectDir}/data/references/sensustricto8.repetitive.intervals"

    // SNP Hard-Filtering Thresholds (GATK Best Practices)
    snp_qd_filter        = 'QD < 5.0'
    snp_qual_filter      = 'QUAL < 30.0'
    snp_sor_filter       = 'SOR > 3.0'
    snp_fs_filter        = 'FS > 60.0'
    snp_mq_filter        = 'MQ < 40.0'
    snp_mqranksum_filter = 'MQRankSum < -12.5'
    snp_readpos_filter   = 'ReadPosRankSum < -8.0'

    // INDEL Hard-Filtering Thresholds (GATK Best Practices)
    indel_qd_filter      = 'QD < 5.0'
    indel_qual_filter    = 'QUAL < 30.0'
    indel_fs_filter      = 'FS > 60.0'
    indel_readpos_filter = 'ReadPosRankSum < -20.0'

```

---

## 7) Output Layout

The output directory responds dynamically to the granular config toggles.

```text
results/
└── Variant_Filtration/
    ├── SNPs/                                 <-- Enabled by params.keep_snp_vcf
    │   ├── cerevisiae_cohort.snpvarfiltered.vcf.gz
    │   └── cerevisiae_cohort.snpvarfiltered.vcf.gz.tbi
    ├── INDELs/                               <-- Enabled by params.keep_indel_vcf
    │   ├── cerevisiae_cohort.indelvarfilteredleftal.vcf.gz
    │   └── cerevisiae_cohort.indelvarfilteredleftal.vcf.gz.tbi
    └── Final_Merged/                         <-- Enabled by params.keep_merged_vcf
        ├── cerevisiae_cohort.final.vcf.gz
        ├── cerevisiae_cohort.final.vcf.gz.tbi
        └── cerevisiae_cohort.variant_eval.txt <-- QC report for MultiQC

```

---

## 8) Tool Reference

| Tool | Version | Purpose |
| --- | --- | --- |
| GATK4 | 4.4.0.0 | SelectVariants, VariantFiltration, LeftAlignAndTrimVariants, MergeVcfs, VariantEval |

## 9) Change Management Rule

Any logic change affecting the GATK hard-filtering thresholds, Java heap memory allocation limits, Tabix tracking, or `.join()` merging logic must update:

1. `docs/MODULE_7_VARIANT_FILTRATION.md` (this file)
2. `modules/local/variant_filtration.nf`
3. The Module 7 routing block of `main.nf`