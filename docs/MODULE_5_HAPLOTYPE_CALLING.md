
---

# 📖 Module 5: Variant Calling (GATK HaplotypeCaller)

* **Repository:** `hanzala-14/saccharomyces-genomics-pipeline`
* **Module Source:** `modules/haplotype_calling.nf`
* **Workflow Source:** `main.nf`
* **Pipeline Configuration:** `nextflow.config`
* **Pipeline Version:** `2.1.5`
* **Last Updated:** `07 August 2026`

---

## 1) Purpose & Scope

Module 5 transitions the pipeline from read alignment into genomic variant detection. Utilizing **GATK HaplotypeCaller**, this module generates per-strain Genomic VCFs (GVCFs). Running in GVCF mode (`-ERC GVCF`) is a mandatory prerequisite for downstream population-level joint genotyping (Module 6), as it preserves confidence scores for both variant and non-variant reference sites.

---

## 2) Design Principles

### 2.1 Scatter-Gather Architecture

To bypass the extreme computational bottleneck of whole-genome variant calling, this module employs a highly parallelized scatter-gather strategy:

* **Scatter:** Each incoming strain BAM is mathematically combined with a dynamically generated list of target chromosomes, spawning independent jobs per chromosome.
* **Gather:** Upon completion, a dedicated aggregation process collects the per-chromosome GVCFs and automatically generates a sample manifest (`gvcf_list.txt`) required for Module 6 `GenomicsDBImport`.

### 2.2 Dynamic Reference Agnosticism

The target genomic intervals are fully decoupled from the pipeline code. Instead of relying on hardcoded lists in the configuration, the pipeline natively parses the FASTA Index (`.fai`) on the fly to auto-detect all contigs. The pipeline effortlessly scales from a standard 16-chromosome *Saccharomyces cerevisiae* run to a 130+ chromosome *Sensu Stricto* complex run without modifying a single line of Nextflow logic.

### 2.3 Biological Ploidy Calibration (Organelle Handling)

Mitochondria (`_M`) and plasmids (`_P`) do not follow diploid Mendelian inheritance. Running GATK with a global `-ploidy 2` introduces mathematical errors for these haploid organelles.

* **Master Switch:** An `include_extrachromosomal` toggle allows researchers to drop organelles from the pipeline entirely before variant calling.
* **Dynamic Ploidy:** For organelles that pass through, a config-driven `ploidy_map` identifies chromosome suffixes and forces GATK into haploid calling mode (`-ploidy 1`), preserving statistical integrity while defaulting to diploid for nuclear chromosomes.

### 2.4 Just-In-Time Reference Formatting

Because Module 3 utilizes a pure-C `samtools` duplicate marking workflow (bypassing Picard to save memory), the required Picard Sequence Dictionary (`.dict`) is absent. Module 5 introduces a lightweight, just-in-time preparation step (`PREPARE_GATK_REF`) that leverages the existing `medium` Conda environment to generate the `.dict` and `.fai` files instantly before execution.

---

## 3) Input Contract

### 3.1 Flow from Module 3 & Config

| Channel | Tuple Shape | Source |
| --- | --- | --- |
| `hc_input` | `(strain_id, mdup.bam, mdup.bam.bai, chrom)` | `MARKDUP` × `chromosomes_ch` |
| `gatk_fasta_ch` | `fasta` | Extracted from `indexed_reference` |
| `ploidy_map` | `val(map)` | Passed directly from `nextflow.config` |

---

## 4) Process Specification

### 4.1 `PREPARE_GATK_REF`

**Role:** Generates the precise reference indices required by GATK.

* **Label:** `medium` (Reuses the Picard/Samtools Conda environment from Module 3)
* **Output:** `(fasta, fasta.fai, fasta.dict)`

### 4.2 `HAPLOTYPE_CALLER`

**Role:** Executes GATK HaplotypeCaller on a single strain-chromosome pair.

* **Label:** `calling`
* **Container:** `biocontainers/gatk4:4.4.0.0`
* **Dynamic Logic:** Parses the chromosome suffix (e.g., `_M`) and looks it up in `ploidy_map` to assign `-ploidy 1` or `-ploidy 2`.
* **Memory Allocation:** Dynamically parses `${task.memory.toGiga()}` to allocate 80% of cluster-provided RAM to the Java Heap (`-Xmx`), preventing container OOM kills.
* **Output:** A block-gzipped `.g.vcf.gz` and its `.tbi` index per chromosome.

### 4.3 `GATHER_STRAIN_GVCFS`

**Role:** Collects all scattered GVCFs for a given strain and authors the downstream manifest.

* **Label:** `tiny`
* **Output:** Aggregated GVCFs/TBIs and `${strain_id}.gvcf_list.txt`.
* **Value Add:** Completely automates sample map creation, eliminating manual text-wrangling before Joint Genotyping.

---

## 5) Wiring in `main.nf`

```groovy
include { 
    PREPARE_GATK_REF
    HAPLOTYPE_CALLER
    GATHER_STRAIN_GVCFS 
} from './modules/haplotype_calling.nf'

// =========================================================================
// MODULE 5: GATK HaplotypeCaller — Per-chromosome GVCF Calling
// =========================================================================

// Extract just the FASTA from Module 3's reference channel
gatk_fasta_ch = indexed_reference.map { fasta, _idx -> fasta }

// Prepare the exact .fai and .dict bundle GATK needs
PREPARE_GATK_REF(gatk_fasta_ch)

// Dynamically extract all chromosome names directly from the FASTA index (.fai)
// and filter based on the extrachromosomal master switch
chromosomes_ch = PREPARE_GATK_REF.out.gatk_ref
    .map { _fasta, fai, _dict -> fai }
    .splitCsv(sep: '\t')
    .map { row -> row[0] }
    .filter { chrom -> 
        if (!params.include_extrachromosomal) {
            // Drop any chromosome ending in _M or _P if the switch is false
            return !(chrom.endsWith('_M') || chrom.endsWith('_P'))
        }
        return true
    }

// Scatter: 1 BAM × N Chromosomes 
hc_input = MARKDUP.out.markdup_bam.combine(chromosomes_ch)

// Execute variant calling utilizing the NEW GATK reference bundle and dynamic ploidy map
HAPLOTYPE_CALLER(
    hc_input, 
    PREPARE_GATK_REF.out.gatk_ref.collect(),
    params.ploidy_map
)

// Gather: Group outputs by strain_id and generate manifest
GATHER_STRAIN_GVCFS(HAPLOTYPE_CALLER.out.gvcf.groupTuple(by: 0))

```

---

## 6) Configuration

In `nextflow.config`:

```groovy
params {
    // --- Module 5: Haplotype Calling ---
    
    // Master switch to include/exclude mitochondria and plasmids
    include_extrachromosomal = true 
    
    ploidy = 2 // Default nuclear ploidy
    
    // Ploidy map: Suffix (M, P) -> Ploidy (1). 
    // Anything not in this list defaults to params.ploidy (2).
    ploidy_map = [
        'M': 1,
        'P': 1
    ]
}

process {
    withLabel: 'calling' {
        cpus   = 2
        memory = 6.GB
        time   = 4.h
    }
    withName: 'HAPLOTYPE_CALLER' {
        maxForks = 2 // Safely throttles concurrency for local development
    }
}

```

---

## 7) Output Layout

```text
results/
└── Haplotype_Calling/
    ├── StrainA/
    │   ├── StrainA.c_I.g.vcf.gz
    │   ├── StrainA.c_I.g.vcf.gz.tbi
    │   ├── ...
    │   ├── StrainA.c_M.g.vcf.gz        <-- Retained and haploid-called if enabled
    │   ├── StrainA.e_XVI.g.vcf.gz
    │   └── StrainA.gvcf_list.txt       <-- Auto-generated manifest
    └── StrainB/
        └── ...

```

---

## 8) Tool Reference

| Tool | Version | Purpose |
| --- | --- | --- |
| Picard | 3.5.0 | Sequence Dictionary Generation |
| GATK4 | 4.4.0.0 | HaplotypeCaller (GVCF mode) |

## 9) Change Management Rule

Any logic change affecting Java memory tuning, dynamic ploidy maps, `.fai` extraction filtering, or GATK container versions must update:

1. `docs/MODULE_5_VARIANT_CALLING.md` (this file)
2. `modules/haplotype_calling.nf`
3. The Module 5 section of `main.nf`
4. The parameters block of `nextflow.config`