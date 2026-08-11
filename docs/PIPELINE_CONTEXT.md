Here is the updated pipeline context document reflecting the transition to samtools markdup:Pipeline PurposeWhole-genome sequencing (WGS) variant calling pipeline for Saccharomyces (yeast) species. Designed for population genomics with dozens-to-hundreds of strains, supporting both ENA/SRA public data and local FASTQ files. Runs locally (Ubuntu, VS Code, conda) with an optional Docker/Singularity profile for production.Architecture Overview
Plaintextsamplesheet.csv
       │
       ▼
┌─────────────────────────┐
│  MODULE 1: download_raw │  → Download/ingest → PE/SE classify → merge by strain
└─────────────────────────┘
       │
       ▼  (strain_id, R1, R2) or (strain_id, SE)
┌─────────────────────────┐
│  MODULE 2: filter_reads │  → fastp quality/adapter trimming
└─────────────────────────┘
       │
       ▼  (strain_id, filt.R1, filt.R2) or (strain_id, filt.SE)
┌─────────────────────────┐
│  MODULE 3: map_reads    │  → BWA alignment → markdup → stats
└─────────────────────────┘
       │
       ▼  (strain_id, .mdup.bam, .bam.bai)
┌─────────────────────────┐
│  MODULE 4: (Next)       │  → Variant calling (GATK / bcftools)
└─────────────────────────┘
File StructurePlaintextsaccharomyces-genomics-pipeline/
├── main.nf                        # Entry workflow (wires all modules)
├── nextflow.config                # Params, resource labels, profiles
├── modules/
│   ├── download_raw.nf            # Module 1: SRA fetch, local ingest, QC classify, merge
│   ├── filter_reads.nf            # Module 2: fastp PE/SE filtering
│   └── map_reads.nf               # Module 3: BWA map, markdup, mapping stats
└── bin/
    └── PIPELINE_CONTEXT.md        # Pipeline context documentation
Completed ModulesModule 1: Data Acquisition (modules/download_raw.nf)Processes: FETCH_SRA, INGEST_LOCAL, CLASSIFY_RUN, MERGE_PE, MERGE_SE, CREATE_STRAIN_LIST, WRITE_RUN_QC, WRITE_STRAIN_QCInput: samplesheet.csv with columns: strain_id, accession (ENA/SRA), r1, r2Logic: Download from ENA or ingest local FASTQs → classify each run as PE_PASS / SE_FALLBACK / SE_INPUT / DROP → group by strain → merge multi-run strains → produce run_qc.tsv + strain_qc.tsvChannel shape out:(strain_id, merged_R1.fastq.gz, merged_R2.fastq.gz) for PE(strain_id, merged_SE.fastq.gz) for SEModule 2: Read Filtration (modules/filter_reads.nf)Processes: FILTER_PE, FILTER_SE, WRITE_FILTER_SUMMARYTool: fastp (adapter trimming, quality filtering, poly-G removal)Params: adapter_fasta, trim_front/tail, min_length, qualified_quality, cut_mean_qualityChannel shape out:FILTER_PE.out.filtered_pe → (strain_id, filt.R1.fastq.gz, filt.R2.fastq.gz)FILTER_SE.out.filtered_se → (strain_id, filt.SE.fastq.gz)Outputs: results/filtered/ (trimmed FASTQs, json reports, filter_summary.tsv)Module 3: Mapping & Duplicate Marking (modules/map_reads.nf)Processes: FETCH_REFERENCE, BWA_INDEX, MAP_PE, MAP_SE, MARKDUP, MAPPING_STATS, WRITE_MAPPING_SUMMARYTools: bwa 0.7.19, samtools 1.24Reference strategy:params.reference = <path> local FASTA path (preferred). If .bwt index exists alongside, skips BWA indexing.params.genome_id = <accession> NCBI Assembly accession (e.g., GCF_000146045.2). Auto-fetches + indexes.Key design decisions:Piped Stream: bwa mem | samtools view -F4 | samtools sort eliminates SAM/unsorted BAM intermediate files on disk (~70% disk savings).Unified Output Stream: PE/SE split only occurs at bwa mem call; both merge into a single BAM stream immediately after sorting.Collision-Free Read Groups: Read groups incorporate file names dynamically (ID:${strain_id}_${r1.simpleName}) while preserving sample identity (SM:${strain_id}).Deduplication via Samtools: Replaced Picard with samtools markdup to bypass fatal Java memory crashes (PairInfoMap) caused by malformed/duplicate QNAMEs commonly found in public SRA datasets. The deduplication logic utilizes coordinate collation and mate-score fixing before marking:Bashsamtools collate -o name_collate.bam ${sorted_bam}
samtools fixmate -m name_collate.bam fixmate.bam
samtools sort -o coord_sorted.bam fixmate.bam
samtools markdup -s coord_sorted.bam ${strain_id}.mdup.bam 2> ${strain_id}.mdup.metrics.txt
Deduplication Policy: samtools markdup identifies and flags duplicates (Bitwise flag 0x400) rather than physically removing them from the file.Flagstat Validation: Runs directly on final duplicate-marked BAM files.Insert Size Metrics: Derived from samtools stats (PE strains record average insert size; SE strains yield "NA").Channel shape out:MARKDUP.out.markdup_bam → (strain_id, strain.mdup.bam, strain.mdup.bam.bai)MARKDUP.out.dup_metrics → (strain_id, strain.mdup.metrics.txt)MAPPING_STATS.out.flagstat → (strain_id, strain.flagstat.txt)Outputs layout:results/mapped/mapping_summary.tsv (cohort summary: strain, total reads, mapped reads, pct mapped, dup%, insert size, status)results/mapped/strains/ (final BAMs + .bai indexes)results/mapped/stats/ (flagstat, idxstats, samtools stats, markdup metrics)Configuration DesignParams (nextflow.config)ParamDefaultPurposesamplesheetnull (required)Path to input CSVoutdirresultsOutput root directoryreferencenullLocal reference FASTA pathgenome_idnullNCBI Assembly accession stringsave_referencefalsePublish indexed FASTA to outdirmin_mapping_pct0Mapping rate threshold warning (%)min_reads10000Minimum reads required to retain runmin_read_length50Minimum initial average read lengthmin_length50fastp post-trimming minimum lengthqualified_quality20fastp base quality thresholdResource LabelsLabelCPUsMemoryTypical Usagetiny1512 MBTSV creation, statistical parsersbase12 GBNCBI/SRA downloads, BWA indexinglow48 GBAlignment (bwa mem + sorting stream)medium416 GBsamtools markdup, heavy indexingProfilesconda — Local development environment (Bioconda channel with pinned package versions).docker — Production containerized deployment.singularity — HPC cluster environment.test — CI test suite with minimal datasets.Conda Environment SpecificationBashconda create -n yeast-wgs -c conda-forge -c bioconda \
    nextflow \
    bwa=0.7.19 \
    samtools=1.24 \
    fastp \
    sra-tools \
    pigz \
    python>=3.9
Execution CommandsBash# Executing with pre-indexed local reference FASTA (bypasses BWA_INDEX):
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --reference /path/to/reference.fasta \
    -profile conda \
    -resume

# Executing with NCBI auto-fetch accession:
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --genome_id GCF_000146045.2 \
    -profile conda \
    -resume
Technical Standards & RulesStream-based processing: Intermediate SAM or unsorted BAM files are piped in-memory.Work directory isolation: No inline script rm statements; cleanup handled via nextflow clean -f. (Exception: Internal cleanup of temporary collate/fixmate files during deduplication is permitted to preserve disk space).Channel discipline: PE/SE separation is restricted strictly to tool commands that demand paired vs single inputs, then unified immediately.Assertions: Every process verifies file completeness (non-empty BAM assertions, read counts).Strict Modularization: Isolated .nf modules in ./modules/ directory.Reproducibility: Pinned tool versions across Conda and Container definitions.Roadmap (Module 4+)The pipeline currently yields indexed, duplicate-marked BAM files (.mdup.bam and .mdup.bam.bai). Next logical steps:Module 4: Variant Calling: Cohort joint variant calling using BCFtools (mpileup + call) or GATK (HaplotypeCaller → GenotypeGVCFs).Module 5: VCF Filtering: Quality score calibration, depth masking, and INDEL/SNP filtering.Module 6: Population Genomics: PCA, phylogenetic tree construction, and admixture analysis.