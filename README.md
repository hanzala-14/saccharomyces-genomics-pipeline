Comparative Phylogenomics and Ancestry Analysis of Saccharomyces

This repository contains the data, scripts, and documentation for the comparative phylogenomics and ancestral lineage analysis of yeast isolates within the Saccharomyces genus. The primary objective of this project is to optimize our existing genomics workflow, transition legacy bash scripts into a robust, reproducible Nextflow pipeline, and try to integrate evolutionary clock analysis to place historic strains(Lagers) in their proper evolutionary context.
📅 Proposed Project Work Plan
Phase 1: Local Script Collection & Environment Setup

    Audit Existing Code: Gather all disparate bash and R scripts currently used to run the workflow steps locally.
    Define Dependencies: Document all required software versions (e.g., fastp, BWA-MEM, GATK, SNPRelate, SplitsTree) to prepare for containerization.
    Repository Structure: Organize raw scripts, reference genomes, and sample metadata schemas within this repository.

Phase 2: Pipeline Optimization & Nextflow Transition

    Workflow Modularization: Port legacy bash command sequences into distinct Nextflow processes (modules).
    Environment Management: Integrate Conda environments or Docker containers for each process to ensure absolute reproducibility.
    Validation: Test the Nextflow pipeline with Saccharomyces control datasets to verify output consistency against previous manual runs.

Phase 3: Phylogenomic Analysis & Tree Generation

    Run the completed pipeline across the full Saccharomyces isolate dataset.
    Generate final distance matrices and reconstruct Neighbor-Joining trees to visualize genomic relationships within the sensu stricto complex.

Phase 4: Evolutionary Clock Analysis (Future Objective)

    Collaboration: Partner with mentors from the One Health Institute.
    Dating Divergence: Use molecular clock models to estimate divergence dates of key Saccharomyces clades (Lagers).
    Historical Synthesis: Correlate evolutionary milestones with historical context to trace lineage domestication and divergence timelines.

🧬 Current Genomics Workflow Reference

An example of our established bioinformatics workflow consists of the following operational steps:

A[1. Collect FTP links from ENA] --> B[2. Download .fastq files via curl]
B --> C[3. Preprocessing & Quality Control via fastp]
C --> D[4. Reference Mapping to sensu stricto via BWA-MEM]
D --> E[5. Variant Calling via GATK HaplotypeCaller]
E --> F[6. Merge .gvcf files via GATK GenomicsDBImport]
F --> G[7. Joint Genotyping via GATK GenotypeGVCFs]
G --> H[8. Select SNPs via GATK SelectVariants]
H --> I[9. Filter SNPs via GATK VariantFiltration]
I --> J[10. Format .vcf to .gds via SNPRelate]
J --> K[11. Calculate IBS State Pairs via snpgdsIBSNum]
K --> L[12. Distance Matrices per Chromosome via R base]
L --> M[13. Calculate Combined Distance Matrix via R base]
M --> N[14. Neighbor-Joining Tree via SplitsTree v4.19]

Ofcourse, we can change them accordingly. Its just an example of a presumable work plan and an image of what i had in mind for this internship. This is the main workflow we currently work on but it changes according to the project we work on . For example, for Lager yeast - we change the reference and everything accordingly