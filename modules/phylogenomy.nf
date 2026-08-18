/*
 * Module 8: IBSx Phylogeny (Part 1 - Allele Preparation)
 * -------------------------------------------------------------------
 *  PREP_ALLELES with:
 *   - Primary Strategy: Dynamic VCF slicing natively via bcftools to isolate chromosomes for granular phylogenetic analysis
 *   - Config-Driven Architecture: Execution modes ('chromosomal' vs 'combined') and monomorphic site filtering mapped directly from nextflow.config (zero hardcoding)
 *   - Biological Accuracy: Transforms standard VCF genotypes into split Allele 1 and Allele 2 matrices required for reticulate evolution network building
 *   - Storage Optimization: Massive intermediate chromosomal VCF slices remain in temporary work directories; only the final analysis-ready TSV tables are published
 *   - QC & Traceability: Automated generation of `.meta.json` checksums tracking sample counts and informative sites per slice for rigorous data provenance
 */
process PREP_ALLELES {
    tag { "${strain_id} - Allele Prep" }
    label 'low'
    
    // Publishes the tables cleanly categorized by strain cohort (if enabled in config)
    publishDir path: { "${params.outdir}/Phylogeny/Allele_Tables/${strain_id}" }, mode: 'copy', enabled: params.keep_allele_tables

    input:
    tuple val(strain_id), path(master_vcf), path(master_tbi)
    
    output:
    tuple val(strain_id), path("tables/*/*.samples.txt"), path("tables/*/*.allele1.tsv"), path("tables/*/*.allele2.tsv"), emit: allele_tables
    path "tables/*/*.meta.json", emit: metadata

    script:
    // Convert boolean config to 1/0 for AWK
    def filter_mono = params.filter_monomorphic ? '1' : '0'

    """
    set -euo pipefail

    mkdir -p sliced_vcfs
    mkdir -p tables

    # ==============================================================================
    # 1. DYNAMIC VCF SLICING (WITH FAIL-SAFES)
    # ==============================================================================
    if [ "${params.ibs_strategy}" == "chromosomal" ]; then
        echo "Slicing master VCF into individual chromosomes..."
        
        # Primary attempt: Extract contigs from header
        bcftools view -h ${master_vcf} | grep "^##contig=" | sed 's/.*ID=//' | cut -d ',' -f 1 | tr -d '>' > contigs.txt
        
        # FAIL-SAFE 1: If header is malformed, extract directly from data column
        if [ ! -s contigs.txt ]; then
            echo "Warning: No contig headers found. Extracting directly from variant data..."
            bcftools query -f '%CHROM\\n' ${master_vcf} | sort -u > contigs.txt
        fi
        
        while read chrom; do
            # FAIL-SAFE 2: Sanitize weird characters out of chromosome names for Linux compatibility
            safe_chrom=\$(echo "\$chrom" | sed 's/[^a-zA-Z0-9_]/_/g')
            bcftools view -r "\$chrom" -O z -o "sliced_vcfs/${strain_id}_\${safe_chrom}.vcf.gz" ${master_vcf} || true
        done < contigs.txt
    else
        echo "Running in combined mode..."
        cp ${master_vcf} sliced_vcfs/${strain_id}_combined.vcf.gz
    fi

    # ==============================================================================
    # 2. VECTORIZED AWK EXTRACTION (WITH DOWNSTREAM SHIELD)
    # ==============================================================================
    
    # FAIL-SAFE 3: Prevent literal '*.vcf.gz' expansion if no files exist
    shopt -s nullglob 
    
    for vcf in sliced_vcfs/*.vcf.gz; do
        
        base_name=\$(basename "\$vcf" .vcf.gz)
        
        # CREATE A DEDICATED SUB-FOLDER FOR THIS CHROMOSOME
        mkdir -p "tables/\${base_name}"
        
        bcftools query -l "\$vcf" > "tables/\${base_name}/\${base_name}.samples.txt"
        N_SAMPLES=\$(wc -l < "tables/\${base_name}/\${base_name}.samples.txt")

        # Pre-create files to prevent 0-variant crashes
        touch "tables/\${base_name}/\${base_name}.allele1.tsv"
        touch "tables/\${base_name}/\${base_name}.allele2.tsv"

        # Core Genotype split logic
        bcftools query -H -f '[\\t%GT]\\n' "\$vcf" \\
        | awk -v FS="\\t" -v OFS="\\t" \\
              -v filter_mono="${filter_mono}" \\
              -v allele1_out="tables/\${base_name}/\${base_name}.allele1.tsv" \\
              -v allele2_out="tables/\${base_name}/\${base_name}.allele2.tsv" \\
        '
        function norm_gt(gt,   a,b,t,n,tmp) {
            gsub(/\\|/, "/", gt)
            if (gt == "." || gt == "./." || gt == ".|.") return ".."
            n = split(gt, t, "/")
            if (n == 1) {
                if (t[1] == ".") return ".."
                return t[1] t[1]
            }
            a = t[1]; b = t[2]
            if (a == "") a = "."
            if (b == "") b = "."
            if (a == "." && b != ".") return "." b
            if (b == "." && a != ".") return "." a
            if (a == "." && b == ".") return ".."
            if (a > b) { tmp=a; a=b; b=tmp }
            return a b
        }
        NR == 1 { next }
        {
            for (i = 2; i <= NF; i++) \$i = norm_gt(\$i)
            if (filter_mono == 1) {
                ref = \$2
                same = 1
                for (i = 3; i <= NF; i++) if (\$i != ref) { same = 0; break }
                if (same == 1) next
            }
            out1 = substr(\$2, 1, 1)
            out2 = substr(\$2, 2, 1)
            for (i = 3; i <= NF; i++) {
                out1 = out1 OFS substr(\$i, 1, 1)
                out2 = out2 OFS substr(\$i, 2, 1)
            }
            print out1 > allele1_out
            print out2 > allele2_out
        }
        '

        # Count the informative sites
        N_SITES_A1=\$(wc -l < "tables/\${base_name}/\${base_name}.allele1.tsv" || echo 0)
        
        # FAIL-SAFE 4: THE DOWNSTREAM R SHIELD
        if [ "\$N_SITES_A1" -eq 0 ]; then
            echo "WARNING: \${base_name} has 0 informative variants. Deleting to prevent R crash."
            rm -rf "tables/\${base_name}"
            continue 
        fi

        # Generate JSON Metadata
        MD5_A1=\$(md5sum "tables/\${base_name}/\${base_name}.allele1.tsv" | awk '{print \$1}')
        
        cat > "tables/\${base_name}/\${base_name}.meta.json" <<EOF
        {
            "slice_id": "\${base_name}",
            "n_samples": \$N_SAMPLES,
            "n_sites_extracted": \$N_SITES_A1,
            "md5_checksum": "\$MD5_A1"
        }
EOF
    done
    """
}

/*
 * Module 8: IBSx Phylogeny (Part 2 - Distance Matrices)
 * -------------------------------------------------------------------
 *  RUN_IBS with:
 *   - CLI-Driven Architecture: Replaced legacy hardcoded R configs with a robust command-line interface.
 *   - Physical Staging (stageInMode): Forces Nextflow to physically copy the scripts, destroying Docker symlink walls.
 *   - Staging Collision Bypass: Renamed the input channel to 'custom_scripts' to bypass Nextflow's strict bin/ executable filter.
 */
process RUN_IBS {
    tag { "${strain_id} - IBS Matrix" }
    label 'phylogeny' 
    
    // THE SILVER BULLET: Copies files instead of symlinking them for Docker
    stageInMode 'copy'

    // Publish only the final matrices
    publishDir path: { "${params.outdir}/Phylogeny/Distance_Matrices/${strain_id}" }, mode: 'copy', pattern: 'ibs_outputs/*'

    input:
    tuple val(strain_id), path(samples), path(allele1), path(allele2)
    path custom_scripts // THE FIX: Safely stages all .cpp and .R files without filtering
    
    output:
    tuple val(strain_id), path("ibs_outputs/*"), emit: distance_matrices

    script:
    """
    set -euo pipefail

    # -------------------------------------------------------------------
    # 1. ENVIRONMENT & STAGING
    # -------------------------------------------------------------------
    # CRITICAL: Force OpenMP C++ engines to respect Nextflow's CPU allocation!
    # Without this, OpenMP will hijack every core on the HPC node.
    export OMP_NUM_THREADS=${task.cpus}

    mkdir -p ibs_outputs
    mkdir -p raw_tables
    
    mv *.allele1.tsv raw_tables/
    mv *.allele2.tsv raw_tables/
    
    MASTER_SAMPLE=\$(ls *.samples.txt | head -n 1)

    # -------------------------------------------------------------------
    # 2. EXECUTION
    # -------------------------------------------------------------------
    Rscript ${custom_scripts}/run_ibs.R \\
        "raw_tables" \\
        "\${MASTER_SAMPLE}" \\
        ${task.cpus} \\
        ${params.ibs_block_rows} \\
        ${params.ibs_report_blocks} \\
        "${custom_scripts}"

    # -------------------------------------------------------------------
    # 3. POST-RUN SANITY CHECKS
    # -------------------------------------------------------------------
    if [ -z "\$(ls -A ibs_outputs)" ]; then
        echo "💥 FATAL: ibs_outputs/ is empty! The C++ matrix engine failed silently."
        exit 1
    fi
    
    find ibs_outputs -type f -name "*matrix*" -size 0 -exec echo "💥 FATAL: Empty matrix detected!" \\; -exec false {} +

    # -------------------------------------------------------------------
    # 4. PROVENANCE (Version Logging)
    # -------------------------------------------------------------------
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | awk '{print \$5}')
    END_VERSIONS
    """
}

