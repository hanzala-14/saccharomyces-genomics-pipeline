/*
 * Module 6: Joint Genotyping
 * Split into three discrete processes for precise debugging, caching, and merging.
 */


// =========================================================================
// STEP 6A: GENOMICSDB IMPORT
// =========================================================================
process GENOMICSDB_IMPORT {
    tag { "${chrom} | DB: ${update_path ? 'INCREMENTAL' : 'FRESH'}" }
    label 'high'
    
    // publishDir with overwrite: true allows Nextflow to safely replace the old DB folder
    publishDir path: { "${params.outdir}/Joint_Genotyping/${chrom}" }, mode: 'symlink', overwrite: true

    errorStrategy { task.exitStatus in [137, 143, 247] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    tuple val(chrom), path(gvcfs), path(tbis)
    val update_path

    output:
    tuple val(chrom), path("genomicsdb_${chrom}"), emit: db
    path("sample_map_${chrom}.txt"), emit: sample_map

    script:
    def avail_mem = (task.memory.toGiga() * 0.8).intValue()
    def db_name = "genomicsdb_${chrom}"
    
    // SOURCE: Where we copy the old database from (e.g., external drive or old results)
    def existing_db = update_path ? "${file(update_path).toAbsolutePath()}/${chrom}/${db_name}" : null  
       
    // DESTINATION: Where Nextflow is going to publish the results
    def target_publish_dir = "${params.outdir}/Joint_Genotyping/${chrom}/${db_name}"

    """
    set -euo pipefail

    echo "=== [GENOMICSDB_IMPORT] Building sample map for ${chrom} ===" >&2
    rm -f sample_map_${chrom}.txt

    for vcf in *.${chrom}.g.vcf.gz; do
        strain=\$(basename "\$vcf" ".${chrom}.g.vcf.gz")
        echo -e "\${strain}\\t\${PWD}/\${vcf}" >> sample_map_${chrom}.txt
    done

    n_samples=\$(wc -l < sample_map_${chrom}.txt)

    if [ "\${n_samples}" -eq 0 ]; then
        echo "ERROR: No GVCFs found for chromosome ${chrom}" >&2
        exit 1
    fi

    # READ-ONLY COPY: Pull the database from the source (external drive or old results)
    if [ -n "${existing_db ?: ''}" ] && [ -d "${existing_db}" ]; then
        echo "=== [GENOMICSDB_IMPORT] MODE: INCREMENTAL ===" >&2
        echo "=== [GENOMICSDB_IMPORT] Existing DB: ${existing_db} ===" >&2
        echo "=== [GENOMICSDB_IMPORT] New GVCFs: \${n_samples} ===" >&2

        cp -r "${existing_db}" "${db_name}"

        IMPORT_FLAG="--genomicsdb-update-workspace-path"

        echo "=== [GENOMICSDB_IMPORT] Using \${IMPORT_FLAG} ===" >&2

    else
        echo "=== [GENOMICSDB_IMPORT] MODE: FRESH ===" >&2
        echo "=== [GENOMICSDB_IMPORT] Existing DB not found: ${existing_db ?: 'none'} ===" >&2
        echo "=== [GENOMICSDB_IMPORT] New GVCFs: \${n_samples} ===" >&2

        IMPORT_FLAG="--genomicsdb-workspace-path"

        echo "=== [GENOMICSDB_IMPORT] Using \${IMPORT_FLAG} ===" >&2
    fi

    echo "=== [GENOMICSDB_IMPORT] Running Import for ${chrom} ===" >&2

    gatk --java-options "-Xmx${avail_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads=2" \\
        GenomicsDBImport \\
        --sample-name-map sample_map_${chrom}.txt \\
        \${IMPORT_FLAG} ${db_name} \\
        --intervals "${chrom}" \\
        --reader-threads ${task.cpus} \\
        --batch-size ${params.batch_size}

    # =========================================================================
    # THE SAFETY FIX: Only delete the database in the LOCAL OUTPUT DIRECTORY 
    # Never touch the source compendium!
    # =========================================================================
    if [ -d "${target_publish_dir}" ]; then
        echo "=== [GENOMICSDB_IMPORT] Clearing local output directory to allow Nextflow overwrite ===" >&2
        rm -rf "${target_publish_dir}"
    fi
    """
}

// =========================================================================
// STEP 6B: GENOTYPE GVCFS
// =========================================================================
process GENOTYPE_GVCFS {
    tag { "Chr: ${chrom}" }
    label 'high'
    
    publishDir path: { "${params.outdir}/Joint_Genotyping/${chrom}" }, mode: 'symlink', pattern: "cohort.*"

    errorStrategy { task.exitStatus in [137, 143, 247] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    tuple val(chrom), path(workspace)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(chrom), path("cohort.${chrom}.vcf.gz"), path("cohort.${chrom}.vcf.gz.tbi"), emit: cohort_vcf

    script:
    def avail_mem = (task.memory.toGiga() * 0.8).intValue()

    """
    set -euo pipefail

    echo "=== [GENOTYPE_GVCFS] Calling variants for ${chrom} ===" >&2

    gatk --java-options "-Xmx${avail_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads=2" \\
        GenotypeGVCFs \\
        -R "${fasta}" \\
        -G StandardAnnotation \\
        -V "gendb://${workspace}" \\
        -O "cohort.${chrom}.vcf.gz"

    if [ ! -s "cohort.${chrom}.vcf.gz" ]; then
        echo "ERROR: Empty VCF output for ${chrom}" >&2
        exit 1
    fi
    """
}

// =========================================================================
// STEP 6C: MERGE VCFS (DYNAMIC)
// =========================================================================
process MERGE_VCFS {
    tag { "Merging: ${prefix}" }
    label 'medium'
    
    publishDir path: "${params.outdir}/Joint_Genotyping/Final_Merged", mode: 'symlink', overwrite: true

    input:
    tuple val(prefix), path(vcfs), path(tbis)

    output:
    tuple val(prefix), path("${prefix}_cohort.vcf.gz"), path("${prefix}_cohort.vcf.gz.tbi"), emit: final_vcf
    script:
    def avail_mem = (task.memory.toGiga() * 0.8).intValue()
    def input_args = vcfs.collect { vcf -> "-I ${vcf}" }.join(" ")

    """
    set -euo pipefail
    echo "=== [MERGE_VCFS] Merging VCFs for ${prefix} ===" >&2
    
    gatk --java-options "-Xmx${avail_mem}g" \\
        MergeVcfs \\
        ${input_args} \\
        -O ${prefix}_cohort.vcf.gz
    """
}