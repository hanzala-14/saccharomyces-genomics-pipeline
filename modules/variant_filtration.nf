/*
 * Module 7: Variant Filtration & Quality Control
 * -------------------------------------------------------------------
 *  FILTER_SNPS & FILTER_INDELS with:
 *   - Primary Strategy: Independent, parallel processing of SNPs and INDELs to cut execution time in half
 *   - Config-Driven Architecture: GATK hard-filtering thresholds mapped directly from nextflow.config (zero hardcoding)
 *   - Biological Accuracy: Early exclusion of repetitive regions and automated left-alignment of INDELs to standardize coordinates
 *   - Storage Optimization: Granular `publishDir` toggles (`keep_snp_vcf`, etc.) for precise storage control
 *   - Resource Efficiency: Dynamic Java heap memory injection (`-Xmx`) mapped to Nextflow task limits to prevent memory leaks
 * 
 *  MERGE_AND_CLEAN with:
 *   - Strict Integrity: Native GATK execution bypassing legacy Java `.jar` dependencies
 *   - Downstream Readiness: Automated Tabix (`.tbi`) index tracking to ensure compatibility with R and SnpEff
 *   - QC Tracking: Automated survivor variant counting via GATK `VariantEval` for MultiQC integration
 */

process FILTER_SNPS {
    tag "${strain_id} - SNPs"
    label 'high'
    
    publishDir "${params.outdir}/Variant_Filtration/SNPs", mode: 'symlink', enabled: params.keep_snp_vcf

    input:
    tuple val(strain_id), path(raw_vcf), path(raw_vcf_index)
    path fasta
    path fasta_fai
    path fasta_dict
    path intervals

    output:
    tuple val(strain_id), path("${strain_id}.snpvarfiltered.vcf.gz"), path("${strain_id}.snpvarfiltered.vcf.gz.tbi"), emit: filtered_snps

    script:
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : "-Xmx4g"
    
    // 0. Config-Driven Prefix Resolution (Groovy)
    // Matches the cohort name (e.g., 'cerevisiae_nuclear') to the config dictionary to get the 'c_' prefix
    def prefix = ""
    params.species_map.each { key, species ->
        if (strain_id.startsWith(species)) {
            prefix = "${key}_"
        }
    }
    
    // Safety override: Extrachromosomal VCFs do not have nuclear repetitive elements.
    // If we pass nuclear intervals to a mitochondrial/plasmid VCF, GATK will crash.
    if (strain_id.contains('mitochondria') || strain_id.contains('plasmid')) {
        prefix = "SKIP_EXTRACHROMOSOMAL_"
    }

    """
    # Extract only the intervals for this specific species using the dynamic prefix
    # If the prefix is empty or not found, it won't crash, it just creates an empty file
    grep "^${prefix}" ${intervals} > cohort_specific.intervals || true
    
    # Safety Check: If the file has contents, use it. If empty, skip the -XL flag.
    if [ -s cohort_specific.intervals ]; then
        XL_FLAG="-XL cohort_specific.intervals"
    else
        XL_FLAG=""
    fi

    # 1. Isolate SNPs and exclude repetitive regions
    gatk --java-options "${avail_mem}" SelectVariants \
        -R ${fasta} \
        -V ${raw_vcf} \
        -select-type SNP \
        \$XL_FLAG \
        -O ${strain_id}.snp.vcf.gz

    # 2. Apply Config-Driven SNP Hard Filters
    gatk --java-options "${avail_mem}" VariantFiltration \
        -R ${fasta} \
        -V ${strain_id}.snp.vcf.gz \
        -filter "${params.snp_qd_filter}" --filter-name "QD5" \
        -filter "${params.snp_qual_filter}" --filter-name "QUAL30" \
        -filter "${params.snp_sor_filter}" --filter-name "SOR3" \
        -filter "${params.snp_fs_filter}" --filter-name "FS60" \
        -filter "${params.snp_mq_filter}" --filter-name "MQ40" \
        -filter "${params.snp_mqranksum_filter}" --filter-name "MQRankSum-12.5" \
        -filter "${params.snp_readpos_filter}" --filter-name "ReadPosRankSum-8" \
        -O ${strain_id}.snpvarfiltered.vcf.gz
    """
}

process FILTER_INDELS {
    tag "${strain_id} - INDELs"
    label 'high'
    
    publishDir "${params.outdir}/Variant_Filtration/INDELs", mode: 'symlink', enabled: params.keep_indel_vcf
    
    input:
    tuple val(strain_id), path(raw_vcf), path(raw_vcf_index)
    path fasta
    path fasta_fai
    path fasta_dict
    path intervals

    output:
    tuple val(strain_id), path("${strain_id}.indelvarfilteredleftal.vcf.gz"), path("${strain_id}.indelvarfilteredleftal.vcf.gz.tbi"), emit: filtered_indels

    script:
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : "-Xmx4g"
    
    // 0. Config-Driven Prefix Resolution (Groovy)
    def prefix = ""
    params.species_map.each { key, species ->
        if (strain_id.startsWith(species)) {
            prefix = "${key}_"
        }
    }
    
    // Safety override: Extrachromosomal VCFs do not have nuclear repetitive elements.
    // If we pass nuclear intervals to a mitochondrial/plasmid VCF, GATK will crash.
    if (strain_id.contains('mitochondria') || strain_id.contains('plasmid')) {
        prefix = "SKIP_EXTRACHROMOSOMAL_"
    }

    """
    # Extract only the intervals for this specific species using the dynamic prefix
    grep "^${prefix}" ${intervals} > cohort_specific.intervals || true
    
    # Safety Check: If the file has contents, use it. If empty, skip the -XL flag.
    if [ -s cohort_specific.intervals ]; then
        XL_FLAG="-XL cohort_specific.intervals"
    else
        XL_FLAG=""
    fi

    # 1. Isolate INDELs and exclude repetitive regions
    gatk --java-options "${avail_mem}" SelectVariants \
        -R ${fasta} \
        -V ${raw_vcf} \
        -select-type INDEL \
        \$XL_FLAG \
        -O ${strain_id}.indel.vcf.gz
 
    # 2. Apply Config-Driven INDEL Hard Filters
    gatk --java-options "${avail_mem}" VariantFiltration \
        -R ${fasta} \
        -V ${strain_id}.indel.vcf.gz \
        -filter "${params.indel_qd_filter}" --filter-name "QD5" \
        -filter "${params.indel_qual_filter}" --filter-name "QUAL30" \
        -filter "${params.indel_fs_filter}" --filter-name "FS60" \
        -filter "${params.indel_readpos_filter}" --filter-name "ReadPosRankSum-20" \
        -O ${strain_id}.indelvarfiltered.vcf.gz
 
    # 3. Left-Align and Trim
    gatk --java-options "${avail_mem}" LeftAlignAndTrimVariants \
        -R ${fasta} \
        -V ${strain_id}.indelvarfiltered.vcf.gz \
        --dont-trim-alleles \
        -O ${strain_id}.indelvarfilteredleftal.vcf.gz
    """
}

process MERGE_AND_CLEAN {
    tag "${strain_id} - Merge & QC"
    label 'high'  
    
    publishDir "${params.outdir}/Variant_Filtration/Final_Merged", mode: 'symlink', enabled: params.keep_merged_vcf

    input:
    tuple val(strain_id), path(snp_vcf), path(snp_tbi), path(indel_vcf), path(indel_tbi)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(strain_id), path("${strain_id}.final.vcf.gz"), path("${strain_id}.final.vcf.gz.tbi"), emit: final_vcf
    path "${strain_id}.variant_eval.txt", emit: stats_report 

    script:
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : "-Xmx4g"
    """
    # 1. Merge filtered SNPs and INDELs back together
    gatk --java-options "${avail_mem}" MergeVcfs \
        -I ${snp_vcf} \
        -I ${indel_vcf} \
        -O merged.vcf.gz
 
    # 2. Drop non-variants and failed variants
    gatk --java-options "${avail_mem}" SelectVariants \
        -R ${fasta} \
        -V merged.vcf.gz \
        --exclude-non-variants \
        --exclude-filtered \
        -O ${strain_id}.final.vcf.gz

    # 3. Generate Survivor QC report (Native GATK) for MultiQC
    gatk --java-options "${avail_mem}" VariantEval \
        -R ${fasta} \
        -eval ${strain_id}.final.vcf.gz \
        -O ${strain_id}.variant_eval.txt
    """
}