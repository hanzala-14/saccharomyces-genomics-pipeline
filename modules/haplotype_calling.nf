/*
 * Module 5: GATK HaplotypeCaller — Per-chromosome GVCF calling
 * -------------------------------------------------------------
 * Scatters variant calling across individual chromosomes.
 * Gathers them per-strain to generate a manifest for Module 6.
 */

process PREPARE_GATK_REF {
    tag { fasta.name }
    label 'medium' // Uses 'medium' because it contains Picard & Samtools

    input:
    path fasta

    output:
    tuple path(fasta), path("*.fai"), path("*.dict"), emit: gatk_ref

    script:
    def dict_name = fasta.baseName + ".dict"
    """
    set -euo pipefail
    
    # 1. Create .fai index using Samtools
    samtools faidx "${fasta}"
    
    # 2. Create .dict using GATK
    gatk CreateSequenceDictionary R="${fasta}" O="${dict_name}"
    """
}

process HAPLOTYPE_CALLER {
    tag { "${strain_id}:${chrom}" }
    label 'calling'
    
    publishDir path: { "${params.outdir}/Haplotype_Calling/${strain_id}" }, mode: 'symlink'

    errorStrategy { task.exitStatus in [137, 143, 247] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    tuple val(strain_id), path(bam), path(bai), val(chrom)
    tuple path(fasta), path(fai), path(dict) 
    // ❌ REMOVED: val ploidy_map channel input

    output:
    tuple val(strain_id), val(chrom), path("${strain_id}.${chrom}.g.vcf.gz"), path("${strain_id}.${chrom}.g.vcf.gz.tbi"), emit: gvcf

    script:
    def avail_mem = (task.memory.toGiga() * 0.8).intValue()
    
    // Dynamic Ploidy: Extract the suffix after the underscore and look it up in the map
    def parts = chrom.tokenize('_')
    def suffix = parts.size() > 1 ? parts[1] : ""
    
    // ✅ FIX: Read directly from the global params map!
    def ploidy = params.ploidy_map.containsKey(suffix) ? params.ploidy_map[suffix] : params.ploidy

    """
    set -euo pipefail

    echo "[HAPLOTYPE_CALLER] Running ${strain_id} on chromosome ${chrom} (Ploidy: ${ploidy})..."

    gatk --java-options "-Xmx${avail_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads=2" \\
        HaplotypeCaller \\
        -R "${fasta}" \\
        -I "${bam}" \\
        -O "${strain_id}.${chrom}.g.vcf.gz" \\
        -G StandardAnnotation \\
        -G StandardHCAnnotation \\
        -L "${chrom}" \\
        -ploidy ${ploidy} \\
        -ERC GVCF

    if [ ! -s "${strain_id}.${chrom}.g.vcf.gz" ]; then
        echo "ERROR: GVCF not produced for ${strain_id} ${chrom}" >&2
        exit 1
    fi
    """
}

process GATHER_STRAIN_GVCFS {
    tag { "${strain_id}" }
    label 'tiny'
    
    // Wrap the publishDir path in a closure by adding "path: { ... }"
    publishDir path: { "${params.outdir}/Haplotype_Calling/${strain_id}" }, mode: 'symlink'

    input:
    tuple val(strain_id), val(chroms), path(gvcfs), path(tbis)

    output:
    tuple val(strain_id), path(gvcfs), path(tbis), emit: strain_gvcfs
    tuple val(strain_id), path("${strain_id}.gvcf_list.txt"), emit: gvcf_list

    script:
    def gvcf_files = gvcfs.collect { f -> f.name }.sort().join('\n')
    """
    cat > "${strain_id}.gvcf_list.txt" << 'EOF'
${gvcf_files}
EOF
    echo "[GATHER] Strain ${strain_id}: ${gvcfs.size()} GVCFs collected and manifest generated."
    """
}