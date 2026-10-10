/*
 * Module 2: Read Filtration (fastp)
 * ----------------------------------
 * Two processes handle the PE and SE streams independently.
 * Each produces filtered FASTQs + JSON/HTML QC reports.
 *
 * Inputs come from Module 1 merge outputs:
 *   MERGE_PE  → (strain_id, PE_R1.fastq.gz, PE_R2.fastq.gz)
 *   MERGE_SE  → (strain_id, SE.fastq.gz)
 *
 * Features:
 *   - Dynamic parameter injection from nextflow.config
 *   - Post-filter gzip integrity check (gzip -t)
 *   - Non-empty output validation (filtered files must have reads)
 *   - Read-count summary logged to stdout for traceability
 *   - Adapter auto-detection enabled (fastp default)
 *   - Overrepresented sequence analysis enabled
 */

process FILTER_PE {
    tag "${strain_id}"
    label 'base'

    publishDir "${params.outdir}/Filtration/filtered", mode: 'symlink', pattern: '*.filt.*.fastq.gz'
    publishDir "${params.outdir}/Filtration/reports",  mode: 'symlink', pattern: '*.{html,json}'

    input:
    tuple val(strain_id), path(r1), path(r2)

    output:
    tuple val(strain_id), path("${strain_id}.filt.R1.fastq.gz"), path("${strain_id}.filt.R2.fastq.gz"), emit: filtered_pe
    tuple val(strain_id), path("${strain_id}.fastp.json"),                                              emit: json_report
    tuple val(strain_id), path("${strain_id}.fastp.html"),                                              emit: html_report

    script:
    def user_args = task.ext.args ?: ''
    
    // Dynamically build arguments from params block
    def opts = []
    if (params.adapter_fasta)         opts << "--adapter_fasta ${params.adapter_fasta}"
    if (params.trim_front1 > 0)       opts << "--trim_front1 ${params.trim_front1}"
    if (params.trim_tail1 > 0)        opts << "--trim_tail1 ${params.trim_tail1}"
    if (params.trim_front2 > 0)       opts << "--trim_front2 ${params.trim_front2}"
    if (params.trim_tail2 > 0)        opts << "--trim_tail2 ${params.trim_tail2}"
    if (params.qualified_quality > 0) opts << "-q ${params.qualified_quality}"
    if (params.cut_mean_quality > 0)  opts << "--cut_right --cut_mean_quality ${params.cut_mean_quality}"
    
    def extra_args = opts.join(' ')

    """
    # --- Run fastp (PE mode) ---
    fastp \\
        -i ${r1} \\
        -I ${r2} \\
        -o ${strain_id}.filt.R1.fastq.gz \\
        -O ${strain_id}.filt.R2.fastq.gz \\
        --dont_overwrite \\
        -G \\
        -l ${params.min_read_length} \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --overrepresentation_analysis \\
        ${extra_args} \\
        -h ${strain_id}.fastp.html \\
        -j ${strain_id}.fastp.json \\
        ${user_args}

    # --- Post-filter integrity checks ---
    gzip -t ${strain_id}.filt.R1.fastq.gz
    gzip -t ${strain_id}.filt.R2.fastq.gz

    # --- Verify outputs are non-empty (contain at least 1 read) ---
    R1_LINES=\$(zcat ${strain_id}.filt.R1.fastq.gz | head -4 | wc -l)
    R2_LINES=\$(zcat ${strain_id}.filt.R2.fastq.gz | head -4 | wc -l)

    if [ "\$R1_LINES" -lt 4 ]; then
        echo "ERROR: ${strain_id} FILTER_PE produced empty R1 after filtering" >&2
        exit 1
    fi
    if [ "\$R2_LINES" -lt 4 ]; then
        echo "ERROR: ${strain_id} FILTER_PE produced empty R2 after filtering" >&2
        exit 1
    fi

    # --- Log read counts for traceability ---
    R1_READS=\$(zcat ${strain_id}.filt.R1.fastq.gz | awk 'END{print NR/4}')
    R2_READS=\$(zcat ${strain_id}.filt.R2.fastq.gz | awk 'END{print NR/4}')
    echo "[FILTER_PE] ${strain_id}: R1=\${R1_READS} reads, R2=\${R2_READS} reads after filtering"
    """
}

process FILTER_SE {
    tag "${strain_id}"
    label 'base'

    publishDir "${params.outdir}/Filtration/filtered", mode: 'symlink', pattern: '*.filt.SE.fastq.gz'
    publishDir "${params.outdir}/Filtration/reports",  mode: 'symlink', pattern: '*.{html,json}'

    input:
    tuple val(strain_id), path(se)

    output:
    tuple val(strain_id), path("${strain_id}.filt.SE.fastq.gz"), emit: filtered_se
    tuple val(strain_id), path("${strain_id}.SE.fastp.json"),    emit: json_report
    tuple val(strain_id), path("${strain_id}.SE.fastp.html"),    emit: html_report

    script:
    def user_args = task.ext.args ?: ''
    
    // Dynamically build arguments (Exclude R2-specific flags for SE)
    def opts = []
    if (params.adapter_fasta)         opts << "--adapter_fasta ${params.adapter_fasta}"
    if (params.trim_front1 > 0)       opts << "--trim_front1 ${params.trim_front1}"
    if (params.trim_tail1 > 0)        opts << "--trim_tail1 ${params.trim_tail1}"
    if (params.qualified_quality > 0) opts << "-q ${params.qualified_quality}"
    if (params.cut_mean_quality > 0)  opts << "--cut_right --cut_mean_quality ${params.cut_mean_quality}"
    
    def extra_args = opts.join(' ')

    """
    # --- Run fastp (SE mode) ---
    fastp \\
        -i ${se} \\
        -o ${strain_id}.filt.SE.fastq.gz \\
        --dont_overwrite \\
        -G \\
        -l ${params.min_read_length} \\
        --thread ${task.cpus} \\
        --overrepresentation_analysis \\
        ${extra_args} \\
        -h ${strain_id}.SE.fastp.html \\
        -j ${strain_id}.SE.fastp.json \\
        ${user_args}

    # --- Post-filter integrity check ---
    gzip -t ${strain_id}.filt.SE.fastq.gz

    # --- Verify output is non-empty ---
    SE_LINES=\$(zcat ${strain_id}.filt.SE.fastq.gz | head -4 | wc -l)

    if [ "\$SE_LINES" -lt 4 ]; then
        echo "ERROR: ${strain_id} FILTER_SE produced empty output after filtering" >&2
        exit 1
    fi

    # --- Log read count for traceability ---
    SE_READS=\$(zcat ${strain_id}.filt.SE.fastq.gz | awk 'END{print NR/4}')
    echo "[FILTER_SE] ${strain_id}: \${SE_READS} reads after filtering"
    """
}

process WRITE_FILTER_SUMMARY {
    tag "filter_summary"
    label 'tiny'

    publishDir "${params.outdir}/Filtration", mode: 'symlink'

    input:
    val(pe_json_list)
    val(se_json_list)

    output:
    path "filtration_summary.tsv", emit: summary

    script:
    """
    #!/usr/bin/env python3
    import json, os, sys

    header = "strain_id\\tmode\\treads_before\\treads_after\\treads_passed_pct\\tq20_rate_before\\tq30_rate_before\\tadapter_trimmed_pct\\tduplicate_rate"
    rows = []

    pe_files = "${pe_json_list}".split()
    se_files = "${se_json_list}".split()

    for f in pe_files + se_files:
        f = f.strip().rstrip(',').strip('[').strip(']')
        if not f or not os.path.isfile(f):
            continue
        with open(f) as fh:
            data = json.load(fh)

        strain = os.path.basename(f).replace('.fastp.json','').replace('.SE.fastp.json','')
        mode = 'SE' if '.SE.fastp' in f else 'PE'

        summary = data.get('summary', {})
        before = summary.get('before_filtering', {})
        after  = summary.get('after_filtering', {})

        reads_before = before.get('total_reads', 0)
        reads_after  = after.get('total_reads', 0)
        pct_passed   = f"{reads_after/reads_before*100:.2f}" if reads_before > 0 else "0.00"
        q20_before   = f"{before.get('q20_rate', 0)*100:.2f}"
        q30_before   = f"{before.get('q30_rate', 0)*100:.2f}"

        adapter = data.get('adapter_cutting', {})
        adapter_pct  = f"{adapter.get('adapter_trimmed_reads', 0)/reads_before*100:.2f}" if reads_before > 0 else "0.00"

        dup_rate = f"{data.get('duplication', {}).get('rate', 0)*100:.2f}"

        rows.append(f"{strain}\\t{mode}\\t{reads_before}\\t{reads_after}\\t{pct_passed}\\t{q20_before}\\t{q30_before}\\t{adapter_pct}\\t{dup_rate}")

    with open('filtration_summary.tsv', 'w') as out:
        out.write(header + '\\n')
        for r in sorted(rows):
            out.write(r + '\\n')

    print(f"[WRITE_FILTER_SUMMARY] Wrote {len(rows)} strain entries")
    """
}