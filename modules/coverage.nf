nextflow.enable.dsl = 2

/*
 * Module 4: Sliding-Window Coverage & Interactive Ploidy Curation
 * -------------------------------------------------------------------
 * - Stage 1 (WRITE_COVERAGE_SUMMARY): Aggregates sliding-window depths into `lager_coverage_master.xlsx`
 *   featuring an Interactive_Plotter sheet (with live embedded subgenome bar charts,
 *   dashed chromosome dividers, G20 Strain_ID auto-lookup, Column E manual paste, and live F:K
 *   chromosome copy-number tables) plus a formula-driven Cohort_Chromosome_Summary sheet.
 * - Stage 2 (PLOT_COVERAGE): Reads the user-curated `lager_coverage_master.xlsx` to batch-render
 *   clean cohort karyograms (`karyograms/*.png` and `cohort_karyograms.pdf`).
 */

process COMPUTE_COVERAGE {
    tag "${sample_id}"
    label 'low'

    publishDir "${params.outdir}/04_coverage/sliding_windows", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome_file
    path windows_bed

    output:
    tuple val(sample_id), path("${sample_id}.slidingwindow.tab"), emit: sliding_window_cov

    script:
    """
    set -euo pipefail
    bedtools coverage -sorted -g ${genome_file} -a ${windows_bed} -b ${bam} -mean > ${sample_id}.slidingwindow.tab
    """
}

process WRITE_COVERAGE_SUMMARY {
    tag "cohort_coverage_master_excel"
    label 'low'

    stageInMode 'copy'

    publishDir "${params.outdir}/04_coverage/summary", mode: 'copy'

    input:
    path window_files
    path bin_dir
    val species_map_str
    val ploidy_mult
    val auto_adjust

    output:
    path "lager_coverage_master.xlsx",  emit: excel
    path "default_ploidy_table.tsv",    emit: ploidy_tsv

    script:
    def auto_flag = auto_adjust ? "--auto-adjust-ploidy" : ""
    """
    set -euo pipefail

    python3 ${bin_dir}/summarize_coverage.py \\
        --inputs ${window_files} \\
        --species-map "${species_map_str}" \\
        --ploidy-multiplier ${ploidy_mult} \\
        ${auto_flag} \\
        --out-excel lager_coverage_master.xlsx
    """
}

process PLOT_COVERAGE {
    tag "cohort_karyograms"
    label 'phylogeny'

    stageInMode 'copy'

    publishDir "${params.outdir}/04_coverage/karyograms", mode: 'copy'

    input:
    path window_files
    path curated_source
    val source_type
    path bin_dir
    val species_map_str
    val ploidy_mult

    output:
    path "karyograms/*.png",            emit: pngs
    path "cohort_karyograms.pdf",       emit: pdf
    path "applied_ploidy_settings.tsv", emit: applied_settings

    script:
    """
    set -euo pipefail

    if [ "${source_type}" == "xlsx" ]; then
        python3 ${bin_dir}/summarize_coverage.py \\
            --extract-curated-excel "${curated_source}" \\
            --out-ploidy-tsv applied_ploidy_settings.tsv
    else
        cp "${curated_source}" applied_ploidy_settings.tsv
    fi

    Rscript ${bin_dir}/plot_coverage.R \\
        --species_map "${species_map_str}" \\
        --ploidy_mult ${ploidy_mult} \\
        --ploidy_tsv applied_ploidy_settings.tsv \\
        ${window_files}
    """
}