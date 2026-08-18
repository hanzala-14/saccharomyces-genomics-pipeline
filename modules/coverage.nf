/*
 * Module 4: Genome Coverage Profiling & Sliding Window Analysis
 * -------------------------------------------------------------
 * Computes per-base depth from deduplicated BAMs, converts it to bedGraph,
 * sorts coordinates, and calculates sliding window median coverage.
 */
process COMPUTE_COVERAGE {
    tag "${strain_id}"
    label 'low'
    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'

    publishDir "${params.outdir}/Coverage", mode: 'copy', pattern: '*.slidingwindow.tab'

    input:
    tuple val(strain_id), path(bam), path(bai)
    path(genome_file)
    path(sliding_windows)

    output:
    tuple val(strain_id), path("${strain_id}.slidingwindow.tab"), emit: sliding_window_cov

    script:
    """
    set -euo pipefail

    echo "[COMPUTE_COVERAGE] Starting coverage profiling for ${strain_id}..."

    # 1. Per-base depth calculation
    bedtools genomecov -ibam ${bam} -g ${genome_file} -d > ${strain_id}.temporary.tab

    # 2. Convert to per-base bedGraph format
    awk '{print \$1"\\t"\$2-1"\\t"\$2"\\t"\$3}' ${strain_id}.temporary.tab > ${strain_id}.per_base.bedGraph

    # 3. Sort bedGraph according to reference genome layout
    bedtools sort -g ${genome_file} -i ${strain_id}.per_base.bedGraph > ${strain_id}.sort_perbase.bedGraph

    # 4. Compute sliding window median coverage
    bedtools map -a ${sliding_windows} \
                 -b ${strain_id}.sort_perbase.bedGraph \
                 -c 4 -o median \
                 -g ${genome_file} \
                 > ${strain_id}.slidingwindow.tab

    # 5. Clean up temporary files
    rm -f ${strain_id}.temporary.tab ${strain_id}.per_base.bedGraph ${strain_id}.sort_perbase.bedGraph

    # 6. Validate output
    if [ ! -s "${strain_id}.slidingwindow.tab" ]; then
        echo "ERROR: ${strain_id} coverage profiling produced an empty output file." >&2
        exit 1
    fi

    echo "[COMPUTE_COVERAGE] Completed sliding window coverage for ${strain_id}."
    """
}