/*
 * Module 1: Data Acquisition, QC Classification & Strain Aggregation
 * -------------------------------------------------------------------
 *  FETCH_SRA with:
 *   - Primary Strategy: Multi-threaded ENA FTP via aria2c to bypass throttling
 *   - Secondary Fallback: NCBI SRA prefetch → fasterq-dump pipeline if ENA fails
 *   - Bandwidth Saturation: Maxes out bandwidth via maxForks parallelism
 *   - Resource Efficiency: Downloads pre-compressed .fastq.gz directly to save local CPU/disk
 *   - Hard Network Rules: Timeouts, connection drops, and exponential retry back-offs
 *   - Strict Integrity: Automated `gzip -t` checks to guarantee non-corrupt reads
 *   - Non-Interactive: Pre-configured SRA Toolkit settings suppressing interactive prompts
 */
process FETCH_SRA {
    tag { "${strain_id} ${accession}" }
    label 'base'

    // Error logic is mapped to nextflow.config parameters
    errorStrategy { task.exitStatus in [1, 137, 143, 255] ? 'retry' : 'finish' }
    maxRetries params.max_Retries

    input:
    tuple val(strain_id), val(accession)

    output:
    tuple val(strain_id), val(accession), path("${accession}_1.fastq.gz"), path("${accession}_2.fastq.gz")

    script:
    if (!(accession ==~ /^[A-Z]{3}\d{6,8}$/)) {
        error "[FETCH_SRA] strain=${strain_id} invalid accession format: ${accession}"
    }

    def len = accession.size()
    def pfx = accession[0..5]
    def subdir = (len == 9) ? null : (len == 10) ? "00${accession[-1]}" : (len == 11) ? "0${accession[-2..-1]}" : (len == 12) ? accession[-3..-1] : null
    
    def base = (len == 9)
        ? "https://ftp.sra.ebi.ac.uk/vol1/fastq/${pfx}/${accession}"
        : "https://ftp.sra.ebi.ac.uk/vol1/fastq/${pfx}/${subdir}/${accession}"

    def retrySleep = 10 * (task.attempt ?: 1)

    // Nextflow Parameters //
    def ncbi_dir         = params.ncbi_dir ?: '~/.ncbi'
    def aria_min_split   = params.aria2c_min_split_size ?: '1M'
    def aria_conn_timeout= params.aria2c_connect_timeout ?: 30
    def aria_timeout     = params.aria2c_timeout ?: 60
    def aria_max_tries   = params.aria2c_max_tries ?: 3
    def aria_retry_wait  = params.aria2c_retry_wait ?: 10
    
    def prefetch_t_out   = params.prefetch_timeout ?: 3600
    def prefetch_size    = params.prefetch_max_size ?: '50G'
    def fasterq_t_out    = params.fasterq_timeout ?: 1800

    """
    set -euo pipefail

    echo "=== FETCH_SRA: strain=${strain_id} accession=${accession} ===" >&2

    # ─── SRA Toolkit Config ────────────────────────────────────────────────
    mkdir -p "${ncbi_dir}"
    cat > "${ncbi_dir}/user-settings.mkfg" <<'MKFG'
/LIBS/IMAGE_GUID = "auto-nextflow-pipeline"
/libs/cloud/report_instance_identity = "false"
/libs/cloud/accept_aws_charges = "false"
/repository/user/main/public/root = "."
MKFG
    export VDB_CONFIG="${ncbi_dir}"
    export NCBI_SETTINGS="${ncbi_dir}/user-settings.mkfg"

    if [ "${task.attempt}" -gt 1 ]; then
        echo "Retry attempt ${task.attempt} — sleeping ${retrySleep}s..." >&2
        sleep ${retrySleep}
    fi

    ena_ok=0

    # ─── Strategy 1: ENA FTP via Aria2c ───────────────────────────────────
    echo "[Strategy 1] Trying ENA via Aria2c..." >&2

    ARIA_OPTS="-x ${params.aria2c_connections} -s ${params.aria2c_connections} -c \\
        --max-connection-per-server=${params.aria2c_connections} \\
        --min-split-size=${aria_min_split} \\
        --connect-timeout=${aria_conn_timeout} \\
        --timeout=${aria_timeout} \\
        --max-tries=${aria_max_tries} \\
        --retry-wait=${aria_retry_wait} \\
        --console-log-level=notice \\
        --summary-interval=10"

    if aria2c \$ARIA_OPTS -o "${accession}_1.fastq.gz" "${base}/${accession}_1.fastq.gz"; then
        if [ -s "${accession}_1.fastq.gz" ] && gzip -t "${accession}_1.fastq.gz" 2>/dev/null; then
            ena_ok=1
            echo "[Strategy 1] R1 downloaded from ENA." >&2

            if aria2c \$ARIA_OPTS -o "${accession}_2.fastq.gz" "${base}/${accession}_2.fastq.gz"; then
                if ! gzip -t "${accession}_2.fastq.gz" 2>/dev/null; then
                    rm -f "${accession}_2.fastq.gz"
                fi
            fi
        else
            rm -f "${accession}_1.fastq.gz" "${accession}_2.fastq.gz"
        fi
    fi

    # ─── Strategy 2: prefetch fallback (If ENA is missing) ────────────────
    if [ "\$ena_ok" -eq 0 ]; then
        echo "[Strategy 2] ENA failed. Using prefetch + fasterq-dump..." >&2
        
        timeout ${prefetch_t_out} prefetch --max-size ${prefetch_size} "${accession}" || exit 1
        timeout ${fasterq_t_out} fasterq-dump --split-files --threads ${task.cpus} --temp . "${accession}" || exit 1
        
        rm -rf "${accession}/" || true

        if command -v pigz &>/dev/null; then ZIPCMD="pigz -p ${task.cpus}"; else ZIPCMD="gzip"; fi
        [ -f "${accession}_1.fastq" ] && \$ZIPCMD "${accession}_1.fastq"
        [ -f "${accession}.fastq" ]   && mv "${accession}.fastq" "${accession}_1.fastq" && \$ZIPCMD "${accession}_1.fastq"
        [ -f "${accession}_2.fastq" ] && \$ZIPCMD "${accession}_2.fastq"
    fi

    # ─── Validation ───────────────────────────────────────────────────────
    if [ ! -s "${accession}_1.fastq.gz" ] || ! gzip -t "${accession}_1.fastq.gz" 2>/dev/null; then
        echo "ERROR: Download failed or corrupt." >&2
        exit 1
    fi
    [ ! -f "${accession}_2.fastq.gz" ] && touch "${accession}_2.fastq.gz"
    echo "=== FETCH_SRA COMPLETE ===" >&2
    """
}

process INGEST_LOCAL {
    tag { "${strain_id}" }
    label 'tiny'

    input:
    tuple val(strain_id), path(r1_file), path(r2_file)

    output:
    tuple val(strain_id), val(run_id), path("${run_id}_R1.fastq.gz"), path("${run_id}_R2.fastq.gz")

    script:
    run_id = "LOCAL_${strain_id}_${r1_file.baseName}".replaceAll(/[^A-Za-z0-9_.-]/, "_")

    """
    set -euo pipefail

    # Copy R1 (already validated by file(checkIfExists:true) in main.nf)
    cp "${r1_file}" "${run_id}_R1.fastq.gz"

    # Copy R2 if it's a real file (not the EMPTY placeholder)
    if [ "${r2_file.name}" != "EMPTY" ] && [ -s "${r2_file}" ]; then
        cp "${r2_file}" "${run_id}_R2.fastq.gz"
    else
        : > "${run_id}_R2.fastq.gz"
    fi
    """
}

process CLASSIFY_RUN {
    tag { "${strain_id} ${run_id}" }
    label 'tiny'

    input:
    tuple val(strain_id), val(run_id), path(r1), path(r2)

    output:
    tuple val(strain_id), val(run_id), stdout, path(r1), path(r2), emit: qc_rows

    script:
    def minR2 = params.min_r2_len as int

    """
    set -euo pipefail

    status="DROP"
    reason="R1_INVALID"
    r1_len=0
    r2_len=0

    r1_ok=0
    r2_ok=0

    if [ -s "${r1}" ] && gzip -t "${r1}" 2>/dev/null; then r1_ok=1; fi
    if [ -s "${r2}" ] && gzip -t "${r2}" 2>/dev/null; then r2_ok=1; fi

    if [ "\$r1_ok" -eq 1 ]; then
      r1_len=\$(python3 - "${r1}" <<'PY'
import gzip, sys
path = sys.argv[1]
count = 0
total = 0
with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as fh:
    for i, line in enumerate(fh, 1):
        if i % 4 == 2:
            total += len(line.rstrip('\\n'))
            count += 1
            if count == 2000:
                break
print(int(total / count) if count else 0)
PY
)
    fi
    if [ "\$r2_ok" -eq 1 ]; then
      r2_len=\$(python3 - "${r2}" <<'PY'
import gzip, sys
path = sys.argv[1]
count = 0
total = 0
with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as fh:
    for i, line in enumerate(fh, 1):
        if i % 4 == 2:
            total += len(line.rstrip('\\n'))
            count += 1
            if count == 2000:
                break
print(int(total / count) if count else 0)
PY
)
    fi

    is_local=0
    echo "${run_id}" | grep -q '^LOCAL_' && is_local=1 || true

    if [ "\$r1_ok" -eq 1 ] && [ "\$r2_ok" -eq 1 ]; then
      md1=\$(md5sum "${r1}" | awk '{print \$1}')
      md2=\$(md5sum "${r2}" | awk '{print \$1}')
      if [ "\$md1" = "\$md2" ]; then
        status="SE_FALLBACK"; reason="R1_R2_IDENTICAL"
      elif [ "\$r2_len" -lt ${minR2} ]; then
        status="SE_FALLBACK"; reason="R2_TOO_SHORT"
      else
        status="PE_PASS"; reason="OK"
      fi
    elif [ "\$r1_ok" -eq 1 ] && [ "\$r2_ok" -eq 0 ]; then
      if [ "\$is_local" -eq 1 ]; then
        status="SE_INPUT"; reason="LOCAL_R1_ONLY_OR_BAD_R2"
      else
        status="SE_FALLBACK"; reason="R2_INVALID_OR_MISSING"
      fi
    else
      status="DROP"; reason="R1_INVALID"
    fi

    printf '%s\t%s\t%s\t%s\n' "\$status" "\$reason" "\$r1_len" "\$r2_len"
    """
}

process MERGE_PE {
    tag { "${strain_id}" }
    label 'base'
    publishDir path: { "${params.outdir}/strains/${strain_id}" }, mode: 'copy'

    input:
    tuple val(strain_id), path(r1_files), path(r2_files)

    output:
    tuple val(strain_id), path("${strain_id}_PE_R1.fastq.gz"), path("${strain_id}_PE_R2.fastq.gz"), emit: merged_pe

    script:
    """
    set -euo pipefail
    cat ${r1_files.join(' ')} > "${strain_id}_PE_R1.fastq.gz"
    cat ${r2_files.join(' ')} > "${strain_id}_PE_R2.fastq.gz"
    gzip -t "${strain_id}_PE_R1.fastq.gz"
    gzip -t "${strain_id}_PE_R2.fastq.gz"
    """
}

process MERGE_SE {
    tag { "${strain_id}" }
    label 'base'
    publishDir path: { "${params.outdir}/strains/${strain_id}" }, mode: 'copy'

    input:
    tuple val(strain_id), path(r1_files)

    output:
    tuple val(strain_id), path("${strain_id}_SE.fastq.gz"), emit: merged_se

    script:
    """
    set -euo pipefail
    cat ${r1_files.join(' ')} > "${strain_id}_SE.fastq.gz"
    gzip -t "${strain_id}_SE.fastq.gz"
    """
}

process CREATE_STRAIN_LIST {
    label 'tiny'
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    val strain_ids

    output:
    path "strains.txt"

    script:
    def list_str = strain_ids.collect { sid -> sid.toString().trim() }.findAll { sid -> sid }.unique().sort().join('\n')
    """
    cat > strains.txt << 'EOF'
${list_str}
EOF
    """
}

process WRITE_RUN_QC {
    label 'tiny'
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    val rows

    output:
    path "run_qc.tsv"

    script:
    def row_list = (rows instanceof java.util.Collection && !rows.isEmpty() && rows[0] instanceof java.util.Collection) ? rows : [rows]
    def body = row_list.collect { r ->
        def values = (r instanceof java.util.Collection) ? r : [r]
        [values[0], values[1], values[2], values[3], values[6], values[7]].join('\t')
    }.join('\n')

    """
    cat > run_qc.tsv << 'EOF'
strain_id\trun_id\tstatus\treason\tr1_mean_len\tr2_mean_len
${body}
EOF
    """
}

process WRITE_STRAIN_QC {
    label 'tiny'
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    val rows

    output:
    path "strain_qc.tsv"

    script:
    def row_list = (rows instanceof java.util.Collection && !rows.isEmpty() && rows[0] instanceof java.util.Collection) ? rows : [rows]
    def body = row_list.collect { r ->
        def values = (r instanceof java.util.Collection) ? r : [r]
        [values[0], values[1], values[2], values[3], values[4]].join('\t')
    }.join('\n')
    """
    cat > strain_qc.tsv << 'EOF'
strain_id\tpe_runs\tse_runs\tdropped_runs\tfinal_mode
${body}
EOF
    """
}