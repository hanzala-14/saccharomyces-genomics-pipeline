/*
 * Module 1: Data Acquisition, QC Classification & Strain Aggregation
 * -------------------------------------------------------------------
 *
 * FETCH_SRA:
 *   1. Query ENA for FASTQ metadata.
 *   2. If ENA metadata is available:
 *        -> download FASTQ directly with aria2c.
 *   3. If ENA metadata is unavailable:
 *        -> construct standard ENA FASTQ path.
 *        -> try aria2c directly.
 *   4. If ENA retrieval fails:
 *        -> NCBI prefetch
 *        -> vdb-validate
 *        -> fasterq-dump (default split-3)
 *   5. Validate PE R1/R2 read counts and read IDs.
 *
 * Notes:
 *   - ENA is always preferred when archive-generated FASTQ is available.
 *   - NCBI is a true fallback.
 *   - gzip -t checks compression integrity.
 *   - PE read-ID validation checks actual pairing integrity.
 */

process FETCH_SRA {

    tag { "${strain_id} ${accession}" }

    label 'base'

    errorStrategy { task.exitStatus in [1, 137, 143, 255] ? 'retry' : 'finish' }

    maxRetries params.max_Retries

    input:
    tuple val(strain_id), val(accession)

    output:
    tuple val(strain_id),
          val(accession),
          path("${accession}_1.fastq.gz"),
          path("${accession}_2.fastq.gz")

    script:

    if (!(accession ==~ /^[A-Z]{3}\d{6,8}$/)) {
        error "[FETCH_SRA] strain=${strain_id} invalid accession format: ${accession}"
    }

    def retrySleep        = params.fetch_retry_sleep * (task.attempt ?: 1)

    def ncbi_dir          = params.ncbi_dir

    def aria_min_split    = params.aria2c_min_split_size
    def aria_conn_timeout = params.aria2c_connect_timeout
    def aria_timeout      = params.aria2c_timeout
    def aria_max_tries    = params.aria2c_max_tries
    def aria_retry_wait   = params.aria2c_retry_wait
    def aria_summary      = params.aria2c_summary_interval

    def prefetch_t_out    = params.prefetch_timeout
    def prefetch_size     = params.prefetch_max_size
    def fasterq_t_out     = params.fasterq_timeout

    def api_max_time      = params.ena_api_max_time
    def api_retries       = params.ena_api_retries
    def api_retry_wait    = params.ena_api_retry_wait

    /*
     * ENA archive-generated FASTQ structure:
     *
     * SRR8799653
     *   prefix = SRR879
     *   subdir = 003   (00 + last digit)
     */
    def ena_prefix = accession.take(6)
    def ena_suffix = "00${accession[-1]}"

    """
    set -euo pipefail

    echo "=== FETCH_SRA: strain=${strain_id} accession=${accession} attempt=${task.attempt} ===" >&2

    # =========================================================================
    # SRA TOOLKIT CONFIGURATION
    # =========================================================================

    mkdir -p "${ncbi_dir}"

    cat > "${ncbi_dir}/user-settings.mkfg" <<'MKFG'
/LIBS/IMAGE_GUID = "auto-nextflow-pipeline"
/libs/cloud/report_instance_identity = "false"
/libs/cloud/accept_aws_charges = "false"
/repository/user/main/public/root = "."
MKFG

    export VDB_CONFIG="${ncbi_dir}"
    export NCBI_SETTINGS="${ncbi_dir}/user-settings.mkfg"

    # =========================================================================
    # RETRY BACK-OFF
    # =========================================================================

    if [ "${task.attempt}" -gt 1 ]; then
        echo "Retry attempt ${task.attempt} — sleeping ${retrySleep}s..." >&2
        sleep "${retrySleep}"
    fi

    # =========================================================================
    # FUNCTION: VALIDATE R1/R2 PAIRING
    # =========================================================================

    validate_pe_pair() {

        local r1="\$1"
        local r2="\$2"
        local accession_name="\$3"

        local r1_reads
        local r2_reads

        echo "[VALIDATE_PE] \${accession_name}: checking read counts..." >&2

        r1_reads=\$(zcat "\${r1}" | awk 'END {print NR/4}')
        r2_reads=\$(zcat "\${r2}" | awk 'END {print NR/4}')

        echo "[VALIDATE_PE] \${accession_name}: R1=\${r1_reads}, R2=\${r2_reads}" >&2

        if [ "\${r1_reads}" -ne "\${r2_reads}" ]; then
            echo "ERROR: \${accession_name}: R1/R2 read counts differ." >&2
            echo "       R1=\${r1_reads}, R2=\${r2_reads}" >&2
            return 1
        fi

        echo "[VALIDATE_PE] \${accession_name}: checking read-name synchronization..." >&2

        if ! paste \
            <(
                zcat "\${r1}" |
                awk 'NR % 4 == 1 {
                    name=\$1
                    sub(/^@/, "", name)
                    if (length(name) >= 2) {
                        suffix=substr(name, length(name)-1, 2)
                        if (suffix == "/1" || suffix == "/2") {
                            name=substr(name, 1, length(name)-2)
                        }
                    }
                    print name
                }'
            ) \
            <(
                zcat "\${r2}" |
                awk 'NR % 4 == 1 {
                    name=\$1
                    sub(/^@/, "", name)
                    if (length(name) >= 2) {
                        suffix=substr(name, length(name)-1, 2)
                        if (suffix == "/1" || suffix == "/2") {
                            name=substr(name, 1, length(name)-2)
                        }
                    }
                    print name
                }'
            ) |
            awk '\$1 != \$2 {
                print "ERROR: R1/R2 mismatch at read", NR > "/dev/stderr"
                print "       R1:", \$1 > "/dev/stderr"
                print "       R2:", \$2 > "/dev/stderr"
                exit 1
            }'
        then
            return 1
        fi

        echo "[VALIDATE_PE] \${accession_name}: R1/R2 pairing validated." >&2

        return 0
    }

    # =========================================================================
    # STEP 0: QUERY ENA METADATA
    # =========================================================================

    layout="UNKNOWN"

    url_r1=""
    url_r2=""
    url_se=""

    API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"

    fastq_field=""

    if command -v curl &>/dev/null; then

        fastq_field=\$(
            curl -fsSL \
                --retry ${api_retries} \
                --retry-delay ${api_retry_wait} \
                --max-time ${api_max_time} \
                "\${API_URL}" 2>/dev/null |
            awk 'NR == 2 {print \$2}' ||
            true
        )

    fi

    if [ -n "\${fastq_field}" ]; then

        IFS=';' read -ra urls <<< "\${fastq_field}"

        for u in "\${urls[@]}"; do

            # Normalize URL scheme.
            u="\${u#ftp://}"
            u="\${u#http://}"
            u="\${u#https://}"

            case "\${u}" in

                *_1.fastq.gz)
                    url_r1="http://\${u}"
                    ;;

                *_2.fastq.gz)
                    url_r2="http://\${u}"
                    ;;

                *.fastq.gz)
                    url_se="http://\${u}"
                    ;;

            esac

        done

        if [ -n "\${url_r1}" ] && [ -n "\${url_r2}" ]; then
            layout="PE"

        elif [ -n "\${url_r1}" ]; then
            layout="SE"
            url_se="\${url_r1}"

        elif [ -n "\${url_se}" ]; then
            layout="SE"
        fi

        echo "[ENA API] ${accession} reported layout=\${layout}" >&2

    else

        echo "[ENA API] No FASTQ metadata returned for ${accession}." >&2
        echo "[ENA API] Will try direct ENA FASTQ retrieval with aria2c." >&2

    fi

    ena_ok=0

    # =========================================================================
    # ARIA2C OPTIONS
    # =========================================================================

    if command -v aria2c &>/dev/null; then

        ARIA_OPTS="-x ${params.aria2c_connections} \
            -s ${params.aria2c_connections} \
            -c \
            --max-connection-per-server=${params.aria2c_connections} \
            --min-split-size=${aria_min_split} \
            --connect-timeout=${aria_conn_timeout} \
            --timeout=${aria_timeout} \
            --max-tries=${aria_max_tries} \
            --retry-wait=${aria_retry_wait} \
            --console-log-level=notice \
            --summary-interval=${aria_summary}"

    fi

    # =========================================================================
    # STRATEGY 1: ENA API URL + aria2c
    # =========================================================================

    if command -v aria2c &>/dev/null && [ "\${layout}" != "UNKNOWN" ]; then

        echo "[Strategy 1] Trying ENA via aria2c (layout=\${layout})..." >&2

        if [ "\${layout}" = "PE" ]; then

            rm -f "${accession}_1.fastq.gz"
            rm -f "${accession}_2.fastq.gz"

            if aria2c \${ARIA_OPTS} \
                -o "${accession}_1.fastq.gz" \
                "\${url_r1}" \
                && [ -s "${accession}_1.fastq.gz" ] \
                && gzip -t "${accession}_1.fastq.gz" 2>/dev/null \
                && aria2c \${ARIA_OPTS} \
                -o "${accession}_2.fastq.gz" \
                "\${url_r2}" \
                && [ -s "${accession}_2.fastq.gz" ] \
                && gzip -t "${accession}_2.fastq.gz" 2>/dev/null
            then

                if validate_pe_pair \
                    "${accession}_1.fastq.gz" \
                    "${accession}_2.fastq.gz" \
                    "${accession}"
                then

                    ena_ok=1
                    echo "[Strategy 1] ENA PE download + pairing validation succeeded." >&2

                else

                    echo "[Strategy 1] ENA PE pairing validation FAILED." >&2

                    rm -f "${accession}_1.fastq.gz"
                    rm -f "${accession}_2.fastq.gz"

                fi

            else

                echo "[Strategy 1] ENA PE download failed/incomplete/corrupt." >&2

                rm -f "${accession}_1.fastq.gz"
                rm -f "${accession}_2.fastq.gz"

            fi

        else

            rm -f "${accession}_1.fastq.gz"

            if aria2c \${ARIA_OPTS} \
                -o "${accession}_1.fastq.gz" \
                "\${url_se}" \
                && [ -s "${accession}_1.fastq.gz" ] \
                && gzip -t "${accession}_1.fastq.gz" 2>/dev/null
            then

                ena_ok=1
                echo "[Strategy 1] ENA SE download validated successfully." >&2

            else

                echo "[Strategy 1] ENA SE download failed/incomplete/corrupt." >&2

                rm -f "${accession}_1.fastq.gz"

            fi

        fi

    else

        echo "[Strategy 1] API-based ENA aria2c retrieval unavailable." >&2

    fi

    # =========================================================================
    # STRATEGY 1B: DIRECT ENA FASTQ PATH + aria2c
    #
    # Used when ENA metadata is unavailable.
    # =========================================================================

    if [ "\${ena_ok}" -eq 0 ] \
        && [ "\${layout}" = "UNKNOWN" ] \
        && command -v aria2c &>/dev/null
    then

        echo "[Strategy 1B] Trying direct ENA FASTQ path with aria2c..." >&2

        direct_base="http://ftp.sra.ebi.ac.uk/vol1/fastq/${ena_prefix}/${ena_suffix}/${accession}"

        rm -f "${accession}_1.fastq.gz"
        rm -f "${accession}_2.fastq.gz"

        # ---------------------------------------------------------------------
        # Try PE
        # ---------------------------------------------------------------------

        if aria2c \${ARIA_OPTS} \
            -o "${accession}_1.fastq.gz" \
            "\${direct_base}/${accession}_1.fastq.gz" \
            && [ -s "${accession}_1.fastq.gz" ] \
            && gzip -t "${accession}_1.fastq.gz" 2>/dev/null
        then

            echo "[Strategy 1B] Direct ENA R1 found." >&2

            if aria2c \${ARIA_OPTS} \
                -o "${accession}_2.fastq.gz" \
                "\${direct_base}/${accession}_2.fastq.gz" \
                && [ -s "${accession}_2.fastq.gz" ] \
                && gzip -t "${accession}_2.fastq.gz" 2>/dev/null
            then

                echo "[Strategy 1B] Direct ENA R2 found." >&2

                if validate_pe_pair \
                    "${accession}_1.fastq.gz" \
                    "${accession}_2.fastq.gz" \
                    "${accession}"
                then

                    layout="PE"
                    ena_ok=1

                    echo "[Strategy 1B] Direct ENA PE download + pairing validation succeeded." >&2

                else

                    echo "[Strategy 1B] Direct ENA PE pairing validation FAILED." >&2

                    rm -f "${accession}_1.fastq.gz"
                    rm -f "${accession}_2.fastq.gz"

                fi

            else

                echo "[Strategy 1B] Direct ENA R2 unavailable; trying SE." >&2

                rm -f "${accession}_1.fastq.gz"
                rm -f "${accession}_2.fastq.gz"

                if aria2c \${ARIA_OPTS} \
                    -o "${accession}_1.fastq.gz" \
                    "\${direct_base}/${accession}.fastq.gz" \
                    && [ -s "${accession}_1.fastq.gz" ] \
                    && gzip -t "${accession}_1.fastq.gz" 2>/dev/null
                then

                    layout="SE"
                    ena_ok=1

                    echo "[Strategy 1B] Direct ENA SE download succeeded." >&2

                else

                    rm -f "${accession}_1.fastq.gz"

                    echo "[Strategy 1B] Direct ENA FASTQ retrieval failed." >&2

                fi

            fi

        else

            rm -f "${accession}_1.fastq.gz"
            rm -f "${accession}_2.fastq.gz"

            echo "[Strategy 1B] Direct ENA R1 unavailable." >&2

        fi

    fi

    # =========================================================================
    # STRATEGY 2: NCBI FALLBACK
    # =========================================================================

    if [ "\${ena_ok}" -eq 0 ]; then

        echo "[Strategy 2] Using NCBI prefetch + vdb-validate + fasterq-dump..." >&2

        rm -f "${accession}_1.fastq.gz"
        rm -f "${accession}_2.fastq.gz"

        if command -v timeout &>/dev/null; then

            timeout ${prefetch_t_out} \
                prefetch \
                --max-size ${prefetch_size} \
                "${accession}" \
                || exit 1

            echo "[NCBI] Validating SRA object..." >&2

            timeout ${prefetch_t_out} \
                vdb-validate \
                "${accession}" \
                || exit 1

            echo "[NCBI] Extracting FASTQ with fasterq-dump..." >&2

            timeout ${fasterq_t_out} \
                fasterq-dump \
                --threads ${task.cpus} \
                --temp . \
                "${accession}" \
                || exit 1

        else

            prefetch \
                --max-size ${prefetch_size} \
                "${accession}" \
                || exit 1

            echo "[NCBI] Validating SRA object..." >&2

            vdb-validate \
                "${accession}" \
                || exit 1

            echo "[NCBI] Extracting FASTQ with fasterq-dump..." >&2

            fasterq-dump \
                --threads ${task.cpus} \
                --temp . \
                "${accession}" \
                || exit 1

        fi

        # ---------------------------------------------------------------------
        # Compress extracted FASTQ files
        # ---------------------------------------------------------------------

        if command -v pigz &>/dev/null; then
            ZIPCMD="pigz -p ${task.cpus}"
        else
            ZIPCMD="gzip"
        fi

        if [ -f "${accession}_1.fastq" ]; then
            \${ZIPCMD} "${accession}_1.fastq"
        fi

        if [ -f "${accession}_2.fastq" ]; then
            \${ZIPCMD} "${accession}_2.fastq"
        fi

        # ---------------------------------------------------------------------
        # fasterq-dump default split-3 singleton output:
        #
        # PE  -> _1.fastq + _2.fastq + optional .fastq singleton file
        # SE  -> .fastq only
        # ---------------------------------------------------------------------

        if [ -f "${accession}.fastq" ]; then

            if [ ! -f "${accession}_1.fastq.gz" ] \
                && [ ! -f "${accession}_2.fastq.gz" ]
            then

                # Genuine SE output
                mv "${accession}.fastq" "${accession}_1.fastq"
                \${ZIPCMD} "${accession}_1.fastq"

            else

                # PE with singleton reads.
                # Singleton reads are excluded from the PE stream.
                \${ZIPCMD} "${accession}.fastq"
                rm -f "${accession}.fastq.gz"

            fi

        fi

        # ---------------------------------------------------------------------
        # Determine layout from extracted files
        # ---------------------------------------------------------------------

        if [ -s "${accession}_1.fastq.gz" ] \
            && [ -s "${accession}_2.fastq.gz" ] \
            && gzip -t "${accession}_1.fastq.gz" 2>/dev/null \
            && gzip -t "${accession}_2.fastq.gz" 2>/dev/null
        then

            layout="PE"

            echo "[Strategy 2] NCBI extraction produced PE reads." >&2

            if ! validate_pe_pair \
                "${accession}_1.fastq.gz" \
                "${accession}_2.fastq.gz" \
                "${accession}"
            then

                echo "ERROR: ${accession}: NCBI-generated R1/R2 failed pairing validation." >&2
                exit 1

            fi

        elif [ -s "${accession}_1.fastq.gz" ] \
            && gzip -t "${accession}_1.fastq.gz" 2>/dev/null
        then

            layout="SE"

            echo "[Strategy 2] NCBI extraction produced SE reads." >&2

        else

            echo "ERROR: ${accession}: NCBI extraction produced no valid FASTQ input." >&2
            exit 1

        fi

    fi

    # =========================================================================
    # FINAL VALIDATION
    # =========================================================================

    if [ ! -s "${accession}_1.fastq.gz" ] \
        || ! gzip -t "${accession}_1.fastq.gz" 2>/dev/null
    then

        echo "ERROR: R1 missing or corrupt for ${accession}." >&2
        exit 1

    fi

    if [ "\${layout}" = "PE" ]; then

        if [ ! -s "${accession}_2.fastq.gz" ] \
            || ! gzip -t "${accession}_2.fastq.gz" 2>/dev/null
        then

            echo "ERROR: ${accession}: PE layout but R2 is missing/corrupt." >&2
            exit 1

        fi

        if ! validate_pe_pair \
            "${accession}_1.fastq.gz" \
            "${accession}_2.fastq.gz" \
            "${accession}"
        then

            echo "ERROR: ${accession}: final R1/R2 pairing validation FAILED." >&2
            exit 1

        fi

        echo "=== FETCH_SRA COMPLETE (PE): ${accession} ===" >&2

    elif [ "\${layout}" = "SE" ]; then

        rm -f "${accession}_2.fastq.gz"
        touch "${accession}_2.fastq.gz"

        echo "=== FETCH_SRA COMPLETE (SE): ${accession} ===" >&2

    else

        echo "ERROR: ${accession}: unable to determine PE/SE layout." >&2
        exit 1

    fi

    echo "=== FETCH_SRA COMPLETE: ${accession} layout=\${layout} ===" >&2

    """
}


process INGEST_LOCAL {

    tag { "${strain_id}" }

    label 'tiny'

    input:
    tuple val(strain_id), path(r1_file), path(r2_file)

    output:
    tuple val(strain_id),
          val(run_id),
          path("${run_id}_R1.fastq.gz"),
          path("${run_id}_R2.fastq.gz")

    script:

    run_id = "LOCAL_${strain_id}_${r1_file.baseName}".replaceAll(/[^A-Za-z0-9_.-]/, "_")

    """
    set -euo pipefail

    cp "${r1_file}" "${run_id}_R1.fastq.gz"

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

    if [ -s "${r1}" ] && gzip -t "${r1}" 2>/dev/null; then
        r1_ok=1
    fi

    if [ -s "${r2}" ] && gzip -t "${r2}" 2>/dev/null; then
        r2_ok=1
    fi

    if [ "\$r1_ok" -eq 1 ]; then

        r1_len=\$(python3 - "${r1}" <<'PY'
import gzip
import sys

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
import gzip
import sys

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

            status="SE_FALLBACK"
            reason="R1_R2_IDENTICAL"

        elif [ "\$r2_len" -lt ${minR2} ]; then

            status="SE_FALLBACK"
            reason="R2_TOO_SHORT"

        else

            status="PE_PASS"
            reason="OK"

        fi

    elif [ "\$r1_ok" -eq 1 ] && [ "\$r2_ok" -eq 0 ]; then

        if [ "\$is_local" -eq 1 ]; then

            status="SE_INPUT"
            reason="LOCAL_R1_ONLY_OR_BAD_R2"

        else

            status="SE_FALLBACK"
            reason="R2_INVALID_OR_MISSING"

        fi

    else

        status="DROP"
        reason="R1_INVALID"

    fi

    printf '%s\\t%s\\t%s\\t%s\\n' \
        "\$status" \
        "\$reason" \
        "\$r1_len" \
        "\$r2_len"

    """
}


process MERGE_PE {

    tag { "${strain_id}" }

    label 'base'

    publishDir path: { "${params.outdir}/Strains/${strain_id}" }, mode: 'symlink'

    input:
    tuple val(strain_id), path(r1_files), path(r2_files)

    output:
    tuple val(strain_id),
          path("${strain_id}_PE_R1.fastq.gz"),
          path("${strain_id}_PE_R2.fastq.gz"),
          emit: merged_pe

    script:

    /*
     * Deterministic ordering.
     *
     * R1/R2 filenames contain the same run_id, so sorting both sides by
     * filename ensures identical run ordering when multiple runs belong
     * to one strain.
     */

    def r1_sorted = r1_files.sort { f -> f.name }
    def r2_sorted = r2_files.sort { f -> f.name }

    if (r1_sorted.size() != r2_sorted.size()) {
        error "[MERGE_PE] ${strain_id}: R1 has ${r1_sorted.size()} files, R2 has ${r2_sorted.size()} files"
    }

    """
    set -euo pipefail

    echo "[MERGE_PE] ${strain_id}: merging ${r1_sorted.size()} PE runs"

    cat ${r1_sorted.join(' ')} > "${strain_id}_PE_R1.fastq.gz"
    cat ${r2_sorted.join(' ')} > "${strain_id}_PE_R2.fastq.gz"

    gzip -t "${strain_id}_PE_R1.fastq.gz"
    gzip -t "${strain_id}_PE_R2.fastq.gz"

    echo "[MERGE_PE] ${strain_id}: merge completed successfully"

    """
}


process MERGE_SE {

    tag { "${strain_id}" }

    label 'base'

    publishDir path: { "${params.outdir}/Strains/${strain_id}" }, mode: 'symlink'

    input:
    tuple val(strain_id), path(r1_files)

    output:
    tuple val(strain_id),
          path("${strain_id}_SE.fastq.gz"),
          emit: merged_se

    script:
    """
    set -euo pipefail

    cat ${r1_files.join(' ')} > "${strain_id}_SE.fastq.gz"

    gzip -t "${strain_id}_SE.fastq.gz"

    """
}


process CREATE_STRAIN_LIST {

    label 'tiny'

    publishDir path: "${params.outdir}/Strains/", mode: 'copy'

    input:
    val strain_ids

    output:
    path "strains.txt"

    script:

    def list_str = strain_ids
        .collect { sid -> sid.toString().trim() }
        .findAll { sid -> sid }
        .unique()
        .sort()
        .join('\n')

    """
    cat > strains.txt << 'EOF'
${list_str}
EOF
    """
}


process WRITE_RUN_QC {

    label 'tiny'

    publishDir path: "${params.outdir}/Strains/", mode: 'copy'

    input:
    val rows

    output:
    path "run_qc.tsv"

    script:

    def row_list =
        (rows instanceof java.util.Collection &&
         !rows.isEmpty() &&
         rows[0] instanceof java.util.Collection)
        ? rows
        : [rows]

    def body = row_list.collect { r ->

        def values = (r instanceof java.util.Collection) ? r : [r]

        [values[0], values[1], values[2], values[3], values[6], values[7]]
            .join('\t')

    }.join('\n')

    """
    cat > run_qc.tsv << 'EOF'
strain_id\\trun_id\\tstatus\\treason\\tr1_mean_len\\tr2_mean_len
${body}
EOF
    """
}


process WRITE_STRAIN_QC {

    label 'tiny'

    publishDir path: "${params.outdir}/Strains/", mode: 'copy'

    input:
    val rows

    output:
    path "strain_qc.tsv"

    script:

    def row_list =
        (rows instanceof java.util.Collection &&
         !rows.isEmpty() &&
         rows[0] instanceof java.util.Collection)
        ? rows
        : [rows]

    def body = row_list.collect { r ->

        def values = (r instanceof java.util.Collection) ? r : [r]

        [values[0], values[1], values[2], values[3], values[4]]
            .join('\t')

    }.join('\n')

    """
    cat > strain_qc.tsv << 'EOF'
strain_id\\tpe_runs\\tse_runs\\tdropped_runs\\tfinal_mode
${body}
EOF
    """
}