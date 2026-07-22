process FETCH_SRA {
    tag { "${strain_id} - ${accession}" }
    cpus 2
    memory 4.GB

    errorStrategy { task.exitStatus in [1, 143, 137, 255] ? 'retry' : 'finish' }
    maxRetries 3

    input:
    tuple val(strain_id), val(accession)

    output:
    tuple val(strain_id), path("${accession}_1.fastq.gz"), path("${accession}_2.fastq.gz")

    script:
    """
    LEN=\$(echo -n "${accession}" | wc -c)
    ACC_PREFIX=\$(echo "${accession}" | cut -c 1-6)

    if [ \$LEN -eq 9 ]; then
        BASE_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/\${ACC_PREFIX}/${accession}"
    elif [ \$LEN -eq 10 ]; then
        SUBDIR="00\$(echo "${accession}" | tail -c 2)"
        BASE_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/\${ACC_PREFIX}/\${SUBDIR}/${accession}"
    elif [ \$LEN -eq 11 ]; then
        SUBDIR="0\$(echo "${accession}" | tail -c 3)"
        BASE_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/\${ACC_PREFIX}/\${SUBDIR}/${accession}"
    else
        SUBDIR="\$(echo "${accession}" | tail -c 4)"
        BASE_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/\${ACC_PREFIX}/\${SUBDIR}/${accession}"
    fi

    CURL_OPTS="-# --retry 5 --retry-delay 2 --retry-max-time 60 -L"

    echo "Downloading ${accession}_1.fastq.gz & ${accession}_2.fastq.gz..."
    curl \$CURL_OPTS "\${BASE_URL}/${accession}_1.fastq.gz" -o "${accession}_1.fastq.gz" &
    PID1=\$!
    curl \$CURL_OPTS "\${BASE_URL}/${accession}_2.fastq.gz" -o "${accession}_2.fastq.gz" &
    PID2=\$!

    wait \$PID1 \$PID2
    """
}

process MERGE_RAW {
    tag { "${strain_id}" }
    cpus 2
    memory 4.GB

    publishDir path: { "results/strains/${strain_id}" }, mode: 'copy'

    input:
    tuple val(strain_id), path(r1_files), path(r2_files)

    output:
    tuple val(strain_id), path("${strain_id}_R1.fastq.gz"), path("${strain_id}_R2.fastq.gz"), emit: merged_reads

    script:
    """
    echo "=== Merging run FASTQs for strain: ${strain_id} ==="
    cat ${r1_files} > "${strain_id}_R1.fastq.gz"
    cat ${r2_files} > "${strain_id}_R2.fastq.gz"

    gzip -t "${strain_id}_R1.fastq.gz"
    gzip -t "${strain_id}_R2.fastq.gz"
    """
}

process CREATE_STRAIN_LIST {
    publishDir path: { "results" }, mode: 'copy'

    input:
    val strain_ids

    output:
    path "strains.txt"

    script:
    def list_str = strain_ids.unique().join("\n")
    """
    echo "${list_str}" > strains.txt
    """
}