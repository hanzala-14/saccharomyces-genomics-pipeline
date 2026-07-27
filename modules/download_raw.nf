process FETCH_SRA {
    tag { "${strain_id} - ${accession}" }
    label 'base'

    errorStrategy { task.exitStatus in [1, 143, 137, 255] ? 'retry' : 'finish' }
    maxRetries 3

    input:
    tuple val(strain_id), val(accession)

    output:
    tuple val(strain_id), val(accession), path("${accession}_1.fastq.gz"), path("${accession}_2.fastq.gz")

    script:
    // Validate in Groovy (safe), not bash regex
    if( !(accession ==~ /^[A-Z]{3}\d{6,8}$/) ) {
        error "Invalid accession format: ${accession}"
    }

    def len = accession.size()
    def accPrefix = accession[0..5]

    def subdir = ''
    if (len == 9) {
        subdir = ''
    } else if (len == 10) {
        subdir = "00${accession[-1]}"
    } else if (len == 11) {
        subdir = "0${accession[-2..-1]}"
    } else if (len == 12) {
        subdir = accession[-3..-1]
    } else {
        error "Unsupported accession length for ${accession}"
    }

    def baseUrl = subdir ?
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/${accPrefix}/${subdir}/${accession}" :
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/${accPrefix}/${accession}"

    """
    set -euo pipefail

    CURL_OPTS="-f -L --retry 5 --retry-delay 2 --retry-max-time 60"

    echo "Downloading ${accession}_1.fastq.gz and ${accession}_2.fastq.gz ..."
    curl \$CURL_OPTS "${baseUrl}/${accession}_1.fastq.gz" -o "${accession}_1.fastq.gz" &
    PID1=\$!
    curl \$CURL_OPTS "${baseUrl}/${accession}_2.fastq.gz" -o "${accession}_2.fastq.gz" &
    PID2=\$!

    wait "\$PID1"; S1=\$?
    wait "\$PID2"; S2=\$?

    if [ "\$S1" -ne 0 ] || [ "\$S2" -ne 0 ]; then
        echo "Download failed for ${accession}: R1 status=\$S1, R2 status=\$S2" >&2
        exit 4
    fi

    test -s "${accession}_1.fastq.gz"
    test -s "${accession}_2.fastq.gz"
    gzip -t "${accession}_1.fastq.gz"
    gzip -t "${accession}_2.fastq.gz"
    """
}

process MERGE_RAW {
    tag { "${strain_id}" }
    label 'base'

    publishDir path: { "${params.outdir}/strains/${strain_id}" }, mode: 'copy'

    input:
    tuple val(strain_id), path(r1_files), path(r2_files)

    output:
    tuple val(strain_id), path("${strain_id}_R1.fastq.gz"), path("${strain_id}_R2.fastq.gz"), emit: merged_reads

    script:
    """
    set -euo pipefail
    cat ${r1_files.join(' ')} > "${strain_id}_R1.fastq.gz"
    cat ${r2_files.join(' ')} > "${strain_id}_R2.fastq.gz"
    gzip -t "${strain_id}_R1.fastq.gz"
    gzip -t "${strain_id}_R2.fastq.gz"
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
    def cleaned = strain_ids.collect { sid -> sid.toString().trim() }
                            .findAll { sid -> sid }
                            .unique()
                            .sort()
    def list_str = cleaned.join('\n')
    """
    cat > strains.txt << 'EOF'
    ${list_str}
    EOF
    """
}