/*
================================================================================
Module 3: BWA Mapping, Strain-Level BAM Merging, Duplicate Marking & Alignment QC
--------------------------------------------------------------------------------
Maps filtered reads to a reference genome, merges PE/SE alignments belonging
to the same biological strain, marks duplicates, and produces
sorted/indexed BAMs with comprehensive QC stats.

Inputs from Module 2:
  FILTER_PE.out.filtered_pe -> (strain_id, filt.R1.fastq.gz, filt.R2.fastq.gz)
  FILTER_SE.out.filtered_se -> (strain_id, filt.SE.fastq.gz)

Reference genome:
  - Local FASTA (params.reference) — preferred
  - Auto-fetch from NCBI (params.genome_id) — fallback if no local ref

Features:
  - Piped bwa->view->sort (no intermediate SAM/BAM on disk)
  - Pre-indexed reference detection (skips BWA index if .bwt exists)
  - Post-mapping non-empty BAM validation
  - Strain-level BAM reconciliation:
      * PE_ONLY  -> one PE BAM passed through
      * SE_ONLY  -> one SE BAM passed through
      * MIXED    -> PE + SE BAMs merged into one strain BAM
  - Duplicate marking performed once per strain
  - Proper read group assignment per mapping branch
  - Cohort-wide mapping summary TSV
================================================================================
*/

process FETCH_REFERENCE {
    tag "${genome_id}"
    label 'base'

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    val(genome_id)

    output:
    path("*.fasta"), emit: fasta

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Resolve NCBI FTP path from assembly accession
    # Example: GCF_000146045.2 → GCF/000/146/045/GCF_000146045.2
    PREFIX=\$(echo "${genome_id}" | cut -c1-3)
    NUM=\$(echo "${genome_id}" | sed 's/^[A-Z]*_//' | sed 's/\\..*//')
    D1=\$(echo "\$NUM" | cut -c1-3)
    D2=\$(echo "\$NUM" | cut -c4-6)
    D3=\$(echo "\$NUM" | cut -c7-9)

    FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/\${PREFIX}/\${D1}/\${D2}/\${D3}"

    # Find the latest assembly directory
    ASSEMBLY_DIR=\$(curl -sL "\${FTP_BASE}/" | grep -oP '${genome_id}[^"<]*' | sort -V | tail -1)

    if [ -z "\$ASSEMBLY_DIR" ]; then
        echo "ERROR: Could not resolve NCBI FTP path for ${genome_id}" >&2
        exit 1
    fi

    FASTA_URL="\${FTP_BASE}/\${ASSEMBLY_DIR}/\${ASSEMBLY_DIR}_genomic.fna.gz"

    echo "[FETCH_REFERENCE] Downloading: \${FASTA_URL}"
    curl -fSL --retry 5 --retry-delay 5 "\${FASTA_URL}" -o genome.fna.gz

    # Decompress and rename
    gunzip genome.fna.gz
    mv genome.fna "${genome_id}.fasta"

    echo "[FETCH_REFERENCE] Downloaded ${genome_id} (\$(grep -c '^>' ${genome_id}.fasta) sequences)"
    """
}

process BWA_INDEX {
    tag "${fasta.simpleName}"
    label 'medium'

    publishDir "${params.outdir}/reference", mode: 'copy', enabled: params.save_reference

    input:
    path(fasta)

    output:
    tuple path(fasta), path("${fasta}.*"), emit: indexed_ref

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    echo "[BWA_INDEX] Indexing ${fasta}"

    # BWA index
    bwa index ${fasta}

    echo "[BWA_INDEX] Indexing complete: \$(ls ${fasta}.*)"
    """
}

process MAP_PE {
    tag "${strain_id}"
    label 'low'

    input:
    tuple val(strain_id), path(r1), path(r2)
    tuple path(fasta), path(index_files)

    output:
    tuple val(strain_id), path("${strain_id}.PE.sorted.bam"), emit: sorted_bam

    script:
    // Unique read-group ID while keeping the biological sample name in SM.
    def rg = "@RG\\tID:${strain_id}_${r1.simpleName}\\tSM:${strain_id}\\tPL:ILLUMINA\\tLB:${strain_id}_lib1\\tPU:${strain_id}_unit1"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    bwa mem \\
        -t ${task.cpus} \\
        -R '${rg}' \\
        ${fasta} \\
        ${r1} ${r2} \\
    | samtools view -Sb -F 4 -@ 2 - \\
    | samtools sort -@ ${task.cpus} -m 2G -o ${strain_id}.PE.sorted.bam -

    READS=\$(samtools view -c ${strain_id}.PE.sorted.bam)

    if [ "\$READS" -eq 0 ]; then
        echo "ERROR: ${strain_id} MAP_PE produced empty BAM (0 mapped reads)" >&2
        exit 1
    fi

    echo "[MAP_PE] ${strain_id}: \${READS} mapped reads in sorted PE BAM"
    """
}

process MAP_SE {
    tag "${strain_id}"
    label 'low'

    input:
    tuple val(strain_id), path(se)
    tuple path(fasta), path(index_files)

    output:
    tuple val(strain_id), path("${strain_id}.SE.sorted.bam"), emit: sorted_bam

    script:
    def rg = "@RG\\tID:${strain_id}_${se.simpleName}\\tSM:${strain_id}\\tPL:ILLUMINA\\tLB:${strain_id}_lib1\\tPU:${strain_id}_unit1"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    bwa mem \\
        -t ${task.cpus} \\
        -R '${rg}' \\
        ${fasta} \\
        ${se} \\
    | samtools view -Sb -F 4 -@ 2 - \\
    | samtools sort -@ ${task.cpus} -m 2G -o ${strain_id}.SE.sorted.bam -

    READS=\$(samtools view -c ${strain_id}.SE.sorted.bam)

    if [ "\$READS" -eq 0 ]; then
        echo "ERROR: ${strain_id} MAP_SE produced empty BAM (0 mapped reads)" >&2
        exit 1
    fi

    echo "[MAP_SE] ${strain_id}: \${READS} mapped reads in sorted SE BAM"
    """
}


/*
================================================================================
NEW PROCESS: STRAIN-LEVEL BAM RECONCILIATION
--------------------------------------------------------------------------------
Receives all mapped BAMs belonging to the same strain.

Examples:
  PE_ONLY:
      strain.PE.sorted.bam
          ↓
      strain.sorted.bam

  SE_ONLY:
      strain.SE.sorted.bam
          ↓
      strain.sorted.bam

  MIXED:
      strain.PE.sorted.bam
      strain.SE.sorted.bam
          ↓
      samtools merge
          ↓
      strain.sorted.bam

This guarantees exactly ONE sorted BAM per biological strain before
duplicate marking.
================================================================================
*/
process MERGE_STRAIN_BAMS {
    tag "${strain_id}"
    label 'medium'

    input:
    tuple val(strain_id), path(sorted_bams)

    output:
    tuple val(strain_id), path("${strain_id}.sorted.bam"), emit: merged_bam

    script:
    def bam_files = sorted_bams instanceof List ? sorted_bams : [sorted_bams]
    def bam_args  = bam_files.collect { bam -> "\"${bam.name}\"" }.join(' ')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    echo "[MERGE_STRAIN_BAMS] ${strain_id}: ${bam_files.size()} mapping BAM(s) received"

    if [ "${bam_files.size()}" -eq 1 ]; then

        # PE_ONLY or SE_ONLY:
        # No actual merge is necessary.
        cp ${bam_args} "${strain_id}.sorted.bam"

        echo "[MERGE_STRAIN_BAMS] ${strain_id}: single BAM passed through"

    else

        # MIXED:
        # Combine the PE and SE alignment files into one strain-level BAM.
        echo "[MERGE_STRAIN_BAMS] ${strain_id}: merging ${bam_files.size()} BAMs"

        samtools merge \\
            -@ ${task.cpus} \\
            -f \\
            "${strain_id}.sorted.bam" \\
            ${bam_args}

        echo "[MERGE_STRAIN_BAMS] ${strain_id}: BAM merge complete"

    fi

    # Validate the resulting BAM.
    if [ ! -s "${strain_id}.sorted.bam" ]; then
        echo "ERROR: ${strain_id} strain-level BAM is empty" >&2
        exit 1
    fi

    samtools quickcheck "${strain_id}.sorted.bam"

    READS=\$(samtools view -c "${strain_id}.sorted.bam")

    if [ "\$READS" -eq 0 ]; then
        echo "ERROR: ${strain_id} strain-level BAM contains 0 reads" >&2
        exit 1
    fi

    echo "[MERGE_STRAIN_BAMS] ${strain_id}: \${READS} reads in final strain-level BAM"
    """
}


process MARKDUP {
    tag "${strain_id}"
    label 'medium'

    publishDir "${params.outdir}/mapped/strains", mode: 'copy', pattern: '*.mdup.bam*'
    publishDir "${params.outdir}/mapped/stats",   mode: 'copy', pattern: '*.metrics.txt'

    input:
    tuple val(strain_id), path(sorted_bam)

    output:
    tuple val(strain_id), path("${strain_id}.mdup.bam"), path("${strain_id}.mdup.bam.bai"), emit: markdup_bam
    tuple val(strain_id), path("${strain_id}.mdup.metrics.txt"), emit: dup_metrics

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # 1. Group by read name.
    samtools collate -o name_collate.bam ${sorted_bam}

    # 2. Add mate score tags.
    samtools fixmate -m name_collate.bam fixmate.bam

    # 3. Coordinate sort again.
    samtools sort -o coord_sorted.bam fixmate.bam

    # 4. Mark duplicates.
    samtools markdup -s coord_sorted.bam ${strain_id}.mdup.bam 2> ${strain_id}.mdup.metrics.txt

    # 5. Index the final BAM.
    samtools index ${strain_id}.mdup.bam

    # Clean up intermediate files to save space.
    rm name_collate.bam fixmate.bam coord_sorted.bam

    READS=\$(samtools view -c -F 1024 ${strain_id}.mdup.bam)

    echo "[MARKDUP] ${strain_id}: \${READS} non-duplicate reads, BAM indexed"
    """
}


process MAPPING_STATS {
    tag "${strain_id}"
    label 'tiny'

    publishDir "${params.outdir}/mapped/stats", mode: 'copy'

    input:
    tuple val(strain_id), path(bam), path(bai)

    output:
    tuple val(strain_id), path("${strain_id}.flagstat.txt"), emit: flagstat
    tuple val(strain_id), path("${strain_id}.idxstats.txt"), emit: idxstats
    tuple val(strain_id), path("${strain_id}.stats.txt"), emit: full_stats

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    samtools flagstat ${bam} > ${strain_id}.flagstat.txt
    samtools idxstats ${bam} > ${strain_id}.idxstats.txt
    samtools stats ${bam} > ${strain_id}.stats.txt

    echo "[MAPPING_STATS] ${strain_id}: stats collected"
    """
}

process WRITE_MAPPING_SUMMARY {
    tag "mapping_summary"
    label 'tiny'

    publishDir "${params.outdir}/mapped", mode: 'copy'

    input:
    val(flagstat_files)
    val(metrics_files)
    val(stats_files)

    output:
    path "mapping_summary.tsv", emit: summary

    script:
    """
    #!/usr/bin/env python3
    import os, re

    header = "strain_id\\ttotal_reads\\tmapped_reads\\tmapped_pct\\tproperly_paired_pct\\tsingletons_pct\\tduplicate_pct\\tinsert_size_mean\\tstatus"
    rows = []

    flagstat_list = "${flagstat_files}".split()
    metrics_list  = "${metrics_files}".split()
    stats_list    = "${stats_files}".split()

    min_mapping_pct = float("${params.min_mapping_pct}")

    for f in flagstat_list:
        f = f.strip().rstrip(',').strip('[').strip(']')
        if not f or not os.path.isfile(f):
            continue

        strain = os.path.basename(f).replace('.flagstat.txt', '')
        total = mapped = paired = singletons = 0

        with open(f) as fh:
            for line in fh:
                if 'in total' in line:
                    total = int(line.split()[0])
                elif 'mapped (' in line and 'primary mapped' not in line:
                    mapped = int(line.split()[0])
                elif 'properly paired' in line:
                    paired = int(line.split()[0])
                elif 'singletons' in line:
                    singletons = int(line.split()[0])

        mapped_pct = (mapped / total * 100) if total > 0 else 0
        paired_pct = (paired / total * 100) if total > 0 else 0
        single_pct = (singletons / total * 100) if total > 0 else 0

        # Parse duplicate rate from metrics file.
        dup_pct = 0.0
        metrics_f = f.replace('.flagstat.txt', '.mdup.metrics.txt')

        for mf in metrics_list:
            mf = mf.strip().rstrip(',').strip('[').strip(']')

            if mf and os.path.isfile(mf) and strain in os.path.basename(mf):
                with open(mf) as mfh:
                    in_metrics = False

                    for line in mfh:
                        if line.startswith('LIBRARY'):
                            in_metrics = True
                            continue

                        if in_metrics and not line.startswith('#') and line.strip():
                            cols = line.strip().split('\\t')

                            if len(cols) >= 9:
                                try:
                                    dup_pct = float(cols[8]) * 100
                                except (ValueError, IndexError):
                                    pass
                            break
                break

        # Parse insert size from samtools stats.
        insert_mean = "NA"

        for sf in stats_list:
            sf = sf.strip().rstrip(',').strip('[').strip(']')

            if sf and os.path.isfile(sf) and strain in os.path.basename(sf):
                with open(sf) as sfh:
                    for line in sfh:
                        if line.startswith('SN\\tinsert size average:'):
                            insert_mean = line.strip().split('\\t')[2]
                            break
                break

        # Status determination.
        status = "PASS"

        if min_mapping_pct > 0 and mapped_pct < min_mapping_pct:
            status = "WARN_LOW_MAPPING"

        rows.append(
            f"{strain}\\t{total}\\t{mapped}\\t{mapped_pct:.2f}\\t"
            f"{paired_pct:.2f}\\t{single_pct:.2f}\\t{dup_pct:.2f}\\t"
            f"{insert_mean}\\t{status}"
        )

    with open('mapping_summary.tsv', 'w') as out:
        out.write(header + '\\n')

        for r in sorted(rows):
            out.write(r + '\\n')

    # Report warnings.
    warn_count = sum(1 for r in rows if 'WARN' in r)

    print(
        f"[WRITE_MAPPING_SUMMARY] "
        f"{len(rows)} strains processed, "
        f"{warn_count} warnings"
    )
    """
}