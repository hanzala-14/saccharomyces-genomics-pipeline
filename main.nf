nextflow.enable.dsl=2

include {
    FETCH_SRA
    INGEST_LOCAL
    CLASSIFY_RUN
    MERGE_PE
    MERGE_SE
    CREATE_STRAIN_LIST
    WRITE_RUN_QC
    WRITE_STRAIN_QC
} from './modules/download_raw.nf'

include {
    FILTER_PE
    FILTER_SE
    WRITE_FILTER_SUMMARY
} from './modules/filter_reads.nf'

include {
    FETCH_REFERENCE
    BWA_INDEX
    MAP_PE
    MAP_SE
    MARKDUP
    MAPPING_STATS
    WRITE_MAPPING_SUMMARY
} from './modules/map_reads.nf'

include {
    COMPUTE_COVERAGE
} from './modules/coverage.nf'

workflow {

    // =========================================================================
    // MODULE 1: Data Acquisition, QC Classification & Strain Aggregation
    // =========================================================================

    rows_ch = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def strain_raw = (row.strain_id ?: '').toString().trim()
            def accession  = (row.accession ?: '').toString().trim()
            def r1         = (row.r1 ?: '').toString().trim()
            def r2         = (row.r2 ?: '').toString().trim()

            if (!strain_raw) {
                error "[samplesheet.csv] Empty strain_id. Row=${row}"
            }

            def strain_id = strain_raw.replaceAll(/[^A-Za-z0-9_.-]/, "_")
            def hasAcc = accession != ''
            def hasR1  = r1 != ''
            def hasR2  = r2 != ''

            if (!hasAcc && !hasR1) {
                error "[samplesheet.csv] Must provide either accession or r1. strain_id=${strain_raw}"
            }

            if (hasAcc && (hasR1 || hasR2)) {
                error "[samplesheet.csv] Cannot mix accession with local r1/r2 in same row. strain_id=${strain_raw}"
            }

            if (hasAcc && !(accession ==~ /^[A-Z]{3}[0-9]{6,8}$/)) {
                error "[samplesheet.csv] Invalid accession='${accession}' (expected ^[A-Z]{3}[0-9]{6,8}$). strain_id=${strain_raw}"
            }

            if (!hasAcc && !hasR1) {
                error "[samplesheet.csv] Local mode requires r1 path. strain_id=${strain_raw}"
            }

            if (hasAcc) {
                tuple('ENA', strain_id, accession, '', '')
            } else {
                tuple('LOCAL', strain_id, '', r1, (hasR2 ? r2 : ''))
            }
        }

    ena_rows = rows_ch
        .filter { source, _strain_id, _accession, _r1, _r2 -> source == 'ENA' }
        .map    { _source, strain_id, accession, _r1, _r2 -> tuple(strain_id, accession) }

    local_rows = rows_ch
        .filter { source, _strain_id, _accession, _r1, _r2 -> source == 'LOCAL' }
        .map    { _source, strain_id, _accession, r1, r2 ->
            def r1_file = file(r1, checkIfExists: true)
            def r2_file = r2 ? file(r2, checkIfExists: true) : file('EMPTY')
            tuple(strain_id, r1_file, r2_file)
        }

    FETCH_SRA(ena_rows)
    INGEST_LOCAL(local_rows)

    unified_runs = FETCH_SRA.out.mix(INGEST_LOCAL.out)
    CLASSIFY_RUN(unified_runs)

    classified_rows = CLASSIFY_RUN.out.qc_rows
        .map { strain, run_id, summary, r1, r2 ->
            def cols = summary.toString().trim().split('\t')
            if (cols.size() != 4) {
                error "Invalid CLASSIFY_RUN summary '${summary}'"
            }
            tuple(strain, run_id, cols[0], cols[1], r1, r2, cols[2] as Integer, cols[3] as Integer)
        }

    // Branch into PE, SE, and DROP streams — each item goes to exactly one branch
    classified_rows.branch { row ->
        pe:   row[2] == 'PE_PASS'
        se:   row[2] == 'SE_FALLBACK' || row[2] == 'SE_INPUT'
        drop: true
    }.set { branched }

    // Collect ALL classified rows for QC reporting
    all_for_qc = branched.pe.mix(branched.se, branched.drop)

    WRITE_RUN_QC(all_for_qc.toList())

    // --- PE merge ---
    pe_merge_in = branched.pe.map { strain, run_id, _status, _reason, r1, r2, _r1len, _r2len ->
        tuple(strain, run_id, r1, r2)
    }

    pe_grouped = pe_merge_in
        .groupTuple(by: 0)
        .map { strain, run_ids, r1s, r2s ->
            def idx = (0..<run_ids.size()).toList().sort { i -> run_ids[i] }
            tuple(strain, idx.collect { index -> r1s[index] }, idx.collect { index -> r2s[index] })
        }

    MERGE_PE(pe_grouped)

    // --- SE merge ---
    se_merge_in = branched.se.map { strain, run_id, _status, _reason, r1, _r2, _r1len, _r2len ->
        tuple(strain, run_id, r1)
    }

    se_grouped = se_merge_in
        .groupTuple(by: 0)
        .map { strain, run_ids, r1s ->
            def idx = (0..<run_ids.size()).toList().sort { i -> run_ids[i] }
            tuple(strain, idx.collect { index -> r1s[index] })
        }

    MERGE_SE(se_grouped)

    // --- Strain list ---
    strain_ids_ch = channel.empty()
        .mix(MERGE_PE.out.merged_pe.map { strain, _r1, _r2 -> strain })
        .mix(MERGE_SE.out.merged_se.map { strain, _r1 -> strain })
        .collect()

    CREATE_STRAIN_LIST(strain_ids_ch)

    // --- Strain QC ---
    strain_qc_rows = all_for_qc
        .map { strain, _run_id, status, _reason, _r1, _r2, _r1len, _r2len ->
            tuple(strain, status)
        }
        .groupTuple(by: 0)
        .map { strain, statuses ->
            def pe = statuses.count { status -> status == 'PE_PASS' }
            def se = statuses.count { status -> status == 'SE_FALLBACK' || status == 'SE_INPUT' }
            def dr = statuses.count { status -> status == 'DROP' }
            def final_mode = (pe > 0 && se > 0) ? 'MIXED' : (pe > 0 ? 'PE_ONLY' : (se > 0 ? 'SE_ONLY' : 'ALL_DROPPED'))
            tuple(strain, pe, se, dr, final_mode)
        }
        .toList()

    WRITE_STRAIN_QC(strain_qc_rows)

    // =========================================================================
    // MODULE 2: Read Filtration (fastp)
    // =========================================================================

    // Both channels are safe to pass even if empty — Nextflow will simply
    // not schedule the process if no items arrive. No hang, no error.
    FILTER_PE(MERGE_PE.out.merged_pe)
    FILTER_SE(MERGE_SE.out.merged_se)

    // --- Filtration Summary Report ---
    pe_jsons = FILTER_PE.out.json_report
        .map { _strain_id, json -> json }
        .collect()
        .ifEmpty([])

    se_jsons = FILTER_SE.out.json_report
        .map { _strain_id, json -> json }
        .collect()
        .ifEmpty([])

    WRITE_FILTER_SUMMARY(pe_jsons, se_jsons)

    // =========================================================================
    // MODULE 3: BWA Mapping, Duplicate Marking & Alignment QC
    // =========================================================================

    // --- Reference genome preparation ---
    // Two modes: local reference (preferred) or auto-fetch from NCBI
    if (params.reference) {
        ref_fasta = channel.fromPath(params.reference, checkIfExists: true)
    } else if (params.genome_id) {
        FETCH_REFERENCE(channel.of(params.genome_id))
        ref_fasta = FETCH_REFERENCE.out.fasta
    } else {
        error "ERROR: Must provide either params.reference (local FASTA path) or params.genome_id (NCBI accession)"
    }

    // --- Index reference (smart: detects if already indexed) ---
    // Check if BWA index files already exist alongside the FASTA
    ref_needs_index = ref_fasta.map { fasta ->
        def bwt = file("${fasta}.bwt")
        if (bwt.exists()) {
            log.info "[MODULE 3] BWA index found for ${fasta.name} — skipping indexing"
        }
        tuple(fasta, bwt.exists())
    }

    // Branch: already indexed vs needs indexing
    ref_needs_index.branch { _fasta, indexed ->
        indexed:     indexed == true
        not_indexed: indexed == false
    }.set { ref_branched }

    // For already-indexed: collect existing index files
    ref_already_indexed = ref_branched.indexed.map { fasta, _flag ->
        def idx_files = files("${fasta}.*")
        tuple(fasta, idx_files)
    }

    // For not-indexed: run BWA_INDEX
    BWA_INDEX(ref_branched.not_indexed.map { fasta, _flag -> fasta })

    // Merge into single reference channel
    indexed_reference = ref_already_indexed.mix(BWA_INDEX.out.indexed_ref)

    // --- Map PE reads ---
    MAP_PE(FILTER_PE.out.filtered_pe, indexed_reference.collect())
    MAP_SE(FILTER_SE.out.filtered_se, indexed_reference.collect())

    // --- Merge PE and SE BAM streams into unified stream ---
    all_sorted_bams = MAP_PE.out.sorted_bam.mix(MAP_SE.out.sorted_bam)

    // --- Mark duplicates ---
    // Reads are already merged by strain in Module 1, so each strain is a single BAM emission
    MARKDUP(all_sorted_bams)

    // --- Collect mapping stats ---
    MAPPING_STATS(MARKDUP.out.markdup_bam)

    // --- Mapping summary report ---
    all_flagstats = MAPPING_STATS.out.flagstat
        .map { _strain_id, f -> f }
        .collect()
        .ifEmpty([])

    all_dup_metrics = MARKDUP.out.dup_metrics
        .map { _strain_id, f -> f }
        .collect()
        .ifEmpty([])

    all_stats = MAPPING_STATS.out.full_stats
        .map { _strain_id, f -> f }
        .collect()
        .ifEmpty([])

    WRITE_MAPPING_SUMMARY(all_flagstats, all_dup_metrics, all_stats)

    // =========================================================================
    // MODULE 4: Coverage Profiling (Conditional Execution)
    // =========================================================================
    if (params.run_coverage) {
        COMPUTE_COVERAGE(
            MARKDUP.out.markdup_bam,
            file(params.genome_file, checkIfExists: true),
            file(params.sliding_windows, checkIfExists: true)
        )
    }
}