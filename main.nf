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
    MERGE_STRAIN_BAMS
    MARKDUP
    MAPPING_STATS
    WRITE_MAPPING_SUMMARY
} from './modules/map_reads.nf'

include {
    COMPUTE_COVERAGE
} from './modules/coverage.nf'

include {
    PREPARE_GATK_REF
    HAPLOTYPE_CALLER
    GATHER_STRAIN_GVCFS
} from './modules/haplotype_calling.nf'

include {
    GENOMICSDB_IMPORT
    GENOTYPE_GVCFS
    MERGE_VCFS
} from './modules/joint_genotyping.nf'

include {
    FILTER_SNPS
    FILTER_INDELS
    MERGE_AND_CLEAN
} from './modules/variant_filtration.nf'

include {
    PREP_ALLELES
    RUN_IBS
} from './modules/phylogenomy.nf'

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
    // MODULE 3: BWA Mapping, Strain-Level BAM Merging & Duplicate Marking
    // =========================================================================

    // --- Reference genome preparation ---
    if (params.reference) {
        ref_fasta = channel.fromPath(params.reference, checkIfExists: true)
    } else if (params.genome_id) {
        FETCH_REFERENCE(channel.of(params.genome_id))
        ref_fasta = FETCH_REFERENCE.out.fasta
    } else {
        error "ERROR: Must provide either params.reference (local FASTA path) or params.genome_id (NCBI accession)"
    }

    // --- Index reference (smart: detects if already indexed) ---
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
    MAP_PE(
        FILTER_PE.out.filtered_pe,
        indexed_reference.collect()
    )

    // --- Map SE reads ---
    MAP_SE(
        FILTER_SE.out.filtered_se,
        indexed_reference.collect()
    )

    // =========================================================================
    // IMPORTANT: Reconcile PE and SE mappings at the strain level
    // =========================================================================
    //
    // PE_ONLY  -> one PE BAM
    // SE_ONLY  -> one SE BAM
    // MIXED    -> one PE BAM + one SE BAM -> merged into ONE strain BAM
    //
    // The grouping is performed by strain_id so that exactly one BAM enters
    // MARKDUP for each biological strain.
    // =========================================================================

    all_sorted_bams = MAP_PE.out.sorted_bam
        .mix(MAP_SE.out.sorted_bam)
        .groupTuple(by: 0)

    MERGE_STRAIN_BAMS(all_sorted_bams)

    // --- Mark duplicates ---
    // Exactly ONE strain-level sorted BAM per strain enters MARKDUP.
    MARKDUP(MERGE_STRAIN_BAMS.out.merged_bam)

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

    // =========================================================================
    // MODULE 5: Variant Calling (GATK HaplotypeCaller)
    // =========================================================================

    // Extract just the FASTA from Module 3's reference channel
    gatk_fasta_ch = indexed_reference.map { fasta, _idx -> fasta }

    // Prepare the exact .fai and .dict bundle GATK needs
    PREPARE_GATK_REF(gatk_fasta_ch)

    // Dynamically extract all chromosome names directly from the FASTA index (.fai)
    // and filter based on the extrachromosomal master switch
    chromosomes_ch = PREPARE_GATK_REF.out.gatk_ref
        .map { _fasta, fai, _dict -> fai }
        .splitCsv(sep: '\t')
        .map { row -> row[0] }
        .filter { chrom ->
            if (!params.include_extrachromosomal) {
                // Drop any chromosome ending in _M or _P if the switch is false
                return !(chrom.endsWith('_M') || chrom.endsWith('_P'))
            }
            return true
        }

    // Scatter: 1 strain-level BAM × N chromosomes
    hc_input = MARKDUP.out.markdup_bam.combine(chromosomes_ch)

    HAPLOTYPE_CALLER(
        hc_input,
        PREPARE_GATK_REF.out.gatk_ref.collect()
    )

    // Gather: Group outputs by strain_id and generate manifest
    GATHER_STRAIN_GVCFS(
        HAPLOTYPE_CALLER.out.gvcf.groupTuple(by: 0)
    )

    // =========================================================================
    // MODULE 6: JOINT GENOTYPING
    // =========================================================================

    // Group Module 5 output by chromosome
    ch_for_genomicsdb = HAPLOTYPE_CALLER.out.gvcf
        .map { _strain_id, chrom, gvcf, tbi -> tuple(chrom, gvcf, tbi) }
        .groupTuple(by: 0)

    // 6A: Build or Update the Database
    GENOMICSDB_IMPORT(
        ch_for_genomicsdb,
        params.genomicsdb_update_path ?: ""
    )

    // 6B: Call the final variants using the Database from 6A
    GENOTYPE_GVCFS(
        GENOMICSDB_IMPORT.out.db,
        PREPARE_GATK_REF.out.gatk_ref.collect()
    )

    // 6C: Intelligent Subgenome Merging (Sensu Stricto Universal)
    if (params.merge_strategy == 'separate') {
        ch_for_merge = GENOTYPE_GVCFS.out.cohort_vcf
            .map { chrom, vcf, tbi ->
                // Split the chromosome name by the underscore
                def parts = chrom.tokenize('_')
                def prefix = parts[0]
                def suffix = parts.size() > 1 ? parts[1] : ""

                // 1. Map the prefix to the species name
                def species_name = params.species_map.get(prefix, "${prefix}_other")
                def final_prefix = species_name

                // 2. Apply the extrachromosomal isolation logic if requested
                if (params.extrachromosomal_strategy == 'isolated') {
                    if (suffix.startsWith('M')) {
                        final_prefix = "${species_name}_mitochondria"
                    } else if (suffix.startsWith('P')) {
                        final_prefix = "${species_name}_plasmid"
                    } else {
                        final_prefix = "${species_name}_nuclear"
                    }
                }

                tuple(final_prefix, vcf, tbi)
            }
            .groupTuple(by: 0)
    } else {
        // Group all chromosomes together under one name
        ch_for_merge = GENOTYPE_GVCFS.out.cohort_vcf
            .map { _chrom, vcf, tbi -> tuple('sensu_stricto_combined', vcf, tbi) }
            .groupTuple(by: 0)
    }

    // Pass the intelligently grouped channel to Picard
    MERGE_VCFS(ch_for_merge)

    // =========================================================================
    // MODULE 7: VARIANT FILTRATION & QC
    // =========================================================================

    // Load the intervals file
    ch_intervals = file(params.repetitive_intervals, checkIfExists: true)

    // Extract the individual reference files from Module 5's prepared GATK bundle
    ch_fasta_clean = PREPARE_GATK_REF.out.gatk_ref
        .map { fasta, _fai, _dict -> fasta }
        .first()

    ch_fai_clean = PREPARE_GATK_REF.out.gatk_ref
        .map { _fasta, fai, _dict -> fai }
        .first()

    ch_dict_clean = PREPARE_GATK_REF.out.gatk_ref
        .map { _fasta, _fai, dict -> dict }
        .first()

    // 1. Launch SNP and INDEL filtering simultaneously on the merged cohort VCFs
    FILTER_SNPS(
        MERGE_VCFS.out,
        ch_fasta_clean,
        ch_fai_clean,
        ch_dict_clean,
        ch_intervals
    )

    FILTER_INDELS(
        MERGE_VCFS.out,
        ch_fasta_clean,
        ch_fai_clean,
        ch_dict_clean,
        ch_intervals
    )

    // 2. Join SNP and INDEL filtering outputs by cohort name
    ch_filtered_joined = FILTER_SNPS.out.filtered_snps
        .join(FILTER_INDELS.out.filtered_indels, by: 0)

    // 3. Merge and clean the surviving variants
    MERGE_AND_CLEAN(
        ch_filtered_joined,
        ch_fasta_clean,
        ch_fai_clean,
        ch_dict_clean
    )

    // =========================================================================
    // MODULE 8: PHYLOGENOMICS (IBSx & Tree Building)
    // =========================================================================

    if (params.run_phylogeny) {
        // Explicitly pass the bin/ directory to avoid Docker symlink crashes
        ch_custom_scripts = file("${projectDir}/bin")

        // 1. Slice the clean VCFs into allele tables
        PREP_ALLELES(MERGE_AND_CLEAN.out.final_vcf)

        // 2. Generate IBS Distance Matrices using the C++ engine
        RUN_IBS(
            PREP_ALLELES.out.allele_tables,
            ch_custom_scripts
        )
    }
}