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
    WRITE_COVERAGE_SUMMARY
    PLOT_COVERAGE
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

include {
    EXTRACT_ALLELE_ROH
    PLOT_ALLELE_ROH
} from './modules/allele_roh.nf'

include {
    PREPARE_SUBGENOME_MATRICES
    COMPARE_SUBGENOME_DISTANCES
    EXTRACT_NUCLEOTIDE_DIVERSITY
    PLOT_NUCLEOTIDE_DIVERSITY
} from './modules/subgenome_diversity.nf'

workflow {

    // Serialize the 8-species Universal Sensu Stricto Dictionary for Python/R scripts
    def species_map_str = params.species_map.collect { k, v -> "${k}:${v}" }.join(",")
    ch_custom_scripts   = file("${projectDir}/bin")

    // =========================================================================
    // SHORTCUT MODE A: Standalone Coverage-Only Execution
    // Triggered when both --coverage_input_dir and --only_coverage are set
    // =========================================================================
     if (params.coverage_input_dir && params.only_coverage) {
    log.info "[STANDALONE MODE] Running Module 4 Coverage Master Excel & Curated Plots from: ${params.coverage_input_dir}"

    ch_all_cov = channel
        .fromPath("${params.coverage_input_dir}/*.slidingwindow.tab", checkIfExists: true)
        .collect()

    if (params.coverage_excel) {
        // Stage 2: User provided their edited Excel file -> Batch generate all 161 R plots!
        PLOT_COVERAGE(
            ch_all_cov,
            file(params.coverage_excel, checkIfExists: true),
            "xlsx",
            ch_custom_scripts,
            species_map_str,
            params.ploidy_multiplier
        )
    } else {
        // Stage 1: Generate the formula-driven Master Excel Workbook for manual G24 curation
        WRITE_COVERAGE_SUMMARY(
            ch_all_cov,
            ch_custom_scripts,
            species_map_str,
            params.ploidy_multiplier,
            params.auto_adjust_ploidy
        )
    }
}

    // =========================================================================
    // SHORTCUT MODE B: Standalone Pairwise Subgenome Distance Comparison
    // Triggered when --comparison_csv or --matrix_input_dir is provided
    // =========================================================================
    else if (params.comparison_csv || params.matrix_input_dir) {
        log.info "[STANDALONE MODE] Running Universal Subgenome Distance Comparison"
        def excl_str = (params.exclude_outlier_strains ?: []).join(",")

        if (params.comparison_csv) {
            ch_comp_in   = channel.fromPath(params.comparison_csv, checkIfExists: true)
            val_comp_mod = "csv"
        } else {
            ch_comp_in   = channel.fromPath(params.matrix_input_dir, checkIfExists: true)
            val_comp_mod = "matrices"
        }

        COMPARE_SUBGENOME_DISTANCES(
            ch_comp_in,
            val_comp_mod,
            ch_custom_scripts,
            species_map_str,
            excl_str
        )
    }
    // =========================================================================
    // SHORTCUT MODE C: Standalone Allele & ROH Plots from VCF Directory
    // =========================================================================
    else if (params.vcf_input_dir && params.only_allele_roh) {
        log.info "[STANDALONE MODE] Running Allele & ROH Plots from VCF directory: ${params.vcf_input_dir}"

        ch_standalone_vcfs = channel
            .fromPath("${params.vcf_input_dir}/*.vcf.gz", checkIfExists: true)
            .map { vcf ->
                def cohort_id = vcf.name.replaceAll(/\.vcf\.gz$/, '')
                def tbi = file("${vcf}.tbi")
                def tbi_arg = tbi.exists() ? tbi : file("NO_TBI")
                tuple(cohort_id, vcf, tbi_arg)
            }

        ch_genome_fallback = params.genome_file ? file(params.genome_file) : file("NO_GENOME")
        ch_target_strains  = (params.target_strains && file(params.target_strains.toString()).exists())
            ? file(params.target_strains.toString())
            : file("NO_STRAIN_FILE")
        val_target_strains = params.target_strains ? params.target_strains.toString() : "ALL"

        EXTRACT_ALLELE_ROH(
            ch_standalone_vcfs,
            ch_genome_fallback,
            ch_target_strains,
            val_target_strains
        )

        PLOT_ALLELE_ROH(
            EXTRACT_ALLELE_ROH.out.extracted_bundle,
            ch_custom_scripts,
            species_map_str
        )
    }

    // =========================================================================
    // SHORTCUT MODE D: Standalone Nucleotide Diversity from VCF Directory
    // =========================================================================
    else if (params.vcf_input_dir && params.only_diversity) {
        log.info "[STANDALONE MODE] Running Nucleotide Diversity (Pi, Theta_W, Tajima's D) from: ${params.vcf_input_dir}"

        ch_standalone_vcfs = channel
            .fromPath("${params.vcf_input_dir}/*.vcf.gz", checkIfExists: true)
            .map { vcf ->
                def cohort_id = vcf.name.replaceAll(/\.vcf\.gz$/, '')
                def tbi = file("${vcf}.tbi")
                def tbi_arg = tbi.exists() ? tbi : file("NO_TBI")
                tuple(cohort_id, vcf, tbi_arg)
            }

        ch_genome_fallback = params.genome_file ? file(params.genome_file) : file("NO_GENOME")
        ch_target_strains  = (params.target_strains && file(params.target_strains.toString()).exists())
            ? file(params.target_strains.toString())
            : file("NO_STRAIN_FILE")
        val_target_strains = params.target_strains ? params.target_strains.toString() : "ALL"

        EXTRACT_NUCLEOTIDE_DIVERSITY(
         ch_standalone_vcfs,
         ch_genome_fallback,
         ch_target_strains,
         val_target_strains,
         params.clade_name ?: "",
         params.pi_window_size,
         params.pi_step_size
     )

        PLOT_NUCLEOTIDE_DIVERSITY(
            EXTRACT_NUCLEOTIDE_DIVERSITY.out.diversity_bundle,
            ch_custom_scripts,
            species_map_str
        )
    }

    // =========================================================================
    // FULL PIPELINE MODE: Modules 1 through 8
    // =========================================================================
    else {

        // =====================================================================
        // MODULE 1: Data Acquisition, QC Classification & Strain Aggregation
        // =====================================================================

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
                    error "[samplesheet.csv] Invalid accession='${accession}' (expected ^[A-Z]{3}[0-9]{6,8}\$). strain_id=${strain_raw}"
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

        // =====================================================================
        // MODULE 2: Read Filtration (fastp)
        // =====================================================================

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

        // =====================================================================
        // MODULE 3: BWA Mapping, Strain-Level BAM Merging & Duplicate Marking
        // =====================================================================

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

        // =====================================================================
        // IMPORTANT: Reconcile PE and SE mappings at the strain level
        // =====================================================================

        all_sorted_bams = MAP_PE.out.sorted_bam
            .mix(MAP_SE.out.sorted_bam)
            .groupTuple(by: 0)

        MERGE_STRAIN_BAMS(all_sorted_bams)

        // --- Mark duplicates ---
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

        // =====================================================================
        // MODULE 4: Coverage Profiling (Conditional Execution)
        // =====================================================================
        if (params.run_coverage) {
    COMPUTE_COVERAGE(
        MARKDUP.out.markdup_bam,
        file(params.genome_file, checkIfExists: true),
        file(params.sliding_windows, checkIfExists: true)
    )

    ch_all_cov = COMPUTE_COVERAGE.out.sliding_window_cov
        .map { _strain_id, tab_file -> tab_file }
        .collect()

    WRITE_COVERAGE_SUMMARY(
        ch_all_cov,
        ch_custom_scripts,
        species_map_str,
        params.ploidy_multiplier,
        params.auto_adjust_ploidy
    )
}

        // =====================================================================
        // MODULE 5: Variant Calling (GATK HaplotypeCaller)
        // =====================================================================

        gatk_fasta_ch = indexed_reference.map { fasta, _idx -> fasta }

        PREPARE_GATK_REF(gatk_fasta_ch)

        chromosomes_ch = PREPARE_GATK_REF.out.gatk_ref
            .map { _fasta, fai, _dict -> fai }
            .splitCsv(sep: '\t')
            .map { row -> row[0] }
            .filter { chrom ->
                if (!params.include_extrachromosomal) {
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

        GATHER_STRAIN_GVCFS(
            HAPLOTYPE_CALLER.out.gvcf.groupTuple(by: 0)
        )

        // =====================================================================
        // MODULE 6: JOINT GENOTYPING
        // =====================================================================

        ch_for_genomicsdb = HAPLOTYPE_CALLER.out.gvcf
            .map { _strain_id, chrom, gvcf, tbi -> tuple(chrom, gvcf, tbi) }
            .groupTuple(by: 0)

        GENOMICSDB_IMPORT(
            ch_for_genomicsdb,
            params.genomicsdb_update_path ?: ""
        )

        GENOTYPE_GVCFS(
            GENOMICSDB_IMPORT.out.db,
            PREPARE_GATK_REF.out.gatk_ref.collect()
        )

        if (params.merge_strategy == 'separate') {
            ch_for_merge = GENOTYPE_GVCFS.out.cohort_vcf
                .map { chrom, vcf, tbi ->
                    def parts = chrom.tokenize('_')
                    def prefix = parts[0]
                    def suffix = parts.size() > 1 ? parts[1] : ""

                    def species_name = params.species_map.get(prefix, "${prefix}_other")
                    def final_prefix = species_name

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
            ch_for_merge = GENOTYPE_GVCFS.out.cohort_vcf
                .map { _chrom, vcf, tbi -> tuple('sensu_stricto_combined', vcf, tbi) }
                .groupTuple(by: 0)
        }

        MERGE_VCFS(ch_for_merge)

        // =====================================================================
        // MODULE 7: VARIANT FILTRATION & QC
        // =====================================================================

        ch_intervals = file(params.repetitive_intervals, checkIfExists: true)

        ch_fasta_clean = PREPARE_GATK_REF.out.gatk_ref
            .map { fasta, _fai, _dict -> fasta }
            .first()

        ch_fai_clean = PREPARE_GATK_REF.out.gatk_ref
            .map { _fasta, fai, _dict -> fai }
            .first()

        ch_dict_clean = PREPARE_GATK_REF.out.gatk_ref
            .map { _fasta, _fai, dict -> dict }
            .first()

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

        ch_filtered_joined = FILTER_SNPS.out.filtered_snps
            .join(FILTER_INDELS.out.filtered_indels, by: 0)

        MERGE_AND_CLEAN(
            ch_filtered_joined,
            ch_fasta_clean,
            ch_fai_clean,
            ch_dict_clean
        )

        // =====================================================================
        // MODULE 7.5: ALLELE RATIO & ROH PROFILING
        // =====================================================================
        if (params.run_allele_roh) {
            ch_genome_fallback = params.genome_file ? file(params.genome_file) : file("NO_GENOME")

            EXTRACT_ALLELE_ROH(
                MERGE_AND_CLEAN.out.final_vcf,
                ch_genome_fallback,
                file("NO_STRAIN_FILE"),
                "ALL"
            )

            PLOT_ALLELE_ROH(
                EXTRACT_ALLELE_ROH.out.extracted_bundle,
                ch_custom_scripts,
                species_map_str
            )
        }

        // =====================================================================
        // MODULE 7.6: NUCLEOTIDE DIVERSITY (Pi, Theta_W, Tajima's D)
        // =====================================================================
        if (params.run_diversity) {
         ch_genome_fallback_div = params.genome_file ? file(params.genome_file) : file("NO_GENOME")
         ch_target_strains_div  = (params.target_strains && file(params.target_strains.toString()).exists())
             ? file(params.target_strains.toString())
             : file("NO_STRAIN_FILE")
         val_target_strains_div = params.target_strains ? params.target_strains.toString() : "ALL"

         EXTRACT_NUCLEOTIDE_DIVERSITY(
             MERGE_AND_CLEAN.out.final_vcf,
             ch_genome_fallback_div,
             ch_target_strains_div,
             val_target_strains_div,
             params.clade_name ?: "",
             params.pi_window_size,
             params.pi_step_size
         )

         PLOT_NUCLEOTIDE_DIVERSITY(
             EXTRACT_NUCLEOTIDE_DIVERSITY.out.diversity_bundle,
             ch_custom_scripts,
             species_map_str
         )
     }

        // =====================================================================
        // MODULE 8: PHYLOGENOMICS (IBSx & Subgenome Distance Comparison)
        // =====================================================================
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

           // 3. Universal Sensu Stricto Subgenome Pairwise Distance Comparison
            if (params.run_subgenome_comp) {
                PREPARE_SUBGENOME_MATRICES(
                    RUN_IBS.out.distance_matrices,
                    species_map_str
                )

                COMPARE_SUBGENOME_DISTANCES(
                    PREPARE_SUBGENOME_MATRICES.out.tagged_matrices.collect(),
                    "matrices",
                    ch_custom_scripts,
                    species_map_str,
                    (params.exclude_outlier_strains ?: []).join(",")
                )
            }
        }
    }
}