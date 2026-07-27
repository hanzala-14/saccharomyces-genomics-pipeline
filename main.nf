nextflow.enable.dsl=2

include { FETCH_SRA; MERGE_RAW; CREATE_STRAIN_LIST } from './modules/download_raw.nf'

workflow {
    samples_ch = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.containsKey('strain_id')) {
                error "CSV Error: Missing required column 'strain_id' in ${params.samplesheet}"
            }
            if (!row.containsKey('accession')) {
                error "CSV Error: Missing required column 'accession' in ${params.samplesheet}"
            }

            def raw_strain = (row.strain_id ?: '').toString().trim()
            def accession  = (row.accession ?: '').toString().trim()

            if (!raw_strain) {
                error "CSV Error: Empty 'strain_id' value found in ${params.samplesheet}"
            }
            if (!accession) {
                error "CSV Error: Empty 'accession' value found in ${params.samplesheet}"
            }

            // keep alnum, underscore, hyphen, dot
            def strain_id = raw_strain.replaceAll(/[^A-Za-z0-9_.-]/, "_")

            tuple(strain_id, accession)
        }

    FETCH_SRA(samples_ch)

    grouped_reads_ch = FETCH_SRA.out
        .map { strain_id, accession, r1, r2 -> tuple(strain_id, accession, r1, r2) }
        .groupTuple(by: 0)
        .map { strain_id, accessions, r1_list, r2_list ->
            def idx = (0..<accessions.size()).toList().sort { i -> accessions[i] }
            def sorted_r1 = idx.collect { i -> r1_list[i] }
            def sorted_r2 = idx.collect { i -> r2_list[i] }
            tuple(strain_id, sorted_r1, sorted_r2)
        }

    MERGE_RAW(grouped_reads_ch)

    MERGE_RAW.out.merged_reads
        .map { strain_id, _r1, _r2 -> strain_id }
        .collect()
        .set { all_strains_ch }

    CREATE_STRAIN_LIST(all_strains_ch)

    MERGE_RAW.out.merged_reads.view { strain_id, file_R1, file_R2 ->
        "MERGE COMPLETE: ${strain_id} -> [ ${file_R1.name}, ${file_R2.name} ]"
    }
}