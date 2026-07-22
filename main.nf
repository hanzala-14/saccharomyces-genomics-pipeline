nextflow.enable.dsl=2

include { FETCH_SRA; MERGE_RAW; CREATE_STRAIN_LIST } from './modules/download_raw.nf'

workflow {
    // 1. Read CSV and validate column names
    samples_ch = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def s_id = row.strain_id ?: row.strain ?: row.Strain ?: row.STRAIN
            def acc  = row.accession ?: row.Accession ?: row.ACCESSION
            
            if (!s_id || !acc) {
                error "CSV Error: Could not find 'strain_id' or 'accession' columns in ${params.samplesheet}"
            }
            return tuple(s_id.trim(), acc.trim())
        }

    // Step 1: Download each accession individually in parallel
    FETCH_SRA(samples_ch)

    // Step 2: Group downloaded reads by strain_id
    grouped_reads_ch = FETCH_SRA.out
        .groupTuple(by: 0) // Groups (strain_id, r1_path, r2_path)

    // Step 3: Merge reads per strain
    MERGE_RAW(grouped_reads_ch)

    // Step 4: Collect strain IDs into strains.txt
    MERGE_RAW.out.merged_reads
        .map { strain_id, _r1, _r2 -> strain_id }
        .collect()
        .set { all_strains_ch }

    CREATE_STRAIN_LIST(all_strains_ch)

    // Console output logger
    MERGE_RAW.out.merged_reads.view { strain_id, file_R1, file_R2 ->
        "🚀 MERGE COMPLETE: ${strain_id} -> [ ${file_R1.name}, ${file_R2.name} ]"
    }
}