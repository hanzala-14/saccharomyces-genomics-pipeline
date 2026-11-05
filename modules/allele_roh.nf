nextflow.enable.dsl = 2

/*
 * Module 7.5: Allele Frequency Ratio & Runs of Homozygosity (ROH) Profiling
 * -------------------------------------------------------------------
 *  EXTRACT_ALLELE_ROH & PLOT_ALLELE_ROH with:
 * - Parallel 8-Chromosome Streaming Engine: Executes concurrent `bcftools query` and `bcftools roh` workers throttled in batches of 8 to maximize throughput while preventing BusyBox/container file-descriptor exhaustion.
 * - Universal Sensu Stricto Subgenome Detection: Automatically partitions multi-species chromosomes (`c_`, `e_`, `u_`, `k_`, `p_`, `m_`, `a_`, `j_`) and sorts nuclear contigs in canonical Roman numeral order (`I`–`XVI`), stripping mitochondrial (`_M`) and plasmid (`_P`) contigs.
 * - Dynamic Sample Subsetting: Supports full-cohort extraction, comma-separated strain lists, or external `.txt` strain files (`--target_strains`) without re-genotyping or splitting the master VCF.
 * - Symlink & Truncation Hardening: Automatically refreshes stale `.tbi` index symlinks inside isolated container mounts and shields downstream R rendering from zero-variant chromosomes.
 * - High-Resolution Karyogram & BAF Rendering: Generates per-strain B-allele frequency (AD ratio) scatter tracks overlaid with HMM/Viterbi-inferred ROH blocks and exports a multi-sheet cohort `.xlsx` summary.
 */
process EXTRACT_ALLELE_ROH {
    tag { "${cohort_id}" }
    label 'medium'

    publishDir path: { "${params.outdir}/08_allele_roh/${cohort_id}/extracted_data" }, mode: 'copy'

    input:
    tuple val(cohort_id), path(vcf), path(tbi)
    path genome_file
    path target_strains_file
    val target_strains_val

    output:
    tuple val(cohort_id),
          path("subgenomes/*/*.chr_length.csv"),
          path("subgenomes/*/allele_csvs/*.alleleplot.csv"),
          path("subgenomes/*/roh_interest/*.roh.interest.txt"), emit: extracted_bundle, optional: true
    path "subgenomes/*/roh_full/*.roh.txt", emit: roh_full, optional: true

    script:
    """
    set -euo pipefail

    # Raise open file descriptor limit if permitted by container
    ulimit -n \$(ulimit -Hn) 2>/dev/null || true

    # 1. Index VCF if .tbi/.csi is missing; if .tbi symlink is older than .vcf.gz, safely copy & touch
    if [ ! -f "${vcf}.tbi" ] && [ ! -f "${vcf}.csi" ]; then
        bcftools index -t --threads ${task.cpus} "${vcf}"
    elif [ -f "${vcf}.tbi" ] && [ "${vcf}" -nt "${vcf}.tbi" ]; then
        REAL_TBI=\$(readlink -f "${vcf}.tbi")
        rm -f "${vcf}.tbi"
        cp "\$REAL_TBI" "${vcf}.tbi"
        touch "${vcf}.tbi"
    fi

    # 2. Extract sample list (strip CRLF) and filter if --target_strains was provided
    bcftools query -l "${vcf}" 2>/dev/null | tr -d '\\r' > vcf_all_samples.txt || true

    if [ "${target_strains_file.name}" != "NO_STRAIN_FILE" ] && [ -f "${target_strains_file}" ]; then
        tr -d '\\r' < "${target_strains_file}" | tr ',' '\\n' | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//' | grep -v '^#' | grep -v '^\$' > requested_samples.txt || true
        awk 'NR==FNR {req[\$1]=1; next} (\$1 in req)' requested_samples.txt vcf_all_samples.txt > samples.txt
    elif [ "${target_strains_val}" != "ALL" ] && [ -n "${target_strains_val}" ]; then
        printf '%s\\n' "${target_strains_val}" | tr -d '\\r' | tr ',' '\\n' | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//' | grep -v '^\$' > requested_samples.txt || true
        awk 'NR==FNR {req[\$1]=1; next} (\$1 in req)' requested_samples.txt vcf_all_samples.txt > samples.txt
    else
        cp vcf_all_samples.txt samples.txt
    fi

    if [ ! -s samples.txt ]; then
        echo "WARN: None of the requested --target_strains were found in ${vcf}. Skipping cleanly." >&2
        exit 0
    fi

    # 3. Extract nuclear contig lengths from VCF header (excluding _M and _P)
    bcftools view -h "${vcf}" 2>/dev/null | tr -d '\\r' | awk -F'[=,>]' '
    /^##contig=/ {
        id=""; len=""
        for (i=1; i<=NF; i++) {
            if (\$i == "ID") id=\$(i+1)
            if (\$i == "length") len=\$(i+1)
        }
        if (id != "" && len != "" && id !~ /_(M|P)\$/) {
            print id "\\t" len
        }
    }' > raw_header_contigs.tsv

    if [ ! -s raw_header_contigs.tsv ] && [ -f "${genome_file}" ] && [ "${genome_file.name}" != "NO_GENOME" ]; then
        tr -d '\\r' < "${genome_file}" | awk '\$1 !~ /_(M|P)\$/ {print \$1 "\\t" \$2}' > raw_header_contigs.tsv
    fi

    # 4. Instant chromosome discovery via tabix -l
    > active_chrs.txt
    if [ -f "${vcf}.tbi" ] && command -v tabix >/dev/null 2>&1; then
        tabix -l "${vcf}" 2>/dev/null | tr -d '\\r' | awk '\$1 !~ /_(M|P)\$/ {print \$1}' > active_chrs.txt || true
    fi

    if [ ! -s active_chrs.txt ]; then
        while IFS=\$'\\t' read -r chr _len; do
            [ -z "\$chr" ] && continue
            first_hit=\$(bcftools view -H -r "\$chr" "${vcf}" 2>/dev/null | head -n 1 || true)
            if [ -n "\$first_hit" ]; then
                echo "\$chr" >> active_chrs.txt
            fi
        done < raw_header_contigs.tsv
    fi

    if [ ! -s active_chrs.txt ]; then
        echo "WARN: No active nuclear chromosomes found in ${vcf}. Skipping cleanly." >&2
        exit 0
    fi

    if [ ! -s raw_header_contigs.tsv ]; then
        { bcftools query -f '%CHROM\\t%POS\\n' "${vcf}" 2>/dev/null || true; } | awk '
        \$1 !~ /_(M|P)\$/ {
            if (\$2 > max_pos[\$1]) max_pos[\$1] = \$2
        }
        END {
            for (c in max_pos) print c "\\t" max_pos[c]
        }' > raw_header_contigs.tsv
    fi

    # Sort active chromosomes in Roman numeral order (I..XVI)
    awk '
    BEGIN {
        split("I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI", rom, " ")
        for (i=1; i<=16; i++) r_ord[rom[i]] = i
    }
    NR==FNR { active[\$1] = 1; next }
    (\$1 in active) {
        chr = \$1
        len = \$2
        prefix = "c"
        label = chr
        if (match(chr, /^[^_]+_/)) {
            prefix = substr(chr, 1, RLENGTH - 1)
            label = substr(chr, RLENGTH + 1)
        }
        ord = (label in r_ord) ? r_ord[label] : 99
        printf "%s\\t%02d\\t%s,%s,%s\\n", prefix, ord, chr, len, label
    }' active_chrs.txt raw_header_contigs.tsv | sort -k1,1 -k2,2n > sorted_contigs.tsv

    cut -f1 sorted_contigs.tsv | sort -u > prefixes.txt

    # 5. Parallel Per-Chromosome Extraction & ROH Calculation (8 chromosomes per batch)
    while IFS= read -r pref; do
        [ -z "\$pref" ] && continue
        mkdir -p "subgenomes/\${pref}/allele_csvs"
        mkdir -p "subgenomes/\${pref}/roh_full"
        mkdir -p "subgenomes/\${pref}/roh_interest"
        mkdir -p "tmp_parallel_\${pref}"

        {
            echo "chr,length,label"
            awk -v p="\$pref" '\$1 == p {print \$3}' sorted_contigs.tsv
        } > "subgenomes/\${pref}/\${pref}.chr_length.csv"

        awk -v p="\$pref" '\$1 == p {split(\$3, a, ","); print a[1]}' sorted_contigs.tsv > "tmp_parallel_\${pref}/chr_list.txt"

        c_idx=0
        while IFS= read -r chr; do
            [ -z "\$chr" ] && continue
            c_idx=\$((c_idx + 1))
            c_tag=\$(printf "%02d" "\$c_idx")
            c_dir="tmp_parallel_\${pref}/chr_\${c_tag}"
            mkdir -p "\$c_dir"

            (
                # A. Fast direct per-chromosome AD query with lightweight 64-line buffer (max 1 open FD!)
                { bcftools query -S samples.txt -r "\$chr" -f '%CHROM;%POS[;%AD]\\n' "${vcf}" 2>/dev/null || true; } \\
                | awk -F';' -v outdir="\$c_dir" '
                BEGIN {
                    idx = 0
                    while ((getline s < "samples.txt") > 0) {
                        gsub(/^[ \\t\\r]+|[ \\t\\r]+\$/, "", s)
                        if (s == "" || s ~ /^#/) continue
                        idx++
                        safe = s
                        gsub(/[^A-Za-z0-9._-]+/, "_", safe)
                        out_files[idx] = outdir "/" safe ".part"
                        buf[idx] = ""
                        cnt[idx] = 0
                    }
                    close("samples.txt")
                    n_samples = idx
                }
                {
                    c = \$1; pos = \$2
                    for (i = 1; i <= n_samples; i++) {
                        ad = \$(i + 2)
                        if (ad != "." && ad != "" && ad != "0,0" && ad !~ /^0(,0)+\$/) {
                            buf[i] = buf[i] c ";" pos ";" ad "\\n"
                            cnt[i]++
                            if (cnt[i] >= 64) {
                                printf "%s", buf[i] >> out_files[i]
                                close(out_files[i])
                                buf[i] = ""
                                cnt[i] = 0
                            }
                        }
                    }
                }
                END {
                    for (i = 1; i <= n_samples; i++) {
                        if (cnt[i] > 0) {
                            printf "%s", buf[i] >> out_files[i]
                            close(out_files[i])
                        }
                    }
                }'

                # B. Per-chromosome ROH HMM calculation (2 threads per chromosome)
                bcftools roh --threads 2 -S samples.txt -r "\$chr" "${vcf}" -o "\$c_dir/roh.part" 2>/dev/null
            ) &

            # Strictly POSIX batching: run 8 chromosomes concurrently, then wait
            if [ \$((c_idx % 8)) -eq 0 ]; then
                wait
            fi
        done < "tmp_parallel_\${pref}/chr_list.txt"
        wait

        # Concatenate per-chromosome AD parts in exact Roman numeral order per sample (POSIX safe)
        while IFS= read -r sample || [ -n "\$sample" ]; do
            sample=\$(printf '%s' "\$sample" | tr -d '\\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//')
            [ -z "\$sample" ] && continue
            case "\$sample" in \\#*) continue ;; esac
            safe=\$(printf '%s' "\$sample" | sed 's/[^A-Za-z0-9._-]/_/g')

            out_csv="subgenomes/\${pref}/allele_csvs/\${safe}.\${pref}.alleleplot.csv"
            echo "chr;pos;ad" > "\$out_csv"
            for part_file in tmp_parallel_\${pref}/chr_*/"\${safe}.part"; do
                if [ -f "\$part_file" ]; then
                    cat "\$part_file" >> "\$out_csv"
                fi
            done
        done < samples.txt

        # Concatenate per-chromosome ROH outputs and split by sample using 64-line buffer
        > "subgenomes/\${pref}/roh_full/cohort_all.roh.txt"
        for roh_part in tmp_parallel_\${pref}/chr_*/roh.part; do
            if [ -f "\$roh_part" ]; then
                cat "\$roh_part" >> "subgenomes/\${pref}/roh_full/cohort_all.roh.txt"
            fi
        done

        awk -F'\\t' -v outdir="subgenomes/\${pref}/roh_interest" -v p="\${pref}" '
        BEGIN {
            hdr = "RG\\tSample\\tChromosome\\tStart\\tEnd\\tLength\\tMarkers\\tQuality"
            while ((getline s < "samples.txt") > 0) {
                gsub(/^[ \\t\\r]+|[ \\t\\r]+\$/, "", s)
                if (s == "" || s ~ /^#/) continue
                safe = s
                gsub(/[^A-Za-z0-9._-]+/, "_", safe)
                f = outdir "/" safe "." p ".roh.interest.txt"
                s_map[s] = f
                print hdr > f
                close(f)
                buf[s] = ""
                cnt[s] = 0
            }
            close("samples.txt")
        }
        /^RG/ {
            s = \$2
            if (s in s_map) {
                buf[s] = buf[s] \$0 "\\n"
                cnt[s]++
                if (cnt[s] >= 64) {
                    printf "%s", buf[s] >> s_map[s]
                    close(s_map[s])
                    buf[s] = ""
                    cnt[s] = 0
                }
            }
        }
        END {
            for (s in s_map) {
                if (cnt[s] > 0) {
                    printf "%s", buf[s] >> s_map[s]
                    close(s_map[s])
                }
            }
        }' "subgenomes/\${pref}/roh_full/cohort_all.roh.txt"

        rm -rf "tmp_parallel_\${pref}"
    done < prefixes.txt
    """
}

process PLOT_ALLELE_ROH {
    tag { "${cohort_id}" }
    label 'phylogeny'

    stageInMode 'copy'

    publishDir path: { "${params.outdir}/08_allele_roh/${cohort_id}" }, mode: 'copy'

    input:
    tuple val(cohort_id),
          path(chr_len_files, stageAs: 'chr_lens/*'),
          path(allele_csvs,   stageAs: 'allele_csvs/*'),
          path(roh_files,     stageAs: 'roh_files/*')
    path custom_scripts
    val species_map_str

    output:
    path "allele_plots/*.png",               emit: allele_pngs,   optional: true
    path "roh_plots/*.png",                  emit: roh_pngs,      optional: true
    path "combined_overlay/*.png",           emit: combined_pngs, optional: true
    path "cohort_allele_plots.pdf",          emit: allele_pdf,    optional: true
    path "cohort_roh_plots.pdf",             emit: roh_pdf,       optional: true
    path "cohort_allele_roh_combined.pdf",   emit: combined_pdf,  optional: true
    path "cohort_allele_roh_summary.tsv",    emit: summary_tsv
    path "cohort_allele_roh_summary.xlsx",   emit: summary_xlsx,  optional: true

    script:
    """
    set -euo pipefail

    Rscript ${custom_scripts}/plot_allele_roh.R \\
        --cohort_id "${cohort_id}" \\
        --species_map "${species_map_str}" \\
        --bin_dir "${custom_scripts}"
    """
}