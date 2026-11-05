nextflow.enable.dsl = 2

/*
 * Module 7.6 & 8.5: Subgenome Divergence Regression & Population Nucleotide Diversity
 * -------------------------------------------------------------------
 *  PREPARE_SUBGENOME_MATRICES, COMPARE_SUBGENOME_DISTANCES, EXTRACT_NUCLEOTIDE_DIVERSITY & PLOT_NUCLEOTIDE_DIVERSITY with:
 * - Combinatorial N-Subgenome Regression: Automatically pairs distance matrices across all present Saccharomyces sensu stricto subgenomes, computing OLS linear regression, Pearson r, R², and a 1,000-permutation empirical Mantel test p-value.
 * - Length-Weighted Chromosomal Aggregation: Combines per-chromosome IBSx distance matrices (`I`–`XVI`) weighted by physical nuclear chromosome lengths and filters incomplete/non-hybrid background noise.
 * - Exact Missing-Data-Aware Diversity Engine: Streams diploid genotypes in 8-chromosome parallel batches to compute exact site-specific allele counts, sliding-window Nucleotide Diversity (π), Watterson's Theta (θw), and Tajima's D normalized by physical window length.
 * - Multi-Clade & Single-Pass Cohort Splitting: Accepts a 2-column clade map (`Strain <tab> Clade`) or 1-column strain list via `--target_strains` to stratify diversity tracks and boxplots across multiple clades (e.g., Saaz, Frohberg, Hybrids, Novalager) in a single pass.
 * - Zero-Dependency Multi-Sheet Reporting: Produces publication-grade ggplot2 chromosomal tracks, cross-clade comparison boxplots, and native multi-sheet `.xlsx` workbooks without external Python/R Excel libraries.
 */

process PREPARE_SUBGENOME_MATRICES {
    tag { "${cohort_id}" }
    label 'low'

    input:
    tuple val(cohort_id), path(ibs_files)
    val species_map_str

    output:
    path "*.subgenome_matrix.tsv", emit: tagged_matrices, optional: true

    script:
    """
    set -euo pipefail

    python3 - "${cohort_id}" "${species_map_str}" << 'PYEOF'
import os, sys, re

cohort_id = sys.argv[1]
species_map_str = sys.argv[2]

sp_map = {}
for item in species_map_str.split(","):
    if ":" in item:
        k, v = item.split(":", 1)
        sp_map[k.strip()] = v.strip()

chr_romans = {"I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"}

def is_chrom_path(p):
    parts = re.split(r"[/\\\\_.-]+", p)
    return any(tok in chr_romans for tok in parts)

candidates = []
for root, dirs, files in os.walk("."):
    dirs[:] = [d for d in dirs if not d.startswith(".")]
    for f in files:
        if f.startswith(".") or f.endswith((".log", ".yml", ".json", ".sh", ".py", ".subgenome_matrix.tsv")):
            continue
        if f.endswith((".tsv", ".csv", ".txt", ".tab")) or any(k in f.lower() for k in ("matrix", "dist", "ibs")):
            rel = os.path.join(root, f)
            # Score: prioritize genome-wide IBSx_all / combined files outside per-chromosome folders
            is_chr = is_chrom_path(rel)
            is_all = any(w in f.lower() for w in ("ibsx_all", "combined", "final", "all", "genome"))
            score = (10 if is_chr else 0) + (0 if is_all else 5)
            candidates.append((score, rel))

if not candidates:
    sys.exit(0)

candidates.sort(key=lambda x: (x[0], x[1]))
chosen = candidates[0][1]

detected_species = cohort_id
for pref, sp_name in sp_map.items():
    if sp_name.lower() in cohort_id.lower() or sp_name.lower() in chosen.lower():
        detected_species = sp_name
        break
    if re.search(rf"(^|_){re.escape(pref)}(_|\$)", cohort_id):
        detected_species = sp_name
        break

safe_sp = re.sub(r"[^A-Za-z0-9._-]+", "_", detected_species)
safe_id = re.sub(r"[^A-Za-z0-9._-]+", "_", cohort_id)
out_name = f"{safe_sp}__{safe_id}.subgenome_matrix.tsv"

with open(chosen, "r", encoding="utf-8", errors="replace") as src, open(out_name, "w", encoding="utf-8") as dst:
    dst.write(src.read())
PYEOF
    """
}

process COMPARE_SUBGENOME_DISTANCES {
    tag { "Sensu Stricto Subgenome Regression" }
    label 'phylogeny'

    stageInMode 'copy'

    publishDir path: { "${params.outdir}/09_subgenome_comparison" }, mode: 'copy'

    input:
    path input_data
    val input_mode
    path custom_scripts
    val species_map_str
    val exclude_strains_str

    output:
    path "*.csv",  emit: csvs,  optional: true
    path "*.xlsx", emit: xlsx,  optional: true
    path "*.png",  emit: pngs,  optional: true
    path "*.pdf",  emit: pdfs,  optional: true

    script:
    def mode_arg = (input_mode == "csv") ? "--input-csv \"${input_data}\"" : "--matrix-dir \".\""
    """
    set -euo pipefail

    python3 ${custom_scripts}/compare_subgenomes.py \\
        ${mode_arg} \\
        --species-map "${species_map_str}" \\
        --exclude-strains "${exclude_strains_str}" \\
        --out-prefix "Subgenome_Distance_Comparison"
    """
}

process EXTRACT_NUCLEOTIDE_DIVERSITY {
    tag { "${cohort_id}" }
    label 'medium'

    publishDir path: { "${params.outdir}/09_nucleotide_diversity/${cohort_id}/extracted_data" }, mode: 'copy'

    input:
    tuple val(cohort_id), path(vcf), path(tbi)
    path genome_file
    path target_strains_file
    val target_strains_val
    val clade_name_val
    val window_size
    val step_size

    output:
    tuple val(cohort_id),
          path("subgenomes/*/*.chr_length.csv"),
          path("subgenomes/*/*.windows.tsv"),
          path("subgenomes/*/*.chrom_summary.tsv"), emit: diversity_bundle, optional: true

    script:
    """
    set -euo pipefail

    ulimit -n \$(ulimit -Hn) 2>/dev/null || true

    # 1. Index VCF if .tbi/.csi is missing; refresh .tbi symlink if older than .vcf.gz
    if [ ! -f "${vcf}.tbi" ] && [ ! -f "${vcf}.csi" ]; then
        bcftools index -t --threads ${task.cpus} "${vcf}"
    elif [ -f "${vcf}.tbi" ] && [ "${vcf}" -nt "${vcf}.tbi" ]; then
        REAL_TBI=\$(readlink -f "${vcf}.tbi")
        rm -f "${vcf}.tbi"
        cp "\$REAL_TBI" "${vcf}.tbi"
        touch "${vcf}.tbi"
    fi

    # 2. Extract full VCF sample list
    bcftools query -l "${vcf}" 2>/dev/null | tr -d '\\r' > vcf_all_samples.txt || true

    # 3. Build Clade -> Sample lists (Supports 2-col Clade Map, 1-col Strain File, comma list, or ALL)
    mkdir -p clade_lists
    > clades_manifest.tsv

    if [ "${target_strains_file.name}" != "NO_STRAIN_FILE" ] && [ -f "${target_strains_file}" ]; then
        # Default clade label from --clade_name or filename (e.g., Saaz_Clade.txt -> Saaz Clade)
        DEFAULT_CLADE="${clade_name_val}"
        if [ -z "\$DEFAULT_CLADE" ] || [ "\$DEFAULT_CLADE" == "null" ]; then
            DEFAULT_CLADE=\$(basename "${target_strains_file}" | sed -E 's/\\.(txt|tsv|csv)\$//I' | tr '_' ' ')
        fi

        tr -d '\\r' < "${target_strains_file}" | grep -v '^[[:space:]]*#' | grep -v '^[[:space:]]*\$' > clean_targets.txt || true

        # Check if 2-column file (Strain <tab/comma/semicolon> Clade)
        HAS_TWO_COLS=\$(awk -F'[\\t,;]+' 'NF >= 2 && tolower(\$1) != "strain" && tolower(\$1) != "sample" && tolower(\$1) != "strain_id" {print 1; exit}' clean_targets.txt || true)

        if [ "\$HAS_TWO_COLS" == "1" ]; then
            awk -F'[\\t,;]+' '
            tolower(\$1) == "strain" || tolower(\$1) == "sample" || tolower(\$1) == "strain_id" { next }
            NF >= 2 {
                s = \$1; c = \$2
                gsub(/^[[:space:]]+|[[:space:]]+\$/, "", s)
                gsub(/^[[:space:]]+|[[:space:]]+\$/, "", c)
                if (s != "" && c != "") print s "\\t" c
            }' clean_targets.txt > strain_clade_pairs.tsv

            cut -f2 strain_clade_pairs.tsv | awk '!seen[\$0]++' > unique_clades.txt
            while IFS= read -r c_label; do
                [ -z "\$c_label" ] && continue
                c_slug=\$(echo "\$c_label" | sed 's/[^A-Za-z0-9._-]/_/g')
                awk -F'\\t' -v cl="\$c_label" '\$2 == cl {print \$1}' strain_clade_pairs.tsv > "clade_lists/\${c_slug}.req"
                awk 'NR==FNR {req[\$1]=1; next} (\$1 in req)' "clade_lists/\${c_slug}.req" vcf_all_samples.txt > "clade_lists/\${c_slug}.samples.txt"
                n_s=\$(wc -l < "clade_lists/\${c_slug}.samples.txt" | tr -d ' ')
                if [ "\$n_s" -ge 2 ]; then
                    echo -e "\${c_slug}\\t\${c_label}\\tclade_lists/\${c_slug}.samples.txt\\t\${n_s}" >> clades_manifest.tsv
                else
                    echo "WARN: Clade '\$c_label' has <2 matching samples (\$n_s) in ${vcf}. Skipping clade." >&2
                fi
            done < unique_clades.txt
        else
            c_slug=\$(echo "\$DEFAULT_CLADE" | sed 's/[^A-Za-z0-9._-]/_/g')
            tr ',' '\\n' < clean_targets.txt | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//' | grep -v '^\$' > "clade_lists/\${c_slug}.req" || true
            awk 'NR==FNR {req[\$1]=1; next} (\$1 in req)' "clade_lists/\${c_slug}.req" vcf_all_samples.txt > "clade_lists/\${c_slug}.samples.txt"
            n_s=\$(wc -l < "clade_lists/\${c_slug}.samples.txt" | tr -d ' ')
            if [ "\$n_s" -ge 2 ]; then
                echo -e "\${c_slug}\\t\${DEFAULT_CLADE}\\tclade_lists/\${c_slug}.samples.txt\\t\${n_s}" >> clades_manifest.tsv
            fi
        fi
    elif [ "${target_strains_val}" != "ALL" ] && [ -n "${target_strains_val}" ]; then
        DEFAULT_CLADE="${clade_name_val}"
        if [ -z "\$DEFAULT_CLADE" ] || [ "\$DEFAULT_CLADE" == "null" ]; then
            DEFAULT_CLADE="Selected Strains"
        fi
        c_slug=\$(echo "\$DEFAULT_CLADE" | sed 's/[^A-Za-z0-9._-]/_/g')
        printf '%s\\n' "${target_strains_val}" | tr -d '\\r' | tr ',' '\\n' | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//' | grep -v '^\$' > "clade_lists/\${c_slug}.req" || true
        awk 'NR==FNR {req[\$1]=1; next} (\$1 in req)' "clade_lists/\${c_slug}.req" vcf_all_samples.txt > "clade_lists/\${c_slug}.samples.txt"
        n_s=\$(wc -l < "clade_lists/\${c_slug}.samples.txt" | tr -d ' ')
        if [ "\$n_s" -ge 2 ]; then
            echo -e "\${c_slug}\\t\${DEFAULT_CLADE}\\tclade_lists/\${c_slug}.samples.txt\\t\${n_s}" >> clades_manifest.tsv
        fi
    else
        DEFAULT_CLADE="${clade_name_val}"
        if [ -z "\$DEFAULT_CLADE" ] || [ "\$DEFAULT_CLADE" == "null" ]; then
            DEFAULT_CLADE="All Strains"
        fi
        c_slug=\$(echo "\$DEFAULT_CLADE" | sed 's/[^A-Za-z0-9._-]/_/g')
        cp vcf_all_samples.txt "clade_lists/\${c_slug}.samples.txt"
        n_s=\$(wc -l < "clade_lists/\${c_slug}.samples.txt" | tr -d ' ')
        if [ "\$n_s" -ge 2 ]; then
            echo -e "\${c_slug}\\t\${DEFAULT_CLADE}\\tclade_lists/\${c_slug}.samples.txt\\t\${n_s}" >> clades_manifest.tsv
        fi
    fi

    if [ ! -s clades_manifest.tsv ]; then
        echo "WARN: No clade/group with at least 2 valid samples found in ${vcf}. Skipping cleanly." >&2
        exit 0
    fi

    # 4. Extract nuclear contig lengths from VCF header or fallback genome_file
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

    # 5. Instant chromosome discovery via tabix -l
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
        printf "%s\\t%02d\\t%s\\t%s\\t%s\\n", prefix, ord, chr, len, label
    }' active_chrs.txt raw_header_contigs.tsv | sort -k1,1 -k2,2n > sorted_contigs.tsv

    cut -f1 sorted_contigs.tsv | sort -u > prefixes.txt

    # 6. Parallel 8-Chromosome Extraction per Subgenome Prefix and per Clade
    while IFS= read -r pref; do
        [ -z "\$pref" ] && continue
        mkdir -p "subgenomes/\${pref}"

        {
            echo "chr,length,label"
            awk -v p="\$pref" '\$1 == p {print \$3 "," \$4 "," \$5}' sorted_contigs.tsv
        } > "subgenomes/\${pref}/\${pref}.chr_length.csv"

        awk -v p="\$pref" '\$1 == p {print \$3 "\\t" \$4 "\\t" \$5}' sorted_contigs.tsv > "chr_info_\${pref}.tsv"

        while IFS=\$'\\t' read -r c_slug c_label s_file n_dip; do
            [ -z "\$c_slug" ] && continue
            tmp_dir="tmp_div_\${pref}_\${c_slug}"
            mkdir -p "\$tmp_dir"

            c_idx=0
            while IFS=\$'\\t' read -r chr chr_len chr_label; do
                [ -z "\$chr" ] && continue
                c_idx=\$((c_idx + 1))
                c_tag=\$(printf "%02d" "\$c_idx")
                c_dir="\$tmp_dir/chr_\${c_tag}"
                mkdir -p "\$c_dir"

                (
                    { bcftools query -S "\$s_file" -r "\$chr" -f '%POS[\\t%GT]\\n' "${vcf}" 2>/dev/null || true; } \\
                    | awk -F'\\t' \\
                        -v clade="\$c_label" \\
                        -v chr="\$chr" \\
                        -v label="\$chr_label" \\
                        -v chr_len="\$chr_len" \\
                        -v win_size="${window_size}" \\
                        -v step_size="${step_size}" \\
                        -v n_dip="\$n_dip" \\
                        -v win_out="\$c_dir/windows.part" \\
                        -v chr_out="\$c_dir/chrom.part" '
                    function calc_tajima_d(pi_sum, S, n_hap,   a1,a2,i,b1,b2,c1,c2,e1,e2,theta_w,var_d) {
                        if (S < 1 || n_hap < 4) return "NA"
                        a1 = 0.0; a2 = 0.0
                        for (i = 1; i < n_hap; i++) {
                            a1 += 1.0 / i
                            a2 += 1.0 / (i * i)
                        }
                        theta_w = S / a1
                        b1 = (n_hap + 1.0) / (3.0 * (n_hap - 1.0))
                        b2 = (2.0 * (n_hap * n_hap + n_hap + 3.0)) / (9.0 * n_hap * (n_hap - 1.0))
                        c1 = b1 - (1.0 / a1)
                        c2 = b2 - ((n_hap + 2.0) / (a1 * n_hap)) + (a2 / (a1 * a1))
                        e1 = c1 / a1
                        e2 = c2 / (a1 * a1 + a2)
                        var_d = e1 * S + e2 * S * (S - 1.0)
                        if (var_d <= 0) return "NA"
                        return sprintf("%.6f", (pi_sum - theta_w) / sqrt(var_d))
                    }
                    function calc_theta_bp(S, n_hap, length_bp,   a1,i) {
                        if (S < 1 || n_hap < 2 || length_bp < 1) return 0.0
                        a1 = 0.0
                        for (i = 1; i < n_hap; i++) a1 += 1.0 / i
                        return (S / a1) / length_bp
                    }
                    BEGIN {
                        if (chr_len <= win_size) {
                            n_win = 1
                        } else {
                            n_win = int((chr_len - 1) / step_size) + 1
                        }
                        for (w = 1; w <= n_win; w++) {
                            ws = (w - 1) * step_size + 1
                            we = ws + win_size - 1
                            if (we > chr_len) we = chr_len
                            w_start[w] = ws
                            w_end[w]   = we
                            w_len[w]   = we - ws + 1
                            w_seg[w]   = 0
                            w_pi[w]    = 0.0
                            w_nhap[w]  = 0
                        }
                        c_seg  = 0
                        c_pi   = 0.0
                        c_nhap = 0
                        span   = int(win_size / step_size) + 2
                    }
                    {
                        pos = \$1 + 0
                        n_called = 0
                        delete ac
                        for (i = 2; i <= NF; i++) {
                            gt = \$i
                            if (gt == "." || gt == "./." || gt == ".|.") continue
                            m = split(gt, alleles, "[/|]")
                            if (m == 1) {
                                a = alleles[1]
                                if (a != ".") {
                                    ac[a] += 2
                                    n_called += 2
                                }
                            } else {
                                for (j = 1; j <= m; j++) {
                                    a = alleles[j]
                                    if (a != ".") {
                                        ac[a]++
                                        n_called++
                                    }
                                }
                            }
                        }
                        if (n_called < 2) next
                        n_alleles = 0
                        sum_sq = 0.0
                        for (a in ac) {
                            n_alleles++
                            freq = ac[a] / n_called
                            sum_sq += freq * freq
                        }
                        if (n_alleles < 2) next

                        site_pi = (n_called / (n_called - 1.0)) * (1.0 - sum_sq)

                        c_seg++
                        c_pi += site_pi
                        c_nhap += n_called

                        w_center = int((pos - 1) / step_size) + 1
                        w_lo = w_center - span
                        if (w_lo < 1) w_lo = 1
                        w_hi = w_center + 1
                        if (w_hi > n_win) w_hi = n_win

                        for (w = w_lo; w <= w_hi; w++) {
                            if (pos >= w_start[w] && pos <= w_end[w]) {
                                w_seg[w]++
                                w_pi[w] += site_pi
                                w_nhap[w] += n_called
                            }
                        }
                    }
                    END {
                        default_hap = n_dip * 2
                        for (w = 1; w <= n_win; w++) {
                            S = w_seg[w]
                            L = w_len[w]
                            nh = (S > 0) ? int((w_nhap[w] / S) + 0.5) : default_hap
                            pi_bp = (L > 0) ? (w_pi[w] / L) : 0.0
                            th_bp = calc_theta_bp(S, nh, L)
                            td    = calc_tajima_d(w_pi[w], S, nh)
                            mid   = int((w_start[w] + w_end[w]) / 2)
                            printf "%s\\t%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%.8f\\t%.8f\\t%s\\n",
                                clade, chr, label, w, w_start[w], w_end[w], mid, L, nh, S, pi_bp, th_bp, td > win_out
                        }
                        close(win_out)

                        nh_c  = (c_seg > 0) ? int((c_nhap / c_seg) + 0.5) : default_hap
                        pi_c  = (chr_len > 0) ? (c_pi / chr_len) : 0.0
                        th_c  = calc_theta_bp(c_seg, nh_c, chr_len)
                        td_c  = calc_tajima_d(c_pi, c_seg, nh_c)
                        printf "%s\\t%s\\t%s\\t%d\\t%d\\t%d\\t%.4f\\t%.8f\\t%.8f\\t%s\\n",
                            clade, chr, label, chr_len, nh_c, c_seg, c_pi, pi_c, th_c, td_c > chr_out
                        close(chr_out)
                    }'
                ) &

                if [ \$((c_idx % 8)) -eq 0 ]; then
                    wait
                fi
            done < "chr_info_\${pref}.tsv"
            wait

            win_file="subgenomes/\${pref}/\${pref}__\${c_slug}.windows.tsv"
            echo -e "Clade\\tChromosome\\tLabel\\tWindow_Idx\\tStart\\tEnd\\tMidpoint\\tWindow_Length\\tCalled_Haplotypes\\tSegregating_Sites\\tPi_per_bp\\tTheta_W_per_bp\\tTajima_D" > "\$win_file"
            for wp in "\$tmp_dir"/chr_*/windows.part; do
                if [ -f "\$wp" ]; then
                    cat "\$wp" >> "\$win_file"
                fi
            done

            chr_file="subgenomes/\${pref}/\${pref}__\${c_slug}.chrom_summary.tsv"
            echo -e "Clade\\tChromosome\\tLabel\\tLength_bp\\tCalled_Haplotypes\\tSegregating_Sites\\tPi_Total_Sum\\tPi_per_bp\\tTheta_W_per_bp\\tTajima_D" > "\$chr_file"
            for cp in "\$tmp_dir"/chr_*/chrom.part; do
                if [ -f "\$cp" ]; then
                    cat "\$cp" >> "\$chr_file"
                fi
            done

            rm -rf "\$tmp_dir"
        done < clades_manifest.tsv
    done < prefixes.txt
    """
}

process PLOT_NUCLEOTIDE_DIVERSITY {
    tag { "${cohort_id}" }
    label 'phylogeny'

    stageInMode 'copy'

    publishDir path: { "${params.outdir}/09_nucleotide_diversity/${cohort_id}" }, mode: 'copy'

    input:
    tuple val(cohort_id),
          path(chr_len_files, stageAs: 'chr_lens/*'),
          path(window_tsvs,   stageAs: 'window_tsvs/*'),
          path(chrom_tsvs,    stageAs: 'chrom_tsvs/*')
    path custom_scripts
    val species_map_str

    output:
    path "plots/*.png",                                emit: pngs,         optional: true
    path "plots/*.pdf",                                emit: pdfs,         optional: true
    path "cohort_nucleotide_diversity_summary.tsv",    emit: summary_tsv
    path "cohort_nucleotide_diversity_summary.xlsx",   emit: summary_xlsx, optional: true

    script:
    """
    set -euo pipefail

    Rscript ${custom_scripts}/plot_nucleotide_diversity.R \\
        --cohort_id "${cohort_id}" \\
        --species_map "${species_map_str}"
    """
}