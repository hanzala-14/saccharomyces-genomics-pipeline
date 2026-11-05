#!/usr/bin/env python3
import os
import re
import csv
import math
import random
import zipfile
import argparse
import subprocess
import xml.sax.saxutils as saxutils

# Standard Saccharomyces sensu stricto nuclear chromosome lengths for weighting I..XVI
CHR_WEIGHTS = {
    "I": 230218, "II": 813184, "III": 316620, "IV": 1531933,
    "V": 576874, "VI": 270161, "VII": 1090940, "VIII": 562643,
    "IX": 439888, "X": 745751, "XI": 666816, "XII": 1078177,
    "XIII": 924431, "XIV": 784333, "XV": 1091291, "XVI": 948066,
}


def write_xlsx_workbook(filename, sheets_dict):
    """Self-contained zero-dependency multi-sheet XLSX writer."""
    with zipfile.ZipFile(filename, "w", zipfile.ZIP_DEFLATED) as zf:
        sheet_names = list(sheets_dict.keys())
        overrides = "\n".join(
            f'  <Override PartName="/xl/worksheets/sheet{i+1}.xml" '
            f'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            for i in range(len(sheet_names))
        )
        zf.writestr(
            "[Content_Types].xml",
            f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">\n'
            f'  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>\n'
            f'  <Default Extension="xml" ContentType="application/xml"/>\n'
            f'  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>\n'
            f"{overrides}\n</Types>",
        )
        zf.writestr(
            "_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
            '  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>\n'
            "</Relationships>",
        )
        sheets_xml = "\n".join(
            f'    <sheet name="{saxutils.escape(name[:31])}" sheetId="{i+1}" r:id="rId{i+1}"/>'
            for i, name in enumerate(sheet_names)
        )
        zf.writestr(
            "xl/workbook.xml",
            f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            f'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">\n'
            f"  <sheets>\n{sheets_xml}\n  </sheets>\n</workbook>",
        )
        rels_xml = "\n".join(
            f'  <Relationship Id="rId{i+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i+1}.xml"/>'
            for i in range(len(sheet_names))
        )
        zf.writestr(
            "xl/_rels/workbook.xml.rels",
            f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
            f"{rels_xml}\n</Relationships>",
        )
        for idx, name in enumerate(sheet_names):
            rows = sheets_dict[name]
            row_xml = []
            for r_i, row in enumerate(rows, start=1):
                cells = []
                for val in row:
                    if val is None or val == "":
                        continue
                    if isinstance(val, (int, float)) and not isinstance(val, bool) and math.isfinite(val):
                        cells.append(f"<c><v>{val}</v></c>")
                    else:
                        esc = saxutils.escape(str(val))
                        cells.append(f'<c t="inlineStr"><is><t>{esc}</t></is></c>')
                row_xml.append(f'    <row r="{r_i}">{"".join(cells)}</row>')
            zf.writestr(
                f"xl/worksheets/sheet{idx+1}.xml",
                '<?xml version="1.0" encoding="UTF-8"?>\n'
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">\n'
                f'  <sheetData>\n{chr(10).join(row_xml)}\n  </sheetData>\n</worksheet>',
            )


def load_matrix_dict(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [l.rstrip("\r\n") for l in f if l.strip()]
    if len(lines) < 2:
        return {}, []
    delim = "\t" if "\t" in lines[0] else (";" if ";" in lines[0] else ",")
    raw_headers = [h.strip().strip('"') for h in lines[0].split(delim)]
    first_row = [p.strip().strip('"') for p in lines[1].split(delim)]

    if len(raw_headers) == len(first_row):
        col_names = raw_headers[1:]
    elif len(raw_headers) == len(first_row) - 1:
        col_names = raw_headers
    elif raw_headers[0] == "" or raw_headers[0].lower() in ("strain", "sample", "id", "row", "taxa"):
        col_names = raw_headers[1:]
    else:
        col_names = raw_headers

    mat = {}
    for line in lines[1:]:
        parts = [p.strip().strip('"') for p in line.split(delim)]
        if len(parts) < 2:
            continue
        row_id = parts[0]
        vals = parts[1:]
        for col_id, v_str in zip(col_names, vals):
            if row_id == col_id:
                continue
            try:
                v = float(v_str.replace(",", "."))
                if math.isfinite(v):
                    mat[(row_id, col_id)] = v
            except ValueError:
                pass
    return mat, col_names


def combine_chrom_matrices(chrom_file_list):
    """Combines 1..16 chromosomal matrices weighted by chromosome length and filters incomplete strains."""
    num_sum = {}
    den_sum = {}
    chr_count = {}
    all_strains = set()
    n_chrs = len(chrom_file_list)

    for roman, path in chrom_file_list:
        w = float(CHR_WEIGHTS.get(roman, 500000))
        mat, cols = load_matrix_dict(path)
        if not mat:
            continue
        all_strains.update(cols)
        for (s1, s2), v in mat.items():
            if s1 < s2:
                key = (s1, s2)
                num_sum[key] = num_sum.get(key, 0.0) + v * w
                den_sum[key] = den_sum.get(key, 0.0) + w
                chr_count[key] = chr_count.get(key, 0) + 1

    # Require pairs to have valid distance calls on at least 75% of chromosomes (removes non-hybrid/missing subgenomes)
    min_chrs = max(1, int(math.ceil(n_chrs * 0.75)))
    combined_mat = {}
    valid_strains = set()
    for (s1, s2), d_tot in den_sum.items():
        if chr_count.get((s1, s2), 0) >= min_chrs and d_tot > 0:
            val = num_sum[(s1, s2)] / d_tot
            combined_mat[(s1, s2)] = val
            combined_mat[(s2, s1)] = val
            valid_strains.add(s1)
            valid_strains.add(s2)

    return combined_mat, sorted(valid_strains)


def mantel_test_pvalue(pairs, permutations=1000):
    strains = sorted({p[0] for p in pairs}.union({p[1] for p in pairs}))
    if len(strains) < 4 or len(pairs) < 3:
        return 1.0
    d_x = {}
    d_y = {}
    for s1, s2, dx, dy in pairs:
        d_x[(s1, s2)] = d_x[(s2, s1)] = dx
        d_y[(s1, s2)] = d_y[(s2, s1)] = dy

    idx_pairs = []
    xs = []
    ys = []
    for i in range(len(strains)):
        for j in range(i + 1, len(strains)):
            a, b = strains[i], strains[j]
            if (a, b) in d_x and (a, b) in d_y:
                idx_pairs.append((i, j))
                xs.append(d_x[(a, b)])
                ys.append(d_y[(a, b)])

    n = len(xs)
    if n < 3:
        return 1.0
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return 1.0
    obs_cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))

    rng = random.Random(42)
    perm_strains = list(strains)
    extreme = 0
    for _ in range(permutations):
        rng.shuffle(perm_strains)
        p_cov = 0.0
        for (i, j), x_val in zip(idx_pairs, xs):
            pa, pb = perm_strains[i], perm_strains[j]
            y_val = d_y.get((pa, pb), my)
            p_cov += (x_val - mx) * (y_val - my)
        if abs(p_cov) >= abs(obs_cov):
            extreme += 1
    return (extreme + 1.0) / (permutations + 1.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-dir", default=".")
    parser.add_argument("--input-csv", default="")
    parser.add_argument("--species-map", default="c:cerevisiae,e:eubayanus,p:paradoxus,m:mikatae,k:kudriavzevii,u:uvarum,a:arboricola,j:jurei")
    parser.add_argument("--exclude-strains", default="")
    parser.add_argument("--out-prefix", default="Subgenome_Distance_Comparison")
    args = parser.parse_args()

    sp_map = {}
    for item in args.species_map.split(","):
        if ":" in item:
            k, v = item.split(":", 1)
            sp_map[k.strip()] = v.strip().capitalize()

    excl_set = {
        s.strip()
        for s in re.split(r"[,;\s\[\]'\"]+", args.exclude_strains or "")
        if s.strip()
    }

    all_comparisons = []

    # Mode A: Pre-existing pairwise CSV (e.g., Lager_Comparison.csv)
    if args.input_csv and os.path.isfile(args.input_csv):
        with open(args.input_csv, "r", encoding="utf-8", errors="replace") as f:
            lines = [l.strip() for l in f if l.strip()]
        if len(lines) >= 2:
            delim = ";" if ";" in lines[0] else ("\t" if "\t" in lines[0] else ",")
            headers = [h.strip() for h in lines[0].split(delim)]
            dist_cols = headers[2:] if len(headers) >= 4 else ["Cerevisiae", "Eubayanus"]
            col_names = [re.sub(r"(?i)_?distance$", "", c).capitalize() for c in dist_cols]
            for c1_i in range(len(col_names)):
                for c2_i in range(c1_i + 1, len(col_names)):
                    sp1, sp2 = col_names[c1_i], col_names[c2_i]
                    pairs = []
                    for line in lines[1:]:
                        parts = [p.strip() for p in line.split(delim)]
                        if len(parts) <= max(2 + c1_i, 2 + c2_i):
                            continue
                        s1, s2 = parts[0], parts[1]
                        if s1 in excl_set or s2 in excl_set:
                            continue
                        try:
                            d1 = float(parts[2 + c1_i].replace(",", "."))
                            d2 = float(parts[2 + c2_i].replace(",", "."))
                            if math.isfinite(d1) and math.isfinite(d2):
                                pairs.append((s1, s2, d1, d2))
                        except ValueError:
                            pass
                    if len(pairs) >= 3:
                        all_comparisons.append((sp1, sp2, pairs))

    # Mode B: Discover N subgenomes (either pre-combined or 16 chromosomal matrices per subgenome)
    if not all_comparisons:
        chr_romans = set(CHR_WEIGHTS.keys())

        def get_chrom_roman(rel_path):
            tokens = re.split(r"[/\\_.-]+", rel_path)
            for t in reversed(tokens):
                if t in chr_romans:
                    return t
            return None

        combined_by_sp = {}
        chroms_by_sp = {}

        for root, dirs, files in os.walk(args.matrix_dir):
            dirs[:] = [d for d in dirs if not d.startswith(".")]
            for f in files:
                if f.startswith(".") or f.endswith((".log", ".yml", ".json", ".sh", ".py", ".png", ".pdf", ".xlsx", ".R")):
                    continue
                if not f.endswith((".tsv", ".csv", ".txt", ".tab")) or f.startswith(args.out_prefix):
                    continue
                # Avoid raw allele tables or sample lists
                if f.endswith((".allele1.tsv", ".allele2.tsv", ".samples.txt")):
                    continue

                rel_path = os.path.join(root, f)
                sp_label = None
                for pref, full_sp in sp_map.items():
                    if re.search(rf"(^|[/\\_.-]){re.escape(full_sp.lower())}([/\\_.-]|$)", rel_path.lower()):
                        sp_label = full_sp
                        break
                if not sp_label:
                    for pref, full_sp in sp_map.items():
                        if re.search(rf"(^|[/\\_.-]){re.escape(pref.lower())}_([ivx]+|[0-9]+)([/\\_.-]|$)", rel_path.lower()):
                            sp_label = full_sp
                            break
                if not sp_label:
                    continue

                roman = get_chrom_roman(rel_path)
                is_all = any(w in f.lower() for w in ("ibsx_all", "subgenome_matrix", "combined", "final", "all", "genome"))
                if roman:
                    if is_all or f.endswith(".tsv"):
                        chroms_by_sp.setdefault(sp_label, {})[roman] = rel_path
                else:
                    score = 0 if is_all else 5
                    combined_by_sp.setdefault(sp_label, []).append((score, rel_path))

        loaded_species = {}
        all_sp = sorted(set(combined_by_sp.keys()) | set(chroms_by_sp.keys()))
        for sp_label in all_sp:
            if sp_label in combined_by_sp and combined_by_sp[sp_label]:
                cands = sorted(combined_by_sp[sp_label], key=lambda x: (x[0], x[1]))
                mat, cols = load_matrix_dict(cands[0][1])
                if mat and len(cols) >= 3:
                    print(f"[INFO] Loaded pre-combined {sp_label} matrix: {cands[0][1]} ({len(cols)} strains)")
                    loaded_species[sp_label] = (mat, cols)
                    continue
            if sp_label in chroms_by_sp and chroms_by_sp[sp_label]:
                chr_list = sorted(chroms_by_sp[sp_label].items())
                mat, cols = combine_chrom_matrices(chr_list)
                if mat and len(cols) >= 3:
                    print(f"[INFO] Combined {len(chr_list)} chromosomes for {sp_label} ({len(cols)} valid strains)")
                    loaded_species[sp_label] = (mat, cols)

        sp_keys = sorted(loaded_species.keys())
        for i in range(len(sp_keys)):
            for j in range(i + 1, len(sp_keys)):
                sp1, sp2 = sp_keys[i], sp_keys[j]
                mat1, cols1 = loaded_species[sp1]
                mat2, cols2 = loaded_species[sp2]
                common = sorted((set(cols1) & set(cols2)) - excl_set)
                pairs = []
                for a_i in range(len(common)):
                    for b_i in range(a_i + 1, len(common)):
                        s1, s2 = common[a_i], common[b_i]
                        d1 = mat1.get((s1, s2), mat1.get((s2, s1)))
                        d2 = mat2.get((s1, s2), mat2.get((s2, s1)))
                        if d1 is not None and d2 is not None and math.isfinite(d1) and math.isfinite(d2):
                            pairs.append((s1, s2, d1, d2))
                if len(pairs) >= 3:
                    all_comparisons.append((sp1, sp2, pairs))

    if not all_comparisons:
        print("WARN: Fewer than 2 subgenomes with >=3 shared strains found. Skipping comparison cleanly.")
        return

    reg_rows = [["Subgenome_X", "Subgenome_Y", "Shared_Strains", "Total_Pairs", "Slope_m", "Intercept_b", "R_Squared", "Pearson_r", "Mantel_p_value"]]
    xlsx_sheets = {"Regression_Summary": reg_rows}
    combined_csv_rows = []

    for sp1, sp2, pairs in all_comparisons:
        n = len(pairs)
        shared_strains = len({p[0] for p in pairs}.union({p[1] for p in pairs}))
        xs = [p[2] for p in pairs]
        ys = [p[3] for p in pairs]
        mx = sum(xs) / n
        my = sum(ys) / n
        sxx = sum((x - mx) ** 2 for x in xs)
        syy = sum((y - my) ** 2 for y in ys)
        sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        slope = sxy / sxx if sxx > 0 else 0.0
        intercept = my - slope * mx
        r_val = sxy / math.sqrt(sxx * syy) if (sxx > 0 and syy > 0) else 0.0
        r_sq = r_val ** 2
        mantel_p = mantel_test_pvalue(pairs, permutations=1000)

        reg_rows.append([
            sp1, sp2, shared_strains, n,
            round(slope, 6), round(intercept, 6),
            round(r_sq, 6), round(r_val, 6), round(mantel_p, 5)
        ])

        pair_tag = f"{sp1}_vs_{sp2}"
        pair_csv = f"{args.out_prefix}_{pair_tag}.csv"
        sheet_rows = [["Strain1", "Strain2", f"{sp1}_distance", f"{sp2}_distance"]]
        with open(pair_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(sheet_rows[0])
            for s1, s2, d1, d2 in pairs:
                row = [s1, s2, round(d1, 8), round(d2, 8)]
                w.writerow(row)
                sheet_rows.append(row)
                combined_csv_rows.append([pair_tag, sp1, sp2, s1, s2, round(d1, 8), round(d2, 8), round(slope, 6), round(intercept, 6), round(r_sq, 6), round(mantel_p, 5)])

        xlsx_sheets[pair_tag[:31]] = sheet_rows

        # Clean publication plot matching the reference Trello style
        r_code = f"""
        suppressPackageStartupMessages(library(tidyverse))
        df <- read_csv("{pair_csv}", show_col_types = FALSE)
        colnames(df) <- c("Strain1", "Strain2", "Dist_X", "Dist_Y")

        max_data_x <- max(df$Dist_X, na.rm = TRUE)
        max_data_y <- max(df$Dist_Y, na.rm = TRUE)
        lim_x <- ifelse(max_data_x <= 0.40, 0.40, max_data_x * 1.08)
        lim_y <- ifelse(max_data_y <= 0.14, 0.14, max_data_y * 1.12)

        x_min_line <- min(df$Dist_X, na.rm = TRUE)
        x_max_line <- max(df$Dist_X, na.rm = TRUE)
        y_at_max   <- {slope} * x_max_line + {intercept}

        sign_chr  <- ifelse({intercept} >= 0, "+", "-")
        label_txt <- sprintf("y = %.4fx %s %.4f\\nR² = %.4f", {slope}, sign_chr, abs({intercept}), {r_sq})

        lbl_x <- min(lim_x * 0.82, max(lim_x * 0.18, x_max_line))
        lbl_y <- min(lim_y * 0.90, max(lim_y * 0.15, y_at_max + lim_y * 0.14))

        p <- ggplot(df, aes(x = Dist_X, y = Dist_Y)) +
          geom_point(color = "#4A7BB7", alpha = 0.80, size = 1.5) +
          geom_segment(
            x = x_min_line, y = max(0, {slope} * x_min_line + {intercept}),
            xend = x_max_line, yend = max(0, y_at_max),
            linetype = "dotted", color = "#4A7BB7", linewidth = 0.8
          ) +
          annotate(
            "text", x = lbl_x, y = lbl_y,
            label = label_txt, size = 3.5, fontface = "bold", color = "#404040"
          ) +
          scale_x_continuous(limits = c(0, lim_x), expand = c(0, 0), breaks = scales::pretty_breaks(n = 8)) +
          scale_y_continuous(limits = c(0, lim_y), expand = c(0, 0), breaks = scales::pretty_breaks(n = 7)) +
          labs(
            title = "Pairwise Genetic Distance of {sp1} vs. {sp2} Subgenomes",
            x = "{sp1} Pairwise Distance",
            y = "{sp2} Pairwise Distance"
          ) +
          theme_minimal(base_size = 11) +
          theme(
            plot.title = element_text(size = 13, color = "#404040", hjust = 0.5, margin = margin(b = 10)),
            axis.title = element_text(size = 10, color = "#595959"),
            axis.text  = element_text(size = 9, color = "#595959"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.4),
            panel.border = element_rect(color = "#D9D9D9", fill = NA, linewidth = 0.6),
            plot.margin = margin(14, 18, 14, 14)
          )

        ggsave("{args.out_prefix}_{pair_tag}.png", plot = p, width = 9.5, height = 5.4, dpi = 300, bg = "white")
        ggsave("{args.out_prefix}_{pair_tag}.pdf", plot = p, width = 9.5, height = 5.4, bg = "white")
        """
        subprocess.run(["Rscript", "-e", r_code], check=True)

    master_csv = f"{args.out_prefix}_All_Pairs.csv"
    with open(master_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Comparison", "Subgenome_X", "Subgenome_Y", "Strain1", "Strain2", "Dist_X", "Dist_Y", "Slope", "Intercept", "R_Squared", "Mantel_p"])
        w.writerows(combined_csv_rows)

    write_xlsx_workbook(f"{args.out_prefix}.xlsx", xlsx_sheets)


if __name__ == "__main__":
    main()