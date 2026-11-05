#!/usr/bin/env python3
import os
import re
import csv
import math
import zipfile
import argparse
import statistics
import xml.etree.ElementTree as ET
import xml.sax.saxutils as saxutils

ROMAN_ORDER = [
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII",
    "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"
]


def col_letter(n):
    string = ""
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        string = chr(65 + remainder) + string
    return string


def col_index_from_ref(cell_ref):
    m = re.match(r"^([A-Z]+)", cell_ref)
    if not m:
        return 1
    letters = m.group(1)
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch) - ord("A") + 1)
    return idx


def parse_species_map(map_str):
    smap = {}
    for item in map_str.split(","):
        if ":" in item:
            k, v = item.strip().split(":", 1)
            smap[k.strip()] = v.strip()
    return smap


def chrom_sort_key(chrom, species_keys):
    if "_" in chrom:
        prefix, suffix = chrom.split("_", 1)
    else:
        prefix, suffix = "", chrom
    p_idx = species_keys.index(prefix) if prefix in species_keys else 99
    if suffix in ROMAN_ORDER:
        s_idx = ROMAN_ORDER.index(suffix)
    elif suffix == "M":
        s_idx = 100
    elif suffix == "P":
        s_idx = 101
    else:
        s_idx = 50
    return (p_idx, s_idx, chrom)


def quantile_p(vals, q):
    if not vals:
        return 0.0
    s = sorted(vals)
    pos = (len(s) - 1) * q
    low = int(math.floor(pos))
    high = int(math.ceil(pos))
    if low == high:
        return float(s[low])
    return float(s[low] + (s[high] - s[low]) * (pos - low))


def compute_dynamic_baseline(chrom_medians):
    nuc_vals = [
        v for c, v in chrom_medians.items()
        if not c.endswith("_M") and not c.endswith("_P") and v > 0
    ]
    if not nuc_vals:
        return 1.0, 1.0, []

    q75 = quantile_p(nuc_vals, 0.75)
    noise_cutoff = max(5.0, 0.20 * q75)

    sc_active = [
        v for c, v in chrom_medians.items()
        if c.startswith("c_") and not c.endswith(("_M", "_P")) and v >= noise_cutoff
    ]
    kept_chroms = [
        c for c, v in chrom_medians.items()
        if not c.endswith(("_M", "_P")) and v >= noise_cutoff
    ]
    all_active = [
        v for c, v in chrom_medians.items()
        if not c.endswith(("_M", "_P")) and v >= noise_cutoff
    ]

    g22_valter = float(statistics.median(sc_active)) if sc_active else (
        float(statistics.median(all_active)) if all_active else 1.0
    )
    g22_balint = float(statistics.median(all_active)) if all_active else 1.0
    return g22_valter, g22_balint, kept_chroms


def rolling_median_5(vals):
    n = len(vals)
    out = [0.0] * n
    for i in range(n):
        w = sorted(vals[max(0, i - 2):min(n, i + 3)])
        m = len(w)
        out[i] = w[m // 2] if m % 2 == 1 else 0.5 * (w[m // 2 - 1] + w[m // 2])
    return out


def detect_breakpoints(chrom_windows, baseline):
    breakpoints = []
    if baseline <= 0:
        return breakpoints
    inv_base = 1.0 / baseline
    for chrom, rows in chrom_windows.items():
        if chrom.endswith(("_M", "_P")) or len(rows) < 10:
            continue
        covs = [r[2] * inv_base for r in rows]
        norm_roll = rolling_median_5(covs)
        n = len(norm_roll)
        for i in range(5, n - 5):
            left_med = sorted(norm_roll[i - 5:i])[2]
            right_med = sorted(norm_roll[i:i + 5])[2]
            diff = abs(right_med - left_med)
            if diff >= 0.75 and (left_med >= 0.4 or right_med >= 0.4):
                breakpoints.append({
                    "chrom": chrom,
                    "Window_Start": rows[i][0],
                    "Window_End": rows[i][1],
                    "Left_Norm_Depth": round(left_med, 3),
                    "Right_Norm_Depth": round(right_med, 3),
                    "Delta": round(right_med - left_med, 3)
                })
    dedup = []
    for bp in breakpoints:
        if not dedup or dedup[-1]["chrom"] != bp["chrom"] or (bp["Window_Start"] - dedup[-1]["Window_Start"]) > 25000:
            dedup.append(bp)
    return dedup


def build_subgenome_chart_xml(title_text, sheet_name, col_let, start_row, end_row, cached_vals, bar_hex_color, has_user_shapes=True):
    pts_xml = []
    for idx, v in enumerate(cached_vals):
        v_num = v if (isinstance(v, (int, float)) and math.isfinite(v)) else 0.0
        pts_xml.append(f'<c:pt idx="{idx}"><c:v>{v_num}</c:v></c:pt>')
    pts_str = "".join(pts_xml)
    pt_count = len(cached_vals)
    user_shapes_tag = '<c:userShapes r:id="rId1"/>' if has_user_shapes else ""

    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<c:chartSpace xmlns:c="http://schemas.openxmlformats.org/drawingml/2006/chart" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <c:date1904 val="0"/>
  <c:lang val="en-US"/>
  <c:roundedCorners val="0"/>
  <c:chart>
    <c:title>
      <c:tx>
        <c:rich>
          <a:bodyPr/>
          <a:lstStyle/>
          <a:p>
            <a:pPr><a:defRPr sz="1150" b="1" i="1"/></a:pPr>
            <a:r><a:rPr lang="en-US" sz="1150" b="1" i="1"/><a:t>{saxutils.escape(title_text)}</a:t></a:r>
          </a:p>
        </c:rich>
      </c:tx>
      <c:layout><c:manualLayout><c:xMode val="edge"/><c:yMode val="edge"/><c:x val="0.02"/><c:y val="0.01"/></c:manualLayout></c:layout>
      <c:overlay val="0"/>
    </c:title>
    <c:autoTitleDeleted val="0"/>
    <c:plotArea>
      <c:layout>
        <c:manualLayout>
          <c:layoutTarget val="inner"/>
          <c:xMode val="edge"/>
          <c:yMode val="edge"/>
          <c:x val="0.025"/>
          <c:y val="0.18"/>
          <c:w val="0.965"/>
          <c:h val="0.74"/>
        </c:manualLayout>
      </c:layout>
      <c:barChart>
        <c:barDir val="col"/>
        <c:grouping val="clustered"/>
        <c:varyColors val="0"/>
        <c:ser>
          <c:idx val="0"/>
          <c:order val="0"/>
          <c:spPr>
            <a:solidFill><a:srgbClr val="{bar_hex_color}"/></a:solidFill>
            <a:ln><a:noFill/></a:ln>
          </c:spPr>
          <c:invertIfNegative val="0"/>
          <c:val>
            <c:numRef>
              <c:f>{sheet_name}!${col_let}${start_row}:${col_let}${end_row}</c:f>
              <c:numCache>
                <c:formatCode>0.0</c:formatCode>
                <c:ptCount val="{pt_count}"/>
                {pts_str}
              </c:numCache>
            </c:numRef>
          </c:val>
        </c:ser>
        <c:gapWidth val="0"/>
        <c:overlap val="-27"/>
        <c:axId val="497148216"/>
        <c:axId val="497148872"/>
      </c:barChart>
      <c:catAx>
        <c:axId val="497148216"/>
        <c:scaling><c:orientation val="minMax"/></c:scaling>
        <c:delete val="1"/>
        <c:axPos val="b"/>
        <c:majorTickMark val="none"/>
        <c:minorTickMark val="none"/>
        <c:tickLblPos val="nextTo"/>
        <c:crossAx val="497148872"/>
        <c:crosses val="autoZero"/>
        <c:auto val="1"/>
        <c:lblAlgn val="ctr"/>
        <c:lblOffset val="100"/>
      </c:catAx>
      <c:valAx>
        <c:axId val="497148872"/>
        <c:scaling>
          <c:orientation val="minMax"/>
          <c:max val="8"/>
          <c:min val="0"/>
        </c:scaling>
        <c:delete val="0"/>
        <c:axPos val="l"/>
        <c:majorGridlines>
          <c:spPr>
            <a:ln w="9525" cap="flat" cmpd="sng" algn="ctr">
              <a:solidFill><a:srgbClr val="D9D9D9"/></a:solidFill>
            </a:ln>
          </c:spPr>
        </c:majorGridlines>
        <c:numFmt formatCode="0" sourceLinked="0"/>
        <c:majorTickMark val="none"/>
        <c:minorTickMark val="none"/>
        <c:tickLblPos val="nextTo"/>
        <c:crossAx val="497148216"/>
        <c:crosses val="autoZero"/>
        <c:crossBetween val="between"/>
        <c:majorUnit val="2"/>
      </c:valAx>
      <c:spPr><a:noFill/><a:ln><a:noFill/></a:ln></c:spPr>
    </c:plotArea>
    <c:plotVisOnly val="0"/>
    <c:dispBlanksAs val="gap"/>
  </c:chart>
  <c:spPr>
    <a:solidFill><a:srgbClr val="FFFFFF"/></a:solidFill>
    <a:ln w="9525"><a:solidFill><a:srgbClr val="BFBFBF"/></a:solidFill></a:ln>
  </c:spPr>
  {user_shapes_tag}
</c:chartSpace>'''


def build_chr_divider_shapes_xml(chr_fractions):
    shapes = ['<c:userShapes xmlns:c="http://schemas.openxmlformats.org/drawingml/2006/chart" xmlns:cdr="http://schemas.openxmlformats.org/drawingml/2006/chartDrawing" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">']
    plot_x0 = 0.025
    plot_w = 0.965
    prev_x = plot_x0

    for idx, (label, end_frac) in enumerate(chr_fractions, start=1):
        x_pos = round(plot_x0 + plot_w * end_frac, 5)
        mid_x = round((prev_x + x_pos) / 2.0, 5)
        prev_x = x_pos

        shapes.append(f'''  <cdr:relSizeAnchor>
    <cdr:from><cdr:x>{x_pos}</cdr:x><cdr:y>0.18000</cdr:y></cdr:from>
    <cdr:to><cdr:x>{x_pos}</cdr:x><cdr:y>0.92000</cdr:y></cdr:to>
    <cdr:cxnSp macro="">
      <cdr:nvCxnSpPr><cdr:cNvPr id="{idx + 10}" name="ChrLine_{idx}"/><cdr:cNvCxnSpPr/></cdr:nvCxnSpPr>
      <cdr:spPr>
        <a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="444500"/></a:xfrm>
        <a:prstGeom prst="line"><a:avLst/></a:prstGeom>
        <a:ln w="9525"><a:solidFill><a:srgbClr val="404040"/></a:solidFill><a:prstDash val="dash"/></a:ln>
      </cdr:spPr>
    </cdr:cxnSp>
  </cdr:relSizeAnchor>''')

        lbl_x1 = max(0.0, round(mid_x - 0.015, 5))
        lbl_x2 = min(1.0, round(mid_x + 0.015, 5))
        shapes.append(f'''  <cdr:relSizeAnchor>
    <cdr:from><cdr:x>{lbl_x1}</cdr:x><cdr:y>0.91500</cdr:y></cdr:from>
    <cdr:to><cdr:x>{lbl_x2}</cdr:x><cdr:y>0.99500</cdr:y></cdr:to>
    <cdr:sp macro="" textlink="">
      <cdr:nvSpPr><cdr:cNvPr id="{idx + 100}" name="ChrLbl_{idx}"/><cdr:cNvSpPr txBox="1"/></cdr:nvSpPr>
      <cdr:spPr>
        <a:xfrm><a:off x="0" y="0"/><a:ext cx="200000" cy="100000"/></a:xfrm>
        <a:prstGeom prst="rect"><a:avLst/></a:prstGeom>
        <a:noFill/>
        <a:ln><a:noFill/></a:ln>
      </cdr:spPr>
      <cdr:txBody>
        <a:bodyPr wrap="none" anchor="ctr"/>
        <a:lstStyle/>
        <a:p>
          <a:pPr algn="ctr"/>
          <a:r><a:rPr lang="en-US" sz="700" b="1"><a:solidFill><a:srgbClr val="333333"/></a:solidFill></a:rPr><a:t>{saxutils.escape(label)}</a:t></a:r>
        </a:p>
      </cdr:txBody>
    </cdr:sp>
  </cdr:relSizeAnchor>''')

    shapes.append("</c:userShapes>")
    return "\n".join(shapes)


def get_cell_style(sheet_name, r_idx, c_idx, val, tier_layout=None):
    """Assigns color-coded section styles, borders, and decimal formats across 1..8 subgenomes."""
    if sheet_name == "Interactive_Plotter":
        if r_idx == 1:
            if c_idx in (1, 2, 3, 4, 5, 13):
                return 1
            elif c_idx in (6, 7, 8):
                return 2
            elif c_idx in (9, 10, 11):
                return 3
            return 0
        if r_idx == 2 and 6 <= c_idx <= 11:
            return 21

        if tier_layout:
            for (hdr_r, chr_r1, chr_r2, org_m1, org_m2, org_c1, org_c2, has_right) in tier_layout:
                if r_idx == hdr_r:
                    if c_idx in (6, 7, 8):
                        return 2
                    elif c_idx in (9, 10, 11) and has_right:
                        return 3
                if chr_r1 <= r_idx <= chr_r2:
                    if c_idx == 6:
                        return 4
                    elif c_idx == 7:
                        return 5
                    elif c_idx == 8:
                        return 6
                    elif has_right:
                        if c_idx == 9:
                            return 7
                        elif c_idx == 10:
                            return 8
                        elif c_idx == 11:
                            return 9
                if org_m1 <= r_idx <= org_m2:
                    if c_idx in (6, 8) or (has_right and c_idx in (9, 11)):
                        return 15
                    elif c_idx == 7 or (has_right and c_idx == 10):
                        return 16
                if org_c1 <= r_idx <= org_c2:
                    if c_idx in (6, 8) or (has_right and c_idx in (9, 11)):
                        return 17
                    elif c_idx == 7 or (has_right and c_idx == 10):
                        return 18

        if r_idx in (20, 22, 24):
            if c_idx == 6:
                return 10
            elif c_idx == 7:
                if r_idx == 20:
                    return 11
                elif r_idx == 22:
                    return 12
                elif r_idx == 24:
                    return 13
            elif c_idx in (8, 9):
                return 14
        if c_idx == 4 and r_idx >= 2:
            return 19
        if c_idx == 5 and 2 <= r_idx <= 35:
            return 21
        if c_idx == 13 and r_idx >= 2:
            return 20
        return 0

    elif sheet_name == "Cohort_Chromosome_Summary":
        if r_idx == 1:
            return 1
        if r_idx == 2:
            return 10 if c_idx == 1 else 12
        if r_idx == 3:
            return 10 if c_idx == 1 else 13
        if c_idx == 1 and r_idx >= 4:
            s_val = str(val)
            if s_val.endswith(("_M", "_P")):
                return 22
            elif s_val.startswith(("c_", "p_", "k_", "a_")):
                return 4
            elif s_val.startswith(("e_", "m_", "u_", "j_")):
                return 7
            return 10
        if c_idx >= 2 and r_idx >= 4:
            return 20
        return 0

    elif sheet_name == "Raw_Chromosome_Medians":
        if r_idx == 1:
            return 1
        if c_idx == 1 and r_idx >= 2:
            s_val = str(val)
            if s_val.endswith(("_M", "_P")):
                return 22
            elif s_val.startswith(("c_", "p_", "k_", "a_")):
                return 4
            elif s_val.startswith(("e_", "m_", "u_", "j_")):
                return 7
            return 10
        if c_idx >= 2 and r_idx >= 2:
            return 19
        return 0

    else:
        if r_idx == 1:
            return 1
        return 0


def build_cols_xml(sheet_name, max_cols):
    if sheet_name == "Interactive_Plotter":
        widths = {
            1: 10,   # A: chrom
            2: 11,   # B: start
            3: 11,   # C: end
            4: 21,   # D: Auto_Lookup_Raw_Depth
            5: 21,   # E: Manual_Paste_Override
            6: 31,   # F: Left Subgenome Panel
            7: 16,   # G: Left Median Depth / Control Value
            8: 16,   # H: Left Copy Number
            9: 29,   # I: Right Subgenome Panel
            10: 15,  # J: Right Median Depth
            11: 15,  # K: Right Copy Number
            12: 4,   # L: Spacer
            13: 23   # M: Normalized_Copy_Number
        }
        cols = ["<cols>"]
        for c_i, w in sorted(widths.items()):
            cols.append(f'<col min="{c_i}" max="{c_i}" width="{w}" customWidth="1"/>')
        cols.append("</cols>")
        return "".join(cols)
    elif sheet_name in ("Cohort_Chromosome_Summary", "Raw_Chromosome_Medians"):
        end_c = max(2, max_cols)
        return (
            f'<cols>'
            f'<col min="1" max="1" width="18" customWidth="1"/>'
            f'<col min="2" max="{end_c}" width="15" customWidth="1"/>'
            f'</cols>'
        )
    elif sheet_name in ("Baseline_QC_Metadata", "Chimeric_Breakpoints"):
        end_c = max(2, max_cols)
        return f'<cols><col min="1" max="{end_c}" width="22" customWidth="1"/></cols>'
    else:
        end_c = max(4, max_cols)
        return (
            f'<cols>'
            f'<col min="1" max="3" width="12" customWidth="1"/>'
            f'<col min="4" max="{end_c}" width="14" customWidth="1"/>'
            f'</cols>'
        )


def write_multi_sheet_xlsx(filename, sheets_dict, color_scale_sheets=None, chart_specs=None, tier_layout=None):
    if color_scale_sheets is None:
        color_scale_sheets = {}
    if chart_specs is None:
        chart_specs = []

    sheet_names = list(sheets_dict.keys())
    has_charts = len(chart_specs) > 0

    with zipfile.ZipFile(filename, "w", zipfile.ZIP_DEFLATED) as zf:
        ct = [
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
            '<Default Extension="xml" ContentType="application/xml"/>',
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        ]
        for idx in range(1, len(sheet_names) + 1):
            ct.append(f'<Override PartName="/xl/worksheets/sheet{idx}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')

        if has_charts:
            ct.append('<Override PartName="/xl/drawings/drawing1.xml" ContentType="application/vnd.openxmlformats-officedocument.drawing+xml"/>')
            ct.append('<Override PartName="/xl/worksheets/sheet_diagram.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
            ct.append('<Override PartName="/xl/drawings/drawing_diag.xml" ContentType="application/vnd.openxmlformats-officedocument.drawing+xml"/>')
            for c_i in range(1, len(chart_specs) * 2 + 1):
                ct.append(f'<Override PartName="/xl/charts/chart{c_i}.xml" ContentType="application/vnd.openxmlformats-officedocument.drawingml.chart+xml"/>')
                ct.append(f'<Override PartName="/xl/drawings/chartshapes{c_i}.xml" ContentType="application/vnd.openxmlformats-officedocument.drawingml.chartshapes+xml"/>')

        ct.append("</Types>")
        zf.writestr("[Content_Types].xml", "".join(ct))

        zf.writestr("_rels/.rels", '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>''')

        wb_rels = [
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        ]
        for idx in range(1, len(sheet_names) + 1):
            wb_rels.append(f'<Relationship Id="rId{idx}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{idx}.xml"/>')

        next_rid = len(sheet_names) + 1
        diag_rid = None
        if has_charts:
            diag_rid = f"rId{next_rid}"
            wb_rels.append(f'<Relationship Id="{diag_rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet_diagram.xml"/>')
            next_rid += 1

        wb_rels.append(f'<Relationship Id="rId{next_rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>')
        wb_rels.append("</Relationships>")
        zf.writestr("xl/_rels/workbook.xml.rels", "".join(wb_rels))

        wb_xml = [
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
            '<calcPr calcId="171027" fullCalcOnLoad="1"/>',
            '<sheets>'
        ]
        sheet_id_counter = 1
        for idx, name in enumerate(sheet_names, 1):
            wb_xml.append(f'<sheet name="{saxutils.escape(name[:31])}" sheetId="{sheet_id_counter}" r:id="rId{idx}"/>')
            sheet_id_counter += 1
            if idx == 1 and has_charts and diag_rid:
                wb_xml.append(f'<sheet name="Diagram2" sheetId="{sheet_id_counter}" r:id="{diag_rid}"/>')
                sheet_id_counter += 1
        wb_xml.append("</sheets></workbook>")
        zf.writestr("xl/workbook.xml", "".join(wb_xml))

        zf.writestr("xl/styles.xml", '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <numFmts count="2">
    <numFmt numFmtId="164" formatCode="0.0"/>
    <numFmt numFmtId="165" formatCode="0.00"/>
  </numFmts>
  <fonts count="7">
    <font><sz val="11"/><color rgb="FF1F2937"/><name val="Calibri"/></font>
    <font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font>
    <font><b/><sz val="11"/><color rgb="FF1F2937"/><name val="Calibri"/></font>
    <font><b/><sz val="11"/><color rgb="FFC00000"/><name val="Calibri"/></font>
    <font><i/><sz val="10"/><color rgb="FF4B5563"/><name val="Calibri"/></font>
    <font><b/><sz val="12"/><color rgb="FFC00000"/><name val="Calibri"/></font>
    <font><sz val="10"/><color rgb="FF1F2937"/><name val="Calibri"/></font>
  </fonts>
  <fills count="14">
    <fill><patternFill patternType="none"/></fill>
    <fill><patternFill patternType="gray125"/></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FF1F4E78"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FF2F5597"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FF375623"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFD9E1F2"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFE2EFDA"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFFFF2CC"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFFFD966"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFFFC000"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFC6E0B4"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFFCE4D6"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFEDE7F6"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFF3F4F6"/><bgColor indexed="64"/></patternFill></fill>
  </fills>
  <borders count="3">
    <border><left/><right/><top/><bottom/><diagonal/></border>
    <border>
      <left style="thin"><color rgb="FFBFBFBF"/></left>
      <right style="thin"><color rgb="FFBFBFBF"/></right>
      <top style="thin"><color rgb="FFBFBFBF"/></top>
      <bottom style="thin"><color rgb="FFBFBFBF"/></bottom>
      <diagonal/>
    </border>
    <border>
      <left style="medium"><color rgb="FF1F2937"/></left>
      <right style="medium"><color rgb="FF1F2937"/></right>
      <top style="medium"><color rgb="FF1F2937"/></top>
      <bottom style="medium"><color rgb="FF1F2937"/></bottom>
      <diagonal/>
    </border>
  </borders>
  <cellStyleXfs count="1">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>
  </cellStyleXfs>
  <cellXfs count="23">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
    <xf numFmtId="0" fontId="1" fillId="2" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="1" fillId="3" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="1" fillId="4" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="2" fillId="5" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="164" fontId="0" fillId="5" borderId="1" xfId="0" applyNumberFormat="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="165" fontId="3" fillId="5" borderId="1" xfId="0" applyNumberFormat="1" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="2" fillId="6" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="164" fontId="0" fillId="6" borderId="1" xfId="0" applyNumberFormat="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="165" fontId="3" fillId="6" borderId="1" xfId="0" applyNumberFormat="1" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="2" fillId="8" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="5" fillId="7" borderId="2" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="164" fontId="5" fillId="7" borderId="2" xfId="0" applyNumberFormat="1" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="164" fontId="5" fillId="9" borderId="2" xfId="0" applyNumberFormat="1" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="4" fillId="7" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="2" fillId="10" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="164" fontId="0" fillId="10" borderId="1" xfId="0" applyNumberFormat="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="2" fillId="7" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="165" fontId="3" fillId="7" borderId="1" xfId="0" applyNumberFormat="1" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>
    <xf numFmtId="165" fontId="0" fillId="13" borderId="0" xfId="0" applyNumberFormat="1" applyFill="1"/>
    <xf numFmtId="0" fontId="2" fillId="11" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
    <xf numFmtId="0" fontId="2" fillId="12" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1"/>
  </cellXfs>
  <cellStyles count="1">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>
  </cellStyles>
</styleSheet>''')

        for idx, name in enumerate(sheet_names, 1):
            rows = sheets_dict[name]
            max_cols = max((len(r) for r in rows), default=1)
            cols_xml = build_cols_xml(name, max_cols)

            ws_xml = [
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
                cols_xml,
                '<sheetData>'
            ]
            for r_idx, row in enumerate(rows, 1):
                ws_xml.append(f'<row r="{r_idx}">')
                for c_idx, val in enumerate(row, 1):
                    style_id = get_cell_style(name, r_idx, c_idx, val, tier_layout=tier_layout)
                    style_attr = f' s="{style_id}"' if style_id > 0 else ""
                    cell_ref = f"{col_letter(c_idx)}{r_idx}"

                    if val is None or val == "":
                        if style_id > 0:
                            ws_xml.append(f'<c r="{cell_ref}"{style_attr}/>')
                        continue

                    if isinstance(val, tuple) and len(val) >= 3 and val[0] == "__FORMULA__":
                        f_str = saxutils.escape(str(val[1]))
                        cached = val[2]
                        if isinstance(cached, (int, float)) and not isinstance(cached, bool) and math.isfinite(cached):
                            ws_xml.append(f'<c r="{cell_ref}"{style_attr}><f>{f_str}</f><v>{cached}</v></c>')
                        else:
                            ws_xml.append(f'<c r="{cell_ref}"{style_attr}><f>{f_str}</f></c>')
                    elif isinstance(val, (int, float)) and not isinstance(val, bool) and math.isfinite(val):
                        ws_xml.append(f'<c r="{cell_ref}"{style_attr}><v>{val}</v></c>')
                    else:
                        ws_xml.append(f'<c r="{cell_ref}" t="inlineStr"{style_attr}><is><t>{saxutils.escape(str(val))}</t></is></c>')
                ws_xml.append("</row>")
            ws_xml.append("</sheetData>")

            if name in color_scale_sheets and len(rows) >= 2 and max_cols >= 2:
                start_r, min_v, mid_v, max_v = color_scale_sheets[name]
                end_col = col_letter(max_cols)
                end_r = len(rows)
                if end_r >= start_r:
                    ws_xml.append(
                        f'<conditionalFormatting sqref="B{start_r}:{end_col}{end_r}">'
                        f'<cfRule type="colorScale" priority="1"><colorScale>'
                        f'<cfvo type="num" val="{min_v}"/><cfvo type="num" val="{mid_v}"/><cfvo type="num" val="{max_v}"/>'
                        f'<color rgb="FFF8696B"/><color rgb="FFFFEB84"/><color rgb="FF63BE7B"/>'
                        f'</colorScale></cfRule></conditionalFormatting>'
                    )

            if idx == 1 and has_charts:
                ws_xml.append('<drawing r:id="rId1"/>')
                zf.writestr(
                    "xl/worksheets/_rels/sheet1.xml.rels",
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/drawing" Target="../drawings/drawing1.xml"/>'
                    '</Relationships>'
                )

            ws_xml.append("</worksheet>")
            zf.writestr(f"xl/worksheets/sheet{idx}.xml", "".join(ws_xml))

        if has_charts:
            zf.writestr(
                "xl/worksheets/sheet_diagram.xml",
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
                '<sheetViews><sheetView showGridLines="0" workbookViewId="0"/></sheetViews>'
                '<sheetData><row r="1"><c r="A1" t="inlineStr" s="1"><is><t>Live Subgenome Coverage Diagram (Linked to Interactive_Plotter G20 / G22 / G24)</t></is></c></row></sheetData>'
                '<drawing r:id="rId1"/>'
                '</worksheet>'
            )
            zf.writestr(
                "xl/worksheets/_rels/sheet_diagram.xml.rels",
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/drawing" Target="../drawings/drawing_diag.xml"/>'
                '</Relationships>'
            )

            d1_anchors = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?><xdr:wsDr xmlns:xdr="http://schemas.openxmlformats.org/drawingml/2006/spreadsheetDrawing" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">']
            d1_rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']

            dd_anchors = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?><xdr:wsDr xmlns:xdr="http://schemas.openxmlformats.org/drawingml/2006/spreadsheetDrawing" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">']
            dd_rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']

            for s_idx, spec in enumerate(chart_specs):
                c_id_1 = s_idx * 2 + 1
                c_id_2 = s_idx * 2 + 2

                for cid in (c_id_1, c_id_2):
                    c_xml = build_subgenome_chart_xml(
                        spec["title"], "Interactive_Plotter", spec["col_let"],
                        spec["start_row"], spec["end_row"], spec["cached_vals"],
                        spec["color"], has_user_shapes=True
                    )
                    zf.writestr(f"xl/charts/chart{cid}.xml", c_xml)
                    zf.writestr(f"xl/drawings/chartshapes{cid}.xml", build_chr_divider_shapes_xml(spec["chr_fractions"]))
                    zf.writestr(
                        f"xl/charts/_rels/chart{cid}.xml.rels",
                        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                        f'<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/chartUserShapes" Target="../drawings/chartshapes{cid}.xml"/>'
                        '</Relationships>'
                    )

                r_top = s_idx * 16 + 1
                r_bot = r_top + 15
                d1_anchors.append(f'''<xdr:twoCellAnchor>
  <xdr:from><xdr:col>14</xdr:col><xdr:colOff>0</xdr:colOff><xdr:row>{r_top}</xdr:row><xdr:rowOff>0</xdr:rowOff></xdr:from>
  <xdr:to><xdr:col>30</xdr:col><xdr:colOff>0</xdr:colOff><xdr:row>{r_bot}</xdr:row><xdr:rowOff>0</xdr:rowOff></xdr:to>
  <xdr:graphicFrame macro="">
    <xdr:nvGraphicFramePr><xdr:cNvPr id="{c_id_1 + 1}" name="SubgenomeChart_{c_id_1}"/><xdr:cNvGraphicFramePr/></xdr:nvGraphicFramePr>
    <xdr:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/></xdr:xfrm>
    <a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/chart"><c:chart xmlns:c="http://schemas.openxmlformats.org/drawingml/2006/chart" r:id="rId{s_idx + 1}"/></a:graphicData></a:graphic>
  </xdr:graphicFrame>
  <xdr:clientData/>
</xdr:twoCellAnchor>''')
                d1_rels.append(f'<Relationship Id="rId{s_idx + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/chart" Target="../charts/chart{c_id_1}.xml"/>')

                dr_top = s_idx * 18 + 2
                dr_bot = dr_top + 17
                dd_anchors.append(f'''<xdr:twoCellAnchor>
  <xdr:from><xdr:col>0</xdr:col><xdr:colOff>0</xdr:colOff><xdr:row>{dr_top}</xdr:row><xdr:rowOff>0</xdr:rowOff></xdr:from>
  <xdr:to><xdr:col>24</xdr:col><xdr:colOff>0</xdr:colOff><xdr:row>{dr_bot}</xdr:row><xdr:rowOff>0</xdr:rowOff></xdr:to>
  <xdr:graphicFrame macro="">
    <xdr:nvGraphicFramePr><xdr:cNvPr id="{c_id_2 + 1}" name="DiagramChart_{c_id_2}"/><xdr:cNvGraphicFramePr/></xdr:nvGraphicFramePr>
    <xdr:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/></xdr:xfrm>
    <a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/chart"><c:chart xmlns:c="http://schemas.openxmlformats.org/drawingml/2006/chart" r:id="rId{s_idx + 1}"/></a:graphicData></a:graphic>
  </xdr:graphicFrame>
  <xdr:clientData/>
</xdr:twoCellAnchor>''')
                dd_rels.append(f'<Relationship Id="rId{s_idx + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/chart" Target="../charts/chart{c_id_2}.xml"/>')

            d1_anchors.append("</xdr:wsDr>")
            d1_rels.append("</Relationships>")
            zf.writestr("xl/drawings/drawing1.xml", "".join(d1_anchors))
            zf.writestr("xl/drawings/_rels/drawing1.xml.rels", "".join(d1_rels))

            dd_anchors.append("</xdr:wsDr>")
            dd_rels.append("</Relationships>")
            zf.writestr("xl/drawings/drawing_diag.xml", "".join(dd_anchors))
            zf.writestr("xl/drawings/_rels/drawing_diag.xml.rels", "".join(dd_rels))


def extract_curated_xlsx_to_tsv(xlsx_path, out_tsv):
    ns = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(xlsx_path, "r") as zf:
        shared_strings = []
        if "xl/sharedStrings.xml" in zf.namelist():
            ss_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in ss_root.findall("main:si", ns):
                texts = [t.text or "" for t in si.findall(".//main:t", ns)]
                shared_strings.append("".join(texts))

        wb_root = ET.fromstring(zf.read("xl/workbook.xml"))
        rels_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_map = {}
        for rel in rels_root:
            rel_map[rel.attrib.get("Id")] = rel.attrib.get("Target")

        target_sheet_path = "xl/worksheets/sheet2.xml"
        for sheet in wb_root.findall(".//main:sheet", ns):
            s_name = sheet.attrib.get("name", "")
            r_id = sheet.attrib.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
            if s_name in ("Cohort_Chromosome_Summary", "Baseline_QC_Metadata"):
                if r_id in rel_map:
                    t_path = rel_map[r_id].lstrip("/")
                    if not t_path.startswith("xl/"):
                        t_path = "xl/" + t_path
                    target_sheet_path = t_path
                    if s_name == "Cohort_Chromosome_Summary":
                        break

        sheet_root = ET.fromstring(zf.read(target_sheet_path))
        row_dict = {}
        for row_el in sheet_root.findall(".//main:row", ns):
            r_num = int(row_el.attrib.get("r", "0"))
            if r_num not in (1, 2, 3):
                continue
            row_dict[r_num] = {}
            for c_el in row_el.findall("main:c", ns):
                c_ref = c_el.attrib.get("r", "")
                c_idx = col_index_from_ref(c_ref)
                c_type = c_el.attrib.get("t", "")
                val_el = c_el.find("main:v", ns)
                is_el = c_el.find(".//main:t", ns)
                val_str = ""
                if c_type == "s" and val_el is not None and val_el.text is not None:
                    s_idx = int(val_el.text)
                    if 0 <= s_idx < len(shared_strings):
                        val_str = shared_strings[s_idx]
                elif c_type == "inlineStr" and is_el is not None and is_el.text is not None:
                    val_str = is_el.text
                elif val_el is not None and val_el.text is not None:
                    val_str = val_el.text
                row_dict[r_num][c_idx] = val_str.strip()

    samples_row = row_dict.get(1, {})
    g22_row = row_dict.get(2, {})
    g24_row = row_dict.get(3, {})

    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Sample_ID", "G22_Baseline", "G24_Ploidy"])
        for c_idx in sorted(samples_row.keys()):
            if c_idx < 2:
                continue
            s_id = samples_row[c_idx]
            if not s_id:
                continue
            try:
                g22_val = float(g22_row.get(c_idx, 1.0))
            except ValueError:
                g22_val = 1.0
            try:
                g24_val = float(g24_row.get(c_idx, 2.0))
            except ValueError:
                g24_val = 2.0
            w.writerow([s_id, round(g22_val, 4), round(g24_val, 4)])


def build_matrix_rows(samples, row_keys, records_by_sample, digits=4):
    rows = [["chrom"] + samples]
    for rk in row_keys:
        rows.append([rk] + [round(records_by_sample[s].get(rk, 0.0), digits) for s in samples])
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+")
    parser.add_argument("--species-map", default="c:cerevisiae,e:eubayanus,p:paradoxus,m:mikatae,k:kudriavzevii,u:uvarum,a:arboricola,j:jurei")
    parser.add_argument("--ploidy-multiplier", type=float, default=2.0)
    parser.add_argument("--auto-adjust-ploidy", action="store_true")
    parser.add_argument("--out-excel", default="lager_coverage_master.xlsx")
    parser.add_argument("--extract-curated-excel", default="")
    parser.add_argument("--out-ploidy-tsv", default="applied_ploidy_settings.tsv")
    args = parser.parse_args()

    if args.extract_curated_excel:
        extract_curated_xlsx_to_tsv(args.extract_curated_excel, args.out_ploidy_tsv)
        return

    if not args.inputs:
        raise ValueError("Must provide --inputs when building the coverage summary.")

    species_map = parse_species_map(args.species_map)
    species_keys = list(species_map.keys())
    os.makedirs("subgenome_rosters", exist_ok=True)

    raw_records = {}
    curated_records = {}
    g22_by_sample = {}
    g24_by_sample = {}
    window_cov_by_sample = {}
    canonical_windows = []

    meta_rows = [["Sample_ID", "G22_Baseline", "Default_G24_Ploidy", "Active_Chromosomes_Count", "Active_Chromosomes_Used"]]
    bp_rows = [["Sample_ID", "chrom", "Window_Start", "Window_End", "Left_Norm_Depth", "Right_Norm_Depth", "Delta"]]
    rosters = {k: [] for k in species_keys}
    all_chrom_set = set()
    samples = []

    for fpath in sorted(args.inputs):
        sample_id = os.path.basename(fpath).replace(".slidingwindow.tab", "")
        samples.append(sample_id)
        chrom_windows = {}
        sample_win_map = {}
        with open(fpath, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                c = parts[0]
                try:
                    s_pos, e_pos, cov = int(float(parts[1])), int(float(parts[2])), float(parts[3])
                except ValueError:
                    s_pos, e_pos, cov = 0, 0, 0.0
                chrom_windows.setdefault(c, []).append((s_pos, e_pos, cov))
                sample_win_map[(c, s_pos, e_pos)] = round(cov, 3)

        if not canonical_windows:
            for c in sorted(chrom_windows.keys(), key=lambda x: chrom_sort_key(x, species_keys)):
                for (s_pos, e_pos, _cov) in chrom_windows[c]:
                    canonical_windows.append((c, s_pos, e_pos))

        window_cov_by_sample[sample_id] = sample_win_map

        chrom_meds = {}
        for c, w_list in chrom_windows.items():
            all_chrom_set.add(c)
            if c.endswith("_P"):
                chrom_meds[c] = 0.0
            else:
                chrom_meds[c] = float(statistics.median([x[2] for x in w_list]))

        g22_base, g22_all_nuc, kept_chrs = compute_dynamic_baseline(chrom_meds)

        eff_g24 = args.ploidy_multiplier
        if args.auto_adjust_ploidy and args.ploidy_multiplier == 2.0:
            se_active = [
                v for c, v in chrom_meds.items()
                if c.startswith("e_") and not c.endswith(("_M", "_P")) and v >= 5.0
            ]
            sc_active = [
                v for c, v in chrom_meds.items()
                if c.startswith("c_") and not c.endswith(("_M", "_P")) and v >= 5.0
            ]
            if se_active and sc_active:
                if (float(statistics.median(se_active)) / g22_base) >= 1.8:
                    eff_g24 = 1.0

        g22_by_sample[sample_id] = round(g22_base, 4)
        g24_by_sample[sample_id] = eff_g24

        for prefix in species_keys:
            active_cnt = sum(
                1 for c, v in chrom_meds.items()
                if c.startswith(f"{prefix}_") and not c.endswith(("_M", "_P")) and v >= 5.0
            )
            if active_cnt >= 3:
                rosters[prefix].append(sample_id)

        c_dict = {}
        for c, raw_val in chrom_meds.items():
            if c.endswith("_P"):
                c_dict[c] = 0.0
            elif c.endswith("_M"):
                c_dict[c] = raw_val / g22_base if g22_base > 0 else 0.0
            else:
                c_dict[c] = (raw_val / g22_base) * eff_g24 if g22_base > 0 else 0.0

        raw_records[sample_id] = chrom_meds
        curated_records[sample_id] = c_dict

        meta_rows.append([
            sample_id, round(g22_base, 2), eff_g24,
            len(kept_chrs), ",".join(kept_chrs)
        ])

        for bp in detect_breakpoints(chrom_windows, g22_all_nuc):
            bp_rows.append([
                sample_id, bp["chrom"], bp["Window_Start"], bp["Window_End"],
                bp["Left_Norm_Depth"], bp["Right_Norm_Depth"], bp["Delta"]
            ])

    with open("default_ploidy_table.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Sample_ID", "G22_Baseline", "G24_Ploidy"])
        for s in samples:
            w.writerow([s, g22_by_sample[s], g24_by_sample[s]])

    for prefix, s_list in rosters.items():
        if s_list:
            with open(f"subgenome_rosters/callable_samples_{prefix}.txt", "w") as out_f:
                out_f.write("\n".join(sorted(s_list)) + "\n")

    all_chroms = sorted(all_chrom_set, key=lambda x: chrom_sort_key(x, species_keys))
    organellar_rows = [c for c in all_chroms if c.endswith(("_M", "_P"))]
    nuclear_rows = [c for c in all_chroms if not c.endswith(("_M", "_P"))]
    gathered_order = organellar_rows + nuclear_rows

    with open("chimeric_breakpoints.tsv", "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows(bp_rows)

    raw_med_rows = build_matrix_rows(samples, gathered_order, raw_records, 3)

    summary_formula_rows = [
        ["chrom"] + samples,
        ["G22_Baseline"] + [g22_by_sample[s] for s in samples],
        ["G24_Ploidy"] + [g24_by_sample[s] for s in samples]
    ]
    for idx_c, rk in enumerate(gathered_order):
        raw_sheet_row = idx_c + 2
        row_cells = [rk]
        for idx_s, s in enumerate(samples):
            col_let = col_letter(idx_s + 2)
            raw_cell = f"Raw_Chromosome_Medians!{col_let}{raw_sheet_row}"
            if rk.endswith("_P"):
                f_expr = "0"
                cached_v = 0.0
            elif rk.endswith("_M"):
                f_expr = f"IF({col_let}$2>0, {raw_cell}/{col_let}$2, 0)"
                cached_v = round(curated_records[s].get(rk, 0.0), 4)
            else:
                f_expr = f"IF({col_let}$2>0, ({raw_cell}/{col_let}$2)*{col_let}$3, 0)"
                cached_v = round(curated_records[s].get(rk, 0.0), 4)
            row_cells.append(("__FORMULA__", f_expr, cached_v))
        summary_formula_rows.append(row_cells)

    raw_win_rows = [["chrom", "start", "end"] + samples]
    for (c, s_pos, e_pos) in canonical_windows:
        raw_win_rows.append([c, s_pos, e_pos] + [window_cov_by_sample[s].get((c, s_pos, e_pos), 0.0) for s in samples])

    chr_row_ranges = {}
    for w_idx, (c, _s, _e) in enumerate(canonical_windows, start=2):
        if c not in chr_row_ranges:
            chr_row_ranges[c] = [w_idx, w_idx]
        else:
            chr_row_ranges[c][1] = w_idx

    # Discover all present species prefixes (1..8 species)
    present_prefixes = [
        pref for pref in species_keys
        if any(f"{pref}_{rom}" in chr_row_ranges for rom in ROMAN_ORDER)
    ]
    if not present_prefixes:
        present_prefixes = ["c"]

    first_sample = samples[0] if samples else "Sample_1"
    first_g22 = g22_by_sample.get(first_sample, 1.0)
    first_g24 = g24_by_sample.get(first_sample, 2.0)
    last_sample_col = col_letter(len(samples) + 3)
    last_summary_col = col_letter(len(samples) + 1)

    p0_name = f"S. {species_map.get(present_prefixes[0], present_prefixes[0])}"
    p1_name = f"S. {species_map.get(present_prefixes[1], present_prefixes[1])}" if len(present_prefixes) > 1 else "Subgenome 2"

    plotter_rows = [[
        "chrom", "start", "end", "Auto_Lookup_Raw_Depth", "Manual_Paste_Override",
        f"{p0_name} Panel", "Median Depth (x)", "Copy Number (n)",
        f"{p1_name} Panel", "Median Depth (x)", "Copy Number (n)",
        "", "Normalized_Copy_Number"
    ]]

    side_panel = {}
    side_panel[2] = [
        "DO NOT EDIT FORMULAS IN F:K TABLES", "", "",
        "LIVE DIAGRAMS UPDATE WHEN G20, G22, OR G24 CHANGES", "", ""
    ]

    # Build dynamic 2-species-per-tier tables in Columns F:K for all present species (1 to 8 species)
    tiers = [present_prefixes[i:i + 2] for i in range(0, len(present_prefixes), 2)]
    tier_layout = []

    for t_idx, pair in enumerate(tiers):
        if t_idx == 0:
            hdr_r = 3
            chr_r1, chr_r2 = 4, 19
            org_m1, org_m2 = 27, 28
            org_c1, org_c2 = 30, 31
        else:
            hdr_r = 34 + (t_idx - 1) * 24
            chr_r1, chr_r2 = hdr_r + 1, hdr_r + 16
            org_m1, org_m2 = hdr_r + 18, hdr_r + 19
            org_c1, org_c2 = hdr_r + 20, hdr_r + 21

        left_p = pair[0]
        right_p = pair[1] if len(pair) > 1 else None
        tier_layout.append((hdr_r, chr_r1, chr_r2, org_m1, org_m2, org_c1, org_c2, right_p is not None))

        left_sp = f"S. {species_map.get(left_p, left_p)}"
        right_sp = f"S. {species_map.get(right_p, right_p)}" if right_p else ""

        side_panel[hdr_r] = [
            f"{left_sp} Chromosome", "Median Depth (x)", "Copy Number (n)",
            f"{right_sp} Chromosome" if right_p else "", "Median Depth (x)" if right_p else "", "Copy Number (n)" if right_p else ""
        ]

        for c_i, rom in enumerate(ROMAN_ORDER, start=1):
            r_num = chr_r1 + c_i - 1
            l_chr = f"{left_p}_{rom}"
            l_r1, l_r2 = chr_row_ranges.get(l_chr, [2, 2])
            l_med = round(raw_records.get(first_sample, {}).get(l_chr, 0.0), 2)
            l_cn = round(curated_records.get(first_sample, {}).get(l_chr, 0.0), 3)

            if right_p:
                r_chr = f"{right_p}_{rom}"
                r_r1, r_r2 = chr_row_ranges.get(r_chr, [2, 2])
                r_med = round(raw_records.get(first_sample, {}).get(r_chr, 0.0), 2)
                r_cn = round(curated_records.get(first_sample, {}).get(r_chr, 0.0), 3)
                right_cells = [
                    f"median of {right_sp}, chr {c_i} ({rom})",
                    ("__FORMULA__", f"MEDIAN(D{r_r1}:D{r_r2})", r_med),
                    ("__FORMULA__", f"J{r_num}/G$22*G$24", r_cn)
                ]
            else:
                right_cells = ["", "", ""]

            side_panel[r_num] = [
                f"median of {left_sp}, chr {c_i} ({rom})",
                ("__FORMULA__", f"MEDIAN(D{l_r1}:D{l_r2})", l_med),
                ("__FORMULA__", f"G{r_num}/G$22*G$24", l_cn)
            ] + right_cells

        # Organellar (Mitochondrium & 2-micron Plasmid) rows for this tier
        lm_r1, lm_r2 = chr_row_ranges.get(f"{left_p}_M", [2, 2])
        lp_r1, lp_r2 = chr_row_ranges.get(f"{left_p}_P", [2, 2])
        l_m_f = ("__FORMULA__", f"MEDIAN(D{lm_r1}:D{lm_r2})", 0.0) if f"{left_p}_M" in chr_row_ranges else 0.0
        l_p_f = ("__FORMULA__", f"MEDIAN(D{lp_r1}:D{lp_r2})", 0.0) if f"{left_p}_P" in chr_row_ranges else 0.0

        if right_p:
            rm_r1, rm_r2 = chr_row_ranges.get(f"{right_p}_M", [2, 2])
            rp_r1, rp_r2 = chr_row_ranges.get(f"{right_p}_P", [2, 2])
            r_m_f = ("__FORMULA__", f"MEDIAN(D{rm_r1}:D{rm_r2})", 0.0) if f"{right_p}_M" in chr_row_ranges else 0.0
            r_p_f = ("__FORMULA__", f"MEDIAN(D{rp_r1}:D{rp_r2})", 0.0) if f"{right_p}_P" in chr_row_ranges else 0.0
            r_m_cells = [f"median of {right_sp}, Mitochondrium", r_m_f, ""]
            r_p_cells = [f"median of {right_sp}, Plasmid", r_p_f, ""]
            r_m_cn = [f"{right_sp} mtDNA copy number / 1n", ("__FORMULA__", f"J{org_m1}/G$22", 0.0), ""]
            r_p_cn = [f"{right_sp} plasmid copy number / 1n", ("__FORMULA__", f"J{org_m2}/G$22", 0.0), ""]
        else:
            r_m_cells, r_p_cells, r_m_cn, r_p_cn = ["", "", ""], ["", "", ""], ["", "", ""], ["", "", ""]

        side_panel[org_m1] = [f"median of {left_sp}, Mitochondrium", l_m_f, ""] + r_m_cells
        side_panel[org_m2] = [f"median of {left_sp}, Plasmid", l_p_f, ""] + r_p_cells
        side_panel[org_c1] = [f"{left_sp} mtDNA copy number / 1n", ("__FORMULA__", f"G{org_m1}/G$22", 0.0), ""] + r_m_cn
        side_panel[org_c2] = [f"{left_sp} plasmid copy number / 1n", ("__FORMULA__", f"G{org_m2}/G$22", 0.0), ""] + r_p_cn

    # 8-Species Safe G22 Formula: Uses Cohort_Chromosome_Summary active baseline when G20 is used, or falls back to active median
    g22_formula = (
        f"IF(ISNUMBER(E2), IF(MEDIAN(G4:G19)>0, MEDIAN(G4:G19), MEDIAN(J4:J19)), "
        f"IFERROR(INDEX(Cohort_Chromosome_Summary!$B$2:${last_summary_col}$2, 1, "
        f"MATCH($G$20, Cohort_Chromosome_Summary!$B$1:${last_summary_col}$1, 0)), MEDIAN(G4:G19)))"
    )

    side_panel[20] = ["Active Strain_ID (G20):", first_sample, "<-- Type Strain_ID here", "Or paste raw depths in Col E", "", ""]
    side_panel[22] = ["1x Baseline Depth (G22):", ("__FORMULA__", g22_formula, first_g22), "<-- Auto-baseline (or overwrite)", "", "", ""]
    side_panel[24] = ["Strain Ploidy Multiplier (G24):", first_g24, "<-- ENTER PLOIDY HERE!", "", "", ""]

    norm_cached_by_row = {}
    for w_idx, (c, s_pos, e_pos) in enumerate(canonical_windows, start=2):
        raw_val = window_cov_by_sample.get(first_sample, {}).get((c, s_pos, e_pos), 0.0)
        norm_val = round((raw_val / first_g22) * first_g24, 4) if first_g22 > 0 else 0.0
        norm_cached_by_row[w_idx] = norm_val

        lookup_formula = (
            f"IF(ISNUMBER(E{w_idx}), E{w_idx}, "
            f"IFERROR(INDEX(Raw_Sliding_Windows!$D{w_idx}:${last_sample_col}{w_idx}, 1, "
            f"MATCH($G$20, Raw_Sliding_Windows!$D$1:${last_sample_col}$1, 0)), 0))"
        )
        norm_formula = f"IF($G$22>0, (D{w_idx}/$G$22)*$G$24, 0)"

        sp_cols = side_panel.get(w_idx, ["", "", "", "", "", ""])
        row = [
            c,
            s_pos,
            e_pos,
            ("__FORMULA__", lookup_formula, raw_val),
            "",
            sp_cols[0], sp_cols[1], sp_cols[2],
            sp_cols[3], sp_cols[4], sp_cols[5],
            "",
            ("__FORMULA__", norm_formula, norm_val)
        ]
        plotter_rows.append(row)

    palette_map = {
        "c": "4F81BD", "e": "00B050", "u": "2A9D8F", "k": "9B5DE5",
        "p": "E07A5F", "m": "81B29A", "a": "F4A261", "j": "E76F51"
    }
    chart_specs = []
    for pref in present_prefixes:
        pref_chrs = [f"{pref}_{rom}" for rom in ROMAN_ORDER if f"{pref}_{rom}" in chr_row_ranges]
        if not pref_chrs:
            continue
        start_r = chr_row_ranges[pref_chrs[0]][0]
        end_r = chr_row_ranges[pref_chrs[-1]][1]
        total_wins = max(1, end_r - start_r + 1)

        chr_fractions = []
        for ch_name in pref_chrs:
            c_end = chr_row_ranges[ch_name][1]
            rom_lbl = ch_name.split("_", 1)[1] if "_" in ch_name else ch_name
            frac = (c_end - start_r + 1) / float(total_wins)
            chr_fractions.append((rom_lbl, frac))

        cached_vals = [norm_cached_by_row.get(r_i, 0.0) for r_i in range(start_r, end_r + 1)]
        sp_full = species_map.get(pref, pref)
        chart_specs.append({
            "title": f"S. {sp_full}",
            "col_let": "M",
            "start_row": start_r,
            "end_row": end_r,
            "cached_vals": cached_vals,
            "color": palette_map.get(pref, "4F81BD"),
            "chr_fractions": chr_fractions
        })

    nuc_start_row = len(organellar_rows) + 4
    sheets = {
        "Interactive_Plotter": plotter_rows,
        "Cohort_Chromosome_Summary": summary_formula_rows,
        "Raw_Chromosome_Medians": raw_med_rows,
        "Raw_Sliding_Windows": raw_win_rows,
        "Baseline_QC_Metadata": meta_rows,
        "Chimeric_Breakpoints": bp_rows
    }
    write_multi_sheet_xlsx(
        args.out_excel,
        sheets,
        color_scale_sheets={
            "Cohort_Chromosome_Summary": (nuc_start_row, 0.0, 2.0, 4.0)
        },
        chart_specs=chart_specs,
        tier_layout=tier_layout
    )


if __name__ == "__main__":
    main()