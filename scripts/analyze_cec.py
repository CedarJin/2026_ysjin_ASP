#!/usr/bin/env python3
import csv
import json
from pathlib import Path
import html
import math
import os
import re
import statistics
import sys
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from xml.etree import ElementTree as ET


NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkgrel": "http://schemas.openxmlformats.org/package/2006/relationships",
}

INPUT = "raw_data/cec_raw_all_plate.xlsx"
OUTDIR = "outputs/cec_analysis"
CONFIG_PATH = Path(__file__).resolve().parents[1] / "private_inputs.json"
LOCAL_CONFIG = json.loads(CONFIG_PATH.read_text()) if CONFIG_PATH.exists() else {}
DEDICATED_SAMPLE_PREFIX = LOCAL_CONFIG.get("dedicated_sample_prefix", "DEDICATED-")
DEDICATED_BASENAME = LOCAL_CONFIG.get("dedicated_basename", "dedicated_cec")


def col_to_num(col):
    out = 0
    for ch in col:
        out = out * 26 + ord(ch.upper()) - 64
    return out


def num_to_col(num):
    out = ""
    while num:
        num, rem = divmod(num - 1, 26)
        out = chr(65 + rem) + out
    return out


def split_cell(ref):
    m = re.match(r"([A-Z]+)(\d+)", ref)
    return col_to_num(m.group(1)), int(m.group(2))


def read_shared_strings(zf):
    try:
        xml = zf.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(xml)
    strings = []
    for si in root.findall("main:si", NS):
        parts = []
        for t in si.findall(".//main:t", NS):
            parts.append(t.text or "")
        strings.append("".join(parts))
    return strings


def read_workbook(zf):
    wb = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rel_map = {
        r.attrib["Id"]: r.attrib["Target"].lstrip("/")
        for r in rels.findall("pkgrel:Relationship", NS)
    }
    sheets = []
    for sheet in wb.findall("main:sheets/main:sheet", NS):
        rid = sheet.attrib[f"{{{NS['rel']}}}id"]
        target = rel_map[rid]
        if not target.startswith("xl/"):
            target = "xl/" + target
        sheets.append((sheet.attrib["name"], target))
    return sheets


def read_sheet(zf, path, shared_strings):
    root = ET.fromstring(zf.read(path))
    cells = {}
    for c in root.findall(".//main:sheetData/main:row/main:c", NS):
        ref = c.attrib["r"]
        ctype = c.attrib.get("t")
        value = None
        if ctype == "inlineStr":
            t = c.find(".//main:t", NS)
            value = t.text if t is not None else ""
        else:
            v = c.find("main:v", NS)
            if v is None:
                continue
            raw = v.text
            if ctype == "s":
                value = shared_strings[int(raw)]
            else:
                try:
                    value = float(raw)
                except (TypeError, ValueError):
                    value = raw
        cells[ref] = value
    return cells


def value_at(cells, col, row):
    return cells.get(f"{num_to_col(col)}{row}")


def clean_label(value):
    if value is None:
        return None
    s = str(value).strip()
    return s if s else None


def is_nc(label):
    return clean_label(label) and clean_label(label).upper() == "NC"


def is_qc(label):
    label = clean_label(label)
    return label and label.upper().startswith("QC")


def is_layout_blank(label):
    label = clean_label(label)
    return label is None or label.upper() in {"NA", "N/A"}


def has_dedicated_sample_prefix(label):
    label = clean_label(label)
    return bool(label and label.startswith(DEDICATED_SAMPLE_PREFIX))


def mean(values):
    vals = [v for v in values if v is not None and not math.isnan(v)]
    return statistics.fmean(vals) if vals else None


def stdev(values):
    vals = [v for v in values if v is not None and not math.isnan(v)]
    return statistics.stdev(vals) if len(vals) > 1 else None


def fmt(value, digits=6):
    if value is None:
        return ""
    if isinstance(value, float):
        return round(value, digits)
    return value


def compute(input_path):
    with zipfile.ZipFile(input_path) as zf:
        shared = read_shared_strings(zf)
        sheets = read_workbook(zf)
        all_wells = []
        sample_rows = []
        plate_qc_means = {}

        for plate_name, path in sheets:
            cells = read_sheet(zf, path, shared)
            plate_wells = []
            # The named blocks are A2:M10, A13:M21, and A24:M32, but
            # their first row/column are plate coordinates. The actual
            # 8 x 12 well data are B3:M10, B14:M21, and B25:M32.
            for r_off, row in enumerate(range(3, 11)):
                lyse_row = row + 11
                super_row = row + 22
                for col in range(2, 14):
                    well = f"{chr(65 + r_off)}{col - 1}"
                    label = clean_label(value_at(cells, col, row))
                    lyse = value_at(cells, col, lyse_row)
                    sup = value_at(cells, col, super_row)
                    if label is None and lyse is None and sup is None:
                        continue
                    if not isinstance(lyse, (int, float)) or not isinstance(sup, (int, float)):
                        raw_pct = None
                    elif lyse + sup == 0:
                        raw_pct = None
                    else:
                        raw_pct = sup / (sup + lyse) * 100
                    plate_wells.append(
                        {
                            "plate": plate_name,
                            "well": well,
                            "sample": label,
                            "lyse": lyse,
                            "supernatant": sup,
                            "raw_pct_cec": raw_pct,
                        }
                    )

            nc_mean = mean([w["raw_pct_cec"] for w in plate_wells if is_nc(w["sample"])])
            for w in plate_wells:
                w["plate_nc_mean_raw_pct_cec"] = nc_mean
                w["nc_corrected_pct_cec"] = (
                    w["raw_pct_cec"] - nc_mean
                    if w["raw_pct_cec"] is not None and nc_mean is not None
                    else None
                )

            grouped = defaultdict(list)
            for w in plate_wells:
                if w["sample"] and not is_layout_blank(w["sample"]):
                    grouped[w["sample"]].append(w)

            plate_samples = []
            for sample, wells in grouped.items():
                corrected = [w["nc_corrected_pct_cec"] for w in wells]
                raw_vals = [w["raw_pct_cec"] for w in wells]
                row = {
                    "plate": plate_name,
                    "sample": sample,
                    "sample_type": "NC" if is_nc(sample) else "QC" if is_qc(sample) else "Sample",
                    "n_wells": len(wells),
                    "wells": ", ".join(w["well"] for w in wells),
                    "raw_pct_cec_mean": mean(raw_vals),
                    "plate_nc_mean_raw_pct_cec": nc_mean,
                    "pct_cec_mean": mean(corrected),
                    "pct_cec_sd": stdev(corrected),
                }
                plate_samples.append(row)

            qc_mean = mean([r["pct_cec_mean"] for r in plate_samples if r["sample_type"] == "QC"])
            plate_qc_means[plate_name] = qc_mean
            for row in plate_samples:
                row["plate_qc_mean_pct_cec"] = qc_mean
            all_wells.extend(plate_wells)
            sample_rows.extend(plate_samples)

        global_qc_mean = mean([r["pct_cec_mean"] for r in sample_rows if r["sample_type"] == "QC"])
        for row in sample_rows:
            row["global_qc_mean_pct_cec"] = global_qc_mean
            row["cec_index_plate_qc"] = (
                row["pct_cec_mean"] / row["plate_qc_mean_pct_cec"]
                if row["pct_cec_mean"] is not None and row["plate_qc_mean_pct_cec"]
                else None
            )
            row["cec_index_global_qc"] = (
                row["pct_cec_mean"] / global_qc_mean
                if row["pct_cec_mean"] is not None and global_qc_mean
                else None
            )

        return all_wells, sample_rows, global_qc_mean, plate_qc_means


def write_csv(path, rows, headers):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({h: fmt(row.get(h)) for h in headers})


def escape_attr(s):
    return html.escape(str(s), quote=True)


def cell_xml(row_idx, col_idx, value, style=0):
    ref = f"{num_to_col(col_idx)}{row_idx}"
    s_attr = f' s="{style}"' if style else ""
    if value is None or value == "":
        return f'<c r="{ref}"{s_attr}/>'
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if math.isnan(value):
            return f'<c r="{ref}"{s_attr}/>'
        return f'<c r="{ref}"{s_attr}><v>{value:.12g}</v></c>'
    return f'<c r="{ref}" t="inlineStr"{s_attr}><is><t>{escape_attr(value)}</t></is></c>'


def sheet_xml(rows, widths=None, freeze=False, autofilter=True):
    max_row = len(rows)
    max_col = max((len(r) for r in rows), default=1)
    dim = f"A1:{num_to_col(max_col)}{max_row}"
    cols = ""
    if widths:
        col_elems = []
        for i, width in enumerate(widths, start=1):
            col_elems.append(f'<col min="{i}" max="{i}" width="{width}" customWidth="1"/>')
        cols = "<cols>" + "".join(col_elems) + "</cols>"
    views = ""
    if freeze:
        views = (
            '<sheetViews><sheetView workbookViewId="0">'
            '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
            "</sheetView></sheetViews>"
        )
    row_elems = []
    for r_idx, row in enumerate(rows, start=1):
        cells = [cell_xml(r_idx, c_idx, value, 1 if r_idx == 1 else 0) for c_idx, value in enumerate(row, start=1)]
        row_elems.append(f'<row r="{r_idx}">{"".join(cells)}</row>')
    filter_xml = f'<autoFilter ref="{dim}"/>' if autofilter and max_row > 1 else ""
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<dimension ref="{dim}"/>{views}{cols}<sheetData>{"".join(row_elems)}</sheetData>{filter_xml}</worksheet>'
    )


def write_xlsx(path, sheets):
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    content_types = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
        '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>',
        '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>',
    ]
    for idx in range(1, len(sheets) + 1):
        content_types.append(
            f'<Override PartName="/xl/worksheets/sheet{idx}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )
    content_types.append("</Types>")
    workbook_sheets = "".join(
        f'<sheet name="{escape_attr(name[:31])}" sheetId="{idx}" r:id="rId{idx}"/>'
        for idx, (name, _, _) in enumerate(sheets, start=1)
    )
    wb_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{workbook_sheets}</sheets></workbook>"
    )
    wb_rels = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    ]
    for idx in range(1, len(sheets) + 1):
        wb_rels.append(
            f'<Relationship Id="rId{idx}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{idx}.xml"/>'
        )
    wb_rels.append(
        f'<Relationship Id="rId{len(sheets)+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
    )
    wb_rels.append("</Relationships>")
    styles = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
        '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
        '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        '</styleSheet>'
    )
    rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>'
        '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>'
        '</Relationships>'
    )
    core = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        '<dc:creator>Codex</dc:creator><cp:lastModifiedBy>Codex</cp:lastModifiedBy>'
        f'<dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>'
        f'<dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified></cp:coreProperties>'
    )
    app = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
        'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
        '<Application>Codex</Application></Properties>'
    )
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", "".join(content_types))
        zf.writestr("_rels/.rels", rels)
        zf.writestr("docProps/core.xml", core)
        zf.writestr("docProps/app.xml", app)
        zf.writestr("xl/workbook.xml", wb_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", "".join(wb_rels))
        zf.writestr("xl/styles.xml", styles)
        for idx, (_, rows, widths) in enumerate(sheets, start=1):
            zf.writestr(f"xl/worksheets/sheet{idx}.xml", sheet_xml(rows, widths=widths, freeze=True))


def main():
    input_path = sys.argv[1] if len(sys.argv) > 1 else INPUT
    outdir = sys.argv[2] if len(sys.argv) > 2 else OUTDIR
    os.makedirs(outdir, exist_ok=True)

    wells, samples, global_qc_mean, plate_qc_means = compute(input_path)

    well_headers = [
        "plate",
        "well",
        "sample",
        "lyse",
        "supernatant",
        "raw_pct_cec",
        "plate_nc_mean_raw_pct_cec",
        "nc_corrected_pct_cec",
    ]
    sample_headers = [
        "plate",
        "sample",
        "sample_type",
        "n_wells",
        "wells",
        "raw_pct_cec_mean",
        "plate_nc_mean_raw_pct_cec",
        "pct_cec_mean",
        "pct_cec_sd",
        "plate_qc_mean_pct_cec",
        "global_qc_mean_pct_cec",
        "cec_index_plate_qc",
        "cec_index_global_qc",
    ]
    dedicated_samples = [r for r in samples if has_dedicated_sample_prefix(r.get("sample", ""))]
    dedicated_wells = [r for r in wells if has_dedicated_sample_prefix(r.get("sample", ""))]
    main_samples = [r for r in samples if not has_dedicated_sample_prefix(r.get("sample", ""))]
    main_wells = [r for r in wells if not has_dedicated_sample_prefix(r.get("sample", ""))]

    summary_rows = [
        ["Metric", "Value"],
        ["Input file", input_path],
        ["Plate count", len(plate_qc_means)],
        ["Well rows", len(main_wells)],
        ["Sample rows incl. NC/QC", len(main_samples)],
        [f"Dedicated {DEDICATED_SAMPLE_PREFIX} sample rows excluded from main outputs", len(dedicated_samples)],
        ["Global QC mean %CEC", fmt(global_qc_mean)],
        ["Formula", "%CEC = supernatant / (supernatant + lyse) * 100 - plate NC mean raw %CEC"],
        ["Plate QC index", "sample %CEC mean / same-plate QC mean %CEC"],
        ["Global QC index", "sample %CEC mean / all-plate QC mean %CEC"],
    ]
    summary_rows.append(["", ""])
    summary_rows.append(["Plate", "Plate QC mean %CEC"])
    for plate, qc_mean in plate_qc_means.items():
        summary_rows.append([plate, fmt(qc_mean)])

    sample_rows = [[h for h in sample_headers]]
    sample_rows.extend([[fmt(row.get(h)) for h in sample_headers] for row in main_samples])
    well_rows = [[h for h in well_headers]]
    well_rows.extend([[fmt(row.get(h)) for h in well_headers] for row in main_wells])

    write_csv(os.path.join(outdir, "cec_sample_summary.csv"), main_samples, sample_headers)
    write_csv(os.path.join(outdir, "cec_well_level.csv"), main_wells, well_headers)
    xlsx_path = os.path.join(outdir, "cec_analyzed.xlsx")
    write_xlsx(
        xlsx_path,
        [
            ("Summary", summary_rows, [26, 72]),
            ("Sample Summary", sample_rows, [14, 18, 12, 10, 18, 18, 24, 16, 14, 20, 22, 18, 18]),
            ("Well Level", well_rows, [14, 10, 18, 12, 14, 16, 24, 20]),
        ],
    )

    dedicated_sample_path = os.path.join(outdir, f"{DEDICATED_BASENAME}_summary.csv")
    dedicated_well_path = os.path.join(outdir, f"{DEDICATED_BASENAME}_well_level.csv")
    dedicated_xlsx_path = os.path.join(outdir, f"{DEDICATED_BASENAME}.xlsx")
    if dedicated_samples:
        dedicated_order = {row["sample"]: idx for idx, row in enumerate(dedicated_samples)}
        dedicated_wells = sorted(
            dedicated_wells,
            key=lambda row: (dedicated_order.get(row["sample"], 999), row["well"]),
        )
        write_csv(dedicated_sample_path, dedicated_samples, sample_headers)
        write_csv(dedicated_well_path, dedicated_wells, well_headers)

        dedicated_summary_rows = [
            ["Metric", "Value"],
            ["Input file", input_path],
            ["Dedicated sample prefix", DEDICATED_SAMPLE_PREFIX],
            ["Dedicated sample rows", len(dedicated_samples)],
            ["Dedicated well rows", len(dedicated_wells)],
            ["Plate7 QC mean %CEC", fmt(plate_qc_means.get("plate7"))],
            ["Global QC mean %CEC", fmt(global_qc_mean)],
            ["Formula", "%CEC = supernatant / (supernatant + lyse) * 100 - plate NC mean raw %CEC"],
            ["Plate QC index", "sample %CEC mean / same-plate QC mean %CEC"],
            ["Global QC index", "sample %CEC mean / all-plate QC mean %CEC"],
        ]
        dedicated_sample_rows = [[h for h in sample_headers]]
        dedicated_sample_rows.extend([[fmt(row.get(h)) for h in sample_headers] for row in dedicated_samples])
        dedicated_well_rows = [[h for h in well_headers]]
        dedicated_well_rows.extend([[fmt(row.get(h)) for h in well_headers] for row in dedicated_wells])
        write_xlsx(
            dedicated_xlsx_path,
            [
                ("Summary", dedicated_summary_rows, [28, 72]),
                ("Sample Summary", dedicated_sample_rows, [14, 22, 12, 10, 18, 18, 24, 16, 14, 20, 22, 18, 18]),
                ("Well Level", dedicated_well_rows, [14, 10, 22, 12, 14, 16, 24, 20]),
            ],
        )

    for old_name in [
        "plate7_new_samples_cec_summary.csv",
        "plate7_new_samples_cec_well_level.csv",
        "plate7_new_samples_cec_results.xlsx",
    ]:
        old_path = os.path.join(outdir, old_name)
        if os.path.exists(old_path):
            os.remove(old_path)

    print(f"Wrote {xlsx_path}")
    if dedicated_samples:
        print(f"Wrote {dedicated_xlsx_path}")
    print(f"Global QC mean %CEC: {global_qc_mean:.6f}" if global_qc_mean is not None else "Global QC mean %CEC: NA")
    for plate, qc_mean in plate_qc_means.items():
        print(f"{plate}: plate QC mean %CEC = {qc_mean:.6f}" if qc_mean is not None else f"{plate}: plate QC mean %CEC = NA")


if __name__ == "__main__":
    main()
