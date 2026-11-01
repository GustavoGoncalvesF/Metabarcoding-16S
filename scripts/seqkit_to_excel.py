
import sys
import re
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from openpyxl.formatting.rule import ColorScaleRule


# Input files from Snakemake
raw_file, cut_file, poly_file, final_file, output_file = sys.argv[1:6]

STAGES = [
    ("01_Raw", raw_file),
    ("02_Cutadapt", cut_file),
    ("03_BBDuk", poly_file),
    ("04_Final", final_file),
]

METRICS = ["num_seqs", "min_len", "avg_len", "max_len"]


###################################################
# 1. READ SEQKIT STATISTICS
###################################################
def read_seqkit(path):
    records = {}

    with open(path, encoding="utf-8") as handle:
        header = handle.readline().split()

        for line in handle:
            if not line.strip():
                continue

            fields = line.split()

            if len(fields) != len(header):
                raise ValueError(f"Invalid row in {path}: {line}")

            row = dict(zip(header, fields))
            filename = Path(row["file"]).name

            match = re.match(r"(.+)_R([12])\.", filename)

            if not match:
                raise ValueError(f"Invalid filename: {filename}")

            sample = match.group(1)
            read = f"R{match.group(2)}"
            key = (sample, read)

            if key in records:
                raise ValueError(f"Duplicate record: {key}")

            records[key] = {
                "num_seqs": int(row["num_seqs"].replace(",", "")),
                "min_len": int(row["min_len"]),
                "avg_len": float(row["avg_len"]),
                "max_len": int(row["max_len"]),
            }

    return records


stage_data = {
    stage: read_seqkit(path)
    for stage, path in STAGES
}

# Sample IDs
samples = sorted(
    {sample for sample, read in stage_data["01_Raw"]},
    key=lambda s: [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", s)
    ]
)

# Validate sample completeness and R1/R2 synchronisation
expected = {(sample, read) for sample in samples for read in ("R1", "R2")}

for stage, data in stage_data.items():
    if set(data) != expected:
        raise ValueError(f"Missing or unexpected samples in {stage}")

    for sample in samples:
        if data[(sample, "R1")]["num_seqs"] != data[(sample, "R2")]["num_seqs"]:
            raise ValueError(f"R1/R2 mismatch: {sample}, {stage}")


###################################################
# 2. CREATE EXCEL WORKBOOK
###################################################
wb = Workbook()
wb.remove(wb.active)

COLORS = {
    "01_Raw": "64748B",
    "02_Cutadapt": "2563EB",
    "03_BBDuk": "0F766E",
    "04_Final": "7C3AED",
    "05_Read_Tracking": "17365D",
}


def format_sheet(ws, color):
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = ws.dimensions
    ws.sheet_view.showGridLines = False
    ws.row_dimensions[1].height = 30

    for cell in ws[1]:
        cell.fill = PatternFill("solid", fgColor=color)
        cell.font = Font(color="FFFFFF", bold=True)
        cell.alignment = Alignment(
            horizontal="center",
            vertical="center",
            wrap_text=True
        )

    for col in range(1, ws.max_column + 1):
        letter = get_column_letter(col)
        ws.column_dimensions[letter].width = 18

    ws.column_dimensions["A"].width = 15

    for row in range(2, ws.max_row + 1):
        if row % 2 == 0:
            for cell in ws[row]:
                cell.fill = PatternFill(
                    "solid", fgColor="F1F5F9"
                )

        for cell in ws[row]:
            if isinstance(cell.value, (int, float)):
                cell.number_format = "#,##0.0" if isinstance(
                    cell.value, float
                ) else "#,##0"


###################################################
# 3. ONE SHEET PER PROCESSING STAGE
###################################################
for stage, path in STAGES:
    ws = wb.create_sheet(stage)
    data = stage_data[stage]

    ws.append([
        "SampleID",
        "R1_num_seqs",
        "R1_min_len",
        "R1_avg_len",
        "R1_max_len",
        "R2_num_seqs",
        "R2_min_len",
        "R2_avg_len",
        "R2_max_len"
    ])

    for sample in samples:
        values = [sample]

        for read in ("R1", "R2"):
            record = data[(sample, read)]
            values.extend(record[m] for m in METRICS)

        ws.append(values)

    format_sheet(ws, COLORS[stage])


###################################################
# 4. READ TRACKING
###################################################
ws = wb.create_sheet("05_Read_Tracking")

ws.append([
    "SampleID",
    "Raw",
    "Cutadapt",
    "BBDuk",
    "Final",
    "Ret_Cutadapt",
    "Ret_BBDuk",
    "Ret_Length",
    "Ret_Total",
    "Removed_Cutadapt",
    "Removed_BBDuk",
    "Removed_Length",
])

for i, sample in enumerate(samples, start=2):
    # Formulas linked to the four processing sheets
    ws.append([
        sample,
        f"='01_Raw'!B{i}",
        f"='02_Cutadapt'!B{i}",
        f"='03_BBDuk'!B{i}",
        f"='04_Final'!B{i}",
        f'=IFERROR(C{i}/B{i},0)',
        f'=IFERROR(D{i}/C{i},0)',
        f'=IFERROR(E{i}/D{i},0)',
        f'=IFERROR(E{i}/B{i},0)',
        f'=B{i}-C{i}',
        f'=C{i}-D{i}',
        f'=D{i}-E{i}',
    ])

format_sheet(ws, COLORS["05_Read_Tracking"])

# Percentage formatting
for row in range(2, ws.max_row + 1):
    for col in range(6, 10):
        ws.cell(row, col).number_format = "0.00%"

# Retention color scale
ws.conditional_formatting.add(
    f"I2:I{ws.max_row}",
    ColorScaleRule(
        start_type="num",
        start_value=0,
        start_color="F8696B",
        mid_type="num",
        mid_value=0.5,
        mid_color="FFEB84",
        end_type="num",
        end_value=1,
        end_color="63BE7B",
    )
)


###################################################
# 5. SAVE FILE
###################################################
Path(output_file).parent.mkdir(parents=True, exist_ok=True)

wb.save(output_file)

print(f"Excel report generated: {output_file}")
print(f"Number of samples: {len(samples)}")
print(f"Worksheets: {', '.join(wb.sheetnames)}")
