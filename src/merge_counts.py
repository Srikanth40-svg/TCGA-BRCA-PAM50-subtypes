# Run from the project root:  python src/merge_counts.py
import pandas as pd, os

sheet = pd.read_csv("manifests/gdc_sample_sheet.2026-10-24.tsv", sep="\t")
print(sheet.columns.tolist())

# Sample ID is the TCGA barcode (use the first one if a file lists several)
sheet["Sample ID"] = sheet["Sample ID"].astype(str).str.split(",").str[0].str.strip()
# Characters 14-15 of the barcode give the sample type (01 = primary tumor, 11 = normal, ...)
code = sheet["Sample ID"].str[13:15]
sheet["Sample Type"] = code.map({
    "01": "Primary Tumor", "02": "Recurrent Tumor", "06": "Metastatic",
    "11": "Solid Tissue Normal"}).fillna("Other")
print(sheet["Sample Type"].value_counts())

counts, gene_info, missing = {}, None, 0
for _, row in sheet.iterrows():
    path = os.path.join("data/raw/gdc_data", row["File ID"], row["File Name"])
    if not os.path.exists(path):
        missing += 1
        continue
    df = pd.read_csv(path, sep="\t", comment="#")
    df = df[~df["gene_id"].str.startswith("N_")]          # drop QC rows
    if gene_info is None:
        gene_info = df[["gene_id", "gene_name", "gene_type"]].set_index("gene_id")
    counts[row["File ID"]] = df.set_index("gene_id")["unstranded"]   # raw counts

print("missing files:", missing)
matrix = pd.DataFrame(counts)
print("matrix shape:", matrix.shape)

info = sheet.set_index("File ID").loc[matrix.columns]
info.index.name = "File ID"
matrix.columns = info["Sample ID"].values
print("duplicate sample IDs:", matrix.columns.duplicated().sum())

matrix.to_csv("data/processed/TCGA_BRCA_raw_counts.csv")
gene_info.to_csv("data/processed/TCGA_BRCA_gene_info.csv")
info.reset_index().to_csv("data/processed/TCGA_BRCA_sample_info.csv", index=False)
print("saved")