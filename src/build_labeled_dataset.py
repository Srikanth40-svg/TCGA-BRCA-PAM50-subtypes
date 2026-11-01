# Run from the project root:  python src/build_labeled_dataset.py
import pandas as pd
from pathlib import Path

# Input files
metadata_file = "data/processed/TCGA_BRCA_sample_info.csv"
counts_file = "data/processed/TCGA_BRCA_raw_counts.csv"
clinical_file = "data/raw/BRCA_clinicalMatrix.tsv"

# Load metadata and clinical subtype labels
meta = pd.read_csv(metadata_file, dtype=str)
clinical = pd.read_csv(clinical_file, sep="\t", low_memory=False)

# Confirm required columns exist
required_meta = {
    "Sample ID", "Tissue Type", "Tumor Descriptor", "File ID"
}
if not required_meta.issubset(meta.columns):
    raise ValueError("Required sample metadata columns are missing.")

if not {"sampleID", "PAM50Call_RNAseq"}.issubset(clinical.columns):
    raise ValueError("Required clinical subtype columns are missing.")

# Normalize identifiers to the TCGA sample-level barcode
meta["sample_key"] = meta["Sample ID"].str[:15]
clinical["sample_key"] = clinical["sampleID"].astype(str).str[:15]

# Restrict to primary tumour samples
primary = meta[
    (meta["Tissue Type"] == "Tumor") &
    (meta["Tumor Descriptor"] == "Primary")
].copy()

# Check count-matrix column names without loading the full matrix
count_columns = pd.read_csv(counts_file, nrows=0).columns
available_columns = set(count_columns)

# Keep one file per sample key; prefer aliquots ending in 01A
primary["priority"] = primary["Sample ID"].apply(
    lambda x: 0 if str(x).endswith("01A") else 1
)
primary = primary.sort_values(
    ["sample_key", "priority", "Sample ID"]
)
primary = primary.drop_duplicates("sample_key", keep="first")

# Get the subtype label for each clinical sample key
labels = clinical[
    ["sample_key", "PAM50Call_RNAseq"]
].drop_duplicates("sample_key")

primary = primary.merge(
    labels, on="sample_key", how="left", validate="one_to_one"
)

# Retain labelled samples whose count columns actually exist
primary = primary[
    primary["PAM50Call_RNAseq"].notna() &
    primary["Sample ID"].isin(available_columns)
].copy()

if primary.empty:
    raise ValueError("No labelled primary tumour samples matched the count matrix.")

# Read only the required count columns to reduce memory use
selected_samples = primary["Sample ID"].tolist()
usecols = ["gene_id"] + selected_samples

counts = pd.read_csv(
    counts_file,
    usecols=usecols,
    low_memory=False
)

# Save genes x samples raw-count matrix
counts.to_csv("data/processed/TCGA_BRCA_labeled_counts.csv", index=False)

# Save matching sample metadata and labels
output_meta = primary[
    [
        "Sample ID", "sample_key", "Case ID", "File ID",
        "Tissue Type", "Tumor Descriptor", "PAM50Call_RNAseq"
    ]
].rename(columns={
    "Sample ID": "count_matrix_sample_id",
    "Case ID": "patient_id",
    "PAM50Call_RNAseq": "subtype"
})

output_meta.to_csv("data/processed/TCGA_BRCA_labeled_sample_info.csv", index=False)

print("\nDataset creation complete.")
print("Unique labelled samples:", len(output_meta))
print("Unique patients:", output_meta["patient_id"].nunique())
print("\nSubtype counts:")
print(output_meta["subtype"].value_counts().to_string())
print("\nCount matrix shape (genes x samples):", counts.shape)
print("\nFiles created:")
print("  TCGA_BRCA_labeled_counts.csv")
print("  TCGA_BRCA_labeled_sample_info.csv")
