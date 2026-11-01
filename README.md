# TCGA-BRCA PAM50 Subtype Analysis

BT3041 – Analysis and Interpretation of Biological Data: group term project.

Team: Naveen, Srikanth, Piyush, Madhavan,Vimal

## Aim
Use TCGA-BRCA RNA-seq and clinical data to (1) explore expression structure and test whether
unsupervised clusters recover the PAM50 subtypes, (2) identify differentially expressed genes and
pathways between subtypes (FDR-corrected), (3) build ML classifiers to predict PAM50 subtype and
find the most informative genes, and (4) relate subtypes to clinical variables and survival,
validating against published studies.

## Project structure
```
data/raw/          GDC STAR counts (gdc_data/), clinical files, PAM50 labels   (not in git)
data/processed/    merged / cleaned matrices produced by src/ and notebooks/  (not in git)
manifests/         GDC manifest + sample sheet (to re-download the data)
src/               data-merging scripts
notebooks/         analysis notebooks (01 ... 08)
results/figures/   saved plots
results/tables/    result tables (DE genes, model metrics, ...)
report/            final report & slides
```

## Reproducing
1. Install dependencies
   ```
   python -m venv .venv
   .venv\Scripts\activate        # Windows  (source .venv/bin/activate on Linux/macOS)
   pip install -r requirements.txt
   ```
2. Download RNA-seq counts from GDC (needs the [gdc-client](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)):
   ```
   gdc-client download -m manifests/gdc_manifest.2026-10-25.011040.txt -d data/raw/gdc_data
   ```
3. Clinical / label files go in `data/raw/`:
   - `BRCA_clinicalMatrix.tsv` – UCSC Xena TCGA-BRCA clinical matrix
   - `brca_tcga_pan_can_atlas_2018_clinical_data.tsv` – cBioPortal PanCanAtlas clinical data
   - `TCGA_BRCA_subtypes.csv` – PAM50 subtype labels
4. Merge the count files (run from the project root):
   ```
   python src/merge_counts.py
   ```
5. Run the notebooks in `notebooks/` in order.

## Data
TCGA-BRCA, GDC harmonised RNA-seq (STAR – Counts, GENCODE v36): 1,231 samples x 60,660 genes.
After cleaning: 1,083 patients with PAM50 labels (LumA 562, LumB 209, Basal 190, Her2 82, Normal-like 40).
