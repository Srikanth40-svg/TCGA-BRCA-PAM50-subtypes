# Breast cancer subtypes in TCGA-BRCA

Course project for **BT3041 – Analysis and Interpretation of Biological Data** (IIT Madras).

**Team:** Naveen, Srikanth, Piyush, Vimal, Madhavan

Breast cancer is usually split into five molecular subtypes (PAM50): Luminal A, Luminal B, HER2-enriched,
Basal-like and Normal-like. In this project we took the raw RNA-seq counts from the TCGA breast cancer
cohort and asked a few simple questions:

- Do the subtypes show up on their own if we just look at the expression data (PCA, UMAP, clustering)?
- Which genes and pathways actually differ between them?
- Can a machine-learning model predict the subtype from gene expression, and does it pick the "right" genes?
- Do the subtypes differ in survival once you account for stage and age?

## What we found (short version)

- Basal-like tumours are clearly their own group. LumA and LumB are not — they sit on a continuum, mainly
  separated by how fast the cells are dividing.
- Proliferation pathways (E2F targets, G2M checkpoint, MYC) are up in LumB, Her2 and Basal; estrogen-response
  genes are down in Her2 and Basal. Her2 tumours show the ERBB2 / GRB7 amplicon on chromosome 17q12.
- Logistic regression predicted the subtype with **92% accuracy (macro-F1 0.90)** on held-out patients, and
  33 of the 49 PAM50 genes present in our data ended up in its top 500 genes.
- Her2 and Basal patients had a higher risk of death than LumA (Cox HR 2.7 and 1.7), independent of stage and age.

The full write-up is in `report/`.

## Data

All data are public and are **not** stored in this repo (the count files alone are several GB).

| What | Where it came from |
|---|---|
| RNA-seq counts (STAR – Counts, GENCODE v36), 1,231 samples | [GDC portal](https://portal.gdc.cancer.gov/), project TCGA-BRCA |
| PAM50 subtype labels | `TCGA_BRCA_subtypes.csv` (TCGA PAM50 calls) |
| Clinical data (age, stage, ER/PR/HER2) | UCSC Xena – TCGA-BRCA clinical matrix |
| Overall survival | cBioPortal – PanCanAtlas 2018 clinical data |

After cleaning we kept **1,083 patients** (one primary tumour each): LumA 562, LumB 209, Basal 190, Her2 82,
Normal-like 40, and **15,959** expressed protein-coding genes.

## Repository layout

```
TCGA_BRCA_project/
│
├── manifests/                      GDC manifest + sample sheet (used to re-download the data)
│
├── src/
│   └── merge_counts.py             merges the 1,231 GDC count files into one matrix
│
├── notebooks/
│   ├── 01_data_cleaning.ipynb      PAM50 labels, one sample per patient, clinical + survival merge
│   ├── 02_preprocessing.ipynb      gene filtering, log2-CPM normalisation, clinical cleanup
│   ├── 03_eda_clinical.ipynb       subtype balance, age, stage, ER/PR/HER2
│   ├── 04_eda_expression.ipynb     marker genes, heatmaps, PCA, UMAP, batch-effect check
│   ├── 05_clustering.ipynb         k-means and hierarchical clustering vs PAM50
│   ├── 06_differential_expression.ipynb   Kruskal–Wallis + BH-FDR, PyDESeq2, volcano plots
│   ├── 07_pathway_enrichment.ipynb ORA (Fisher's exact test) and GSEA on Hallmark gene sets
│   ├── 08_classification.ipynb     5 ML models, cross-validation, feature importance
│   └── 09_survival_regression.ipynb  Kaplan–Meier, Cox model, linear regression
│
├── results/
│   ├── figures/                    every plot used in the report
│   └── tables/                     DE gene lists, GSEA tables, model scores, Cox results ...
│
├── report/                         final report, supplementary material, slides, LaTeX sources
│
├── data/                           (not in git) raw downloads and processed matrices
├── requirements.txt
└── README.md
```

## Running it yourself

You need Python 3.10+ and roughly 4 GB of free disk space.

**1. Set up the environment**

```bash
python -m venv .venv
.venv\Scripts\activate            # Windows
# source .venv/bin/activate       # macOS / Linux
pip install -r requirements.txt
```

**2. Download the RNA-seq counts** with the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool):

```bash
gdc-client download -m manifests/gdc_manifest.2026-10-25.011040.txt -d data/raw/gdc_data
```

**3. Put the clinical and label files in `data/raw/`:**

- `BRCA_clinicalMatrix.tsv` (UCSC Xena)
- `brca_tcga_pan_can_atlas_2018_clinical_data.tsv` (cBioPortal)
- `TCGA_BRCA_subtypes.csv` (PAM50 labels)

**4. Merge the count files** (run from the project root):

```bash
python src/merge_counts.py
```

This writes `TCGA_BRCA_raw_counts.csv`, `TCGA_BRCA_gene_info.csv` and `TCGA_BRCA_sample_info.csv` into `data/processed/`.

**5. Run the notebooks in order (01 → 09).** Each one reads the outputs of the previous ones from
`data/processed/` and saves its figures and tables to `results/`. Notebook 06 (DESeq2) takes a few minutes;
the rest run in under a minute each. Notebook 07 downloads the MSigDB Hallmark gene sets the first time it runs,
so it needs an internet connection.

## Main tools

pandas, NumPy, SciPy, statsmodels, scikit-learn, XGBoost, LightGBM, PyDESeq2, GSEApy, lifelines,
UMAP, matplotlib and seaborn.




