# KongLab-ImmuneMLSurvivalAnalysis
Machine learning research project studying immune gene signatures in breast cancer prognosis using relapse-free survival (RFS) as the clinical endpoint.

## Project Overview

This project identifies immune gene signatures predictive of breast cancer relapse-free survival using publicly available microarray datasets from GEO. The pipeline begins with raw Affymetrix CEL files and progresses through normalization, differential expression analysis, Cox regression feature selection, and ensemble survival modeling.

**Research Goal:** Validate and extend a previously published 8-gene immune signature (TSLP, BIRC5, S100B, MDK, S100P, RARRES3, BLNK, ACO1) for breast cancer prognosis.

---

## Datasets

| Dataset | Role | Samples (tumor) | Source |
|---------|------|-----------------|--------|
| GSE42568 | Training | 101 | GEO |
| GSE21653 | Test 1 | 266 | GEO |
| GSE20711 | Test 2 | 90 | GEO |
| GSE88770 | Test 3 | 117 | GEO |

All datasets use Affymetrix microarray expression data, normalized with RMA. Normal tissue samples are excluded prior to modeling.

---

## Pipeline

### 1. Data Acquisition — `load_data.ipynb` (R)
Downloads raw `.CEL` files from GEO using `GEOquery`.

### 2. Preprocessing — `preprocess_train.ipynb`, `preprocess_test_*.ipynb` (R)
- RMA normalization via `affy`
- Probe-to-gene mapping using the GPL570 platform annotation
- Multi-probe genes averaged; non-tumor samples filtered out
- Outputs expression matrices and clinical metadata CSVs for each dataset

### 3. Differential Expression Analysis — `deg_analysis.ipynb` (R)
- Linear model (tumor vs. normal) using `limma`
- Empirical Bayes shrinkage; Benjamini-Hochberg FDR correction
- Filters: p < 0.05, adj. p < 0.05, |log2FC| > 1
- **1,557 DEGs** identified → intersected with 3,120 immune genes → **263 immune DEGs**
- All 8 signature genes present in the final gene set

### 4. Univariate Cox Regression — `univariate_cox.ipynb` (Python)
- Per-gene Cox PH model for each of the 263 immune DEGs
- **55 significant genes** (p < 0.05): 18 danger (HR > 1), 37 protective (HR < 1)
- 7 of 8 signature genes validated (TSLP borderline at p = 0.069)

### 5. Penalized Cox Regression — `penalized_cox.ipynb` (Python)
- LASSO and Elastic Net via `scikit-survival` on 55 significant genes
- 5-fold GridSearchCV to select optimal alpha
- Both methods select **22 genes** (includes MDK, S100P, BLNK, ACO1)

### 6. Multivariate Cox Regression — `mutivariate_cox.ipynb` (Python)
- Simultaneous Cox model with all 22 penalized genes
- **4 independently significant genes**: S100P (danger), FBXL16, MICB, KLHL5 (protective)
- Risk score computed as a weighted sum of univariate betas × gene expression

### 7. Random Survival Forest — `random_survival_forests.ipynb` (Python)
- 1,000 trees; features: all 263 immune DEGs
- Parameters: `min_samples_split=20`, `min_samples_leaf=30`, `max_features="sqrt"`

| Split | C-index |
|-------|---------|
| Train | 0.799 |
| OOB | 0.634 |
| Test 1 (GSE21653) | 0.638 |
| Test 2 (GSE20711) | 0.579 |
| Test 3 (GSE88770) | 0.699 |

### 8. Gradient Boosting Survival Analysis — `gradient_boosting.ipynb` (Python)
- Cox PH loss with regression tree base learners via `scikit-survival`
- Early stopping: 10-iteration window, stops after 5 iterations without improvement
- Stopped at **62 base learners**

| Split | C-index |
|-------|---------|
| Train | 0.914 |
| Test 1 (GSE21653) | 0.548 |
| Test 2 (GSE20711) | 0.439 |
| Test 3 (GSE88770) | 0.586 |

The large train/test gap indicates overfitting; further regularization and cross-validation are needed.

---

## Key Results

**Signature Gene Validation:**

| Gene | Direction | Univariate p | Paper Match |
|------|-----------|--------------|-------------|
| S100P | Danger | 0.004 | ✓ |
| BIRC5 | Danger | 0.012 | ✓ |
| MDK | Danger | 0.026 | ✓ |
| S100B | Danger | 0.045 | ✓ |
| BLNK | Protective | 0.001 | ✓ |
| RARRES3 | Protective | 0.001 | ✓ |
| ACO1 | Protective | 0.033 | ✓ |
| TSLP | Danger | 0.069 | ~ |

**Best model:** Random Survival Forest (OOB C-index: 0.634, Test 3 C-index: 0.699)

---

## Repository Structure

```
kong-lab-ml-project/
├── data/
│   ├── cel_files/              # Raw Affymetrix .CEL files
│   ├── tar_files/              # Downloaded tar archives
│   ├── immune_genes.json       # 3,120 immune-related genes
│   └── transcription_factors.csv
├── datasets/
│   ├── csv_files/              # Preprocessed expression matrices & clinical metadata
│   └── rds_files/              # R-format datasets
├── models/
│   ├── rsf.joblib              # Fitted Random Survival Forest
│   └── est_cph_tree.joblib     # Fitted Gradient Boosting model
├── src/                        # Analysis notebooks
│   ├── load_data.ipynb
│   ├── preprocess_train.ipynb
│   ├── preprocess_test_one.ipynb
│   ├── preprocess_test_two.ipynb
│   ├── preprocess_test_three.ipynb
│   ├── clean_csv.ipynb
│   ├── deg_analysis.ipynb
│   ├── univariate_cox.ipynb
│   ├── penalized_cox.ipynb
│   ├── mutivariate_cox.ipynb
│   ├── random_survival_forests.ipynb
│   ├── gradient_boosting.ipynb
│   ├── survival_support_vector_machines.ipynb  # in progress
│   ├── cox_learning.ipynb          # learning only
│   ├── deg_learning.ipynb          # learning only
│   └── ml_extension_learning.ipynb # learning only
└── visuals/
    └── forest_plot.png         # Univariate Cox hazard ratio plot
```

---

## Installation

**R Libraries**
```
tidyverse, GEOquery, affy, limma, survival
```

**Python Libraries**
```
pandas, numpy, matplotlib, lifelines, scikit-learn, scikit-survival, joblib
```

---

## Notes

- Notebooks labeled `*_learning.ipynb` are not part of the primary analysis pipeline — they are used to independently learn and practice concepts.
- `survival_support_vector_machines.ipynb` is currently in progress.
- All models use only tumor samples (is_tumor == 1) filtered from each dataset.
