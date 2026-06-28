# Analysis Pipeline

## Overview

The pipeline takes raw microarray CEL files from GEO, normalizes and filters them down to immune-relevant differentially expressed genes, then feeds those into survival ML models evaluated by C-index.

```
raw CEL files
    ↓  ReadAffy + RMA          load_data.r / preprocess.R
expression matrix (~11K genes per dataset)
    ↓  limma DEG               deg.R
differentially expressed genes (tumor vs normal)
    ↓  immune gene filter      deg.R
immune DEGs only
    ↓  transpose + join RFS    deg.R
expr_with_rfs.rds
    ↓
RSF / CoxBoost / SSVM          ml/
    ↓
C-index evaluation
```

---

## Stage 1 — Preprocessing (`load_data.r`, `preprocess.R`)

### Noise sources and how each is addressed

| Noise source | Step | Method |
|---|---|---|
| Array-level technical noise | RMA normalization | Background correction → quantile normalization → log2 transform |
| Probe-level noise | Probe collapse | Exclude ambiguous probes (`///`), mean of remaining probes per gene |
| Batch effects (multi-dataset) | ComBat | Empirical Bayes batch correction via `sva::ComBat()` |

### RMA (Robust Multi-array Average)

Performed via `affy::rma()`. Three sub-steps run internally:

1. **Background correction** — subtracts non-specific signal estimated from probe intensities
2. **Quantile normalization** — forces every array to the same intensity distribution, removing chip-to-chip scaling differences
3. **Log2 transform** — compresses dynamic range so high-expressing genes do not dominate

### Probe annotation and collapse

Probes are annotated using the GPL570 platform feature data (column 11 of `featureData`). For genes represented by multiple probes:
- Probes mapping to multiple genes (containing `///`) are excluded
- Remaining probes are averaged per gene (`mean`)

### ComBat batch correction

When merging multiple datasets, systematic inter-study shifts (different reagent lots, scanners, lab protocols) are removed using `sva::ComBat()`. It models each batch as having its own additive and multiplicative effect per gene, then removes both. The empirical Bayes shrinkage makes estimates stable even for genes with low sample counts.

**Key assumption:** batch and biology must not be perfectly confounded. If one dataset is all-normal and another all-tumor, ComBat cannot separate batch from signal.

PCA before/after ComBat is used to validate that samples no longer cluster by dataset.

### Reference implementations

- `sva` package vignette — canonical ComBat documentation with worked examples
- `inSilicoMerging` Bioconductor package — built specifically for merging GPL570 microarray studies with ComBat; vignette is a direct reference for the merge step

---

## Stage 2 — DEG Analysis (`deg.R`)

Differential expression is run on the training dataset (GSE42568) only, using `limma`.

### Steps

1. Build design matrix: tumor vs normal (one-hot encoded)
2. Fit per-gene linear model (`lmFit`)
3. Apply tumor-vs-normal contrast (`contrasts.fit`)
4. Empirical Bayes moderation (`eBayes`) — shrinks per-gene variance estimates toward a common prior, stabilising statistics for small cohorts
5. Extract DEGs: `adj.P.Val < 0.05`, `P.Value < 0.05`, `|logFC| > 1` (≥2-fold change)

### Immune gene filtering

DEGs are filtered against two gene catalogs:

| Catalog | Source | Output file |
|---|---|---|
| ImmPort (GO/Reactome) | `immune_genes.rds` | `_deg_expression_matrix_with_rfs.rds` |
| Legacy ImmPort | `legacy_immune_genes_ImmPort.txt` | `_deg_expression_matrix_with_rfs_legacy_immune.rds` |

The final matrix is transposed (samples × genes) and joined with RFS clinical columns. Non-tumor rows and rows with missing/zero RFS are dropped.

---

## Stage 3 — Survival ML Models

The ML models receive `expr_with_rfs.rds` — a ready-made matrix with tumor samples as rows, immune DEG expression values as features, and `rfs_time` / `rfs_event` columns. They do not perform normalization or filtering themselves.

### Models and R packages

| Model | Package | Notes |
|---|---|---|
| Random Survival Forest (RSF) | `randomForestSRC` | Ensemble of survival trees; handles high-dimensional data and interactions; provides variable importance |
| CoxBoost | `CoxBoost` | Gradient boosting for Cox model; built-in covariate selection via penalized likelihood |
| Survival SVM (SSVM) | `survivalsvm` | Supports regression, ranking, and hybrid kernel approaches |

### Vignettes for reference

```r
vignette("randomForestSRC")
vignette("CoxBoost")
vignette("survivalsvm")
```

The `randomForestSRC` package author (Hemant Ishwaran) also has tutorial PDFs with breast cancer examples. The **CRAN Task View: Survival Analysis** lists every relevant survival package grouped by method.

---

## Stage 4 — Evaluation

| Metric | Package | Function |
|---|---|---|
| Harrell C-index | `Hmisc` | `rcorr.cens()` |
| C-index with confidence intervals | `survcomp` | `concordance.index()` |
| Simple concordance | `survival` | `concordance()` |

C-index = 0.5 is random, 1.0 is perfect discrimination. Values above 0.7 are generally considered good for survival models.

---

## Data flow summary

| File | Location | Produced by | Consumed by |
|---|---|---|---|
| `{dataset}_raw_affy.rds` | `data/raw_affy_rds_files/` | `load_data.r` | `preprocess.R` |
| `{dataset}_metadata.rds` | `data/geo_metadata_rds_files/` | `load_data.r` | `preprocess.R` |
| `{dataset}_expression_matrix.rds` | `data/preprocessed_rds_files/` | `preprocess.R` | `deg.R` |
| `{dataset}_clinical_metadata.rds` | `data/preprocessed_rds_files/` | `preprocess.R` | `deg.R` |
| `{dataset}_deg_stats_all.rds` | `data/deg_rds_files/` | `deg.R` | reference |
| `{dataset}_deg_expression_matrix_with_rfs.rds` | `data/deg_rds_files/` | `deg.R` | ML models |
| `{dataset}_deg_expression_matrix_with_rfs_legacy_immune.rds` | `data/deg_rds_files/` | `deg.R` | ML models |
