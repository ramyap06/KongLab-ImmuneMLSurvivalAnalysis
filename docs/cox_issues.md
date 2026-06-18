# Cox Pipeline — Reviewer Issues

Review of `src/scripts/cox.R`. Issues are grouped by severity.

---

## Critical Bugs (Invalidate Results)

### 1. `max(cv_fit$cvm)` selects the *worst* penalized model, not the best (line 86)

`cv.glmnet` with `family = "cox"` defaults to `type.measure = "deviance"` (partial-likelihood deviance). Deviance is a *loss* — lower is better, and `lambda.min` is chosen to *minimize* it. `max(cv_fit$cvm)` therefore returns the score at the most-over-penalized (worst-fitting) lambda in the path. The downstream comparison on line 251:

```r
if (enet_result$best_cindex >= lasso_result$best_cindex) enet_result else lasso_result
```

then selects the model with *higher* deviance, i.e., the *worse* model. The y-axis label "Concordance index (CV)" in `plot_cv_lambda` and the `geom_hline(yintercept = 0.5)` reference line confirm the intent was Harrell's C (where 0.5 = chance), but `type.measure = "C"` was never passed to `cv.glmnet`. The metric, the extraction logic, the comparison direction, and the plot label are all simultaneously wrong.

**Fix:** Pass `type.measure = "C"` to `cv.glmnet`, and extract the C-index as `cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.min]`.

---

### 2. Univariate Cox uses raw p < 0.05 across all genes — no multiple testing correction (lines 69, 229–230)

The entire downstream signature is built from genes that pass `p < 0.05` in single-gene Cox tests across potentially hundreds of genes. With 200 genes tested, you expect ~10 false positives by chance alone. These inflate the pool fed into the penalized step. No Bonferroni, Benjamini-Hochberg, or any other correction is applied at any stage. This is the single most common reason survival signatures fail independent validation.

**Fix:** Apply BH/FDR correction and filter on `adj_p < 0.05` in `run_univariate_cox`.

---

### 3. No train/test split — all performance metrics are in-sample (lines 228–310)

The full pipeline — univariate Cox → penalized Cox → multivariate Cox → risk score computation → KM stratification — is applied to the same ~100 patients throughout. The reported log-rank p-value and any C-index figures describe how well the model fits its own training data. For this sample size, overfitting-induced optimism can easily shift the apparent log-rank p-value by an order of magnitude.

**Fix:** Use bootstrap resampling or cross-validation for KM evaluation. At minimum, acknowledge the limitation explicitly.

---

### 4. Betas from a joint clinical+gene model are applied as a gene-only score (lines 287–299)

The multivariate Cox model is fit jointly on `penal_genes + CLINICAL_COVARS`. Gene betas are *partial effects* — already adjusted for ER status, grade, nodal status, and T-stage. `compute_risk_scores` then applies only the gene betas (clinical variables are absent from `penal_expr`) to produce a risk score. The linear predictor is incomplete, and the effect sizes are not interpretable as standalone gene contributions. Any claimed prognostic independence of the signature is circular.

**Fix:** Either include clinical variables in the scoring data, or fit a separate gene-only Cox model for scoring.

---

## Major Methodological Issues

### 5. p-values rounded before filtering (lines 57, 69)

`round(..., 4)` is applied before `p_value < P_CUTOFF`. A gene with true p = 0.049999 rounds to 0.0500 and is dropped. Genes with p < 0.0001 are stored as 0.0000, so forest plots print "p=0.000" — not acceptable in a publication without a "p < 0.0001" note.

**Fix:** Store full precision in data frames; apply rounding only at the display/label stage.

---

### 6. Median dichotomisation computed on training data (lines 300–304)

`train_median` is the median risk score over the same patients used to train all three models. Splitting by the in-sample median guarantees the groups will appear separated; it is not an independent evaluation.

**Fix:** If no external validation cohort is available, report this as a descriptive result and not a prospective performance claim.

---

### 7. Silent fallback to all penalized genes when multivariate model finds nothing (line 288)

```r
multi_sig_genes <- if (nrow(multi$sig_results) > 0) multi$sig_results$gene else penal_genes
```

If multivariate Cox finds no significant genes, the script falls back to scoring on *all* penalized genes using adjusted betas from a model in which they were non-significant. A KM plot is then produced and saved with no indication this is a null result.

**Fix:** Stop the pipeline and report a null result explicitly rather than silently producing misleading output.

---

### 8. Proportional hazards assumption never tested

There is no call to `cox.zph()` anywhere. For time-to-relapse endpoints in breast cancer, the PH assumption is frequently violated (early vs late relapse dynamics differ). Violating PH inflates or deflates HRs in a time-dependent way invisible in the model output.

**Fix:** Add `cox.zph()` checks after fitting each Cox model.

---

### 9. Clinical covariates silently coerced to numeric (lines 282–284)

```r
train_data[CLINICAL_COVARS] <- lapply(train_data[CLINICAL_COVARS], function(x) as.numeric(as.character(x)))
```

If any covariate is stored as character labels (e.g., "positive"/"negative" for ER status), `as.numeric(as.character(...))` produces all NAs. `complete.cases()` on the next line then silently discards those patients. Depending on clinical metadata encoding, a large fraction of the cohort could be dropped without warning.

**Fix:** Validate and pre-encode clinical variables explicitly before Cox modelling.

---

## Code Quality / Reproducibility

### 10. `sys.frame(1)$ofile` is NULL when run via `Rscript cox.R` (line 6)

`sys.frame(1)$ofile` is only set when a script is `source()`d interactively. Running via `Rscript` leaves it NULL, and `dirname(NULL)` silently returns `"."`, giving the wrong directory. This breaks the pipeline outside of an interactive session.

**Fix:** Use `commandArgs(trailingOnly = FALSE)` parsing or `here::here()` for portable path resolution.

---

### 11. `set.seed` is too far from the random operation (lines 22, 248–249)

The seed is set at load time, 200+ lines before `cv.glmnet`. Any upstream RNG call consumes random state and makes fold assignments de facto irreproducible across different execution contexts.

**Fix:** Call `set.seed(42)` immediately before the `cv.glmnet` calls.

---

### 12. 5-fold CV is borderline for n ≈ 100 (line 20)

With ~100 tumor patients, 5-fold CV produces ~20-patient test folds. Harrell's C estimated on 20 survival observations has very wide confidence intervals. 10-fold CV is standard at this sample size.

**Fix:** Change `N_FOLDS <- 10`.

---

## Why There Is Small Overlap With the Baseline Signature

The 8-gene baseline signature in `deg.R` is:

```r
SIGNATURE_GENES <- c("TSLP", "BIRC5", "S100B", "MDK", "S100P", "RARRES3", "BLNK", "ACO1")
```

`cox.R` reads `_deg_expression_matrix_with_rfs_legacy_immune.rds`, which at this point in the `deg.R` pipeline already contains only genes that are **both** DEGs **and** in the ImmPort legacy immune gene list.

Of the 8 baseline genes:

| Gene | Biology | In ImmPort? |
|---|---|---|
| BIRC5 (survivin) | Antiapoptotic | No |
| S100B | Calcium-binding protein | No |
| MDK (Midkine) | Growth factor | No |
| S100P | Calcium-binding protein | No |
| ACO1 (Aconitase) | Metabolic enzyme | No |
| RARRES3 | Retinoic acid responder | Unlikely |
| TSLP | Cytokine | Yes |
| BLNK | B cell signalling | Yes |

Most baseline signature genes are tumor biology or metabolic genes, not immune genes. They are filtered out in `deg.R` before the file `cox.R` reads is even written. The two analyses are operating on structurally different gene universes — the small overlap is by design, not a numerical coincidence.

If the intent is to evaluate whether the baseline 8-gene signature is prognostic, `cox.R` needs to either use the non-immune-filtered expression matrix (`_deg_expression_matrix_with_rfs.rds`) or hard-code those 8 genes as the input feature set.

---

## Summary Table

| Severity | Issue | Line(s) |
|---|---|---|
| Critical | `max(cv_fit$cvm)` inverts model selection | 86, 251 |
| Critical | No multiple testing correction in univariate Cox | 69 |
| Critical | No train/test split — all metrics are in-sample | 228–310 |
| Critical | Risk scores use partial-effect betas without clinical terms | 287–299 |
| Major | Median cutoff derived from training data | 300 |
| Major | Silent fallback to null-result genes for KM plotting | 288 |
| Major | PH assumption untested | — |
| Major | Clinical covariates silently NA-coerced, patients dropped | 282–284 |
| Minor | p-values rounded before threshold filtering | 57, 69 |
| Minor | `set.seed` too distant from random call | 22 |
| Minor | `sys.frame(1)$ofile` breaks `Rscript` execution | 6 |
| Minor | 5-fold CV underpowered for n ≈ 100 | 20 |
| Conceptual | Small baseline overlap is structural: ImmPort filters out most signature genes | deg.R:281 |
