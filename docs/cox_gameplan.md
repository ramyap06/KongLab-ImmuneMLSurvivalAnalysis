# Cox Pipeline — Fix Game Plan

This document maps each problem from `cox_issues.md` to a concrete code change. Work through the steps in order — each one builds on the previous.

---

## Step 1 — Fix the CV metric (lines 81–88)

**What's wrong:** `cv.glmnet` is using deviance (lower = better) but the code calls `max()` on it, picks the worst model, and labels it a "C-index".

**What to change in `run_penalized_cox`:**

```r
# BEFORE
cv_fit     <- cv.glmnet(X_train, y_train, family = "cox",
                        alpha = alpha, nfolds = nfolds, standardize = TRUE)
best_coefs <- as.vector(coef(cv_fit, s = "lambda.min"))
cindex     <- max(cv_fit$cvm)   # WRONG — this is max deviance, not C-index

# AFTER
cv_fit     <- cv.glmnet(X_train, y_train, family = "cox",
                        alpha = alpha, nfolds = nfolds, standardize = TRUE,
                        type.measure = "C")   # explicitly request concordance index
best_coefs <- as.vector(coef(cv_fit, s = "lambda.min"))
cindex     <- cv_fit$cvm[which(cv_fit$lambda == cv_fit$lambda.min)]  # C-index at best lambda
```

The `plot_cv_lambda` function already has the right labels ("Concordance index", 0.5 reference line) — it will be correct once the metric itself is fixed.

---

## Step 2 — Add FDR correction to univariate Cox (lines 66–74)

**What's wrong:** Filtering on raw p < 0.05 across hundreds of genes produces many false positives that contaminate the penalized step.

**What to change in `run_univariate_cox`:**

```r
# Add this line after building all_results, before filtering:
all_results$adj_p_value <- p.adjust(all_results$p_value, method = "BH")

# Then filter on adjusted p:
sig_results <- all_results[all_results$adj_p_value < P_CUTOFF, ]
```

Also add `adj_p_value` to the data.frame returned per gene (or compute it post-hoc on the combined table as shown above). Update `P_CUTOFF` documentation comments to clarify it now applies to FDR-adjusted p-values.

---

## Step 3 — Store full-precision p-values; round only at display (lines 51–58)

**What's wrong:** `round(..., 4)` baked into the data frame before filtering means borderline genes get dropped and very small p-values appear as 0.0000 in plots.

**What to change:**

Remove `round()` from every `data.frame(...)` in `run_univariate_cox` and `run_multivariate_cox`. Store raw values. Then apply rounding only in the forest plot label:

```r
# In plot_forest, the geom_text line already formats to 3 decimal places.
# For very small p-values, add a conditional label:
label = ifelse(p_value < 0.001, "p<0.001", sprintf("p=%.3f", p_value))
```

---

## Step 4 — Fix the `set.seed` placement (line 22 → move to ~line 247)

**What's wrong:** The seed is set at module load, not next to the random operation.

**What to change:** Delete `set.seed(42)` from line 22 and add it right before the `cv.glmnet` calls in the main script body:

```r
# Before the penalized cox block:
set.seed(42)
lasso_result <- run_penalized_cox(X_train, y_train, LASSO_ALPHA, N_FOLDS, "LASSO")
enet_result  <- run_penalized_cox(X_train, y_train, ENET_ALPHA,  N_FOLDS, "Elastic Net")
```

---

## Step 5 — Increase CV folds to 10 (line 20)

One-line change:

```r
N_FOLDS <- 10
```

---

## Step 6 — Validate and encode clinical variables explicitly (lines 280–285)

**What's wrong:** `as.numeric(as.character(...))` silently converts non-numeric strings to NA and drops patients via `complete.cases()` without any warning.

**What to change:** Before the merge, add explicit encoding with a loud warning:

```r
# After merging penal_expr with clinical_sel:
for (col in CLINICAL_COVARS) {
    n_before <- sum(!is.na(train_data[[col]]))
    train_data[[col]] <- suppressWarnings(as.numeric(as.character(train_data[[col]])))
    n_after  <- sum(!is.na(train_data[[col]]))
    if (n_after < n_before)
        warning("Column '", col, "': ", n_before - n_after,
                " values coerced to NA — check encoding in clinical metadata.")
}
train_data <- train_data[complete.cases(train_data), ]
message("  Patients with complete clinical data: ", nrow(train_data))
```

---

## Step 7 — Fix risk scoring to use a consistent patient set (lines 287–299)

**What's wrong:** The model is fit on `train_data` (patients with complete clinical data) but scores are computed on `penal_expr` (all patients). This means risk scores for patients with missing clinical data are extrapolations from a model they never influenced.

**Two options — pick one:**

**Option A (simpler):** Score only the patients who were in `train_data`:

```r
scored <- compute_risk_scores(train_data, multi$all_results, multi_sig_genes)
```

**Option B (more principled):** Fit a separate gene-only Cox model (no clinical covariates) for scoring, so the gene betas are not adjusted for clinical variables:

```r
gene_only_cox  <- coxph(as.formula(sprintf("Surv(rfs_time, rfs_event) ~ %s",
                         paste(sprintf("`%s`", penal_genes), collapse = " + "))),
                        data = train_data)
gene_only_coefs <- data.frame(gene = names(coef(gene_only_cox)),
                               beta = coef(gene_only_cox))
scored <- compute_risk_scores(penal_expr, gene_only_coefs, penal_genes)
```

Option B is preferred if you want to claim the gene score has independent prognostic value.

---

## Step 8 — Remove the silent null-result fallback (line 288)

**What's wrong:** When no genes pass multivariate significance, the pipeline silently scores on all penalized genes and produces a KM plot.

**What to change:**

```r
# BEFORE
multi_sig_genes <- if (nrow(multi$sig_results) > 0) multi$sig_results$gene else penal_genes

# AFTER
if (nrow(multi$sig_results) == 0) {
    message("No genes reached p < ", P_CUTOFF, " in multivariate Cox. ",
            "Skipping risk scoring and KM plot — report as null result.")
    # save the model results but stop here
    saveRDS(multi, file.path(COX_DIR, paste0(DATASET, "_multivariate_cox_null.rds")))
    quit(save = "no", status = 0)
}
multi_sig_genes <- multi$sig_results$gene
```

---

## Step 9 — Add PH assumption check after each Cox model

Add a `cox.zph()` call after fitting any `coxph` object. The simplest place is at the end of `run_multivariate_cox`:

```r
ph_test <- cox.zph(cph)
message("  PH assumption test (global p-value): ", round(ph_test$table["GLOBAL", "p"], 4))
if (ph_test$table["GLOBAL", "p"] < 0.05)
    warning("PH assumption may be violated (global p < 0.05). Check cox.zph() plots.")
```

Do the same at the end of `run_univariate_cox` (on a per-gene basis it's noisier, but a GLOBAL check on the top few genes is informative).

---

## Step 10 — Fix the root path resolution (line 6)

**What's wrong:** `sys.frame(1)$ofile` is NULL when the script is run via `Rscript`.

**What to change:**

```r
# BEFORE
ROOT_DIR <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))

# AFTER — works for both source() and Rscript
script_path <- tryCatch(
    normalizePath(sys.frame(1)$ofile),
    error = function(e) normalizePath(commandArgs(trailingOnly = FALSE)[
        startsWith(commandArgs(trailingOnly = FALSE), "--file=")
    ][1] |> sub("--file=", "", x = _))
)
ROOT_DIR <- normalizePath(file.path(dirname(script_path), "../.."))
```

Or just use the `here` package if it is already in the project dependencies:

```r
library(here)
ROOT_DIR <- here::here()
```

---

## Summary Checklist

| Step | Change | Severity Addressed |
|---|---|---|
| 1 | Add `type.measure = "C"` to `cv.glmnet`; fix C-index extraction | Critical |
| 2 | Apply BH correction in `run_univariate_cox` | Critical |
| 3 | Remove `round()` from data frames; round at display only | Minor |
| 4 | Move `set.seed(42)` next to `cv.glmnet` calls | Minor |
| 5 | `N_FOLDS <- 10` | Minor |
| 6 | Warn explicitly when clinical variable coercion produces NAs | Major |
| 7 | Score on consistent patient set or use gene-only Cox for scoring | Critical |
| 8 | Stop pipeline explicitly on null multivariate result | Major |
| 9 | Add `cox.zph()` check after Cox fits | Major |
| 10 | Fix `ROOT_DIR` for `Rscript` execution | Minor |

Steps 1, 2, and 7 are the highest priority — they directly affect which genes end up in the signature and how its performance is reported.

---

## What to Do About the Baseline Signature Overlap

The small overlap is not a bug in `cox.R` — it is a consequence of using a different input file. The legacy-immune filtered matrix excludes most of the 8-gene baseline signature (BIRC5, S100B, MDK, S100P, ACO1, RARRES3) because those genes are not in ImmPort.

**Decision to make before writing code:**

- If the goal is to *validate* the published 8-gene signature prognostically: hard-code those genes as the feature set and skip the DEG-filtering step entirely. Use `_deg_expression_matrix_with_rfs.rds` as input.
- If the goal is to *discover a new* immune survival signature from scratch: the current input file is correct, but acknowledge explicitly in any write-up that this signature is entirely independent of the original paper's signature.
- If the goal is to *compare* the two: run both in parallel and report the overlap.

These are scientifically different questions and the choice should be made deliberately, not by accident of which `.rds` file happens to be read on line 28.
