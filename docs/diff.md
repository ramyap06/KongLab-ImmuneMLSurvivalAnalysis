# Pipeline Decisions & Changes

## preprocess.R

### Drop non-tumor samples
`process_clinical_metadata` now filters to `is_tumor == 1` only before any other cleaning. Non-tumor (normal) samples are excluded entirely. Expression data is naturally synced via `sync_expression_to_clinical`, which intersects on `geo_accession`.

### prepare_combined_data — ComBat-ready format
Rewrote to avoid character coercion and produce the correct shape for ComBat:
- Load RDS as a data frame, separate `gene_symbol` before calling `as.matrix` so the matrix is numeric from the start (`storage.mode(ex) <- "numeric"`)
- No transpose — keeps **genes × samples** throughout (what ComBat expects)
- Use `cbind` to combine datasets instead of transpose + `rbind`
- Returns `list(expr, batch, datasets)` — `expr` is genes × samples, `batch` is a vector of dataset labels per sample

### perform_batch_correction
- Fixed variable name bug: was referencing `result$` instead of `data$` (the actual parameter name)
- Removed unused `datasets` parameter
- ComBat is called with defaults: `par.prior = TRUE` (parametric empirical Bayes)

### export_batch_corrected_datasets
- Fixed stray `\` before closing `}` that caused a syntax error
- Exports expression matrix only to `BATCH_DIR`; clinical metadata remains in `PRE_DIR` (batch correction does not change clinical data)

---

## pca_batch_check.ipynb

### prepare_preprocessed_data — coercion fix
Same fix as `prepare_combined_data`: separate `gene_symbol` before `as.matrix` to avoid the whole matrix being coerced to character. Numeric conversion with `storage.mode(ex) <- "numeric"` happens immediately, before the transpose and `rbind`.

### prepare_batch_corrected_data — separate function
Batch-corrected RDS files are already clean numeric matrices with gene rownames (no `gene_symbol` column). Rather than adding if/else logic to `prepare_preprocessed_data`, a dedicated function handles this cleaner format directly.

### LISI bar plot
- Used `lisi::compute_lisi` on PC1/PC2 coordinates to score batch mixing
- Score range: 1 (no mixing, full batch effect) → 4 (perfect mixing across 4 datasets)
- Observed score: **2–3**, indicating partial batch correction
- Expected outcome: perfect mixing (score = 4) is unlikely across 4 independent cohorts due to real biological variation between patient populations and protocols

### Batch correction assessment
LISI score of 2–3 is acceptable. ComBat with parametric empirical Bayes meaningfully reduced batch effects without overcorrecting biological signal. PCA confirmed improved distribution after correction.

---

## Key design decisions

| Decision | Rationale |
|---|---|
| Drop normal samples early in clinical processing | Normal samples have no RFS data and are irrelevant to survival analysis |
| genes × samples orientation in `prepare_combined_data` | ComBat requires this orientation; PCA notebook transposes internally |
| Separate loader for batch-corrected data | Batch-corrected files have a different structure; one function doing both adds fragile if/else logic |
| Clinical metadata stays in PRE_DIR | Batch correction only changes expression values; clinical is unchanged and does not need to be duplicated |
| `par.prior = TRUE` (default) in ComBat | Parametric empirical Bayes is recommended; non-parametric is rarely necessary unless prior fit is very poor |
