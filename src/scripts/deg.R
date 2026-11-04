library(survival)

ROOT_DIR <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))

DEG_DIR <- file.path(ROOT_DIR, "data/deg_rds_files")
BATCH_DIR <- file.path(ROOT_DIR, "data/batch_corrected_rds_files")
PRE_DIR <- file.path(ROOT_DIR, "data/preprocessed_rds_files")

# defining dataset, only using training dataset
DATASET <- "GSE42568"

# thresholds that define "differentially expressed"
ADJ_P_CUTOFF <- 0.05

# ── data loading ──────────────────────────────────────────────────────────────

load_data <- function(dataset_name) {
    message("[", dataset_name, "] Loading expression matrix & clinical metadata...")
    expr     <- readRDS(file.path(BATCH_DIR, paste0(dataset_name, "_expression_matrix.rds")))
    clinical <- readRDS(file.path(PRE_DIR, paste0(dataset_name, "_clinical_metadata.rds")))
    expr <- expr[, clinical$geo_accession, drop = FALSE]
    list(expr = expr, clinical = clinical)
}

# ── preparing dataset ──────────────────────────────────────────────────────────────

build_transposed_expression <- function(expr) {
    # transpose so genes become columns / samples become rows, ready to be
    # joined with clinical data for downstream Cox/survival modelling
    message("Transposing expression matrix (genes -> columns, samples -> rows)...")
    as.data.frame(t(expr))
}

build_expression_with_rfs <- function(expr_transposed, clinical) {
    message("Joining expression matrix with RFS clinical data...")

    tumor_clinical <- clinical[
        clinical$is_tumor == 1 &
            !is.na(clinical$rfs_time) &
            !is.na(clinical$rfs_event),
        c("geo_accession", "rfs_time", "rfs_event")
    ]
    tumor_clinical$rfs_event <- as.integer(tumor_clinical$rfs_event)

    expr_df <- data.frame(
        geo_accession = rownames(expr_transposed),
        expr_transposed,
        check.names = FALSE
    )

    # t() coerces a mixed data frame to character — restore numerics
    gene_cols <- setdiff(colnames(expr_df), "geo_accession")
    expr_df[gene_cols] <- lapply(expr_df[gene_cols], as.numeric)

    merged <- merge(tumor_clinical, expr_df, by = "geo_accession")

    merged
}

# ── univariate cox model + DEG extract ───────────────────────────────────────────────────────────────

extract_cox_degs <- function(expr_rfs, dataset_name, adj_p_cutoff) {
    message("[", dataset_name, "] Fitting univariate Cox models for ",
            ncol(expr_rfs) - 3, " genes...")

    gene_cols <- setdiff(colnames(expr_rfs), c("geo_accession", "rfs_time", "rfs_event"))
    cox_results <- numeric(length(gene_cols))
    names(cox_results) <- gene_cols

    for (gene in gene_cols) {
        fit <- coxph(Surv(expr_rfs$rfs_time, expr_rfs$rfs_event) ~ expr_rfs[[gene]])
        cox_results[gene] <- summary(fit)$coefficients[1, "Pr(>|z|)"]
    }

    adj_p <- p.adjust(cox_results, method = "BH")

    degs <- data.frame(gene_symbol = gene_cols, p_value = cox_results, adj_p_value = adj_p) |>
        subset(adj_p_value < adj_p_cutoff)

    message("[", dataset_name, "] ", nrow(degs), " / ", length(gene_cols),
            " genes pass adjusted p < ", adj_p_cutoff, ".")

    degs
}

# ── export ───────────────────────────────────────────────────────────────

export_deg_outputs <- function(expr_with_rfs, dataset_name) {
    message("[", dataset_name, "] Exporting DEG outputs...")
    dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)

    saveRDS(expr_with_rfs, file.path(DEG_DIR, paste0(dataset_name, "_deg_expression_matrix_with_rfs.rds")))
}

# ── main ──────────────────────────────────────────────────────────────────────

data     <- load_data(DATASET)
expr     <- data$expr
clinical <- data$clinical

# transpose the full expression matrix for downstream survival modelling
# (samples as rows, genes as columns)
expr_transposed <- build_transposed_expression(expr)

# merge with clinical RFS columns and restrict to tumor samples with valid RFS
expr_rfs <- build_expression_with_rfs(expr_transposed, clinical)

# fit a univariate Cox model per gene and keep only genes passing the
# BH-adjusted p-value cutoff
degs <- extract_cox_degs(expr_rfs, DATASET, ADJ_P_CUTOFF)

# restrict the merged data down to sample id + RFS columns + DEG genes only
deg_expr_rfs <- expr_rfs[, c("geo_accession", "rfs_time", "rfs_event", degs$gene_symbol)]

export_deg_outputs(deg_expr_rfs, DATASET)
message("[", DATASET, "] DEG analysis done.")