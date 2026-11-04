# libraries
library(survival)
library(randomForestSRC)
library(CoxBoost)
library(survivalsvm)

# directories & constants
ROOT_DIR <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))

DEG_DIR <- file.path(ROOT_DIR, "data/deg_rds_files")
BATCH_DIR <- file.path(ROOT_DIR, "data/batch_corrected_rds_files")
PRE_DIR <- file.path(ROOT_DIR, "data/preprocessed_rds_files")
MODEL_DIR <- file.path(ROOT_DIR, "data/ml_rds_files")
RESULTS_DIR <- file.path(ROOT_DIR, "data/ml_results")

TRAIN_DATASET <- "GSE42568"
TEST_DATASETS <- c("GSE21653", "GSE20711", "GSE88770")

# ── data loading ──────────────────────────────────────────────────────────────

load_deg_data <- function(deg_dir, dataset_name) {
    message("[", dataset_name, "] Loading DEG-selected expression + RFS data...")
    readRDS(file.path(deg_dir, paste0(dataset_name, "_deg_expression_matrix_with_rfs.rds")))
}

load_batch_data <- function(batch_dir, pre_dir, dataset_name) {
    message("[", dataset_name, "] Loading batch-corrected expression & clinical data...")
    expr     <- readRDS(file.path(batch_dir, paste0(dataset_name, "_expression_matrix.rds")))
    clinical <- readRDS(file.path(pre_dir, paste0(dataset_name, "_clinical_metadata.rds")))
    list(expr = expr, clinical = clinical)
}

# ── reshape for ml ────────────────────────────────────────────────────────────

build_expression_with_rfs <- function(expr, clinical) {
    expr_t <- as.data.frame(t(expr))
    expr_t$geo_accession <- rownames(expr_t)

    clinical_rfs <- clinical[, c("geo_accession", "rfs_time", "rfs_event")]
    clinical_rfs$rfs_event <- as.integer(clinical_rfs$rfs_event)

    merge(clinical_rfs, expr_t, by = "geo_accession")
}

split_x_y <- function(expr_with_rfs, gene_cols) {
    list(
        x = expr_with_rfs[, gene_cols, drop = FALSE],
        y = list(time = expr_with_rfs$rfs_time, status = expr_with_rfs$rfs_event)
    )
}

# Loads the DEG-selected training set, then restricts each test set (from the
# full batch-corrected matrix) to that same gene panel so every model sees
# identical features at train and test time.
prepare_ml_data <- function(deg_dir, batch_dir, pre_dir, train_dataset, test_datasets) {
    result <- list()

    train_expr_with_rfs <- load_deg_data(deg_dir, train_dataset)
    deg_genes <- setdiff(colnames(train_expr_with_rfs), c("geo_accession", "rfs_time", "rfs_event"))
    result[[train_dataset]] <- split_x_y(train_expr_with_rfs, deg_genes)
    message("[", train_dataset, "] Prepared - ", nrow(train_expr_with_rfs), " samples, ",
            length(deg_genes), " DEGs")

    for (dataset in test_datasets) {
        data          <- load_batch_data(batch_dir, pre_dir, dataset)
        expr_with_rfs <- build_expression_with_rfs(data$expr, data$clinical)

        missing_genes <- setdiff(deg_genes, colnames(expr_with_rfs))
        if (length(missing_genes) > 0) {
            stop("[", dataset, "] Missing ", length(missing_genes),
                 " DEG gene(s) from batch-corrected matrix: ",
                 paste(missing_genes, collapse = ", "))
        }

        result[[dataset]] <- split_x_y(expr_with_rfs, deg_genes)
        message("[", dataset, "] Prepared - ", nrow(expr_with_rfs), " samples, ",
                length(deg_genes), " genes (restricted to training DEGs)")
    }

    result
}

# ── ml training ───────────────────────────────────────────────────────────────

run_rsf <- function(x, y) {
    model_path <- file.path(MODEL_DIR, "rsf.rds")

    message("Training RSF (this may take a few minutes)...")
    train_data <- data.frame(time = y$time, status = y$status, x, check.names = FALSE)

    rsf <- rfsrc(
        Surv(time, status) ~ .,
        data = train_data
    )

    dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(rsf, model_path)
    message("Trained and saved.")

    rsf
}

run_coxboost <- function(x, y) {
    model_path <- file.path(MODEL_DIR, "coxboost.rds")

    message("Training CoxBoost (stepwise likelihood-based boosting)...")
    x_mat <- as.matrix(x)

    coxboost <- CoxBoost(
        time   = y$time,
        status = y$status,
        x      = x_mat
    )

    dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(coxboost, model_path)
    message("Trained and saved.")

    coxboost
}

run_ssvm <- function(x, y) {
    model_path <- file.path(MODEL_DIR, "ssvm.rds")

    message("Training SSVM...")
    train_data <- data.frame(time = y$time, status = y$status, x, check.names = FALSE)
    ssvm <- survivalsvm(Surv(time, status) ~ ., data = train_data, type = "regression", diff.meth = NULL, gamma.mu = 1)

    dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(ssvm, model_path)
    message("Trained and saved.")

    ssvm
}

# ── ml testing ────────────────────────────────────────────────────────────────

load_ml_models <- function(x, y) {
    rsf_path      <- file.path(MODEL_DIR, "rsf.rds")
    coxboost_path <- file.path(MODEL_DIR, "coxboost.rds")
    ssvm_path     <- file.path(MODEL_DIR, "ssvm.rds")

    if (file.exists(rsf_path)) {
        message("Loading saved RSF model...")
        rsf <- readRDS(rsf_path)
        message("Loaded.")
    } else {
        rsf <- run_rsf(x, y)
    }

    if (file.exists(coxboost_path)) {
        message("Loading saved CoxBoost model...")
        coxboost <- readRDS(coxboost_path)
        message("Loaded.")
    } else {
        coxboost <- run_coxboost(x, y)
    }

    if (file.exists(ssvm_path)) {
        message("Loading saved SSVM model...")
        ssvm <- readRDS(ssvm_path)
        message("Loaded.")
    } else {
        ssvm <- run_ssvm(x, y)
    }

    list(rsf = rsf, coxboost = coxboost, ssvm = ssvm)
}

predict_risk <- function(model, model_type, x) {
    switch(model_type,
        rsf = predict(model, newdata = x)$predicted,
        coxboost = as.numeric(predict(model, newdata = as.matrix(x), type = "lp")),
        ssvm     = as.numeric(predict(model, newdata = x)$predicted)
    )
}

test_ml_model <- function(model, model_type, datasets) {
    results <- list()

    if (model_type == "rsf") {
        # OOB C-index only exists for RSF since it's the only one bagged
        oob_cindex <- 1 - model$err.rate[length(model$err.rate)]
        message(sprintf("OOB C-index: %.5f", oob_cindex))
        results[["OOB"]] <- data.frame(model = model_type, dataset = "OOB", cindex = oob_cindex)
    }

    # rsf ("mortality") and coxboost ("lp") predict Cox-style risk scores,
    # where higher = shorter survival; concordance() only reverses that
    # convention automatically for coxph objects, so it must be told here.
    reverse <- model_type %in% c("rsf", "coxboost")

    for (name in names(datasets)) {
        risk <- predict_risk(model, model_type, datasets[[name]]$x)
        y <- datasets[[name]]$y
        cindex <- survival::concordance(Surv(y$time, y$status) ~ risk, reverse = reverse)$concordance
        label <- paste0(name, if (name == TRAIN_DATASET) " (train)" else " (test)")
        message(sprintf("%s C-index: %.5f", label, cindex))
        results[[name]] <- data.frame(model = model_type, dataset = label, cindex = cindex)
    }

    do.call(rbind, results)
}

# ── results export ────────────────────────────────────────────────────────────

export_ml_results <- function(results_df) {
    dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
    csv_path <- file.path(RESULTS_DIR, "ml_cindex_results.csv")
    write.csv(results_df, csv_path, row.names = FALSE)
    message("Saved: ", csv_path)

    dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
    timed_csv_path <- file.path(MODEL_DIR, "timed_validation_results.csv")
    write.csv(results_df, timed_csv_path, row.names = FALSE)
    message("Saved: ", timed_csv_path)
}

# ── main ──────────────────────────────────────────────────────────────────────

datasets  <- prepare_ml_data(DEG_DIR, BATCH_DIR, PRE_DIR, TRAIN_DATASET, TEST_DATASETS)
train_set <- datasets[[TRAIN_DATASET]]

models <- load_ml_models(train_set$x, train_set$y)

results_list <- list()
for (model_type in names(models)) {
    message("\n=== ", toupper(model_type), " ===")
    results_list[[model_type]] <- test_ml_model(models[[model_type]], model_type, datasets)
}

export_ml_results(do.call(rbind, results_list))