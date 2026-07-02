# libraries
library(survival)
library(randomForestSRC)
library(CoxBoost)
library(survivalsvm)

# directories & constants
ROOT_DIR <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))

BATCH_DIR <- file.path(ROOT_DIR, "data/batch_corrected_rds_files")
MODEL_DIR <- file.path(ROOT_DIR, "data/ml_rds_files")

DATASETS <- c("GSE42568", "GSE21653", "GSE20711", "GSE88770")

# ── data loading ──────────────────────────────────────────────────────────────

load_batch_data <- function(datapath, dataset_name) {
    message("[", dataset_name, "] Loading batch-corrected expression & clinical data...")
    expr     <- readRDS(file.path(datapath, paste0(dataset_name, "_expression_matrix.rds")))
    clinical <- readRDS(file.path(datapath, paste0(dataset_name, "_clinical_metadata.rds")))
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

split_x_y <- function(expr_with_rfs) {
    gene_cols <- setdiff(colnames(expr_with_rfs), c("geo_accession", "rfs_time", "rfs_event"))
    list(
        x = expr_with_rfs[, gene_cols],
        y = Surv(expr_with_rfs$rfs_time, expr_with_rfs$rfs_event)
    )
}

prepare_ml_data <- function(datapath) {
    result <- list()

    for (dataset in DATASETS) {
        data          <- load_batch_data(datapath, dataset)
        expr_with_rfs <- build_expression_with_rfs(data$expr, data$clinical)
        result[[dataset]] <- split_x_y(expr_with_rfs)
        message("[", dataset, "] Prepared - ", nrow(expr_with_rfs), " samples, ",
                ncol(result[[dataset]]$x), " genes")
    }

    result
}

# ── ml training ───────────────────────────────────────────────────────────────

run_rsf <- function(x, y) {
    model_path <- file.path(MODEL_DIR, "rsf.rds")

    message("Training RSF (this may take a few minutes)...")
    train_data <- data.frame(time = y[, 1], status = y[, 2], x)

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
        time   = y[, 1],
        status = y[, 2],
        x      = x_mat
    )

    dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(coxboost, model_path)
    message("Trained and saved.")

    coxboost
}

run_ssvm <- function(x, y) {
    model_path <- file.path(MODEL_DIR, "ssvm.rds")

    # type="vanbelle1": pure ranking-objective SVM, the R equivalent of
    # rank_ratio=1.0 — higher predicted value = higher risk
    message("Training SSVM...")
    train_data <- data.frame(time = y[, 1], status = y[, 2], x)
    ssvm       <- survivalsvm(Surv(time, status) ~ ., data = train_data, type = "vanbelle1")

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
        rsf      = predict(model, newdata = x)$predicted,
        coxboost = as.numeric(predict(model, newdata = as.matrix(x), type = "lp")),
        # type="vanbelle1" is a ranking SVM, so like rsf/coxboost higher
        # predicted value = higher risk — no sign flip needed
        ssvm     = as.numeric(predict(model, newdata = x)$predicted)
    )
}

test_ml_model <- function(model, model_type, datasets) {
    if (model_type == "rsf") {
        # OOB C-index only exists for RSF since it's the only one bagged
        oob_cindex <- 1 - model$err.rate[length(model$err.rate)]
        message(sprintf("OOB C-index: %.5f", oob_cindex))
    }

    for (name in names(datasets)) {
        risk   <- predict_risk(model, model_type, datasets[[name]]$x)
        cindex <- survival::concordance(datasets[[name]]$y ~ risk)$concordance
        message(sprintf("%s C-index: %.5f", name, cindex))
    }
}

# ── main ──────────────────────────────────────────────────────────────────────
