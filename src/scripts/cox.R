library(survival)
library(glmnet)
library(ggplot2)
library(survminer)

ROOT_DIR <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))

DEG_DIR <- file.path(ROOT_DIR, "data/deg_rds_files")
PRE_DIR <- file.path(ROOT_DIR, "data/preprocessed_rds_files")
COX_DIR <- file.path(ROOT_DIR, "data/cox_rds_files")
VIS_DIR <- file.path(ROOT_DIR, "results/cox")

DATASET         <- "GSE42568"
CLINICAL_COVARS <- c("er_status", "tumor_grade", "lymph_node_status", "t_stage")
SURV_COLS       <- c("geo_accession", "rfs_time", "rfs_event")

P_CUTOFF    <- 0.05
LASSO_ALPHA <- 1.0
ENET_ALPHA  <- 0.9
N_FOLDS     <- 5

set.seed(42)

# ── data loading ──────────────────────────────────────────────────────────────

load_cox_data <- function(dataset_name) {
    message("[", dataset_name, "] Loading DEG expression matrix & clinical metadata...")
    expr_rfs <- readRDS(file.path(DEG_DIR, paste0(dataset_name, "_deg_expression_matrix_with_rfs.rds")))
    clinical <- readRDS(file.path(PRE_DIR, paste0(dataset_name, "_clinical_metadata.rds")))
    list(expr_rfs = expr_rfs, clinical = clinical)
}

# ── univariate cox ────────────────────────────────────────────────────────────

run_univariate_cox <- function(data, gene_cols,
                               time_col  = "rfs_time",
                               event_col = "rfs_event") {
    message("Running univariate Cox regression for ", length(gene_cols), " genes...")
    results <- vector("list", length(gene_cols))

    for (i in seq_along(gene_cols)) {
        gene <- gene_cols[i]
        results[[i]] <- tryCatch({
            sub <- data[complete.cases(data[, c(time_col, event_col, gene)]),
                        c(time_col, event_col, gene)]
            fit <- coxph(
                as.formula(sprintf("Surv(%s, %s) ~ `%s`", time_col, event_col, gene)),
                data = sub
            )
            s <- summary(fit)
            data.frame(
                gene        = gene,
                beta        = round(s$coefficients[1, "coef"],     4),
                HR          = round(s$conf.int[1,   "exp(coef)"],  4),
                HR_lower_95 = round(s$conf.int[1,   "lower .95"],  4),
                HR_upper_95 = round(s$conf.int[1,   "upper .95"],  4),
                p_value     = round(s$coefficients[1, "Pr(>|z|)"], 4),
                stringsAsFactors = FALSE
            )
        }, error = function(e) {
            message("  Skipped ", gene, ": ", conditionMessage(e))
            NULL
        })
    }

    all_results <- do.call(rbind, Filter(Negate(is.null), results))
    all_results <- all_results[order(all_results$p_value), ]

    sig_results <- all_results[all_results$p_value < P_CUTOFF, ]
    sig_results$role <- ifelse(sig_results$HR > 1, "danger", "protective")

    message("  Tested: ", nrow(all_results),
            " | Significant (p < ", P_CUTOFF, "): ", nrow(sig_results))
    list(all_results = all_results, sig_results = sig_results)
}

# ── penalized cox ─────────────────────────────────────────────────────────────

run_penalized_cox <- function(X_train, y_train, alpha, nfolds, label) {
    message("Fitting ", label, " Cox (alpha = ", alpha, ", ", nfolds, "-fold CV)...")
    cv_fit     <- cv.glmnet(X_train, y_train, family = "cox",
                            alpha = alpha, nfolds = nfolds, standardize = TRUE)
    best_coefs <- as.vector(coef(cv_fit, s = "lambda.min"))
    names(best_coefs) <- colnames(X_train)
    selected   <- names(best_coefs[best_coefs != 0])
    cindex     <- max(cv_fit$cvm)
    message("  Non-zero genes: ", length(selected), " | CV C-index: ", round(cindex, 4))
    list(cv_fit = cv_fit, fit = cv_fit$glmnet.fit, selected_genes = selected, best_cindex = cindex)
}

# ── multivariate cox ──────────────────────────────────────────────────────────

run_multivariate_cox <- function(train_data, gene_cols, clinical_cols) {
    message("Running multivariate Cox (",
            length(gene_cols), " genes + ", length(clinical_cols), " clinical variables)...")

    all_covars  <- c(gene_cols, clinical_cols)
    formula_str <- sprintf(
        "Surv(rfs_time, rfs_event) ~ %s",
        paste(sprintf("`%s`", all_covars), collapse = " + ")
    )
    cph <- coxph(as.formula(formula_str), data = train_data)
    s   <- summary(cph)

    all_results <- data.frame(
        gene        = rownames(s$coefficients),
        beta        = round(s$coefficients[, "coef"],     4),
        HR          = round(s$conf.int[,   "exp(coef)"],  4),
        HR_lower_95 = round(s$conf.int[,   "lower .95"],  4),
        HR_upper_95 = round(s$conf.int[,   "upper .95"],  4),
        p_value     = round(s$coefficients[, "Pr(>|z|)"], 4),
        stringsAsFactors = FALSE, row.names = NULL
    )
    all_results <- all_results[order(all_results$p_value), ]

    sig_results <- all_results[all_results$p_value < P_CUTOFF & all_results$gene %in% gene_cols, ]
    sig_results$role <- ifelse(sig_results$HR > 1, "danger", "protective")

    message("  Significant genes (p < ", P_CUTOFF, "): ", nrow(sig_results))
    list(all_results = all_results, sig_results = sig_results)
}

# ── risk scoring ──────────────────────────────────────────────────────────────

compute_risk_scores <- function(df, cox_results, sig_genes) {
    beta_vec <- setNames(cox_results$beta, cox_results$gene)
    common   <- intersect(sig_genes, intersect(names(beta_vec), names(df)))
    risk     <- as.vector(as.matrix(df[, common, drop = FALSE]) %*% beta_vec[common])
    out      <- df[, c(SURV_COLS, common)]
    out$risk_score <- risk
    message("  Risk scores: ", nrow(out), " patients | range [",
            round(min(risk), 4), ", ", round(max(risk), 4), "]")
    out
}

# ── visualisations ────────────────────────────────────────────────────────────

plot_forest <- function(df, title) {
    df        <- df[order(df$HR), ]
    df$gene   <- factor(df$gene, levels = df$gene)
    pal       <- c(danger = "#E64B35", protective = "#8491B4")
    ggplot(df, aes(x = HR, y = gene, colour = role)) +
        geom_vline(xintercept = 1, linetype = "dashed", colour = "black", alpha = 0.6) +
        geom_errorbarh(aes(xmin = HR_lower_95, xmax = HR_upper_95),
                       height = 0.35, linewidth = 1.1) +
        geom_point(size = 3) +
        geom_text(aes(x = HR_upper_95, label = sprintf("p=%.3f", p_value)),
                  hjust = -0.1, size = 2.7, colour = "grey45") +
        scale_colour_manual(values = pal, name = "",
            labels = c(danger = "Danger gene (HR > 1)", protective = "Protective gene (HR < 1)")) +
        labs(title = title, x = "Hazard Ratio (95% CI)", y = NULL) +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
}

plot_cv_lambda <- function(cv_fit, title) {
    df <- data.frame(
        log_lambda = log(cv_fit$lambda),
        mean_cv    = cv_fit$cvm,
        upper      = cv_fit$cvup,
        lower      = cv_fit$cvlo
    )
    ggplot(df, aes(x = log_lambda, y = mean_cv)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, fill = "#4DBBD5") +
        geom_line(colour = "#4DBBD5") +
        geom_vline(xintercept = log(cv_fit$lambda.min),
                   colour = "#E64B35", linetype = "dashed") +
        geom_hline(yintercept = 0.5, colour = "grey50", linetype = "dashed") +
        labs(title = title, x = expression(log(lambda)), y = "Concordance index (CV)") +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(face = "bold"))
}

plot_coef_path <- function(fit, n_highlight = 10, title) {
    mat <- as.matrix(fit$beta)
    df  <- data.frame(
        gene       = rep(rownames(mat), ncol(mat)),
        lambda_idx = rep(seq_len(ncol(mat)), each = nrow(mat)),
        coef       = as.vector(mat)
    )
    df$log_lambda <- log(fit$lambda[df$lambda_idx])
    top_genes <- names(sort(abs(mat[, ncol(mat)]), decreasing = TRUE))[
        seq_len(min(n_highlight, nrow(mat)))]
    df_top  <- df[df$gene %in% top_genes, ]
    df_rest <- df[!df$gene %in% top_genes, ]
    ggplot() +
        geom_line(data = df_rest, aes(x = log_lambda, y = coef, group = gene),
                  colour = "grey80", alpha = 0.5) +
        geom_line(data = df_top, aes(x = log_lambda, y = coef, group = gene, colour = gene)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
        labs(title = title, x = expression(log(lambda)), y = "Coefficient", colour = "Gene") +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

plot_km <- function(scores_df, title, cutoff) {
    df_km <- data.frame(
        time  = as.numeric(scores_df$rfs_time),
        event = as.integer(scores_df$rfs_event),
        group = factor(ifelse(scores_df$risk_score >= cutoff, "High Risk", "Low Risk"),
                       levels = c("Low Risk", "High Risk"))
    )
    df_km <- df_km[complete.cases(df_km), ]
    fit   <- survfit(Surv(time, event) ~ group, data = df_km)
    ggsurvplot(
        fit, data = df_km, pval = TRUE, conf.int = TRUE,
        palette     = c("Low Risk" = "#4DBBD5", "High Risk" = "#E64B35"),
        title       = title, xlab = "Time (days)", ylab = "Relapse-Free Survival",
        legend.labs = c("Low Risk", "High Risk"),
        ggtheme     = theme_bw(base_size = 9)
    )$plot + theme(plot.title = element_text(face = "bold", size = 10))
}

# ── main ──────────────────────────────────────────────────────────────────────

dir.create(COX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(VIS_DIR, recursive = TRUE, showWarnings = FALSE)

data     <- load_cox_data(DATASET)
expr_rfs <- data$expr_rfs
clinical <- data$clinical

gene_cols <- setdiff(names(expr_rfs), SURV_COLS)
message(sprintf("[%s] %d patients | %d gene columns", DATASET, nrow(expr_rfs), length(gene_cols)))

# ── step 1: univariate cox ────────────────────────────────────────────────────

univ      <- run_univariate_cox(expr_rfs, gene_cols)
univ_sig  <- univ$sig_results$gene
univ_expr <- expr_rfs[, c(SURV_COLS, univ_sig)]

saveRDS(univ_expr, file.path(COX_DIR, paste0(DATASET, "_univariate_cox.rds")))
message("Saved: ", DATASET, "_univariate_cox.rds (",
        nrow(univ_expr), " rows, ", ncol(univ_expr), " cols)")

if (nrow(univ$sig_results) > 0) {
    ggsave(file.path(VIS_DIR, "univariate_forest_plot.png"),
           plot_forest(univ$sig_results, "Univariate Cox: Hazard Ratios (p < 0.05)"),
           dpi = 150, width = 9, height = max(4, nrow(univ$sig_results) * 0.3))
    message("Saved: univariate_forest_plot.png")
}

# ── step 2: penalized cox ─────────────────────────────────────────────────────

X_train <- as.matrix(univ_expr[, univ_sig])
y_train <- Surv(as.numeric(univ_expr$rfs_time), as.integer(univ_expr$rfs_event))

lasso_result <- run_penalized_cox(X_train, y_train, LASSO_ALPHA, N_FOLDS, "LASSO")
enet_result  <- run_penalized_cox(X_train, y_train, ENET_ALPHA,  N_FOLDS, "Elastic Net")

best_penalized <- if (enet_result$best_cindex >= lasso_result$best_cindex) enet_result else lasso_result
best_label     <- if (identical(best_penalized, enet_result)) "Elastic Net" else "LASSO"
message("Best penalized model: ", best_label,
        " (C-index = ", round(best_penalized$best_cindex, 4), ")")

ggsave(file.path(VIS_DIR, "lasso_coef_path.png"),
       plot_coef_path(lasso_result$fit, title = "LASSO Cox — Coefficient Path"),
       dpi = 150, width = 9, height = 6)
ggsave(file.path(VIS_DIR, "lasso_cv_lambda.png"),
       plot_cv_lambda(lasso_result$cv_fit, "LASSO Cox — CV Concordance vs Lambda"),
       dpi = 150, width = 9, height = 6)
ggsave(file.path(VIS_DIR, "enet_coef_path.png"),
       plot_coef_path(enet_result$fit, title = "Elastic Net Cox — Coefficient Path"),
       dpi = 150, width = 9, height = 6)
ggsave(file.path(VIS_DIR, "enet_cv_lambda.png"),
       plot_cv_lambda(enet_result$cv_fit, "Elastic Net Cox — CV Concordance vs Lambda"),
       dpi = 150, width = 9, height = 6)
message("Saved: penalized Cox plots")

penal_genes <- best_penalized$selected_genes
penal_expr  <- expr_rfs[, c(SURV_COLS, penal_genes)]

saveRDS(penal_expr, file.path(COX_DIR, paste0(DATASET, "_penalized_cox.rds")))
message("Saved: ", DATASET, "_penalized_cox.rds (",
        nrow(penal_expr), " rows, ", ncol(penal_expr), " cols)")

# ── step 3: multivariate cox ──────────────────────────────────────────────────

# merge penalized expression with clinical covariates for joint model
clinical_sel               <- clinical[, c("geo_accession", CLINICAL_COVARS)]
train_data                 <- merge(penal_expr, clinical_sel, by = "geo_accession")
train_data[CLINICAL_COVARS] <- lapply(train_data[CLINICAL_COVARS],
                                      function(x) as.numeric(as.character(x)))
train_data                 <- train_data[complete.cases(train_data), ]
train_data$rfs_event       <- as.integer(train_data$rfs_event)

multi           <- run_multivariate_cox(train_data, penal_genes, CLINICAL_COVARS)
multi_sig_genes <- if (nrow(multi$sig_results) > 0) multi$sig_results$gene else penal_genes

if (nrow(multi$sig_results) > 0) {
    ggsave(file.path(VIS_DIR, "multivariate_forest_plot.png"),
           plot_forest(multi$sig_results, "Multivariate Cox: Significant Genes (p < 0.05)"),
           dpi = 150, width = 9, height = max(4, nrow(multi$sig_results) * 0.3))
    message("Saved: multivariate_forest_plot.png")
}

# score all patients from penal_expr using betas from the multivariate model
# (not restricted to those with complete clinical data)
scored            <- compute_risk_scores(penal_expr, multi$all_results, multi_sig_genes)
train_median      <- median(scored$risk_score)
scored$risk_group <- ifelse(scored$risk_score >= train_median, "High Risk", "Low Risk")

ggsave(file.path(VIS_DIR, "multivariate_km_stratification.png"),
       plot_km(scored, sprintf("%s: High vs Low Risk", DATASET), train_median),
       dpi = 150, width = 6, height = 5)
message("Saved: multivariate_km_stratification.png")

saveRDS(scored, file.path(COX_DIR, paste0(DATASET, "_multivariate_cox.rds")))
message("Saved: ", DATASET, "_multivariate_cox.rds (",
        nrow(scored), " rows, ", ncol(scored), " cols)")

message("[", DATASET, "] Cox analysis done.")
