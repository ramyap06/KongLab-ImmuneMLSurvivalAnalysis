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
N_FOLDS     <- 10

# ── data loading ──────────────────────────────────────────────────────────────

load_cox_data <- function(dataset_name) {}

# ── univariate cox ────────────────────────────────────────────────────────────

run_univariate_cox <- function(data, gene_cols,
                               time_col  = "rfs_time",
                               event_col = "rfs_event") {}

# ── penalized cox ─────────────────────────────────────────────────────────────

run_penalized_cox <- function(X_train, y_train, alpha, nfolds, label) {}

# ── multivariate cox ──────────────────────────────────────────────────────────

run_multivariate_cox <- function(train_data, gene_cols, clinical_cols) {}

# ── risk scoring ──────────────────────────────────────────────────────────────

compute_risk_scores <- function(df, cox_results, sig_genes) {}

# ── main ──────────────────────────────────────────────────────────────────────

dir.create(COX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(VIS_DIR, recursive = TRUE, showWarnings = FALSE)