library(limma)

ROOT_DIR <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))

PRE_DIR          <- file.path(ROOT_DIR, "data/preprocessed_rds_files")
DEG_DIR          <- file.path(ROOT_DIR, "data/deg_rds_files")
IMMUNE_GENES_DIR <- file.path(ROOT_DIR, "data/immune_genes_and_tf")
TF_PATH          <- file.path(IMMUNE_GENES_DIR, "transcription_factors.csv")
LEGACY_IMMUNE_PATH <- file.path(IMMUNE_GENES_DIR, "legacy_immune_genes_ImmPort.txt")

# defining dataset, only using training dataset
DATASET <- "GSE42568"

# curated immune-signature genes
SIGNATURE_GENES <- c("TSLP", "BIRC5", "S100B", "MDK", "S100P", "RARRES3", "BLNK", "ACO1")

# thresholds that define "differentially expressed"
P_CUTOFF     <- 0.05
ADJ_P_CUTOFF <- 0.05
LFC_CUTOFF   <- 1

# ── data loading ──────────────────────────────────────────────────────────────

load_expression_data <- function(dataset_name) {
    message("[", dataset_name, "] Loading expression matrix & clinical metadata...")
    expr     <- readRDS(file.path(PRE_DIR, paste0(dataset_name, "_expression_matrix.rds")))
    clinical <- readRDS(file.path(PRE_DIR, paste0(dataset_name, "_clinical_metadata.rds")))
    expr <- expr[, c("gene_symbol", clinical$geo_accession)]
    list(expr = expr, clinical = clinical)
}

# ── limma model ───────────────────────────────────────────────────────────────

build_design_and_contrast <- function(clinical) {
    message("Building design & contrast matrices (tumor vs normal)...")

    # one-hot encode tissue type & label columns from the factor's own levels
    groups <- factor(clinical$is_tumor, levels = c(0, 1), labels = c("normal", "tumor"))
    design <- model.matrix(~ 0 + groups)
    colnames(design) <- levels(groups)

    # differences in expression from tumor to normal
    contrast <- makeContrasts(tumor_vs_normal = tumor - normal, levels = design)

    list(design = design, contrast = contrast)
}

fit_limma_model <- function(expr, design, contrast) {
    message("Fitting linear model and applying contrast...")

    # per-gene linear model (generalisation of the t-test)
    fit <- lmFit(expr, design)

    # collapse the fitted coefficients down to the tumor-vs-normal difference
    fit_contrast <- contrasts.fit(fit, contrast)

    # eBayes: empirical-Bayes moderation, shrinks each gene's variance
    # estimate toward a common prior, which stabilises statistics for a
    # cohort this small and yields the moderated t / B statistics
    eBayes(fit_contrast)
}

# ── DEG extraction ────────────────────────────────────────────────────────────

extract_degs <- function(fit_contrast, p_cutoff, adj_p_cutoff, lfc_cutoff, dataset_name) {
    message("Extracting differentially expressed genes...")

    # rank every gene with Benjamini-Hochberg FDR correction applied
    all_genes <- topTable(fit_contrast, number = Inf, adjust = "BH")

    # must have nominal significance, FDR-corrected significance,
    # & biologically meaningful effect size (>= 2-fold change) to be DEG
    deg_genes <- all_genes[
        all_genes$P.Value < p_cutoff &
            all_genes$adj.P.Val < adj_p_cutoff &
            abs(all_genes$logFC) > lfc_cutoff,
    ]

    dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(deg_genes, file.path(DEG_DIR, paste0(dataset_name, "_deg_stats_all.rds")))
    message("Saved: ", dataset_name, "_deg_stats.rds (", nrow(deg_genes), " DEGs before gene-set filtering)")

    list(all_genes = all_genes, deg_genes = deg_genes, deg_gene_symbols = deg_genes$gene_symbol)
}

# ── gene-set annotation ───────────────────────────────────────────────────────

# separate function because .csv
load_transcription_factors <- function(tf_path) {
    message("Loading tumor-related transcription factors...")
    tf_df <- read.csv(tf_path)
    tf_df$Transcription_Factor
}

# separate function because line by line read
load_legacy_immune_genes <- function(path) {
    message("Loading legacy ImmPort immune gene list...")
    readLines(path)
}

filter_degs_by_gene_set <- function(deg_genes, gene_set) {
    deg_genes[deg_genes$gene_symbol %in% intersect(deg_genes$gene_symbol, gene_set), ]
}

# ── export prep ───────────────────────────────────────────────────────────────

build_transposed_expression <- function(expr, deg_gene_symbols) {
    # filter to DEG genes then transpose so genes become columns / samples
    # become rows, ready to be joined with clinical data for downstream
    # Cox/survival modelling
    filtered <- expr[expr$gene_symbol %in% deg_gene_symbols, ]
    expr_t <- as.data.frame(t(filtered))
    colnames(expr_t) <- filtered$gene_symbol
    expr_t[-1, ]
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

export_deg_outputs <- function(tumor_tf_deg_genes, expr_with_rfs, expr_with_rfs_legacy_immune, dataset_name) {
    message("[", dataset_name, "] Exporting DEG outputs...")
    dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)

    saveRDS(tumor_tf_deg_genes,         file.path(DEG_DIR, paste0(dataset_name, "_deg_tf_filtered.rds")))
    saveRDS(expr_with_rfs,              file.path(DEG_DIR, paste0(dataset_name, "_deg_expression_matrix_with_rfs.rds")))
    saveRDS(expr_with_rfs_legacy_immune, file.path(DEG_DIR, paste0(dataset_name, "_deg_expression_matrix_with_rfs_legacy_immune.rds")))
}

# ── main ──────────────────────────────────────────────────────────────────────

data     <- load_expression_data(DATASET)
expr     <- data$expr
clinical <- data$clinical

design_contrast <- build_design_and_contrast(clinical)
fit_contrast    <- fit_limma_model(expr, design_contrast$design, design_contrast$contrast)

degs <- extract_degs(fit_contrast, P_CUTOFF, ADJ_P_CUTOFF, LFC_CUTOFF, DATASET)

# narrow the DEG list down to biologically relevant subsets: genes with known
# immune roles, and genes that are tumor-related transcription factors
immune_genes        <- readRDS(file.path(IMMUNE_GENES_DIR, "immune_genes.rds"))
tumor_tf            <- load_transcription_factors(TF_PATH)
legacy_immune_genes <- load_legacy_immune_genes(LEGACY_IMMUNE_PATH)

immune_deg_genes        <- filter_degs_by_gene_set(degs$deg_genes, immune_genes)
tumor_tf_deg_genes      <- filter_degs_by_gene_set(degs$deg_genes, tumor_tf)
legacy_immune_deg_genes <- filter_degs_by_gene_set(degs$deg_genes, legacy_immune_genes)

# build expression matrices for each gene subset and transpose for downstream
# survival modelling (samples as rows, genes as columns)
expr_transposed                <- build_transposed_expression(expr, immune_deg_genes$gene_symbol)
expr_transposed_legacy_immune  <- build_transposed_expression(expr, legacy_immune_deg_genes$gene_symbol)

# merge with clinical RFS columns and restrict to tumor samples with valid RFS
expr_with_rfs                <- build_expression_with_rfs(expr_transposed, clinical)
expr_with_rfs_legacy_immune  <- build_expression_with_rfs(expr_transposed_legacy_immune, clinical)

export_deg_outputs(tumor_tf_deg_genes, expr_with_rfs, expr_with_rfs_legacy_immune, DATASET)
message("[", DATASET, "] DEG analysis done.")