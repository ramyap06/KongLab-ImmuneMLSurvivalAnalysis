library(limma)
library(ggplot2)
library(pheatmap)

ROOT_DIR <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))

PRE_DIR          <- file.path(ROOT_DIR, "data/preprocessed_rds_files")
DEG_DIR          <- file.path(ROOT_DIR, "data/deg_rds_files")
IMMUNE_GENES_DIR <- file.path(ROOT_DIR, "data/immune_genes_and_tf")
VIS_DIR          <- file.path(ROOT_DIR, "results/deg")
TF_PATH          <- file.path(IMMUNE_GENES_DIR, "transcription_factors.csv")

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
    saveRDS(deg_genes, file.path(DEG_DIR, paste0(dataset_name, "_deg_stats.rds")))
    message("Saved: ", dataset_name, "_deg_stats.rds (", nrow(deg_genes), " DEGs before gene-set filtering)")

    list(all_genes = all_genes, deg_genes = deg_genes, deg_gene_symbols = deg_genes$gene_symbol)
}

# ── visualisation ─────────────────────────────────────────────────────────────

plot_volcano <- function(all_genes, signature_genes, lfc_cutoff, p_cutoff, out_path) {
    message("Plotting volcano plot...")

    volcano_df <- all_genes
    volcano_df$neg_log_p <- -log10(volcano_df$P.Value)

    # colour by direction of change, but only once a gene is FDR-significant
    volcano_df$category <- "Not Significant"
    volcano_df$category[volcano_df$adj.P.Val < p_cutoff & volcano_df$logFC >  lfc_cutoff] <- "Up-regulated"
    volcano_df$category[volcano_df$adj.P.Val < p_cutoff & volcano_df$logFC < -lfc_cutoff] <- "Down-regulated"
    volcano_df$category <- factor(volcano_df$category,
                                  levels = c("Up-regulated", "Down-regulated", "Not Significant"))

    p <- ggplot(volcano_df, aes(x = logFC, y = neg_log_p, color = category)) +
        geom_point(alpha = 0.45, size = 1.2) +
        geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
        geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
        scale_color_manual(values = c(
            "Up-regulated"    = "#E64B35",
            "Down-regulated"  = "#4DBBD5",
            "Not Significant" = "grey75"
        )) +
        labs(
            title    = "Volcano Plot: Tumor vs Normal",
            subtitle = sprintf("DEGs: %d up  |  %d down  |  signature genes highlighted",
                               sum(volcano_df$category == "Up-regulated"),
                               sum(volcano_df$category == "Down-regulated")),
            x        = "log₂ Fold Change",
            y        = "-log₁₀(P-value)",
            color    = NULL
        ) +
        theme_classic(base_size = 13) +
        theme(
            plot.title    = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
            legend.position = "top"
        )

    ggsave(out_path, p, width = 8, height = 6, dpi = 150)
    message("Saved: ", out_path)
    invisible(p)
}

plot_deg_heatmap <- function(expr, deg_genes, clinical, signature_genes, n_top, out_path) {
    message("Plotting top ", n_top, " DEG heatmap...")

    # take the most significant DEGs by adjusted p-value, then make sure the
    # curated immune signature genes are always represented on the heatmap
    top_genes     <- head(deg_genes[order(deg_genes$adj.P.Val), ], n_top)$gene_symbol
    heatmap_genes <- unique(c(top_genes, signature_genes))
    heatmap_genes <- heatmap_genes[heatmap_genes %in% expr$gene_symbol]

    # subset and reshape to a numeric matrix: genes as rows, samples as columns
    heatmap_expr <- expr[expr$gene_symbol %in% heatmap_genes, ]
    heatmap_mat  <- as.matrix(heatmap_expr[, -1])
    rownames(heatmap_mat) <- heatmap_expr$gene_symbol

    # column annotation: which samples are tumor vs normal
    col_annot <- data.frame(
        Tissue = ifelse(clinical$is_tumor == 1, "Tumor", "Normal"),
        row.names = clinical$geo_accession
    )

    # row annotation: flag which genes belong to the curated immune signature
    row_annot <- data.frame(
        Type = ifelse(rownames(heatmap_mat) %in% signature_genes, "Immune Signature", "Top DEG"),
        row.names = rownames(heatmap_mat)
    )

    annot_colors <- list(
        Tissue = c(Tumor = "#E64B35", Normal = "#4DBBD5"),
        Type   = c("Immune Signature" = "#FF9500", "Top DEG" = "grey80")
    )

    # silent = TRUE: build the gtable without popping up an interactive device,
    # since this runs as a script rather than a notebook
    p_heat <- pheatmap(
        heatmap_mat,
        scale            = "row",
        color            = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
        annotation_col   = col_annot,
        annotation_row   = row_annot,
        annotation_colors = annot_colors,
        show_colnames    = FALSE,
        fontsize_row     = 7,
        border_color     = NA,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        main             = "Top DEGs Heatmap — Tumor vs Normal\n(immune signature genes highlighted)",
        silent           = TRUE
    )

    png(out_path, width = 10 * 150, height = 10 * 150, res = 150)
    grid::grid.newpage()
    grid::grid.draw(p_heat$gtable)
    dev.off()
    message("Saved: ", out_path)
    invisible(p_heat)
}

# ── gene-set annotation ───────────────────────────────────────────────────────

load_transcription_factors <- function(tf_path) {
    message("Loading tumor-related transcription factors...")
    tf_df <- read.csv(tf_path)
    tf_df$Target_Gene
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

export_deg_outputs <- function(tumor_tf_deg_genes, expr_with_rfs, dataset_name) {
    message("[", dataset_name, "] Exporting DEG outputs...")
    dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)

    saveRDS(tumor_tf_deg_genes, file.path(DEG_DIR, paste0(dataset_name, "_deg_tf_filtered.rds")))
    saveRDS(expr_with_rfs,      file.path(DEG_DIR, paste0(dataset_name, "_deg_expression_matrix_with_rfs.rds")))
}

# ── main ──────────────────────────────────────────────────────────────────────

dir.create(VIS_DIR, recursive = TRUE, showWarnings = FALSE)

data     <- load_expression_data(DATASET)
expr     <- data$expr
clinical <- data$clinical

design_contrast <- build_design_and_contrast(clinical)
fit_contrast    <- fit_limma_model(expr, design_contrast$design, design_contrast$contrast)

degs <- extract_degs(fit_contrast, P_CUTOFF, ADJ_P_CUTOFF, LFC_CUTOFF, DATASET)

plot_volcano(degs$all_genes, SIGNATURE_GENES, LFC_CUTOFF, P_CUTOFF,
             file.path(VIS_DIR, "volcano_plot.png"))
plot_deg_heatmap(expr, degs$deg_genes, clinical, SIGNATURE_GENES, n_top = 40,
                 file.path(VIS_DIR, "deg_heatmap.png"))

# narrow the DEG list down to biologically relevant subsets: genes with known
# immune roles, and genes that are tumor-related transcription factors
immune_genes <- readRDS(file.path(IMMUNE_GENES_DIR, "immune_genes.rds"))
tumor_tf     <- load_transcription_factors(TF_PATH)

immune_deg_genes   <- filter_degs_by_gene_set(degs$deg_genes, immune_genes)
tumor_tf_deg_genes <- filter_degs_by_gene_set(degs$deg_genes, tumor_tf)

# build the expression matrix for the immune DEGs and transpose it for
# downstream survival modelling (samples as rows, genes as columns)
expr_transposed <- build_transposed_expression(expr, immune_deg_genes$gene_symbol)

# merge with clinical RFS columns and restrict to tumor samples with valid RFS
expr_with_rfs <- build_expression_with_rfs(expr_transposed, clinical)

export_deg_outputs(tumor_tf_deg_genes, expr_with_rfs, DATASET)
message("[", DATASET, "] DEG analysis done.")