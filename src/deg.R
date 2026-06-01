library(tidyverse)
library(limma)
library(jsonlite)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(grid)

PRE_DIR     <- "../data/preprocessed_rds_files"
DEG_DIR     <- "../data/deg_rds_files"
RESULTS_DIR <- "../results/deg"
DATA_DIR    <- "../data"

datasets        <- c("GSE42568", "GSE21653", "GSE20711", "GSE88770")
SIGNATURE_GENES <- c("TSLP", "BIRC5", "S100B", "MDK", "S100P", "RARRES3", "BLNK", "ACO1")


merge_datasets <- function(datasets) {
    message("Merging datasets...")
    expr_list     <- lapply(datasets, function(d) readRDS(file.path(PRE_DIR, paste0(d, "_expression_matrix.rds"))))
    clinical_list <- lapply(datasets, function(d) readRDS(file.path(PRE_DIR, paste0(d, "_clinical_metadata.rds"))))

    merged_expr     <- Reduce(function(a, b) inner_join(a, b, by = "gene_symbol"), expr_list)
    merged_clinical <- bind_rows(clinical_list)

    message("  Merged: ", nrow(merged_expr), " genes x ", ncol(merged_expr) - 1, " samples")
    list(expr = merged_expr, clinical = merged_clinical)
}


build_design_matrix <- function(clinical) {
    groups <- factor(clinical$is_tumor)
    design <- model.matrix(~ 0 + groups)
    colnames(design) <- c("Normal", "Tumor")
    contrast <- makeContrasts(tumor_vs_normal = Tumor - Normal, levels = design)
    list(design = design, contrast = contrast)
}


run_limma_deg <- function(expr, design, contrast) {
    message("Fitting linear model...")
    fit          <- lmFit(expr, design)
    fit_contrast <- contrasts.fit(fit, contrast)
    fit_contrast <- eBayes(fit_contrast)
    topTable(fit_contrast, number = Inf, adjust = "BH")
}


filter_deg <- function(all_genes, pval_cutoff = 0.05, adjp_cutoff = 0.05, logfc_cutoff = 1) {
    message("Filtering DEGs...")
    deg <- all_genes[
        all_genes$P.Value    < pval_cutoff &
        all_genes$adj.P.Val  < adjp_cutoff &
        abs(all_genes$logFC) > logfc_cutoff,
    ]
    message("  ", nrow(deg), " significant DEGs identified")
    deg
}


load_immune_genes <- function(signature_genes) {
    message("Fetching immune gene sets from ImmPort...")
    immune_data  <- fromJSON(file.path(DATA_DIR, "immune_genes.json"))
    base_url     <- "https://s3.immport.org/release/genelists/current/"
    immune_genes <- list()

    for (i in seq_len(nrow(immune_data))) {
        id  <- immune_data$id[i]
        url <- paste0(base_url, id, ".json")
        try({
            json_data <- fromJSON(url)
            syms      <- json_data$genes$Symbol
            if (!is.null(syms)) immune_genes[[id]] <- syms
        }, silent = TRUE)
    }

    result <- unique(c(unlist(immune_genes), signature_genes))
    message("  ", length(result), " unique immune genes loaded")
    result
}


load_tumor_tfs <- function() {
    message("Loading tumor transcription factors...")
    tf_df <- read.csv(file.path(DATA_DIR, "transcription_factors.csv"))
    tf_df$Target_Gene
}


filter_immune_tf_degs <- function(deg_genes, immune_genes, tumor_tf) {
    deg_list    <- deg_genes$gene_symbol
    immune_syms <- intersect(deg_list, immune_genes)
    tf_syms     <- intersect(deg_list, tumor_tf)

    message("  Immune DEGs: ", length(immune_syms), " | Tumor TF DEGs: ", length(tf_syms))
    list(
        immune_deg   = deg_genes[deg_genes$gene_symbol %in% immune_syms, ],
        tumor_tf_deg = deg_genes[deg_genes$gene_symbol %in% tf_syms, ]
    )
}


filter_expression_matrix <- function(expr, immune_deg_genes) {
    filtered <- expr[expr$gene_symbol %in% immune_deg_genes$gene_symbol, ]
    message("Filtered expression matrix: ", nrow(filtered), " immune DEG genes")
    filtered
}


transpose_expression_matrix <- function(filtered_expr) {
    expr_t           <- as.data.frame(t(filtered_expr))
    colnames(expr_t) <- filtered_expr$gene_symbol
    expr_t           <- expr_t[-1, ]
    expr_t[]         <- lapply(expr_t, as.numeric)
    expr_t
}


plot_volcano <- function(all_genes, signature_genes, results_dir) {
    message("Plotting volcano...")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

    df           <- all_genes
    df$neg_log_p <- -log10(df$P.Value)

    df$category <- "Not Significant"
    df$category[df$adj.P.Val < 0.05 & df$logFC >  1] <- "Up-regulated"
    df$category[df$adj.P.Val < 0.05 & df$logFC < -1] <- "Down-regulated"
    df$category <- factor(df$category, levels = c("Up-regulated", "Down-regulated", "Not Significant"))

    sig_points <- df[df$gene_symbol %in% signature_genes, ]

    p <- ggplot(df, aes(x = logFC, y = neg_log_p, color = category)) +
        geom_point(alpha = 0.45, size = 1.2) +
        geom_point(data = sig_points, shape = 18, size = 3, color = "#FF9500") +
        geom_label_repel(
            data = sig_points, aes(label = gene_symbol),
            color = "#FF9500", size = 2.8, max.overlaps = 20, show.legend = FALSE
        ) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
        scale_color_manual(values = c(
            "Up-regulated"    = "#E64B35",
            "Down-regulated"  = "#4DBBD5",
            "Not Significant" = "grey75"
        )) +
        labs(
            title    = "Volcano Plot: Tumor vs Normal",
            subtitle = sprintf("DEGs: %d up  |  %d down  |  signature genes highlighted (◆)",
                               sum(df$category == "Up-regulated"),
                               sum(df$category == "Down-regulated")),
            x     = "log₂ Fold Change",
            y     = "-log₁₀(P-value)",
            color = NULL
        ) +
        theme_classic(base_size = 13) +
        theme(
            plot.title      = element_text(face = "bold", hjust = 0.5),
            plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey40"),
            legend.position = "top"
        )

    out_path <- file.path(results_dir, "volcano_plot.png")
    ggsave(out_path, p, width = 8, height = 6, dpi = 150)
    message("Saved: ", out_path)
    invisible(p)
}


plot_heatmap <- function(expr, deg_genes, clinical, signature_genes, results_dir) {
    message("Plotting heatmap...")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

    top40_genes   <- head(deg_genes[order(deg_genes$adj.P.Val), ], 40)$gene_symbol
    heatmap_genes <- unique(c(top40_genes, signature_genes))
    heatmap_genes <- heatmap_genes[heatmap_genes %in% expr$gene_symbol]

    heatmap_expr          <- expr[expr$gene_symbol %in% heatmap_genes, ]
    heatmap_mat           <- as.matrix(heatmap_expr[, -1])
    rownames(heatmap_mat) <- heatmap_expr$gene_symbol

    sample_ids <- colnames(heatmap_expr)[-1]
    col_annot  <- data.frame(
        Tissue    = ifelse(clinical$is_tumor == 1, "Tumor", "Normal"),
        row.names = sample_ids
    )

    row_annot <- data.frame(
        Type      = ifelse(rownames(heatmap_mat) %in% signature_genes, "Immune Signature", "Top DEG"),
        row.names = rownames(heatmap_mat)
    )

    annot_colors <- list(
        Tissue = c(Tumor = "#E64B35", Normal = "#4DBBD5"),
        Type   = c("Immune Signature" = "#FF9500", "Top DEG" = "grey80")
    )

    p_heat <- pheatmap(
        heatmap_mat,
        scale                    = "row",
        color                    = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
        annotation_col           = col_annot,
        annotation_row           = row_annot,
        annotation_colors        = annot_colors,
        show_colnames            = FALSE,
        fontsize_row             = 7,
        border_color             = NA,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        main                     = "Top DEGs Heatmap — Tumor vs Normal\n(immune signature genes highlighted)"
    )

    out_path <- file.path(results_dir, "deg_heatmap.png")
    png(out_path, width = 10 * 150, height = 10 * 150, res = 150)
    grid::grid.newpage()
    grid::grid.draw(p_heat$gtable)
    dev.off()
    message("Saved: ", out_path)
    invisible(p_heat)
}


save_deg_results <- function(all_genes, deg_genes, immune_deg, tumor_tf_deg,
                              filtered_expr, expr_transposed, clinical) {
    message("Saving DEG results to ", DEG_DIR, "...")
    dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)

    saveRDS(all_genes,       file.path(DEG_DIR, "all_genes.rds"))
    saveRDS(deg_genes,       file.path(DEG_DIR, "deg_genes.rds"))
    saveRDS(immune_deg,      file.path(DEG_DIR, "immune_deg_genes.rds"))
    saveRDS(tumor_tf_deg,    file.path(DEG_DIR, "tumor_tf_deg_genes.rds"))
    saveRDS(filtered_expr,   file.path(DEG_DIR, "filtered_expression_matrix.rds"))
    saveRDS(expr_transposed, file.path(DEG_DIR, "filtered_expression_matrix_transposed.rds"))
    saveRDS(clinical,        file.path(DEG_DIR, "filtered_clinical.rds"))

    message("All DEG results saved.")
}


# ONLY RE-RUN IF STARTING FROM SCRATCH DO NOT RE-RUN REGULARLY

merged   <- merge_datasets(datasets)
expr     <- merged$expr
clinical <- merged$clinical

design_obj <- build_design_matrix(clinical)
all_genes  <- run_limma_deg(expr, design_obj$design, design_obj$contrast)
deg_genes  <- filter_deg(all_genes)

immune_genes <- load_immune_genes(SIGNATURE_GENES)
tumor_tf     <- load_tumor_tfs()
filtered     <- filter_immune_tf_degs(deg_genes, immune_genes, tumor_tf)

immune_deg   <- filtered$immune_deg
tumor_tf_deg <- filtered$tumor_tf_deg

filtered_expr   <- filter_expression_matrix(expr, immune_deg)
expr_transposed <- transpose_expression_matrix(filtered_expr)

plot_volcano(all_genes, SIGNATURE_GENES, RESULTS_DIR)
plot_heatmap(expr, deg_genes, clinical, SIGNATURE_GENES, RESULTS_DIR)

save_deg_results(all_genes, deg_genes, immune_deg, tumor_tf_deg,
                 filtered_expr, expr_transposed, clinical)
