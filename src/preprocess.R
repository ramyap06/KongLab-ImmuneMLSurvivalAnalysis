library(tidyverse)
library(GEOquery)
library(affy)

# defining datapath
CEL_DIR <- "../data/cel_files"
TAR_DIR <- "../data/tar_files"
AFFY_DIR <- "../data/raw_affy_rds_files"
META_DIR <- "../data/geo_metadata_rds_files"
PRE_DIR <- "../data/preprocessed_rds_files"

# defining datasets
datasets <- c("GSE42568", "GSE21653", "GSE20711", "GSE88770")

normalize_affy_dataset <- function(raw_affy, dataset_name) {
    message("[", dataset_name, "] Normalizing...")
    normalized_data <- rma(raw_affy)
    normalized_expr <- as.data.frame(exprs(normalized_data))
    normalized_expr <- tibble::rownames_to_column(normalized_expr, var = "ID")
    colnames(normalized_expr) <- gsub("_.*", "", colnames(normalized_expr))
    return(normalized_expr)
}

annotate_expression_matrix <- function(normalized_expr, gse, dataset_name) {
    message("[", dataset_name, "] Annotating...")
    feature_data <- gse[[1]]@featureData@data
    feature_data <- feature_data[, c(1, 11)]
    colnames(feature_data) <- c("ID", "Gene Symbol")
    annotated_expr <- normalized_expr |>
        inner_join(feature_data, by = "ID")
    return(annotated_expr)
}

collapse_probe_expression <- function(annotated_expr, dataset_name) {
    message("[", dataset_name, "] Collapsing probes...")
    collapsed_expr <- annotated_expr |>
        filter(!grepl("///", `Gene Symbol`)) |>
        filter(`Gene Symbol` != "") |>
        group_by(`Gene Symbol`) |>
        summarise(across(where(is.numeric), mean)) |>
        ungroup() |>
        rename(gene_symbol = `Gene Symbol`)
    return(collapsed_expr)
}

# revisit this -> need to safely change data types
process_clinical_metadata <- function(gse, dataset_name) {
    message("[", dataset_name, "] Processing clinical metadata...")
    clinical <- pData(gse[[1]])
    clinical[] <- lapply(clinical, as.character)
    clinical$is_tumor <- as.integer(
        ifelse(grepl("tumor", clinical$source_name_ch1, ignore.case = TRUE), 1, 0)
    )
    rownames(clinical) <- NULL
    return(clinical)
}

export_processed_datasets <- function(expression_matrix, clinical_metadata, dataset_name) {
    message("[", dataset_name, "] Exporting...")
    dir.create(PRE_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(expression_matrix, file.path(PRE_DIR, paste0(dataset_name, "_expression_matrix.rds")))
    saveRDS(clinical_metadata, file.path(PRE_DIR, paste0(dataset_name, "_clinical_metadata.rds")))
}

# ONLY RE-RUN IF STARTING FROM SCRATCH DO NOT RE-RUN REGULARLY

for (dataset in datasets) {
    raw_affy <- readRDS(file.path(AFFY_DIR, paste0(dataset, "_raw_affy.rds")))
    gse <- readRDS(file.path(META_DIR,  paste0(dataset, "_metadata.rds")))

    normalized_expr <- normalize_affy_dataset(raw_affy, dataset)
    annotated_expr <- annotate_expression_matrix(normalized_expr, gse, dataset)
    final_expr <- collapse_probe_expression(annotated_expr, dataset)
    clinical <- process_clinical_metadata(gse, dataset)

    export_processed_datasets(final_expr, clinical, dataset)
    message("[", dataset, "] Done.\n")
}