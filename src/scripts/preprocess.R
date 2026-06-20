library(tidyverse)
library(GEOquery)
library(affy)

ROOT_DIR <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))

CEL_DIR  <- file.path(ROOT_DIR, "data/cel_files")
TAR_DIR  <- file.path(ROOT_DIR, "data/tar_files")
AFFY_DIR <- file.path(ROOT_DIR, "data/raw_affy_rds_files")
META_DIR <- file.path(ROOT_DIR, "data/geo_metadata_rds_files")
CLIN_DIR <- file.path(ROOT_DIR, "data/raw_clinical_rds_files")
PRE_DIR  <- file.path(ROOT_DIR, "data/preprocessed_rds_files")

DATASET_CONFIGS <- list(
    GSE42568 = list(
        role = "train",
        rename = c(
            is_tumor = "tissue:ch1",
            rfs_event         = "relapse free survival event:ch1",
            rfs_time          = "relapse free survival time_days:ch1",
            er_status         = "er_status:ch1",
            tumor_grade       = "grade:ch1",
            lymph_node_status = "lymph node status:ch1",
            t_stage           = "size:ch1"
        ),
        keep = c(
            "geo_accession", "is_tumor", "rfs_event", "rfs_time",
            "er_status", "tumor_grade", "lymph_node_status", "t_stage"
        ),
        tumor_label = "breast cancer",
        time_unit   = "days"
    ),
    GSE21653 = list(
        role = "test",
        rename = c(
            is_tumor  = "tissue:ch1",
            rfs_event = "dfs evt:ch1",
            rfs_time  = "dfs time (months):ch1"
        ),
        keep        = c("geo_accession", "is_tumor", "rfs_event", "rfs_time"),
        tumor_label = "breast cancer tumor",
        time_unit   = "months"
    ),
    GSE20711 = list(
        role = "test",
        rename = c(
            is_tumor  = "source_name_ch1",
            rfs_event = "e.rfs:ch1",
            rfs_time  = "t.rfs:ch1"
        ),
        keep        = c("geo_accession", "is_tumor", "rfs_event", "rfs_time"),
        tumor_label = "Breast tumor",
        time_unit   = "years"
    ),
    GSE88770 = list(
        role = "test",
        rename = c(
            is_tumor  = "tissue:ch1",
            rfs_event = "drfs_event:ch1",
            rfs_time  = "drfs_or_last_contact_years:ch1"
        ),
        keep        = c("geo_accession", "is_tumor", "rfs_event", "rfs_time"),
        tumor_label = "Breast cancer tumor",
        time_unit   = "years"
    )
)

# ── expression helpers ────────────────────────────────────────────────────────

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
    normalized_expr |> inner_join(feature_data, by = "ID")
}

collapse_probe_expression <- function(annotated_expr, dataset_name) {
    message("[", dataset_name, "] Collapsing probes...")
    annotated_expr |>
        filter(!grepl("///", `Gene Symbol`), `Gene Symbol` != "") |>
        group_by(`Gene Symbol`) |>
        summarise(across(where(is.numeric), mean)) |>
        ungroup() |>
        rename(gene_symbol = `Gene Symbol`)
}

# ── clinical helpers ──────────────────────────────────────────────────────────

rename_columns <- function(df, rename_map) {
    if (length(rename_map) == 0) return(df)
    dplyr::rename(df, !!!rlang::syms(rename_map))
}

select_features <- function(df, keep) {
    df[, intersect(keep, colnames(df)), drop = FALSE]
}

# Converts rfs_time to days regardless of original unit
to_days <- function(time, unit) {
    switch(unit,
        days   = time,
        months = time * 30.44,
        years  = time * 365.25,
        stop("Unknown time unit: ", unit)
    )
}

# Normalises rfs_event to 1 (relapse) / 0 (no relapse) / NA
parse_rfs_event <- function(x) {
    x <- na_if(trimws(tolower(x)), "na")
    if (all(x[!is.na(x)] %in% c("0", "1"))) return(as.integer(x))
    as.integer(case_when(
        x == "yes" ~ 1L,
        x == "no"  ~ 0L
    ))
}

# Drops rows that are unusable for survival analysis.
drop_invalid <- function(df) {
    df[!is.na(df$rfs_event) & !is.na(df$rfs_time) & df$rfs_time > 0, ]
}

process_clinical_metadata <- function(gse, cfg, dataset_name) {
    message("[", dataset_name, "] Processing clinical metadata...")
    clinical <- pData(gse[[1]])
    clinical[] <- lapply(clinical, as.character)
    rownames(clinical) <- NULL

    df <- clinical |>
        rename_columns(cfg$rename) |>
        mutate(
            is_tumor  = as.integer(is_tumor == cfg$tumor_label),
            rfs_event = parse_rfs_event(rfs_event),
            rfs_time  = to_days(as.numeric(na_if(rfs_time, "NA")), cfg$time_unit)
        )

    tumor_rows     <- df |> filter(is_tumor == 1) |> drop_invalid() |> tidyr::drop_na()
    non_tumor_rows <- df |> filter(is_tumor != 1 | is.na(is_tumor))

    bind_rows(tumor_rows, non_tumor_rows) |>
        select_features(cfg$keep)
}

sync_expression_to_clinical <- function(expr, clinical) {
    retained_samples <- intersect(clinical$geo_accession, colnames(expr))
    expr[, c("gene_symbol", retained_samples)]
}

export_processed_datasets <- function(expression_matrix, clinical_metadata, dataset_name) {
    message("[", dataset_name, "] Exporting...")
    dir.create(PRE_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(expression_matrix, file.path(PRE_DIR, paste0(dataset_name, "_expression_matrix.rds")))
    saveRDS(clinical_metadata, file.path(PRE_DIR, paste0(dataset_name, "_clinical_metadata.rds")))
}

# ── main loop ─────────────────────────────────────────────────────────────────

for (dataset in names(DATASET_CONFIGS)) {
    cfg      <- DATASET_CONFIGS[[dataset]]
    raw_affy <- readRDS(file.path(AFFY_DIR, paste0(dataset, "_raw_affy.rds")))
    gse      <- readRDS(file.path(META_DIR,  paste0(dataset, "_metadata.rds")))

    dir.create(CLIN_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(pData(gse[[1]]), file.path(CLIN_DIR, paste0(dataset, "_raw_clinical.rds")))

    normalized_expr <- normalize_affy_dataset(raw_affy, dataset)
    annotated_expr  <- annotate_expression_matrix(normalized_expr, gse, dataset)
    final_expr      <- collapse_probe_expression(annotated_expr, dataset)
    clinical        <- process_clinical_metadata(gse, cfg, dataset)

    final_expr <- sync_expression_to_clinical(final_expr, clinical)

    export_processed_datasets(final_expr, clinical, dataset)
    message("[", dataset, "] Done.\n")
}
