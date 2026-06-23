library(tidyverse)
library(GEOquery)
library(affy)
library(jsonlite)

ROOT_DIR <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))

CEL_DIR          <- file.path(ROOT_DIR, "data/cel_files")
TAR_DIR          <- file.path(ROOT_DIR, "data/tar_files")
AFFY_DIR         <- file.path(ROOT_DIR, "data/raw_affy_rds_files")
META_DIR         <- file.path(ROOT_DIR, "data/geo_metadata_rds_files")
IMMUNE_GENES_DIR <- file.path(ROOT_DIR, "data/immune_genes_and_tf")
IMMUNE_PATH      <- file.path(IMMUNE_GENES_DIR, "all_immune_genes_ImmPort.json")
LEGACY_DIR       <- file.path(IMMUNE_GENES_DIR, "legacy")
LEGACY_OUT_PATH  <- file.path(IMMUNE_GENES_DIR, "legacy_immune_genes_ImmPort.txt")

IMMPORT_BASE_URL <- "https://s3.immport.org/release/genelists/current/"

SIGNATURE_GENES <- c("TSLP", "BIRC5", "S100B", "MDK", "S100P", "RARRES3", "BLNK", "ACO1")

# defining datasets
datasets <- c("GSE42568", "GSE21653", "GSE20711", "GSE88770")

load_geo_data <- function(name) {
    message("[", name, "] Downloading...")
    dir.create(TAR_DIR, recursive = TRUE, showWarnings = FALSE)
    dir.create(CEL_DIR, recursive = TRUE, showWarnings = FALSE)
    dir.create(META_DIR, recursive = TRUE, showWarnings = FALSE)

    getGEOSuppFiles(name, baseDir = TAR_DIR)

    tar_path <- file.path(TAR_DIR, name, paste0(name, "_RAW.tar"))
    if (!file.exists(tar_path)) stop("RAW tar not found: ", name)

    tar_size_mb <- round(file.info(tar_path)$size / 1024^2, 2)
    message("[", name, "] Archive: ", tar_size_mb, " MB")

    gse <- getGEO(name, GSEMatrix = TRUE)
    saveRDS(gse, file.path(META_DIR, paste0(name, "_metadata.rds")))

    message("[", name, "] Done.")
}

prepare_raw_affy_dataset <- function(name) {
    message("[", name, "] Preparing Affy dataset...")
    dir.create(CEL_DIR, recursive = TRUE, showWarnings = FALSE)
    dir.create(AFFY_DIR, recursive = TRUE, showWarnings = FALSE)

    tar_path <- file.path(TAR_DIR, name, paste0(name, "_RAW.tar"))
    if (!file.exists(tar_path)) stop("RAW tar not found: ", name)

    extract_path <- file.path(CEL_DIR, name)
    untar(tarfile = tar_path, exdir = extract_path)

    raw_affy <- ReadAffy(celfile.path = extract_path)
    saveRDS(raw_affy, file.path(AFFY_DIR, paste0(name, "_raw_affy.rds")))

    message("[", name, "] Done.")
}

# ── immune gene catalog ───────────────────────────────────────────────────────

fetch_immune_gene_catalog <- function(catalog_path, base_url, signature_genes) {
    immune_data  <- fromJSON(catalog_path)
    immune_genes <- list()

    for (i in seq_len(nrow(immune_data))) {
        id  <- immune_data$id[i]
        url <- paste0(base_url, id, ".json")

        # if the GO/Reactome term has no published gene list skip silently
        try({
            json_data   <- fromJSON(url)
            gene_symbol <- json_data$genes$Symbol
            if (!is.null(gene_symbol)) immune_genes[[id]] <- gene_symbol
        }, silent = TRUE)
    }

    immune_genes <- unique(unlist(immune_genes))
    unique(c(immune_genes, signature_genes))
}

build_immune_gene_catalog <- function(out_dir, catalog_path, base_url, signature_genes) {
    message("Building ImmPort immune gene catalog (GO/Reactome term lists)...")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    immune_genes <- fetch_immune_gene_catalog(catalog_path, base_url, signature_genes)
    out_path <- file.path(out_dir, "immune_genes.rds")
    saveRDS(immune_genes, out_path)
    message("Saved ", length(immune_genes), " unique genes: ", out_path)
    invisible(immune_genes)
}

# ── legacy immune gene merge ─────────────────────────────────────────────────

build_legacy_immune_gene_list <- function(legacy_dir, out_path) {
    message("Building legacy immune gene list from ", legacy_dir, "...")
    txt_files   <- list.files(legacy_dir, pattern = "\\.txt$", full.names = TRUE)
    all_symbols <- sort(unique(unlist(lapply(txt_files, function(f) {
        read.delim(f, stringsAsFactors = FALSE)$Symbol
    }))))
    writeLines(all_symbols, out_path)
    message("Saved ", length(all_symbols), " unique genes: ", out_path)
    invisible(all_symbols)
}

# ── main ──────────────────────────────────────────────────────────────────────

# ONLY RE-RUN IF STARTING FROM SCRATCH DO NOT RE-RUN REGULARLY

# for (dataset in datasets) {
#     load_geo_data(dataset)
# }

# for (dataset in datasets) {
#     prepare_raw_affy_dataset(dataset)
# }

# build_immune_gene_catalog(IMMUNE_GENES_DIR, IMMUNE_PATH, IMMPORT_BASE_URL, SIGNATURE_GENES)
# build_legacy_immune_gene_list(LEGACY_DIR, LEGACY_OUT_PATH)