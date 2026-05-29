library(tidyverse)
library(GEOquery)
library(affy)

# defining datapath
CEL_DIR <- "../data/cel_files"
TAR_DIR <- "../data/tar_files"
AFFY_DIR <- "../data/raw_affy_rds_files"
META_DIR <- "../data/geo_metadata_rds_files"

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

# ONLY RE-RUN IF STARTING FROM SCRATCH DO NOT RE-RUN REGULARLY

for (dataset in datasets) {
    load_geo_data(dataset)
}

for (dataset in datasets) {
    prepare_raw_affy_dataset(dataset)
}