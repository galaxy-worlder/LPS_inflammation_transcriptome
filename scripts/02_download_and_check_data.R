# scripts/02_download_and_check_data.R

# ================================
# Download and inspect GEO datasets
# ================================

library(GEOquery)
library(tidyverse)
library(data.table)
library(janitor)

# ----------------
# 1. Create folders
# ----------------

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("data/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# ----------------
# 2. GEO datasets
# ----------------

gse_ids <- c(
  "GSE193336",  # LPS-stimulated human macrophages, RNA-seq
  "GSE5504",    # LPS time-series human monocytes, microarray
  "GSE154918"  # sepsis blood RNA-seq
      
)

# ----------------
# 3. Helper function
# ----------------

download_geo_matrix <- function(gse_id) {
  message("Downloading GEO matrix: ", gse_id)
  
  gse <- getGEO(
    GEO = gse_id,
    GSEMatrix = TRUE,
    AnnotGPL = FALSE,
    getGPL = FALSE
  )
  
  if (length(gse) > 1) {
    message(gse_id, " has multiple platforms. Using the first ExpressionSet for now.")
  }
  
  eset <- gse[[1]]
  
  saveRDS(eset, file = file.path("data/raw", paste0(gse_id, "_eset.rds")))
  
  pheno <- pData(eset) %>%
    rownames_to_column("sample_id") %>%
    clean_names()
  
  write.csv(
    pheno,
    file = file.path("data/metadata", paste0(gse_id, "_metadata.csv")),
    row.names = FALSE
  )
  
  expr_dim <- dim(exprs(eset))
  
  summary_row <- tibble(
    accession = gse_id,
    title = experimentData(eset)@title,
    platform = annotation(eset),
    n_features = expr_dim[1],
    n_samples = expr_dim[2],
    metadata_file = paste0("data/metadata/", gse_id, "_metadata.csv"),
    eset_file = paste0("data/raw/", gse_id, "_eset.rds")
  )
  
  return(summary_row)
}

# ----------------
# 4. Download GEO matrix files
# ----------------

dataset_summary <- map_dfr(gse_ids, download_geo_matrix)

write.csv(
  dataset_summary,
  file = "data/metadata/dataset_summary.csv",
  row.names = FALSE
)

print(dataset_summary)

# ----------------
# 5. Download supplementary files for GSE193336
# ----------------
# GSE193336 provides processed raw count CSV, which is useful for DESeq2.

message("Downloading supplementary files for GSE193336...")

getGEOSuppFiles(
  GEO = "GSE193336",
  baseDir = "data/supplementary",
  makeDirectory = TRUE
)

supp_files <- list.files(
  "data/supplementary/GSE193336",
  recursive = TRUE,
  full.names = TRUE
)

writeLines(
  supp_files,
  con = "data/metadata/GSE193336_supplementary_files.txt"
)

print(supp_files)

# ----------------
# 6. Preview GSE193336 raw count file
# ----------------

count_file <- supp_files[
  grepl("Raw_Count|raw_count|count|Count", supp_files, ignore.case = TRUE) &
    grepl("\\.csv\\.gz$|\\.csv$|\\.txt\\.gz$|\\.txt$", supp_files, ignore.case = TRUE)
]

if (length(count_file) > 0) {
  message("Possible count file found:")
  print(count_file)
  
  count_preview <- fread(count_file[1], nrows = 10)
  
  write.csv(
    count_preview,
    file = "data/metadata/GSE193336_count_preview.csv",
    row.names = FALSE
  )
  
  message("Preview saved to data/metadata/GSE193336_count_preview.csv")
} else {
  warning("No obvious count file found for GSE193336. Check supplementary file list manually.")
}

message("Data download and initial check completed.")