# scripts/06_GSE154918_metadata_check.R

# ============================================================
# Inspect metadata of GSE154918
# Sepsis relevance analysis preparation
# ============================================================

library(tidyverse)
library(GEOquery)
library(Biobase)
library(stringr)
library(janitor)

# Avoid function conflicts
select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
mutate <- dplyr::mutate

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/sepsis", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Load GSE154918 ExpressionSet
# -----------------------------

gse154918_rds <- "data/raw/GSE154918_eset.rds"

if (file.exists(gse154918_rds)) {
  message("Loading existing GSE154918 ExpressionSet...")
  eset <- readRDS(gse154918_rds)
} else {
  message("GSE154918 ExpressionSet not found. Downloading...")
  gse <- getGEO(
    "GSE154918",
    GSEMatrix = TRUE,
    AnnotGPL = FALSE,
    getGPL = FALSE
  )
  
  if (length(gse) > 1) {
    message("GSE154918 has multiple platforms. Using the first ExpressionSet.")
  }
  
  eset <- gse[[1]]
  saveRDS(eset, gse154918_rds)
}

expr <- exprs(eset)
pheno <- pData(eset)

cat("Expression matrix dimension:\n")
print(dim(expr))

cat("Phenotype columns:\n")
print(colnames(pheno))

# -----------------------------
# 2. Save raw metadata
# -----------------------------

pheno_clean <- pheno %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  clean_names()

write.csv(
  pheno_clean,
  "data/processed/GSE154918_metadata_clean.csv",
  row.names = FALSE
)

# -----------------------------
# 3. Combine all metadata text
# -----------------------------

pheno_text <- pheno_clean %>%
  mutate(
    combined_text = apply(
      dplyr::select(., -sample_id),
      1,
      function(x) paste(x, collapse = " | ")
    ),
    combined_text_lower = tolower(combined_text)
  )

write.csv(
  pheno_text,
  "data/processed/GSE154918_metadata_combined_text.csv",
  row.names = FALSE
)

# -----------------------------
# 4. Try to infer group labels
# -----------------------------

sample_group_guess <- pheno_text %>%
  mutate(
    group_guess = case_when(
      str_detect(combined_text_lower, "healthy|control|normal") ~ "Healthy",
      str_detect(combined_text_lower, "septic shock") ~ "Septic_shock",
      str_detect(combined_text_lower, "sepsis") ~ "Sepsis",
      str_detect(combined_text_lower, "infection") ~ "Infection",
      str_detect(combined_text_lower, "follow") ~ "Follow_up",
      TRUE ~ "Unknown"
    )
  )

write.csv(
  sample_group_guess,
  "data/processed/GSE154918_sample_group_guess.csv",
  row.names = FALSE
)

group_summary <- sample_group_guess %>%
  count(group_guess, name = "n_samples") %>%
  arrange(desc(n_samples))

write.csv(
  group_summary,
  "results/sepsis/GSE154918_group_guess_summary.csv",
  row.names = FALSE
)

print(group_summary)

# -----------------------------
# 5. Print potentially useful columns
# -----------------------------

cat("\nColumns that may contain group information:\n")

possible_group_cols <- colnames(pheno_clean)[
  str_detect(
    colnames(pheno_clean),
    regex("group|disease|diagnosis|status|condition|phenotype|characteristics|title|source", ignore_case = TRUE)
  )
]

print(possible_group_cols)

if (length(possible_group_cols) > 0) {
  for (col in possible_group_cols) {
    cat("\n==============================\n")
    cat("Column:", col, "\n")
    cat("==============================\n")
    print(table(pheno_clean[[col]], useNA = "ifany"))
  }
}

message("GSE154918 metadata check completed.")