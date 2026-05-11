# scripts/00_project_setup.R

# ================================
# Project setup
# LPS inflammation transcriptome
# ================================

dirs <- c(
  "data",
  "data/raw",
  "data/processed",
  "data/metadata",
  "data/supplementary",
  "scripts",
  "results",
  "results/DEG",
  "results/enrichment",
  "results/validation",
  "figures",
  "manuscript"
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

cat("Project folders created successfully.\n")
cat("Working directory:\n")
cat(getwd(), "\n")