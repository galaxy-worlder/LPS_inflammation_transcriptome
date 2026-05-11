# scripts/05_GSE5504_time_series_validation.R

# ============================================================
# Time-series validation using GSE5504
# LPS-stimulated human monocytes
# Validate core LPS-responsive genes from GSE193336
# ============================================================

library(tidyverse)
library(GEOquery)
library(Biobase)
library(limma)
library(pheatmap)
library(stringr)

# Avoid function conflicts
select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
mutate <- dplyr::mutate

# -----------------------------
# 1. Create folders
# -----------------------------

dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/validation", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Load or download GSE5504
# -----------------------------

gse5504_rds <- "data/raw/GSE5504_eset.rds"

if (file.exists(gse5504_rds)) {
  message("Loading existing GSE5504 ExpressionSet...")
  eset <- readRDS(gse5504_rds)
} else {
  message("GSE5504 ExpressionSet not found. Downloading with AnnotGPL = TRUE...")
  gse <- getGEO(
    "GSE5504",
    GSEMatrix = TRUE,
    AnnotGPL = TRUE,
    getGPL = TRUE
  )
  
  if (length(gse) > 1) {
    message("GSE5504 has multiple platforms. Using the first ExpressionSet.")
  }
  
  eset <- gse[[1]]
  saveRDS(eset, gse5504_rds)
}

expr <- exprs(eset)
pheno <- pData(eset)
fdata <- fData(eset)

cat("Expression matrix dimension:\n")
print(dim(expr))

cat("Phenotype data columns:\n")
print(colnames(pheno))

cat("Feature data columns:\n")
print(colnames(fdata))

# -----------------------------
# 3. Prepare sample metadata
# -----------------------------

# Combine phenotype columns into one text field for robust time extraction.
pheno_text <- pheno %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(
    combined_text = apply(
      dplyr::select(., -sample_id),
      1,
      function(x) paste(x, collapse = " | ")
    )
  )

sample_info <- tibble(
  sample_id = colnames(expr)
) %>%
  left_join(pheno_text, by = "sample_id") %>%
  mutate(
    combined_text_lower = tolower(combined_text),
    time_raw = str_extract(
      combined_text_lower,
      regex("\\b(0|2|4|8|24)\\s*(h|hr|hrs|hour|hours)\\b", ignore_case = TRUE)
    ),
    time_hour = suppressWarnings(as.numeric(str_extract(time_raw, "\\d+")))
  )

# Rescue baseline samples described as no stimulation / control.
sample_info <- sample_info %>%
  mutate(
    time_hour = case_when(
      is.na(time_hour) & str_detect(
        combined_text_lower,
        "no stimulation|no-stimulation|no_stimulation|control|untreated|unstimulated|baseline"
      ) ~ 0,
      TRUE ~ time_hour
    )
  )

# Save metadata preview for troubleshooting.
write.csv(
  sample_info,
  "data/processed/GSE5504_sample_info_raw_parsed.csv",
  row.names = FALSE
)

if (any(is.na(sample_info$time_hour))) {
  warning("Some sample time points could not be parsed automatically.")
  warning("Please inspect data/processed/GSE5504_sample_info_raw_parsed.csv")
  
  print(sample_info %>% dplyr::select(sample_id, title, description, time_raw, time_hour))
  
  stop("Time point parsing failed for some samples. Send me the parsed metadata file if this happens.")
}

sample_info <- sample_info %>%
  mutate(
    time_hour = as.numeric(time_hour),
    time_label = paste0(time_hour, "h"),
    time_label = factor(time_label, levels = c("0h", "2h", "4h", "8h", "24h"))
  ) %>%
  arrange(time_hour, sample_id)

write.csv(
  sample_info,
  "data/processed/GSE5504_sample_info_clean.csv",
  row.names = FALSE
)

print(sample_info %>% dplyr::select(sample_id, title, time_hour, time_label))

# Reorder expression matrix
expr <- expr[, sample_info$sample_id]

# If expression values are not log2-transformed, transform them.
if (max(expr, na.rm = TRUE) > 50) {
  message("Expression values appear not log2-transformed. Applying log2(x + 1).")
  expr <- log2(expr + 1)
}

# -----------------------------
# 4. Probe annotation: map probes to gene symbols
# -----------------------------

# Try to find a gene symbol column automatically.
symbol_col_candidates <- colnames(fdata)[
  str_detect(colnames(fdata), regex("gene.?symbol|symbol", ignore_case = TRUE))
]

if (length(symbol_col_candidates) == 0) {
  message("No gene symbol column found in existing fData.")
  message("Trying to re-download GSE5504 with AnnotGPL = TRUE...")
  
  gse <- getGEO(
    "GSE5504",
    GSEMatrix = TRUE,
    AnnotGPL = TRUE,
    getGPL = TRUE
  )
  
  eset <- gse[[1]]
  expr <- exprs(eset)
  pheno <- pData(eset)
  fdata <- fData(eset)
  
  symbol_col_candidates <- colnames(fdata)[
    str_detect(colnames(fdata), regex("gene.?symbol|symbol", ignore_case = TRUE))
  ]
}

if (length(symbol_col_candidates) == 0) {
  stop("Still no gene symbol column found. Please send me colnames(fData(eset)).")
}

symbol_col <- symbol_col_candidates[1]
message("Using gene symbol column: ", symbol_col)

probe_anno <- fdata %>%
  as.data.frame() %>%
  rownames_to_column("probe_id") %>%
  transmute(
    probe_id = probe_id,
    symbol_raw = .data[[symbol_col]]
  ) %>%
  mutate(
    SYMBOL = str_split(symbol_raw, "///|//|;|,") %>% map_chr(~ .x[1]),
    SYMBOL = str_trim(SYMBOL),
    SYMBOL = toupper(SYMBOL),
    SYMBOL = ifelse(SYMBOL == "" | SYMBOL == "---" | is.na(SYMBOL), NA, SYMBOL)
  ) %>%
  filter(!is.na(SYMBOL))

write.csv(
  probe_anno,
  "data/processed/GSE5504_probe_annotation_clean.csv",
  row.names = FALSE
)

# -----------------------------
# 5. Collapse probes to gene level
# -----------------------------

expr_df <- as.data.frame(expr) %>%
  rownames_to_column("probe_id") %>%
  inner_join(probe_anno, by = "probe_id")

sample_cols <- sample_info$sample_id

expr_gene <- expr_df %>%
  mutate(mean_expr = rowMeans(dplyr::select(., all_of(sample_cols)), na.rm = TRUE)) %>%
  arrange(SYMBOL, desc(mean_expr)) %>%
  group_by(SYMBOL) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(SYMBOL, all_of(sample_cols)) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

write.csv(
  expr_gene,
  "data/processed/GSE5504_gene_expression_matrix.csv"
)

cat("Gene-level expression matrix dimension:\n")
print(dim(expr_gene))

# -----------------------------
# 6. Select validation genes from GSE193336
# -----------------------------

deg_all <- read.csv(
  "results/DEG/GSE193336_limma_voom_DEG_all.csv",
  check.names = FALSE
)

deg_all <- deg_all %>%
  mutate(
    SYMBOL = toupper(SYMBOL),
    regulation = case_when(
      adj.P.Val < 0.05 & logFC >= 1 ~ "Up",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

manual_core_genes <- c(
  "IL1B", "PTGS2", "CCL4", "CCL5", "CXCL3",
  "TNFAIP6", "MAP3K8", "NFKBIA", "CXCL8", "TNF",
  "IL6", "CCL2", "CXCL10", "STAT1", "IRF1",
  "NOD2", "JAK3", "LRRK2", "SOCS3", "ICAM1"
)

top_up_genes <- deg_all %>%
  filter(regulation == "Up", !is.na(SYMBOL), SYMBOL != "") %>%
  arrange(adj.P.Val) %>%
  pull(SYMBOL) %>%
  unique()

candidate_genes <- unique(c(manual_core_genes, top_up_genes))

present_genes <- intersect(candidate_genes, rownames(expr_gene))

if (length(present_genes) < 5) {
  stop("Too few validation genes found in GSE5504 platform. Check gene annotation.")
}

# Prefer manual inflammatory genes if present, then fill with top DE genes.
manual_present <- intersect(manual_core_genes, rownames(expr_gene))
top_present <- setdiff(present_genes, manual_present)

selected_genes <- unique(c(
  manual_present,
  head(top_present, max(0, 20 - length(manual_present)))
))

selected_genes <- selected_genes[1:min(20, length(selected_genes))]

write.csv(
  tibble(selected_genes = selected_genes),
  "results/validation/GSE5504_selected_validation_genes.csv",
  row.names = FALSE
)

cat("Selected validation genes:\n")
print(selected_genes)

# -----------------------------
# 7. Long-format expression data
# -----------------------------

expr_selected <- expr_gene[selected_genes, , drop = FALSE]

expr_long <- as.data.frame(expr_selected) %>%
  rownames_to_column("SYMBOL") %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample_id",
    values_to = "expression"
  ) %>%
  left_join(sample_info %>% dplyr::select(sample_id, time_hour, time_label), by = "sample_id")

# Gene-wise z-score
expr_z_long <- expr_long %>%
  group_by(SYMBOL) %>%
  mutate(z_expression = as.numeric(scale(expression))) %>%
  ungroup()

expr_summary <- expr_z_long %>%
  group_by(SYMBOL, time_hour, time_label) %>%
  summarise(
    mean_z = mean(z_expression, na.rm = TRUE),
    sd_z = sd(z_expression, na.rm = TRUE),
    n = n(),
    se_z = sd_z / sqrt(n),
    .groups = "drop"
  )

write.csv(
  expr_long,
  "results/validation/GSE5504_selected_gene_expression_long.csv",
  row.names = FALSE
)

write.csv(
  expr_summary,
  "results/validation/GSE5504_selected_gene_timecourse_summary.csv",
  row.names = FALSE
)

# -----------------------------
# 8. Validation summary table
# -----------------------------

baseline <- expr_summary %>%
  filter(time_hour == 0) %>%
  dplyr::select(SYMBOL, baseline_mean_z = mean_z)

validation_summary <- expr_summary %>%
  left_join(baseline, by = "SYMBOL") %>%
  mutate(delta_from_0 = mean_z - baseline_mean_z) %>%
  group_by(SYMBOL) %>%
  summarise(
    max_delta = max(delta_from_0, na.rm = TRUE),
    min_delta = min(delta_from_0, na.rm = TRUE),
    peak_time = time_hour[which.max(delta_from_0)],
    trough_time = time_hour[which.min(delta_from_0)],
    .groups = "drop"
  ) %>%
  left_join(
    deg_all %>%
      dplyr::select(SYMBOL, GSE193336_logFC = logFC, GSE193336_adjP = adj.P.Val, GSE193336_regulation = regulation) %>%
      distinct(SYMBOL, .keep_all = TRUE),
    by = "SYMBOL"
  ) %>%
  arrange(desc(max_delta))

write.csv(
  validation_summary,
  "results/validation/GSE5504_validation_summary.csv",
  row.names = FALSE
)

print(validation_summary)

# -----------------------------
# 9. Time-course line plot
# -----------------------------

p_timecourse <- ggplot(
  expr_summary,
  aes(x = time_hour, y = mean_z, group = SYMBOL)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
    width = 0.3,
    linewidth = 0.4
  ) +
  facet_wrap(~ SYMBOL, scales = "free_y", ncol = 5) +
  scale_x_continuous(breaks = c(0, 2, 4, 8, 24)) +
  theme_bw(base_size = 12) +
  labs(
    title = "Time-series validation of LPS-responsive genes in GSE5504",
    x = "Time after LPS stimulation (hours)",
    y = "Gene-wise z-scored expression"
  )

ggsave(
  "figures/validation/Figure6_GSE5504_timecourse_selected_genes.pdf",
  p_timecourse,
  width = 13,
  height = 9
)

ggsave(
  "figures/validation/Figure6_GSE5504_timecourse_selected_genes.png",
  p_timecourse,
  width = 13,
  height = 9,
  dpi = 300
)

# -----------------------------
# 10. Heatmap across time points
# -----------------------------

expr_z_mat <- t(scale(t(expr_selected)))
expr_z_mat <- expr_z_mat[, sample_info$sample_id]

annotation_col <- sample_info %>%
  dplyr::select(sample_id, time_label) %>%
  column_to_rownames("sample_id")

pdf(
  "figures/validation/Figure6_GSE5504_validation_heatmap.pdf",
  width = 8,
  height = 8
)

pheatmap(
  expr_z_mat,
  annotation_col = annotation_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 8,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "Validation of LPS-responsive genes in GSE5504"
)

dev.off()

png(
  "figures/validation/Figure6_GSE5504_validation_heatmap.png",
  width = 2400,
  height = 2400,
  res = 300
)

pheatmap(
  expr_z_mat,
  annotation_col = annotation_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 8,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "Validation of LPS-responsive genes in GSE5504"
)

dev.off()

# -----------------------------
# 11. Save session information
# -----------------------------

sink("results/validation/GSE5504_validation_sessionInfo.txt")
sessionInfo()
sink()

message("GSE5504 time-series validation completed.")