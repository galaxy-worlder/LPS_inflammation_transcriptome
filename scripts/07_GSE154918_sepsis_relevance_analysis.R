# scripts/07_GSE154918_sepsis_relevance_analysis.R

# ============================================================
# Clinical sepsis relevance analysis using GSE154918
# Fixed version: read expression matrix from GEO supplementary file
# ============================================================

library(tidyverse)
library(GEOquery)
library(Biobase)
library(limma)
library(pheatmap)
library(ggrepel)
library(stringr)
library(janitor)

select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
mutate <- dplyr::mutate

options(timeout = 600)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/sepsis", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/sepsis", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Load metadata from GEO
# -----------------------------

gse154918_rds <- "data/raw/GSE154918_eset_metadata.rds"

if (file.exists(gse154918_rds)) {
  message("Loading existing GSE154918 metadata ExpressionSet...")
  eset <- readRDS(gse154918_rds)
} else {
  message("Downloading GSE154918 metadata...")
  gse <- getGEO(
    "GSE154918",
    GSEMatrix = TRUE,
    AnnotGPL = FALSE,
    getGPL = FALSE
  )
  
  eset <- gse[[1]]
  saveRDS(eset, gse154918_rds)
}

pheno <- pData(eset)

metadata <- pheno %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  clean_names()

write.csv(
  metadata,
  "data/processed/GSE154918_metadata_clean.csv",
  row.names = FALSE
)

if (!"status_ch1" %in% colnames(metadata)) {
  stop("Column status_ch1 not found. Please check metadata.")
}

# -----------------------------
# 2. Download and read supplementary expression matrix
# -----------------------------

supp_file <- "data/raw/GSE154918_Schughart_Sepsis_200320.txt.gz"

if (!file.exists(supp_file)) {
  message("Downloading GSE154918 supplementary expression matrix...")
  
  download.file(
    url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE154918&file=GSE154918_Schughart_Sepsis_200320.txt.gz&format=file",
    destfile = supp_file,
    mode = "wb",
    method = "libcurl"
  )
}

message("Reading supplementary expression matrix...")

expr_raw <- read.table(
  gzfile(supp_file),
  header = TRUE,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\"",
  comment.char = ""
)

cat("Raw supplementary expression table dimension:\n")
print(dim(expr_raw))

cat("First 20 column names:\n")
print(colnames(expr_raw)[1:min(20, ncol(expr_raw))])

# -----------------------------
# 3. Match expression columns to metadata
# -----------------------------

all_expr_cols <- colnames(expr_raw)

sample_cols_available <- intersect(
  all_expr_cols,
  c(metadata$sample_id, metadata$geo_accession, metadata$title)
)

if (length(sample_cols_available) == 0) {
  stop("No sample columns in expression matrix matched metadata sample_id / geo_accession / title.")
}

non_sample_cols <- setdiff(all_expr_cols, sample_cols_available)

cat("Number of matched sample columns:\n")
print(length(sample_cols_available))

cat("Non-sample columns:\n")
print(non_sample_cols)

# Choose gene symbol column.
symbol_candidates <- non_sample_cols[
  str_detect(
    non_sample_cols,
    regex("symbol|gene.?symbol|gene.?name|external.?gene.?name|hgnc", ignore_case = TRUE)
  )
]

if (length(symbol_candidates) > 0) {
  gene_col <- symbol_candidates[1]
} else {
  gene_col <- non_sample_cols[1]
}

message("Using gene column: ", gene_col)

# -----------------------------
# 4. Prepare sample information
# -----------------------------

sample_info <- metadata %>%
  transmute(
    sample_id = sample_id,
    geo_accession = geo_accession,
    title = title,
    status_ch1 = status_ch1,
    sex = if ("sex_ch1" %in% colnames(metadata)) {
      str_replace(sex_ch1, regex("^sex:\\s*", ignore_case = TRUE), "")
    } else {
      NA_character_
    },
    clinical_group = case_when(
      status_ch1 == "Hlty" ~ "Healthy",
      status_ch1 %in% c("Seps_P", "Shock_P") ~ "Acute_sepsis",
      status_ch1 == "Inf1_P" ~ "Infection",
      status_ch1 %in% c("Seps_FU", "Shock_FU") ~ "Follow_up",
      TRUE ~ "Other"
    ),
    sepsis_subgroup = case_when(
      status_ch1 == "Hlty" ~ "Healthy",
      status_ch1 == "Seps_P" ~ "Sepsis",
      status_ch1 == "Shock_P" ~ "Septic_shock",
      status_ch1 == "Inf1_P" ~ "Infection",
      status_ch1 == "Seps_FU" ~ "Sepsis_follow_up",
      status_ch1 == "Shock_FU" ~ "Shock_follow_up",
      TRUE ~ "Other"
    ),
    expr_col = case_when(
      title %in% all_expr_cols ~ title,
      geo_accession %in% all_expr_cols ~ geo_accession,
      sample_id %in% all_expr_cols ~ sample_id,
      TRUE ~ NA_character_
    )
  )

write.csv(
  sample_info,
  "data/processed/GSE154918_sample_info_clean.csv",
  row.names = FALSE
)

group_summary <- sample_info %>%
  count(clinical_group, sepsis_subgroup, name = "n_samples") %>%
  arrange(clinical_group, sepsis_subgroup)

write.csv(
  group_summary,
  "results/sepsis/GSE154918_clean_group_summary.csv",
  row.names = FALSE
)

print(group_summary)

analysis_info <- sample_info %>%
  filter(clinical_group %in% c("Healthy", "Acute_sepsis")) %>%
  filter(!is.na(expr_col)) %>%
  mutate(
    clinical_group = factor(clinical_group, levels = c("Healthy", "Acute_sepsis")),
    sex = factor(sex)
  )

cat("Analysis group sample numbers:\n")
print(table(analysis_info$clinical_group, analysis_info$sepsis_subgroup))

# -----------------------------
# 5. Build gene-level expression matrix
# -----------------------------

expr_df <- expr_raw %>%
  dplyr::select(all_of(c(gene_col, analysis_info$expr_col))) %>%
  rename(SYMBOL = all_of(gene_col)) %>%
  mutate(
    SYMBOL = as.character(SYMBOL),
    SYMBOL = str_split(SYMBOL, "///|//|;|,|\\|") %>% map_chr(~ .x[1]),
    SYMBOL = str_trim(SYMBOL),
    SYMBOL = str_replace_all(SYMBOL, "\"", ""),
    SYMBOL = toupper(SYMBOL),
    SYMBOL = ifelse(SYMBOL == "" | SYMBOL == "---" | is.na(SYMBOL), NA, SYMBOL)
  ) %>%
  filter(!is.na(SYMBOL))

expr_numeric <- expr_df %>%
  mutate(
    across(
      all_of(analysis_info$expr_col),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

expr_gene <- expr_numeric %>%
  mutate(mean_expr = rowMeans(dplyr::select(., all_of(analysis_info$expr_col)), na.rm = TRUE)) %>%
  arrange(SYMBOL, desc(mean_expr)) %>%
  group_by(SYMBOL) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(SYMBOL, all_of(analysis_info$expr_col)) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

# Remove genes with too many missing values
expr_gene <- expr_gene[rowSums(is.na(expr_gene)) == 0, , drop = FALSE]

# Transform only if values are obviously not log-scale.
if (max(expr_gene, na.rm = TRUE) > 100) {
  message("Expression values appear not log2-transformed. Applying log2(x + 1).")
  expr_gene <- log2(expr_gene + 1)
}

cat("Final gene-level expression matrix dimension:\n")
print(dim(expr_gene))

write.csv(
  expr_gene,
  "data/processed/GSE154918_gene_expression_healthy_acute_sepsis.csv"
)

# -----------------------------
# 6. Differential expression analysis
# -----------------------------

use_sex <- all(!is.na(analysis_info$sex)) &&
  length(unique(analysis_info$sex)) > 1

if (use_sex) {
  design <- model.matrix(~ sex + clinical_group, data = analysis_info)
} else {
  design <- model.matrix(~ clinical_group, data = analysis_info)
}

colnames(design) <- make.names(colnames(design))

print(colnames(design))

coef_name <- "clinical_groupAcute_sepsis"

if (!coef_name %in% colnames(design)) {
  stop("Coefficient clinical_groupAcute_sepsis not found.")
}

fit <- lmFit(expr_gene, design)
fit <- eBayes(fit)

sepsis_deg <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("SYMBOL") %>%
  mutate(
    regulation = case_when(
      adj.P.Val < 0.05 & logFC >= 0.5 ~ "Up",
      adj.P.Val < 0.05 & logFC <= -0.5 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

write.csv(
  sepsis_deg,
  "results/sepsis/GSE154918_acute_sepsis_vs_healthy_limma_all.csv",
  row.names = FALSE
)

sepsis_sig <- sepsis_deg %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 0.5)

write.csv(
  sepsis_sig,
  "results/sepsis/GSE154918_acute_sepsis_vs_healthy_sig_FDR005_logFC05.csv",
  row.names = FALSE
)

sepsis_deg_summary <- tibble(
  total_tested_genes = nrow(sepsis_deg),
  significant_FDR_005 = sum(sepsis_deg$adj.P.Val < 0.05, na.rm = TRUE),
  significant_FDR_005_logFC05 = nrow(sepsis_sig),
  up_FDR_005_logFC05 = sum(sepsis_sig$regulation == "Up", na.rm = TRUE),
  down_FDR_005_logFC05 = sum(sepsis_sig$regulation == "Down", na.rm = TRUE),
  healthy_n = sum(analysis_info$clinical_group == "Healthy"),
  acute_sepsis_n = sum(analysis_info$clinical_group == "Acute_sepsis")
)

write.csv(
  sepsis_deg_summary,
  "results/sepsis/GSE154918_DEG_summary.csv",
  row.names = FALSE
)

print(sepsis_deg_summary)

# -----------------------------
# 7. PCA plot
# -----------------------------

expr_for_pca <- expr_gene[apply(expr_gene, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

pca <- prcomp(t(expr_for_pca), scale. = TRUE)

pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("expr_col") %>%
  left_join(analysis_info, by = "expr_col")

percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = clinical_group, shape = sepsis_subgroup)
) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 13) +
  labs(
    title = "PCA of GSE154918 blood transcriptomes",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Clinical group",
    shape = "Subgroup"
  )

ggsave(
  "figures/sepsis/Figure7_GSE154918_PCA_healthy_vs_acute_sepsis.png",
  p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 8. Volcano plot
# -----------------------------

volcano_df <- sepsis_deg %>%
  mutate(
    neg_log10_fdr = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    volcano_group = regulation
  )

label_genes <- volcano_df %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 15)

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_fdr)) +
  geom_point(aes(color = volcano_group), alpha = 0.7, size = 1.4) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = label_genes,
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = 30
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Differential expression in acute sepsis blood samples",
    x = "log2 fold change: acute sepsis vs healthy",
    y = "-log10 adjusted P value",
    color = "Regulation"
  )

ggsave(
  "figures/sepsis/Figure7_GSE154918_volcano_acute_sepsis_vs_healthy.png",
  p_volcano,
  width = 7,
  height = 6,
  dpi = 300
)

# -----------------------------
# 9. Compare with GSE193336 LPS-response genes
# -----------------------------

lps_file_candidates <- c(
  "results/DEG/GSE193336_limma_voom_DEG_all.csv",
  "results/GSE193336_limma_voom_DEG_all.csv",
  "GSE193336_limma_voom_DEG_all.csv"
)

lps_file <- lps_file_candidates[file.exists(lps_file_candidates)][1]

if (is.na(lps_file)) {
  stop("Cannot find GSE193336_limma_voom_DEG_all.csv. Please check the DEG result path.")
}

message("Using LPS DEG file: ", lps_file)

lps_deg <- read.csv(lps_file, check.names = FALSE) %>%
  mutate(
    SYMBOL = toupper(SYMBOL),
    LPS_regulation = case_when(
      adj.P.Val < 0.05 & logFC >= 1 ~ "Up",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  ) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::select(
    SYMBOL,
    LPS_logFC = logFC,
    LPS_adjP = adj.P.Val,
    LPS_regulation
  ) %>%
  distinct(SYMBOL, .keep_all = TRUE)

sepsis_deg_simple <- sepsis_deg %>%
  mutate(SYMBOL = toupper(SYMBOL)) %>%
  dplyr::select(
    SYMBOL,
    Sepsis_logFC = logFC,
    Sepsis_adjP = adj.P.Val,
    Sepsis_regulation = regulation
  ) %>%
  distinct(SYMBOL, .keep_all = TRUE)

lps_sepsis_common <- inner_join(
  lps_deg,
  sepsis_deg_simple,
  by = "SYMBOL"
) %>%
  mutate(
    LPS_sig = LPS_regulation %in% c("Up", "Down"),
    Sepsis_sig = Sepsis_regulation %in% c("Up", "Down"),
    both_sig = LPS_sig & Sepsis_sig,
    same_direction = sign(LPS_logFC) == sign(Sepsis_logFC),
    overlap_class = case_when(
      both_sig & same_direction ~ "Both significant, same direction",
      both_sig & !same_direction ~ "Both significant, opposite direction",
      LPS_sig & !Sepsis_sig ~ "LPS significant only",
      !LPS_sig & Sepsis_sig ~ "Sepsis significant only",
      TRUE ~ "Not significant in either"
    )
  )

write.csv(
  lps_sepsis_common,
  "results/sepsis/GSE193336_LPS_vs_GSE154918_sepsis_common_genes.csv",
  row.names = FALSE
)

overlap_sig <- lps_sepsis_common %>%
  filter(both_sig)

write.csv(
  overlap_sig,
  "results/sepsis/GSE193336_LPS_GSE154918_sepsis_significant_overlap.csv",
  row.names = FALSE
)

overlap_summary <- tibble(
  common_tested_genes = nrow(lps_sepsis_common),
  LPS_significant_genes_in_common = sum(lps_sepsis_common$LPS_sig),
  sepsis_significant_genes_in_common = sum(lps_sepsis_common$Sepsis_sig),
  significant_overlap_genes = nrow(overlap_sig),
  same_direction_overlap = sum(overlap_sig$same_direction, na.rm = TRUE),
  opposite_direction_overlap = sum(!overlap_sig$same_direction, na.rm = TRUE)
)

write.csv(
  overlap_summary,
  "results/sepsis/GSE193336_LPS_GSE154918_sepsis_overlap_summary.csv",
  row.names = FALSE
)

print(overlap_summary)

cor_res <- cor.test(
  lps_sepsis_common$LPS_logFC,
  lps_sepsis_common$Sepsis_logFC,
  method = "spearman",
  exact = FALSE
)

cor_summary <- tibble(
  spearman_rho = unname(cor_res$estimate),
  p_value = cor_res$p.value,
  n_common_genes = nrow(lps_sepsis_common)
)

write.csv(
  cor_summary,
  "results/sepsis/GSE193336_LPS_GSE154918_sepsis_logFC_correlation.csv",
  row.names = FALSE
)

print(cor_summary)

# -----------------------------
# 10. Scatter plot
# -----------------------------

scatter_label <- lps_sepsis_common %>%
  filter(both_sig, same_direction) %>%
  arrange(Sepsis_adjP) %>%
  slice_head(n = 15)

p_scatter <- ggplot(
  lps_sepsis_common,
  aes(x = LPS_logFC, y = Sepsis_logFC)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(color = overlap_class), alpha = 0.65, size = 1.5) +
  geom_text_repel(
    data = scatter_label,
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = 30
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Concordance between LPS response and sepsis-associated expression changes",
    subtitle = paste0(
      "Spearman rho = ",
      round(cor_summary$spearman_rho, 3),
      ", P = ",
      signif(cor_summary$p_value, 3)
    ),
    x = "GSE193336 log2FC: LPS vs unstimulated",
    y = "GSE154918 log2FC: acute sepsis vs healthy",
    color = "Gene category"
  )

ggsave(
  "figures/sepsis/Figure8_LPS_sepsis_logFC_concordance_scatter.png",
  p_scatter,
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------
# 11. Heatmap of representative overlapping genes
# -----------------------------

manual_genes <- c(
  "IL1B", "IL6", "TNF", "CXCL8", "CCL2", "CCL4", "CCL5",
  "CXCL10", "PTGS2", "NFKBIA", "SOCS3", "STAT1", "IRF1",
  "MAP3K8", "TNFAIP6", "ICAM1", "NOD2", "JAK3", "LRRK2"
)

top_overlap_genes <- overlap_sig %>%
  filter(same_direction) %>%
  arrange(Sepsis_adjP, LPS_adjP) %>%
  pull(SYMBOL) %>%
  unique()

heatmap_genes <- unique(c(
  intersect(manual_genes, rownames(expr_gene)),
  intersect(top_overlap_genes, rownames(expr_gene))
))

heatmap_genes <- heatmap_genes[1:min(40, length(heatmap_genes))]

if (length(heatmap_genes) >= 5) {
  
  heatmap_mat <- expr_gene[heatmap_genes, analysis_info$expr_col, drop = FALSE]
  heatmap_mat <- heatmap_mat[apply(heatmap_mat, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
  heatmap_mat_z <- t(scale(t(heatmap_mat)))
  
  annotation_col <- analysis_info %>%
    dplyr::select(expr_col, clinical_group, sepsis_subgroup, sex) %>%
    column_to_rownames("expr_col")
  
  png(
    "figures/sepsis/Figure8_GSE154918_overlap_gene_heatmap.png",
    width = 3000,
    height = 2700,
    res = 300
  )
  
  pheatmap(
    heatmap_mat_z,
    annotation_col = annotation_col,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 8,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "LPS-sepsis overlapping genes in GSE154918"
  )
  
  dev.off()
  
} else {
  warning("Too few overlapping genes for heatmap.")
}

sink("results/sepsis/GSE154918_sepsis_relevance_sessionInfo.txt")
sessionInfo()
sink()

message("GSE154918 sepsis relevance analysis completed.")