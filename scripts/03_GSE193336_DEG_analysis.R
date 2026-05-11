# scripts/03_GSE193336_DEG_analysis.R

# ============================================================
# Differential expression analysis for GSE193336
# LPS-stimulated human macrophages vs unstimulated controls
# Method: limma-voom with donor-paired design
# ============================================================

library(tidyverse)
library(data.table)
library(janitor)
library(edgeR)
library(limma)
library(pheatmap)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)

select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
mutate <- dplyr::mutate

# -----------------------------
# 1. Create output folders
# -----------------------------

dir.create("results/DEG", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Load metadata
# -----------------------------

metadata <- read.csv("data/metadata/GSE193336_metadata.csv", check.names = FALSE)

# Check key columns
required_cols <- c("title", "geo_accession", "donor_ch1", "treatment_ch1")
missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

sample_info <- metadata %>%
  transmute(
    sample_name = title,
    geo_accession = geo_accession,
    donor = donor_ch1,
    treatment_raw = treatment_ch1,
    condition = ifelse(grepl("LPS", treatment_ch1, ignore.case = TRUE),
                       "LPS", "Unstimulated")
  ) %>%
  mutate(
    donor = factor(donor),
    condition = factor(condition, levels = c("Unstimulated", "LPS"))
  )

write.csv(
  sample_info,
  "data/processed/GSE193336_sample_info_clean.csv",
  row.names = FALSE
)

print(sample_info)

# -----------------------------
# 3. Locate and load count matrix
# -----------------------------

supp_dir <- "data/supplementary/GSE193336"

count_files <- list.files(
  supp_dir,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.csv$|\\.csv\\.gz$|\\.txt$|\\.txt\\.gz$"
)

count_files <- count_files[
  grepl("count|Count|Raw|raw", basename(count_files), ignore.case = TRUE)
]

if (length(count_files) == 0) {
  stop("No count file found in data/supplementary/GSE193336. Please check supplementary files.")
}

count_file <- count_files[1]

message("Using count file: ", count_file)

count_df <- fread(count_file)

# Preview
print(dim(count_df))
print(head(count_df))

# -----------------------------
# 4. Prepare count matrix
# -----------------------------

# The first column should be gene ID
gene_col <- colnames(count_df)[1]

count_mat <- count_df %>%
  as.data.frame()

rownames(count_mat) <- count_mat[[gene_col]]
count_mat[[gene_col]] <- NULL

# Make sure sample names match metadata titles
count_sample_names <- colnames(count_mat)
meta_sample_names <- sample_info$sample_name

if (!all(meta_sample_names %in% count_sample_names)) {
  cat("Metadata sample names:\n")
  print(meta_sample_names)
  cat("Count matrix sample names:\n")
  print(count_sample_names)
  stop("Some metadata sample names are not found in count matrix.")
}

# Reorder count matrix columns according to metadata
count_mat <- count_mat[, meta_sample_names]

# Convert to numeric matrix
count_mat <- as.matrix(count_mat)
mode(count_mat) <- "numeric"

# Remove Ensembl version if present
rownames(count_mat) <- sub("\\..*$", "", rownames(count_mat))

# Remove duplicated gene IDs if any
dup_genes <- duplicated(rownames(count_mat))
if (any(dup_genes)) {
  message("Duplicated gene IDs detected. Keeping the row with the highest mean expression.")
  
  count_mat <- as.data.frame(count_mat) %>%
    rownames_to_column("ensembl_id") %>%
    mutate(mean_expr = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>%
    arrange(ensembl_id, desc(mean_expr)) %>%
    group_by(ensembl_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(-mean_expr) %>%
    column_to_rownames("ensembl_id") %>%
    as.matrix()
}

# Save cleaned count matrix
write.csv(
  count_mat,
  "data/processed/GSE193336_count_matrix_clean.csv"
)

# -----------------------------
# 5. Create DGEList and filter genes
# -----------------------------

dge <- DGEList(counts = count_mat)

design <- model.matrix(~ donor + condition, data = sample_info)
colnames(design) <- make.names(colnames(design))

print(design)

keep <- filterByExpr(dge, design = design)

cat("Number of genes before filtering:", nrow(dge), "\n")
cat("Number of genes after filtering:", sum(keep), "\n")

dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge, method = "TMM")

# -----------------------------
# 6. limma-voom differential expression
# -----------------------------

v <- voom(dge, design = design, plot = FALSE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

# The LPS effect coefficient should be conditionLPS
print(colnames(design))

coef_name <- "conditionLPS"

if (!coef_name %in% colnames(design)) {
  stop("Coefficient conditionLPS not found. Check design matrix column names.")
}

deg <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  sort.by = "P"
)

deg <- deg %>%
  rownames_to_column("ensembl_id")

# -----------------------------
# 7. Gene annotation
# -----------------------------

anno <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = deg$ensembl_id,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
)

anno <- anno %>%
  distinct(ENSEMBL, .keep_all = TRUE)

deg_anno <- deg %>%
  left_join(anno, by = c("ensembl_id" = "ENSEMBL")) %>%
  mutate(
    gene_label = ifelse(is.na(SYMBOL) | SYMBOL == "", ensembl_id, SYMBOL),
    regulation = case_when(
      adj.P.Val < 0.05 & logFC >= 1  ~ "Up",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

write.csv(
  deg_anno,
  "results/DEG/GSE193336_limma_voom_DEG_all.csv",
  row.names = FALSE
)

sig_deg <- deg_anno %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)

write.csv(
  sig_deg,
  "results/DEG/GSE193336_limma_voom_DEG_sig_FDR005_logFC1.csv",
  row.names = FALSE
)

deg_summary <- tibble(
  total_tested_genes = nrow(deg_anno),
  significant_FDR_005 = sum(deg_anno$adj.P.Val < 0.05, na.rm = TRUE),
  significant_FDR_005_logFC1 = nrow(sig_deg),
  up_FDR_005_logFC1 = sum(sig_deg$regulation == "Up", na.rm = TRUE),
  down_FDR_005_logFC1 = sum(sig_deg$regulation == "Down", na.rm = TRUE)
)

write.csv(
  deg_summary,
  "results/DEG/GSE193336_DEG_summary.csv",
  row.names = FALSE
)

print(deg_summary)

# -----------------------------
# 8. PCA plot
# -----------------------------

log_cpm <- cpm(dge, log = TRUE, prior.count = 1)

pca <- prcomp(t(log_cpm), scale. = TRUE)

pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample_name") %>%
  left_join(sample_info, by = "sample_name")

percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = donor)) +
  geom_point(size = 4) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA of GSE193336 macrophage RNA-seq samples",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition",
    shape = "Donor"
  )

ggsave(
  "figures/Figure2_GSE193336_PCA.pdf",
  p_pca,
  width = 7,
  height = 5
)

ggsave(
  "figures/Figure2_GSE193336_PCA.png",
  p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 9. Volcano plot
# -----------------------------

volcano_df <- deg_anno %>%
  mutate(
    neg_log10_fdr = -log10(adj.P.Val),
    volcano_group = case_when(
      adj.P.Val < 0.05 & logFC >= 1  ~ "Up",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

label_genes <- volcano_df %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 12)

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_fdr)) +
  geom_point(aes(color = volcano_group), alpha = 0.7, size = 1.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = label_genes,
    aes(label = gene_label),
    size = 3,
    max.overlaps = 30
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Differential expression induced by LPS in macrophages",
    x = "log2 fold change: LPS vs unstimulated",
    y = "-log10 adjusted P value",
    color = "Regulation"
  )

ggsave(
  "figures/Figure3_GSE193336_volcano.pdf",
  p_volcano,
  width = 7,
  height = 6
)

ggsave(
  "figures/Figure3_GSE193336_volcano.png",
  p_volcano,
  width = 7,
  height = 6,
  dpi = 300
)

# -----------------------------
# 10. MA plot
# -----------------------------

ma_df <- deg_anno %>%
  mutate(
    group = case_when(
      adj.P.Val < 0.05 & logFC >= 1  ~ "Up",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

p_ma <- ggplot(ma_df, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = group), alpha = 0.7, size = 1.4) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    title = "MA plot of LPS-induced differential expression",
    x = "Average expression",
    y = "log2 fold change: LPS vs unstimulated",
    color = "Regulation"
  )

ggsave(
  "figures/Figure3_GSE193336_MAplot.pdf",
  p_ma,
  width = 7,
  height = 5
)

ggsave(
  "figures/Figure3_GSE193336_MAplot.png",
  p_ma,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 11. Heatmap of top DE genes
# -----------------------------

top_heatmap_genes <- deg_anno %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  pull(ensembl_id)

heatmap_mat <- log_cpm[top_heatmap_genes, ]

# Replace rownames with gene symbols where available
gene_labels <- deg_anno %>%
  filter(ensembl_id %in% top_heatmap_genes) %>%
  select(ensembl_id, gene_label)

rownames(heatmap_mat) <- gene_labels$gene_label[
  match(rownames(heatmap_mat), gene_labels$ensembl_id)
]

annotation_col <- sample_info %>%
  select(sample_name, condition, donor) %>%
  column_to_rownames("sample_name")

pdf("figures/Figure4_GSE193336_top50_heatmap.pdf", width = 8, height = 10)

pheatmap(
  heatmap_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7,
  main = "Top 50 LPS-responsive genes in macrophages"
)

dev.off()

png("figures/Figure4_GSE193336_top50_heatmap.png", width = 2400, height = 3000, res = 300)

pheatmap(
  heatmap_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7,
  main = "Top 50 LPS-responsive genes in macrophages"
)

dev.off()

# -----------------------------
# 12. Save normalized expression
# -----------------------------

log_cpm_df <- as.data.frame(log_cpm) %>%
  rownames_to_column("ensembl_id")

write.csv(
  log_cpm_df,
  "data/processed/GSE193336_logCPM_TMM.csv",
  row.names = FALSE
)

# -----------------------------
# 13. Save session information
# -----------------------------

sink("results/DEG/GSE193336_sessionInfo.txt")
sessionInfo()
sink()

message("GSE193336 limma-voom differential expression analysis completed.")