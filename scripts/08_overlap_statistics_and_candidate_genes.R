# scripts/08_overlap_statistics_and_candidate_genes.R

library(tidyverse)
library(ggrepel)

dir.create("results/sepsis", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/sepsis", recursive = TRUE, showWarnings = FALSE)

common_file <- "results/sepsis/GSE193336_LPS_vs_GSE154918_sepsis_common_genes.csv"

if (!file.exists(common_file)) {
  stop("Cannot find common gene file: results/sepsis/GSE193336_LPS_vs_GSE154918_sepsis_common_genes.csv")
}

common <- read.csv(common_file, check.names = FALSE)

common <- common %>%
  mutate(
    LPS_sig = as.logical(LPS_sig),
    Sepsis_sig = as.logical(Sepsis_sig),
    both_sig = as.logical(both_sig),
    same_direction = as.logical(same_direction),
    direction_pair = case_when(
      LPS_logFC > 0 & Sepsis_logFC > 0 ~ "Up in both",
      LPS_logFC < 0 & Sepsis_logFC < 0 ~ "Down in both",
      LPS_logFC > 0 & Sepsis_logFC < 0 ~ "LPS up, sepsis down",
      LPS_logFC < 0 & Sepsis_logFC > 0 ~ "LPS down, sepsis up",
      TRUE ~ "Other"
    )
  )

# -----------------------------
# 1. Fisher exact test for overlap enrichment
# -----------------------------

n_both <- sum(common$LPS_sig & common$Sepsis_sig, na.rm = TRUE)
n_lps_only <- sum(common$LPS_sig & !common$Sepsis_sig, na.rm = TRUE)
n_sepsis_only <- sum(!common$LPS_sig & common$Sepsis_sig, na.rm = TRUE)
n_neither <- sum(!common$LPS_sig & !common$Sepsis_sig, na.rm = TRUE)

fisher_mat <- matrix(
  c(n_both, n_lps_only, n_sepsis_only, n_neither),
  nrow = 2,
  byrow = TRUE
)

rownames(fisher_mat) <- c("LPS_sig", "LPS_not_sig")
colnames(fisher_mat) <- c("Sepsis_sig", "Sepsis_not_sig")

fisher_res <- fisher.test(fisher_mat, alternative = "greater")

fisher_summary <- tibble(
  common_tested_genes = nrow(common),
  lps_significant = sum(common$LPS_sig, na.rm = TRUE),
  sepsis_significant = sum(common$Sepsis_sig, na.rm = TRUE),
  observed_overlap = n_both,
  expected_overlap = sum(common$LPS_sig, na.rm = TRUE) *
    sum(common$Sepsis_sig, na.rm = TRUE) / nrow(common),
  odds_ratio = unname(fisher_res$estimate),
  fisher_p_value = fisher_res$p.value
)

write.csv(
  fisher_summary,
  "results/sepsis/Figure8_overlap_fisher_enrichment_summary.csv",
  row.names = FALSE
)

print(fisher_summary)

# -----------------------------
# 2. Direction summary
# -----------------------------

direction_summary <- common %>%
  filter(both_sig) %>%
  count(direction_pair, name = "n_genes") %>%
  mutate(
    percentage = round(100 * n_genes / sum(n_genes), 2)
  ) %>%
  arrange(desc(n_genes))

write.csv(
  direction_summary,
  "results/sepsis/Figure8_overlap_direction_summary.csv",
  row.names = FALSE
)

print(direction_summary)

p_direction <- ggplot(
  direction_summary,
  aes(x = reorder(direction_pair, n_genes), y = n_genes)
) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 13) +
  labs(
    title = "Direction of shared LPS-sepsis genes",
    x = NULL,
    y = "Number of genes"
  )

ggsave(
  "figures/sepsis/Figure8B_overlap_direction_barplot.png",
  p_direction,
  width = 7,
  height = 4.5,
  dpi = 300
)

# -----------------------------
# 3. Candidate gene ranking
# -----------------------------

candidate_genes <- common %>%
  filter(both_sig, same_direction) %>%
  mutate(
    LPS_neglog10FDR = -log10(pmax(LPS_adjP, .Machine$double.xmin)),
    Sepsis_neglog10FDR = -log10(pmax(Sepsis_adjP, .Machine$double.xmin)),
    candidate_score =
      abs(LPS_logFC) +
      abs(Sepsis_logFC) +
      0.1 * LPS_neglog10FDR +
      0.1 * Sepsis_neglog10FDR
  ) %>%
  arrange(desc(candidate_score)) %>%
  dplyr::select(
    SYMBOL,
    LPS_logFC,
    LPS_adjP,
    Sepsis_logFC,
    Sepsis_adjP,
    direction_pair,
    candidate_score
  )

write.csv(
  candidate_genes,
  "results/sepsis/Figure8_candidate_same_direction_genes_ranked.csv",
  row.names = FALSE
)

top30_candidates <- candidate_genes %>%
  slice_head(n = 30)

write.csv(
  top30_candidates,
  "results/sepsis/Figure8_top30_candidate_genes.csv",
  row.names = FALSE
)

print(top30_candidates)

# -----------------------------
# 4. Cleaner scatter plot
# -----------------------------

plot_df <- common %>%
  mutate(
    plot_class = case_when(
      both_sig & same_direction ~ "Both significant, same direction",
      both_sig & !same_direction ~ "Both significant, opposite direction",
      LPS_sig & !Sepsis_sig ~ "LPS significant only",
      !LPS_sig & Sepsis_sig ~ "Sepsis significant only",
      TRUE ~ "Not significant in either"
    )
  )

label_df <- candidate_genes %>%
  slice_head(n = 15)

p_scatter_clean <- ggplot(
  plot_df,
  aes(x = LPS_logFC, y = Sepsis_logFC)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(color = plot_class), alpha = 0.55, size = 1.3) +
  geom_text_repel(
    data = label_df,
    aes(x = LPS_logFC, y = Sepsis_logFC, label = SYMBOL),
    size = 3,
    max.overlaps = 30
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "LPS-sepsis transcriptional concordance",
    subtitle = "Shared genes between GSE193336 and GSE154918",
    x = "LPS response log2FC",
    y = "Acute sepsis log2FC",
    color = "Gene category"
  )

ggsave(
  "figures/sepsis/Figure8A_LPS_sepsis_concordance_scatter_clean.png",
  p_scatter_clean,
  width = 8.5,
  height = 6,
  dpi = 300
)

# -----------------------------
# 5. Top candidate gene barplot
# -----------------------------

p_candidate <- top30_candidates %>%
  mutate(SYMBOL = factor(SYMBOL, levels = rev(SYMBOL))) %>%
  ggplot(aes(x = SYMBOL, y = candidate_score)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top candidate LPS-sepsis concordant genes",
    x = NULL,
    y = "Candidate score"
  )

ggsave(
  "figures/sepsis/Figure8C_top30_candidate_genes_barplot.png",
  p_candidate,
  width = 7,
  height = 7,
  dpi = 300
)

message("Overlap statistics and candidate gene prioritization completed.")