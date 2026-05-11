# scripts/04b_clean_enrichment_figures.R

# ============================================================
# Clean enrichment figures for manuscript
# ============================================================

library(tidyverse)
library(ggplot2)
library(stringr)

dir.create("figures/enrichment_clean", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Load GSEA KEGG result
# -----------------------------

gsea_kegg <- read.csv(
  "results/enrichment/GSE193336_GSEA_KEGG.csv",
  check.names = FALSE
)

# -----------------------------
# 2. Select top activated and suppressed pathways
# -----------------------------

top_kegg_activated <- gsea_kegg %>%
  filter(NES > 0) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(direction = "Activated")

top_kegg_suppressed <- gsea_kegg %>%
  filter(NES < 0) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(direction = "Suppressed")

top_kegg_plot <- bind_rows(top_kegg_activated, top_kegg_suppressed) %>%
  mutate(
    Description_wrapped = str_wrap(Description, width = 35),
    Description_wrapped = factor(
      Description_wrapped,
      levels = rev(unique(Description_wrapped[order(NES)]))
    )
  )

write.csv(
  top_kegg_plot,
  "results/enrichment/GSE193336_GSEA_KEGG_top10_activated_suppressed.csv",
  row.names = FALSE
)

# -----------------------------
# 3. Plot clean KEGG GSEA result
# -----------------------------

p_kegg_clean <- ggplot(
  top_kegg_plot,
  aes(x = NES, y = Description_wrapped)
) +
  geom_point(aes(size = setSize, color = p.adjust)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ direction, scales = "free_y") +
  scale_color_continuous(trans = "reverse") +
  theme_bw(base_size = 13) +
  labs(
    title = "GSEA of KEGG pathways in LPS-stimulated macrophages",
    x = "Normalized enrichment score",
    y = NULL,
    size = "Gene set size",
    color = "Adjusted P value"
  )

ggsave(
  "figures/enrichment_clean/Figure5_GSEA_KEGG_clean_top10.pdf",
  p_kegg_clean,
  width = 11,
  height = 7
)

ggsave(
  "figures/enrichment_clean/Figure5_GSEA_KEGG_clean_top10.png",
  p_kegg_clean,
  width = 11,
  height = 7,
  dpi = 300
)

# -----------------------------
# 4. Load GSEA GO BP result
# -----------------------------

gsea_go <- read.csv(
  "results/enrichment/GSE193336_GSEA_GO_BP.csv",
  check.names = FALSE
)

top_go_activated <- gsea_go %>%
  filter(NES > 0) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(direction = "Activated")

top_go_suppressed <- gsea_go %>%
  filter(NES < 0) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(direction = "Suppressed")

top_go_plot <- bind_rows(top_go_activated, top_go_suppressed) %>%
  mutate(
    Description_wrapped = str_wrap(Description, width = 35),
    Description_wrapped = factor(
      Description_wrapped,
      levels = rev(unique(Description_wrapped[order(NES)]))
    )
  )

write.csv(
  top_go_plot,
  "results/enrichment/GSE193336_GSEA_GO_BP_top10_activated_suppressed.csv",
  row.names = FALSE
)

p_go_clean <- ggplot(
  top_go_plot,
  aes(x = NES, y = Description_wrapped)
) +
  geom_point(aes(size = setSize, color = p.adjust)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ direction, scales = "free_y") +
  scale_color_continuous(trans = "reverse") +
  theme_bw(base_size = 13) +
  labs(
    title = "GSEA of GO Biological Process terms",
    x = "Normalized enrichment score",
    y = NULL,
    size = "Gene set size",
    color = "Adjusted P value"
  )

ggsave(
  "figures/enrichment_clean/Figure5_GSEA_GO_BP_clean_top10.pdf",
  p_go_clean,
  width = 11,
  height = 7
)

ggsave(
  "figures/enrichment_clean/Figure5_GSEA_GO_BP_clean_top10.png",
  p_go_clean,
  width = 11,
  height = 7,
  dpi = 300
)

message("Clean enrichment figures generated.")