# ============================================================
# 99_combine_main_figures.R
# Combine existing PNG figures into publication-style multi-panel figures
# ============================================================

rm(list = ls())

# -----------------------------
# 0. Install and load packages
# -----------------------------

packages <- c("ggplot2", "cowplot", "patchwork")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(cowplot)
library(patchwork)

# -----------------------------
# 1. Output directory
# -----------------------------

dir.create("figures/combined", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Helper functions
# -----------------------------

find_one <- function(folder, pattern) {
  files <- list.files(
    folder,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  
  if (length(files) == 0) {
    stop(
      "No file found in folder: ", folder,
      "\nPattern: ", pattern,
      "\nPlease check the filename."
    )
  }
  
  if (length(files) > 1) {
    message("More than one file matched pattern: ", pattern)
    message("Using the first one:")
    message(files[1])
  }
  
  return(files[1])
}

image_panel <- function(path) {
  ggdraw() +
    draw_image(path, scale = 0.98) +
    theme_void()
}

# -----------------------------
# 3. Define figure paths
# -----------------------------

fig_root <- "figures"

# Figure 2: discovery dataset GSE193336
pca_file <- find_one(fig_root, "Figure2_GSE193336_PCA\\.png$")
ma_file <- find_one(fig_root, "Figure3_GSE193336_MAplot\\.png$")
volcano_file <- find_one(fig_root, "Figure3_GSE193336_volcano\\.png$")
heatmap_file <- find_one(fig_root, "Figure4_GSE193336_top50_heatmap\\.png$")

# Figure 3: enrichment
go_gsea_file <- find_one("figures/enrichment_clean", "Figure5_GSEA_GO_BP_clean_top10\\.png$")
kegg_gsea_file <- find_one("figures/enrichment_clean", "Figure5_GSEA_KEGG_clean_top10\\.png$")

# Figure 4: GSE5504 validation
timecourse_file <- find_one("figures/validation", "Figure6_GSE5504_timecourse_selected_genes\\.png$")
validation_heatmap_file <- find_one("figures/validation", "Figure6_GSE5504_validation_heatmap\\.png$")

# Figure 5: GSE154918 sepsis validation
sepsis_pca_file <- find_one("figures/sepsis", "Figure7_GSE154918_PCA.*\\.png$")
sepsis_volcano_file <- find_one("figures/sepsis", "Figure7_GSE154918_volcano.*\\.png$")
concordance_file <- find_one("figures/sepsis", "Figure8A_LPS_sepsis_concordance_scatter_clean\\.png$")
direction_file <- find_one("figures/sepsis", "Figure8B_overlap_direction_barplot\\.png$")
candidate_file <- find_one("figures/sepsis", "Figure8C_top30_candidate_genes_barplot\\.png$")
sepsis_heatmap_file <- find_one("figures/sepsis", "Figure8_GSE154918_overlap_gene_heatmap\\.png$")

# -----------------------------
# 4. Convert image files to panels
# -----------------------------

p_pca <- image_panel(pca_file)
p_ma <- image_panel(ma_file)
p_volcano <- image_panel(volcano_file)
p_heatmap <- image_panel(heatmap_file)

p_gsea_go <- image_panel(go_gsea_file)
p_gsea_kegg <- image_panel(kegg_gsea_file)

p_timecourse <- image_panel(timecourse_file)
p_validation_heatmap <- image_panel(validation_heatmap_file)

p_sepsis_pca <- image_panel(sepsis_pca_file)
p_sepsis_volcano <- image_panel(sepsis_volcano_file)
p_concordance <- image_panel(concordance_file)
p_direction <- image_panel(direction_file)
p_candidate <- image_panel(candidate_file)
p_sepsis_heatmap <- image_panel(sepsis_heatmap_file)

# ============================================================
# Figure 2
# Discovery analysis of GSE193336
# ============================================================

fig2 <- (
  (p_pca | p_ma) /
    (p_volcano | p_heatmap)
) +
  plot_layout(heights = c(1, 1.25)) +
  plot_annotation(
    tag_levels = "A",
    title = "LPS-induced transcriptional remodeling in human macrophages"
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.tag = element_text(face = "bold", size = 18)
  )

ggsave(
  "figures/combined/Figure2_GSE193336_discovery_combined.png",
  fig2,
  width = 15,
  height = 12,
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  "figures/combined/Figure2_GSE193336_discovery_combined.pdf",
  fig2,
  width = 15,
  height = 12,
  limitsize = FALSE
)

# ============================================================
# Figure 3
# Functional enrichment analysis
# ============================================================

fig3 <- (
  p_gsea_go | p_gsea_kegg
) +
  plot_annotation(
    tag_levels = "A",
    title = "Functional enrichment of LPS-responsive transcriptional programs"
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.tag = element_text(face = "bold", size = 18)
  )

ggsave(
  "figures/combined/Figure3_GSEA_enrichment_combined.png",
  fig3,
  width = 15,
  height = 7,
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  "figures/combined/Figure3_GSEA_enrichment_combined.pdf",
  fig3,
  width = 15,
  height = 7,
  limitsize = FALSE
)

# ============================================================
# Figure 4
# GSE5504 time-series validation
# ============================================================

fig4 <- (
  p_timecourse /
    p_validation_heatmap
) +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    tag_levels = "A",
    title = "Independent time-series validation of LPS-responsive genes"
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.tag = element_text(face = "bold", size = 18)
  )

ggsave(
  "figures/combined/Figure4_GSE5504_validation_combined.png",
  fig4,
  width = 15,
  height = 13,
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  "figures/combined/Figure4_GSE5504_validation_combined.pdf",
  fig4,
  width = 15,
  height = 13,
  limitsize = FALSE
)

# ============================================================
# Figure 5
# GSE154918 clinical sepsis validation
# 推荐主文版：不放超大热图，避免太挤
# ============================================================

fig5_main <- (
  p_concordance |
    (p_direction / p_candidate)
) +
  plot_layout(widths = c(1.25, 1), heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    title = "Clinical sepsis validation of LPS-responsive macrophage genes"
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.tag = element_text(face = "bold", size = 18)
  )

ggsave(
  "figures/combined/Figure5_GSE154918_sepsis_validation_main_combined.png",
  fig5_main,
  width = 15,
  height = 9,
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  "figures/combined/Figure5_GSE154918_sepsis_validation_main_combined.pdf",
  fig5_main,
  width = 15,
  height = 9,
  limitsize = FALSE
)

# ============================================================
# Figure S1 or optional Figure 5D
# GSE154918 overlap heatmap
# 建议作为补充图，或者单独主文图
# ============================================================

figS1 <- p_sepsis_heatmap +
  plot_annotation(
    title = "Expression pattern of LPS-sepsis overlapping genes in GSE154918"
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 16)
  )

ggsave(
  "figures/combined/FigureS1_GSE154918_overlap_heatmap.png",
  figS1,
  width = 14,
  height = 12,
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  "figures/combined/FigureS1_GSE154918_overlap_heatmap.pdf",
  figS1,
  width = 14,
  height = 12,
  limitsize = FALSE
)

message("All combined figures saved to figures/combined/")