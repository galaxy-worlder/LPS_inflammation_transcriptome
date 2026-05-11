# ============================================================
# 09_make_highres_combined_figures_from_files.R
# Combine existing figure files into high-resolution manuscript figures
# ============================================================

rm(list = ls())

# ============================================================
# 0. Install and load packages
# ============================================================

options(timeout = 600)

pkgs <- c("magick")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", type = "binary")
  }
}

library(magick)

message("magick configuration:")
print(magick::magick_config())

# ============================================================
# 1. Detect project root
# ============================================================

find_project_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  
  repeat {
    if (dir.exists(file.path(current, "figures"))) {
      return(current)
    }
    
    parent <- dirname(current)
    
    if (identical(parent, current)) {
      stop(
        "Cannot find figures folder.\n",
        "Current working directory is:\n",
        current, "\n\n",
        "Please manually set project_dir in the script, for example:\n",
        "project_dir <- 'D:/流星的文件/科研实践/LPS_inflammation_transcriptome'"
      )
    }
    
    current <- parent
  }
}

# 自动寻找项目根目录
project_dir <- find_project_root()

# 如果自动寻找失败，可以手动取消下一行注释并改成你的项目路径：
# project_dir <- "D:/流星的文件/科研实践/你的项目文件夹名"

message("Project directory detected:")
message(project_dir)

fig_dir <- file.path(project_dir, "figures")

if (!dir.exists(fig_dir)) {
  stop("Cannot find figures folder: ", fig_dir)
}

out_dir <- file.path(fig_dir, "combined_highres")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Figure directory:")
message(fig_dir)

message("Output directory:")
message(out_dir)

# ============================================================
# 2. Define figure source folders
# ============================================================

source_dirs <- list.dirs(fig_dir, recursive = TRUE, full.names = TRUE)
source_dirs <- unique(c(fig_dir, source_dirs))
source_dirs <- source_dirs[dir.exists(source_dirs)]

message("Searching figures in:")
print(source_dirs)

# ============================================================
# 3. Helper functions
# ============================================================

find_figure_file <- function(candidates, source_dirs) {
  for (cand in candidates) {
    hits <- unlist(lapply(source_dirs, function(d) {
      list.files(
        path = d,
        pattern = paste0("^", gsub("\\.", "\\\\.", cand), "$"),
        full.names = TRUE,
        recursive = FALSE,
        ignore.case = FALSE
      )
    }))
    
    if (length(hits) > 0) {
      return(normalizePath(hits[1], winslash = "/", mustWork = TRUE))
    }
  }
  
  all_png_pdf <- unlist(lapply(source_dirs, function(d) {
    list.files(
      path = d,
      pattern = "\\.(png|pdf|jpg|jpeg|tiff)$",
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    )
  }))
  
  stop(
    "Cannot find figure file. Tried candidates:\n",
    paste(candidates, collapse = "\n"),
    "\n\nAvailable figure files include:\n",
    paste(basename(all_png_pdf), collapse = "\n")
  )
}

read_figure_safely <- function(candidates,
                               source_dirs,
                               pdf_density = 600) {
  path <- find_figure_file(candidates, source_dirs)
  
  message("Reading: ", path)
  
  img <- tryCatch(
    {
      if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
        magick::image_read_pdf(path, density = pdf_density)
      } else {
        magick::image_read(path)
      }
    },
    error = function(e) {
      stop(
        "Failed to read figure file:\n",
        path,
        "\n\nError message:\n",
        e$message,
        "\n\nIf this is a PDF-reading issue, install Ghostscript or use the PNG version."
      )
    }
  )
  
  img <- img[1]
  img <- magick::image_background(img, color = "white", flatten = TRUE)
  
  return(img)
}

make_panel <- function(candidates,
                       label,
                       width,
                       height,
                       source_dirs,
                       pdf_density = 600,
                       trim = TRUE,
                       label_size = 90) {
  img <- read_figure_safely(
    candidates = candidates,
    source_dirs = source_dirs,
    pdf_density = pdf_density
  )
  
  if (trim) {
    img <- magick::image_trim(img, fuzz = 8)
  }
  
  img <- magick::image_resize(img, paste0(width, "x", height))
  
  img <- magick::image_extent(
    img,
    geometry = paste0(width, "x", height),
    gravity = "center",
    color = "white"
  )
  
  img <- magick::image_annotate(
    img,
    text = label,
    size = label_size,
    gravity = "northwest",
    location = "+45+35",
    weight = 700,
    color = "black"
  )
  
  return(img)
}

make_row <- function(img_list) {
  magick::image_append(do.call(c, img_list), stack = FALSE)
}

make_column <- function(img_list) {
  magick::image_append(do.call(c, img_list), stack = TRUE)
}

add_main_title <- function(img,
                           title,
                           subtitle = NULL,
                           title_height = 260,
                           title_size = 95,
                           subtitle_size = 58) {
  info <- magick::image_info(img)
  
  title_img <- magick::image_blank(
    width = info$width[1],
    height = title_height,
    color = "white"
  )
  
  title_img <- magick::image_annotate(
    title_img,
    text = title,
    size = title_size,
    gravity = "northwest",
    location = "+45+25",
    weight = 700,
    color = "black"
  )
  
  if (!is.null(subtitle)) {
    title_img <- magick::image_annotate(
      title_img,
      text = subtitle,
      size = subtitle_size,
      gravity = "northwest",
      location = "+48+145",
      color = "black"
    )
  }
  
  magick::image_append(c(title_img, img), stack = TRUE)
}

write_combined_figure <- function(img,
                                  filename_prefix,
                                  out_dir) {
  png_file <- file.path(out_dir, paste0(filename_prefix, ".png"))
  pdf_file <- file.path(out_dir, paste0(filename_prefix, ".pdf"))
  
  magick::image_write(img, path = png_file, format = "png")
  magick::image_write(img, path = pdf_file, format = "pdf")
  
  message("Saved PNG: ", png_file)
  message("Saved PDF: ", pdf_file)
}

# ============================================================
# 4. Candidate file names
# ============================================================
# 如果你的某张图文件名不同，只需要修改这里。
# 默认优先找 PDF，其次找 PNG。

fig_files <- list(
  
  # GSE193336 discovery
  gse193336_pca = c(
    "Figure2_GSE193336_PCA.pdf",
    "Figure2_GSE193336_PCA.png"
  ),
  
  gse193336_ma = c(
    "Figure3_GSE193336_MAplot.pdf",
    "Figure3_GSE193336_MAplot.png"
  ),
  
  gse193336_volcano = c(
    "Figure3_GSE193336_volcano.pdf",
    "Figure3_GSE193336_volcano.png"
  ),
  
  gse193336_heatmap = c(
    "Figure4_GSE193336_top50_heatmap.pdf",
    "Figure4_GSE193336_top50_heatmap.png"
  ),
  
  # Enrichment clean figures
  gsea_go_bp = c(
    "Figure5_GSEA_GO_BP_clean_top10.pdf",
    "Figure5_GSEA_GO_BP_clean_top10.png",
    "Figure5_GSEA_GO_BP_dotplot.pdf",
    "Figure5_GSEA_GO_BP_dotplot.png"
  ),
  
  gsea_kegg = c(
    "Figure5_GSEA_KEGG_clean_top10.pdf",
    "Figure5_GSEA_KEGG_clean_top10.png",
    "Figure5_GSEA_KEGG_dotplot.pdf",
    "Figure5_GSEA_KEGG_dotplot.png"
  ),
  
  # GSE5504 validation
  gse5504_timecourse = c(
    "Figure6_GSE5504_timecourse_selected_genes.pdf",
    "Figure6_GSE5504_timecourse_selected_genes.png"
  ),
  
  gse5504_heatmap = c(
    "Figure6_GSE5504_validation_heatmap.pdf",
    "Figure6_GSE5504_validation_heatmap.png"
  ),
  
  # GSE154918 sepsis
  gse154918_pca = c(
    "Figure7_GSE154918_PCA_healthy_vs_acute_sepsis.pdf",
    "Figure7_GSE154918_PCA_healthy_vs_acute_sepsis.png",
    "Figure7_GSE154918_PCA.png"
  ),
  
  gse154918_volcano = c(
    "Figure7_GSE154918_volcano_acute_sepsis_vs_healthy.pdf",
    "Figure7_GSE154918_volcano_acute_sepsis_vs_healthy.png",
    "Figure7_GSE154918_volcano.png"
  ),
  
  concordance_scatter = c(
    "Figure8A_LPS_sepsis_concordance_scatter_clean.pdf",
    "Figure8A_LPS_sepsis_concordance_scatter_clean.png",
    "Figure8_LPS_sepsis_logFC_concordance_scatter.pdf",
    "Figure8_LPS_sepsis_logFC_concordance_scatter.png"
  ),
  
  direction_barplot = c(
    "Figure8B_overlap_direction_barplot.pdf",
    "Figure8B_overlap_direction_barplot.png"
  ),
  
  candidate_barplot = c(
    "Figure8C_top30_candidate_genes_barplot.pdf",
    "Figure8C_top30_candidate_genes_barplot.png"
  ),
  
  sepsis_overlap_heatmap = c(
    "Figure8_GSE154918_overlap_gene_heatmap.pdf",
    "Figure8_GSE154918_overlap_gene_heatmap.png",
    "FigureS1_GSE154918_overlap_heatmap.pdf",
    "FigureS1_GSE154918_overlap_heatmap.png"
  )
)

# ============================================================
# 5. Figure 2
# GSE193336 discovery analysis
# ============================================================

message("\nMaking Figure 2...")

p2_A <- make_panel(
  fig_files$gse193336_pca,
  label = "A",
  width = 2200,
  height = 1700,
  source_dirs = source_dirs
)

p2_B <- make_panel(
  fig_files$gse193336_ma,
  label = "B",
  width = 2200,
  height = 1700,
  source_dirs = source_dirs
)

p2_C <- make_panel(
  fig_files$gse193336_volcano,
  label = "C",
  width = 2200,
  height = 1700,
  source_dirs = source_dirs
)

p2_D <- make_panel(
  fig_files$gse193336_heatmap,
  label = "D",
  width = 6600,
  height = 4300,
  source_dirs = source_dirs
)

fig2_top <- make_row(list(p2_A, p2_B, p2_C))
fig2_body <- make_column(list(fig2_top, p2_D))

fig2 <- add_main_title(
  fig2_body,
  title = "LPS-induced transcriptional remodeling in human macrophages",
  subtitle = "Discovery analysis in GSE193336"
)

write_combined_figure(
  fig2,
  filename_prefix = "Figure2_GSE193336_discovery_combined_highres",
  out_dir = out_dir
)

# ============================================================
# 6. Figure 3
# Functional enrichment
# ============================================================

message("\nMaking Figure 3...")

p3_A <- make_panel(
  fig_files$gsea_go_bp,
  label = "A",
  width = 3300,
  height = 4200,
  source_dirs = source_dirs
)

p3_B <- make_panel(
  fig_files$gsea_kegg,
  label = "B",
  width = 3300,
  height = 4200,
  source_dirs = source_dirs
)

fig3_body <- make_row(list(p3_A, p3_B))

fig3 <- add_main_title(
  fig3_body,
  title = "Functional enrichment of LPS-responsive transcriptional programs",
  subtitle = "Gene set enrichment analysis of GO Biological Process and KEGG pathways"
)

write_combined_figure(
  fig3,
  filename_prefix = "Figure3_GSEA_enrichment_combined_highres",
  out_dir = out_dir
)

# ============================================================
# 7. Figure 4
# GSE5504 validation
# ============================================================

message("\nMaking Figure 4...")

p4_A <- make_panel(
  fig_files$gse5504_timecourse,
  label = "A",
  width = 6600,
  height = 3600,
  source_dirs = source_dirs
)

p4_B <- make_panel(
  fig_files$gse5504_heatmap,
  label = "B",
  width = 6600,
  height = 4300,
  source_dirs = source_dirs
)

fig4_body <- make_column(list(p4_A, p4_B))

fig4 <- add_main_title(
  fig4_body,
  title = "Independent time-series validation of LPS-responsive genes",
  subtitle = "Temporal expression patterns in GSE5504"
)

write_combined_figure(
  fig4,
  filename_prefix = "Figure4_GSE5504_validation_combined_highres",
  out_dir = out_dir
)

# ============================================================
# 8. Figure 5
# Sepsis validation and LPS-sepsis concordance
# ============================================================

message("\nMaking Figure 5...")

p5_A <- make_panel(
  fig_files$gse154918_pca,
  label = "A",
  width = 3300,
  height = 2500,
  source_dirs = source_dirs
)

p5_B <- make_panel(
  fig_files$gse154918_volcano,
  label = "B",
  width = 3300,
  height = 2500,
  source_dirs = source_dirs
)

p5_C <- make_panel(
  fig_files$concordance_scatter,
  label = "C",
  width = 3300,
  height = 3000,
  source_dirs = source_dirs
)

p5_D <- make_panel(
  fig_files$direction_barplot,
  label = "D",
  width = 3300,
  height = 3000,
  source_dirs = source_dirs
)

p5_E <- make_panel(
  fig_files$candidate_barplot,
  label = "E",
  width = 6600,
  height = 3600,
  source_dirs = source_dirs
)

fig5_row1 <- make_row(list(p5_A, p5_B))
fig5_row2 <- make_row(list(p5_C, p5_D))
fig5_body <- make_column(list(fig5_row1, fig5_row2, p5_E))

fig5 <- add_main_title(
  fig5_body,
  title = "LPS-sepsis transcriptional concordance",
  subtitle = "Cross-dataset comparison between LPS-stimulated macrophages and acute sepsis"
)

write_combined_figure(
  fig5,
  filename_prefix = "Figure5_GSE154918_sepsis_concordance_combined_highres",
  out_dir = out_dir
)

# ============================================================
# 9. Supplementary Figure S1
# Full overlap heatmap
# ============================================================

message("\nMaking Supplementary Figure S1...")

figS1_panel <- make_panel(
  fig_files$sepsis_overlap_heatmap,
  label = "A",
  width = 5200,
  height = 7200,
  source_dirs = source_dirs
)

figS1 <- add_main_title(
  figS1_panel,
  title = "Expression pattern of LPS-sepsis overlapping genes in GSE154918",
  subtitle = "Heatmap of representative overlapping genes across clinical groups"
)

write_combined_figure(
  figS1,
  filename_prefix = "FigureS1_GSE154918_overlap_heatmap_highres",
  out_dir = out_dir
)

# ============================================================
# 10. Finish
# ============================================================

message("\nAll combined high-resolution figures finished.")
message("Please check this folder:")
message(out_dir)

sessionInfo()