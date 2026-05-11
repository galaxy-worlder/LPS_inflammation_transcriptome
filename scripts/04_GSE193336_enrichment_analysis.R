# scripts/04_GSE193336_enrichment_analysis.R

# ============================================================
# Functional enrichment analysis for GSE193336 DEGs
# GO, KEGG, and GSEA
# ============================================================

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Avoid function conflicts
select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
mutate <- dplyr::mutate

# -----------------------------
# 1. Create output folders
# -----------------------------

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/enrichment", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Load DEG results
# -----------------------------

deg_all <- read.csv(
  "results/DEG/GSE193336_limma_voom_DEG_all.csv",
  check.names = FALSE
)

required_cols <- c("ensembl_id", "logFC", "adj.P.Val", "SYMBOL", "ENTREZID")
missing_cols <- setdiff(required_cols, colnames(deg_all))

if (length(missing_cols) > 0) {
  stop("Missing columns in DEG result: ", paste(missing_cols, collapse = ", "))
}

# Remove genes without ENTREZ ID
deg_entrez <- deg_all %>%
  filter(!is.na(ENTREZID), ENTREZID != "") %>%
  mutate(
    ENTREZID = as.character(ENTREZID),
    regulation = case_when(
      adj.P.Val < 0.05 & logFC >= 1  ~ "Up",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

# If multiple Ensembl IDs map to the same Entrez ID,
# keep the gene with the strongest absolute logFC.
deg_entrez_unique <- deg_entrez %>%
  arrange(ENTREZID, desc(abs(logFC))) %>%
  group_by(ENTREZID) %>%
  slice(1) %>%
  ungroup()

background_genes <- unique(deg_entrez_unique$ENTREZID)

up_genes <- deg_entrez_unique %>%
  filter(regulation == "Up") %>%
  pull(ENTREZID) %>%
  unique()

down_genes <- deg_entrez_unique %>%
  filter(regulation == "Down") %>%
  pull(ENTREZID) %>%
  unique()

all_sig_genes <- c(up_genes, down_genes) %>% unique()

gene_number_summary <- tibble(
  background_genes = length(background_genes),
  significant_genes = length(all_sig_genes),
  upregulated_genes = length(up_genes),
  downregulated_genes = length(down_genes)
)

write.csv(
  gene_number_summary,
  "results/enrichment/GSE193336_enrichment_gene_number_summary.csv",
  row.names = FALSE
)

print(gene_number_summary)

# -----------------------------
# 3. GO Biological Process enrichment
# -----------------------------

ego_up_bp <- enrichGO(
  gene = up_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_down_bp <- enrichGO(
  gene = down_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

write.csv(
  as.data.frame(ego_up_bp),
  "results/enrichment/GSE193336_GO_BP_up.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_down_bp),
  "results/enrichment/GSE193336_GO_BP_down.csv",
  row.names = FALSE
)

# GO BP dotplot for upregulated genes
p_go_up <- dotplot(ego_up_bp, showCategory = 20) +
  ggtitle("GO Biological Process enrichment of LPS-upregulated genes")

ggsave(
  "figures/enrichment/Figure5_GO_BP_up_dotplot.pdf",
  p_go_up,
  width = 9,
  height = 7
)

ggsave(
  "figures/enrichment/Figure5_GO_BP_up_dotplot.png",
  p_go_up,
  width = 9,
  height = 7,
  dpi = 300
)

# GO BP dotplot for downregulated genes
p_go_down <- dotplot(ego_down_bp, showCategory = 20) +
  ggtitle("GO Biological Process enrichment of LPS-downregulated genes")

ggsave(
  "figures/enrichment/Figure5_GO_BP_down_dotplot.pdf",
  p_go_down,
  width = 9,
  height = 7
)

ggsave(
  "figures/enrichment/Figure5_GO_BP_down_dotplot.png",
  p_go_down,
  width = 9,
  height = 7,
  dpi = 300
)

# -----------------------------
# 4. KEGG pathway enrichment
# -----------------------------

# KEGG sometimes requires internet access.
# Use tryCatch so that GO/GSEA results are not affected if KEGG fails.

ekegg_up <- tryCatch(
  {
    enrichKEGG(
      gene = up_genes,
      universe = background_genes,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2
    )
  },
  error = function(e) {
    message("KEGG enrichment for upregulated genes failed: ", e$message)
    NULL
  }
)

ekegg_down <- tryCatch(
  {
    enrichKEGG(
      gene = down_genes,
      universe = background_genes,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2
    )
  },
  error = function(e) {
    message("KEGG enrichment for downregulated genes failed: ", e$message)
    NULL
  }
)

if (!is.null(ekegg_up) && nrow(as.data.frame(ekegg_up)) > 0) {
  ekegg_up_readable <- setReadable(
    ekegg_up,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )
  
  write.csv(
    as.data.frame(ekegg_up_readable),
    "results/enrichment/GSE193336_KEGG_up.csv",
    row.names = FALSE
  )
  
  p_kegg_up <- dotplot(ekegg_up_readable, showCategory = 20) +
    ggtitle("KEGG enrichment of LPS-upregulated genes")
  
  ggsave(
    "figures/enrichment/Figure5_KEGG_up_dotplot.pdf",
    p_kegg_up,
    width = 9,
    height = 7
  )
  
  ggsave(
    "figures/enrichment/Figure5_KEGG_up_dotplot.png",
    p_kegg_up,
    width = 9,
    height = 7,
    dpi = 300
  )
}

if (!is.null(ekegg_down) && nrow(as.data.frame(ekegg_down)) > 0) {
  ekegg_down_readable <- setReadable(
    ekegg_down,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )
  
  write.csv(
    as.data.frame(ekegg_down_readable),
    "results/enrichment/GSE193336_KEGG_down.csv",
    row.names = FALSE
  )
  
  p_kegg_down <- dotplot(ekegg_down_readable, showCategory = 20) +
    ggtitle("KEGG enrichment of LPS-downregulated genes")
  
  ggsave(
    "figures/enrichment/Figure5_KEGG_down_dotplot.pdf",
    p_kegg_down,
    width = 9,
    height = 7
  )
  
  ggsave(
    "figures/enrichment/Figure5_KEGG_down_dotplot.png",
    p_kegg_down,
    width = 9,
    height = 7,
    dpi = 300
  )
}

# -----------------------------
# 5. GSEA using GO Biological Process
# -----------------------------

# Create ranked gene list by logFC
gene_list <- deg_entrez_unique$logFC
names(gene_list) <- deg_entrez_unique$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_go_bp <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_go_bp_readable <- setReadable(
  gsea_go_bp,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

write.csv(
  as.data.frame(gsea_go_bp_readable),
  "results/enrichment/GSE193336_GSEA_GO_BP.csv",
  row.names = FALSE
)

# GSEA dotplot
p_gsea_dot <- dotplot(gsea_go_bp_readable, showCategory = 20, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA of GO Biological Process terms")

ggsave(
  "figures/enrichment/Figure5_GSEA_GO_BP_dotplot.pdf",
  p_gsea_dot,
  width = 12,
  height = 7
)

ggsave(
  "figures/enrichment/Figure5_GSEA_GO_BP_dotplot.png",
  p_gsea_dot,
  width = 12,
  height = 7,
  dpi = 300
)

# -----------------------------
# 6. GSEA using KEGG
# -----------------------------

gsea_kegg <- tryCatch(
  {
    gseKEGG(
      geneList = gene_list,
      organism = "hsa",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
  },
  error = function(e) {
    message("GSEA KEGG failed: ", e$message)
    NULL
  }
)

if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
  gsea_kegg_readable <- setReadable(
    gsea_kegg,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )
  
  write.csv(
    as.data.frame(gsea_kegg_readable),
    "results/enrichment/GSE193336_GSEA_KEGG.csv",
    row.names = FALSE
  )
  
  p_gsea_kegg_dot <- dotplot(gsea_kegg_readable, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle("GSEA of KEGG pathways")
  
  ggsave(
    "figures/enrichment/Figure5_GSEA_KEGG_dotplot.pdf",
    p_gsea_kegg_dot,
    width = 12,
    height = 7
  )
  
  ggsave(
    "figures/enrichment/Figure5_GSEA_KEGG_dotplot.png",
    p_gsea_kegg_dot,
    width = 12,
    height = 7,
    dpi = 300
  )
}

# -----------------------------
# 7. Save top enrichment summary
# -----------------------------

top_go_up <- as.data.frame(ego_up_bp) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  mutate(category = "GO_BP_up")

top_go_down <- as.data.frame(ego_down_bp) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  mutate(category = "GO_BP_down")

top_gsea_go <- as.data.frame(gsea_go_bp_readable) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  mutate(category = "GSEA_GO_BP")

top_summary <- bind_rows(
  top_go_up,
  top_go_down,
  top_gsea_go
)

write.csv(
  top_summary,
  "results/enrichment/GSE193336_top_enrichment_summary.csv",
  row.names = FALSE
)

# -----------------------------
# 8. Save session information
# -----------------------------

sink("results/enrichment/GSE193336_enrichment_sessionInfo.txt")
sessionInfo()
sink()

message("GSE193336 enrichment analysis completed.")