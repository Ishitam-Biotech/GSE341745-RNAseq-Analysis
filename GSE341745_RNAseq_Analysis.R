# RNA-seq Analysis of GSE341745
# Mouse neutrophil transcriptomic analysis

# 1. Load required packages
library(readxl)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)

# 2. Check project structure and import count data

input_file <- "data/GSE341745_merged_gene_counts.xlsx"
if (!file.exists(input_file)) {
 stop(
 "Input file not found. Please run this script from the GitHub project root."
 )
}
counts <- read_excel(input_file)
# Basic input validation
if (ncol(counts) < 3) {
 stop("Input file must contain a Gene ID column and at least two sample columns.")
}
if (!"Gene ID" %in% colnames(counts)) {
 stop("The first column must be named 'Gene ID'.")
}

# 3. Prepare expression matrix

expr <- as.matrix(counts[, -1])
rownames(expr) <- counts$`Gene ID`
# Ensure expression values are numeric
mode(expr) <- "numeric"
if (anyNA(expr)) {
 stop("Missing or non-numeric values detected in the expression matrix.")
}
# Check for duplicated gene IDs
if (anyDuplicated(rownames(expr))) {
 stop("Duplicated Gene IDs detected. Resolve duplicates before analysis.")
}


# 4. Create sample metadata

sample_names <- colnames(counts)[-1]

condition <-  ifelse(
  grepl("12h", sample_names, ignore.case = TRUE), 
  "12h",
  ifelse(
     grepl("24h", sample_names, ignore.case = TRUE),
     "24h",
      NA
    )
)
if (any(is.na(condition))) {
  stop(
    "Could not determine condition for sample(s): ",
    paste(sample_names[is.na(condition)], collapse = ", ")
  )
}

condition <- factor(
  condition,
  levels = c("12h", "24h")
)
sample_info <- data.frame(
  sample = sample_names,
  condition = condition
)
if (any(table(sample_info$condition) < 2)) {
 stop("Each condition must contain at least two samples.")
}

# 5. Create DGEList and design matrix

dge <- DGEList(
 counts = expr,
 group = sample_info$condition
)
design <- model.matrix(
 ~ condition,
 data = sample_info
)

# 6. Filter low-expression genes

keep <- filterByExpr(
 dge,
 design = design
)
dge <- dge[
 keep,
 ,
 keep.lib.sizes = FALSE
]


# 7. TMM normalization
dge <- calcNormFactors(
 dge,
 method = "TMM"
)


# 8. PCA - Principal Component Analysis

logCPM <- cpm(
  dge,
  log = TRUE,
  prior.count = 2
)

pca <- prcomp(
  t(logCPM),
  scale. = FALSE
)

percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
pca_df <- data.frame(
 PC1 = pca$x[, 1],
 PC2 = pca$x[, 2],
 condition = sample_info$condition,
 sample = sample_info$sample
)
pca_plot <- ggplot(
 pca_df,
 aes(x = PC1, y = PC2, label = sample)
) +
 geom_point(size = 4) +
 geom_text(vjust = -1) +
 labs(
 title = "PCA of RNA-seq Samples",
 x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
 y = paste0("PC2 (", round(percent_var[2], 1), "%)")
 ) +
 theme_minimal()


# 9. Differential Expression Analysis
#    Comparison: 24 h vs 12 h

# Estimate dispersion
dge <- estimateDisp(
  dge,
  design
)

# Fit quasi-likelihood model
fit <- glmQLFit(
  dge,
  design
)

# Test 24 h vs 12 h
qlf <- glmQLFTest(
  fit,
  coef = "condition24h"
)

# Extract all DEGs
deg_results <- topTags(
  qlf,
  n = Inf
)$table


# 10. Identify Significant DEGs
#    Criteria: FDR < 0.05 and |log2FC| >= 1

deg_sig <- deg_results[
  deg_results$FDR < 0.05 &
    abs(deg_results$logFC) >= 1,
]


# ============================================================
# 11. Volcano Plot
# ============================================================

deg_results$significance <- "Not significant"

deg_results$significance[
  deg_results$FDR < 0.05 &
    deg_results$logFC >= 1
] <- "Upregulated"

deg_results$significance[
  deg_results$FDR < 0.05 &
    deg_results$logFC <= -1
] <- "Downregulated"

# Protect against -log10(0)
deg_results$plot_FDR <- pmax(
 deg_results$FDR,
 .Machine$double.xmin
)
volcano_plot <- ggplot(
 deg_results,
 aes(
 x = logFC,
 y = -log10(plot_FDR),
 color = significance
 )
) +
 geom_point(alpha = 0.6, size = 1.5) +
 geom_vline(
 xintercept = c(-1, 1),
 linetype = "dashed"
 ) +
 geom_hline(
 yintercept = -log10(0.05),
 linetype = "dashed"
 ) +
 labs(
 title = "Differential Gene Expression: 24 h vs 12 h",
 x = "Log2 Fold Change",
 y = "-Log10 FDR"
 ) +
 theme_minimal()

# 12. Heatmap of Top 50 Differentially Expressed Genes

sig_for_heatmap <- deg_results[
 deg_results$FDR < 0.05 &
 abs(deg_results$logFC) >= 1,
]
n_top <- min(50, nrow(sig_for_heatmap))
if (n_top > 0) {
 top_genes <- rownames(
 sig_for_heatmap[
 order(sig_for_heatmap$FDR),
 ][seq_len(n_top), ]
 )
 heatmap_data <- logCPM[
 top_genes,
 ,
 drop = FALSE
 ]
 pheatmap(
 heatmap_data,
 scale = "row",
 clustering_distance_rows = "euclidean",
 clustering_distance_cols = "euclidean",
 clustering_method = "complete",
 main = paste0("Top ", n_top, " Significant DEGs")
 )
}

# 13. Annotate Ensembl Gene IDs with Gene Symbols

ensembl_ids <- sub(
 "\\..*$",
 "",
 rownames(deg_sig)
)
id_map <- mapIds(
 org.Mm.eg.db,
 keys = ensembl_ids,
 keytype = "ENSEMBL",
 column = "SYMBOL",
 multiVals = "first"
)
names(id_map) <- rownames(deg_sig)
deg_sig$GeneSymbol <- id_map[
 rownames(deg_sig)
]
deg_annotated <- deg_sig[
 !is.na(deg_sig$GeneSymbol) &
 deg_sig$GeneSymbol != "",
]
# Remove duplicate symbols for downstream gene-list analyses
deg_annotated_unique <- deg_annotated[
 !duplicated(deg_annotated$GeneSymbol),
]
# Separate upregulated and downregulated genes
up_genes <- unique(
 deg_annotated_unique[
 deg_annotated_unique$logFC >= 1,
 "GeneSymbol"
 ]
)
down_genes <- unique(
 deg_annotated_unique[
 deg_annotated_unique$logFC <= -1,
 "GeneSymbol"
 ]
)

# 14. Prepare enrichment background/universe
# Background = genes retained after expression filtering

filtered_ensembl <- sub(
 "\\..*$",
 "",
 rownames(dge)
)
background_symbols <- mapIds(
 org.Mm.eg.db,
 keys = filtered_ensembl,
 keytype = "ENSEMBL",
 column = "SYMBOL",
 multiVals = "first"
)
background_symbols <- unique(
 na.omit(background_symbols)
)

# 15. GO Biological Process Enrichment

ego_up <- enrichGO(
  gene = up_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_down <- enrichGO(
  gene = down_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# View top enriched terms
head(as.data.frame(ego_up), 10)
head(as.data.frame(ego_down), 10)


# ============================================================
# 14. GO Enrichment Plots
# ============================================================

dotplot(
  ego_up,
  showCategory = 15,
  title = "GO Biological Process Enrichment - Upregulated Genes"
)

dotplot(
  ego_down,
  showCategory = 15,
  title = "GO Biological Process Enrichment - Downregulated Genes"
)

# 16. KEGG Pathway Enrichment - Upregulated Genes
# Analysis retained for upregulated genes as in original workflow

kegg_up <- bitr(
  up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

background_entrez <- bitr(
 filtered_ensembl,
 fromType = "ENSEMBL",
 toType = "ENTREZID",
 OrgDb = org.Mm.eg.db
)$ENTREZID

background_entrez <- unique(
 na.omit(background_entrez)
)

ekegg_up <- enrichKEGG(
  gene = kegg_up$ENTREZID,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 17. Prepare Ranked Gene List for GSEA
#     Ranking metric: log2 fold change

gene_rank <- data.frame(
  GeneID = rownames(deg_results),
  logFC = deg_results$logFC
)

gene_rank <- gene_rank[
  !is.na(gene_rank$logFC),
]

gene_rank$GeneID <- sub(
 "\\..*$",
 "",
 gene_rank$GeneID
)

gene_rank_entrez <- bitr(
  gene_rank$GeneID,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

gene_rank_entrez <- merge(
  gene_rank_entrez,
  gene_rank,
  by.x = "ENSEMBL",
  by.y = "GeneID"
)

# Retain the strongest absolute logFC when multiple
# Ensembl IDs map to the same Entrez ID
gene_rank_entrez <- gene_rank_entrez[
 order(
 abs(gene_rank_entrez$logFC),
 decreasing = TRUE
 ),
]
gene_rank_entrez <- gene_rank_entrez[
 !duplicated(gene_rank_entrez$ENTREZID),
]
gene_list <- gene_rank_entrez$logFC
names(gene_list) <- gene_rank_entrez$ENTREZID
gene_list <- sort(
 gene_list,
 decreasing = TRUE
)

# 18. GSEA - GO Biological Process

gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View top enriched gene sets
head(as.data.frame(gsea_go), 10)


# 19. Prepare Top 50 Upregulated Genes for PPI Analysis
#     This script prepares gene lists; STRING network analysis
#     is performed externally unless separately documented.

top50_up <- deg_annotated_unique[
 deg_annotated_unique$logFC >= 1,
]
top50_up <- top50_up[
 order(top50_up$FDR),
]
top50_up <- head(top50_up, 50)
ppi_genes <- unique(
 top50_up$GeneSymbol
)
top50_down <- deg_annotated_unique[
 deg_annotated_unique$logFC <= -1,
]
top50_down <- top50_down[
 order(top50_down$FDR),
]
top50_down <- head(top50_down, 50)
ppi_genes_down <- unique(
 top50_down$GeneSymbol
)


# 20. Create Results Directory


dir.create(
  "RNAseq_results",
  showWarnings = FALSE
)


# 21. Save Differential Expression Results

write.csv(
  deg_results,
  "RNAseq_results/All_DEG_results.csv",
  row.names = TRUE
)

write.csv(
  deg_annotated,
  "RNAseq_results/Significant_DEGs_annotated.csv",
  row.names = TRUE
)


# ============================================================
# 22. Save GO Enrichment Results
# ============================================================

write.csv(
  as.data.frame(ego_up),
  "RNAseq_results/GO_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_down),
  "RNAseq_results/GO_Downregulated.csv",
  row.names = FALSE
)


# 23. Save KEGG Enrichment Results

write.csv(
  as.data.frame(ekegg_up),
  "RNAseq_results/KEGG_Upregulated.csv",
  row.names = FALSE
)


# 24. Save GSEA Results

write.csv(
  as.data.frame(gsea_go),
  "RNAseq_results/GSEA_GO_Biological_Process.csv",
  row.names = FALSE
)


# ============================================================
# 25. Save PPI Gene Lists
# ============================================================

write.table(
  ppi_genes,
  "RNAseq_results/PPI_Top50_Upregulated.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  ppi_genes_down,
  "RNAseq_results/PPI_Top50_Downregulated.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# 28. Save Important Figures

#PCA plot
ggsave(
  "RNAseq_results/PCA_plot.pdf",
  pca_plot,
  width = 8,
  height = 6
)

# Volcano plot

ggsave(
  "RNAseq_results/Volcano_plot.pdf",
  volcano_plot,
  width = 8,
  height = 6
)

# GO plots
if (nrow(as.data.frame(ego_up)) > 0) {
pdf(
  "RNAseq_results/GO_Upregulated_dotplot.pdf",
  width = 9,
  height = 7
)

print(
  dotplot(
    ego_up,
    showCategory = 15,
    title = "GO Biological Process Enrichment - Upregulated Genes"
  )
)
dev.off()
}

if (nrow(as.data.frame(ego_down)) > 0) {
pdf(
  "RNAseq_results/GO_Downregulated_dotplot.pdf",
  width = 9,
  height = 7
)

print(
  dotplot(
    ego_down,
    showCategory = 15,
    title = "GO Biological Process Enrichment - Downregulated Genes"
  )
)

dev.off()
}

# KEGG plot
if (nrow(as.data.frame(ekegg_up)) > 0) {
pdf(
  "RNAseq_results/KEGG_Upregulated_dotplot.pdf",
  width = 9,
  height = 7
)

print(
  dotplot(
    ekegg_up,
    showCategory = 15,
    title = "KEGG Pathway Enrichment - Upregulated Genes"
  )
)

dev.off()
}

# GSEA plot
if (nrow(as.data.frame(gsea_go)) > 0) {
pdf(
  "RNAseq_results/GSEA_GO_dotplot.pdf",
  width = 9,
  height = 7
)

print(
  dotplot(
    gsea_go,
    showCategory = 15,
    title = "GSEA - GO Biological Processes"
  )
)

dev.off()
  }

# 27. Final Analysis Summary

cat("RNA-seq Analysis Summary\n")
cat("========================\n")

cat("Total genes:", nrow(expr), "\n")
cat("Genes after filtering:", nrow(expr_filtered), "\n")
cat("Significant DEGs:", nrow(deg_sig), "\n")
cat("Upregulated genes:", sum(deg_sig$logFC >= 1), "\n")
cat("Downregulated genes:", sum(deg_sig$logFC <= -1), "\n")
cat("Annotated significant DEGs:", nrow(deg_annotated), "\n")

cat("\nTop upregulated genes:\n")
print(head(
  top50_up$GeneSymbol,
  10
))

cat("\nTop downregulated genes:\n")
print(head(
  top50_down$GeneSymbol,
  10
))

cat("\nAnalysis completed successfully.\n")

# Save reproducibility information
writeLines(
 capture.output(sessionInfo()),
 "RNAseq_results/sessionInfo.txt"
)
