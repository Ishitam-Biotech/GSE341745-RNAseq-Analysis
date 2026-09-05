# ============================================================
# RNA-seq Analysis of GSE341745
# Mouse neutrophil transcriptomic analysis
# ============================================================

# 1. Load required packages
library(readxl)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)

# 2. Import count data
counts <- read_excel(
  "data/GSE341745_merged_gene_counts.xlsx"
)

# Check dimensions and column names
dim(counts)
colnames(counts)
# ============================================================
# 3. Prepare expression matrix
# ============================================================

expr <- as.matrix(counts[, -1])

rownames(expr) <- counts$`Gene ID`

# Check dimensions
dim(expr)


# ============================================================
# 4. Create sample metadata
# ============================================================

sample_names <- colnames(counts)[-1]

condition <- factor(
  ifelse(grepl("12h", sample_names), "12h", "24h"),
  levels = c("12h", "24h")
)

sample_info <- data.frame(
  sample = sample_names,
  condition = condition
)

# View sample information
sample_info
# ============================================================
# 5. Filter low-expression genes
# ============================================================

cpm_expr <- cpm(expr)

keep <- rowSums(cpm_expr >= 1) >= 4

expr_filtered <- expr[keep, ]

# Check how many genes remain
dim(expr_filtered)


# ============================================================
# 6. TMM normalization
# ============================================================

dge <- DGEList(
  counts = expr_filtered,
  group = sample_info$condition
)

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

# View normalization information
dge$samples
# ============================================================
# 7. PCA - Principal Component Analysis
# ============================================================

logCPM <- cpm(
  dge,
  log = TRUE,
  prior.count = 2
)

pca <- prcomp(
  t(logCPM),
  scale. = FALSE
)

# PCA plot
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = sample_info$condition,
  sample = sample_info$sample
)

ggplot(
  pca_df,
  aes(x = PC1, y = PC2, label = sample)
) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  labs(
    title = "PCA of RNA-seq Samples",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()
# ============================================================
# 8. Differential Expression Analysis
#    Comparison: 24 h vs 12 h
# ============================================================

design <- model.matrix(
  ~ condition,
  data = sample_info
)

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

# View top DEGs
head(deg_results)
# ============================================================
# 9. Identify Significant DEGs
#    Criteria: FDR < 0.05 and |log2FC| >= 1
# ============================================================

deg_sig <- deg_results[
  deg_results$FDR < 0.05 &
    abs(deg_results$logFC) >= 1,
]

# Count significant DEGs
nrow(deg_sig)

# Count upregulated and downregulated genes
sum(deg_sig$logFC >= 1)
sum(deg_sig$logFC <= -1)


# ============================================================
# 10. Volcano Plot
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

ggplot(
  deg_results,
  aes(x = logFC, y = -log10(FDR), color = significance)
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
# ============================================================
# 11. Heatmap of Top 50 Differentially Expressed Genes
# ============================================================

top50_genes <- rownames(
  deg_results[
    order(deg_results$FDR),
  ][1:50, ]
)

heatmap_data <- logCPM[top50_genes, ]

pheatmap(
  heatmap_data,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)
# ============================================================
# 12. Annotate Ensembl Gene IDs with Gene Symbols
# ============================================================

id_map <- mapIds(
  org.Mm.eg.db,
  keys = rownames(deg_sig),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

deg_sig$GeneSymbol <- id_map[rownames(deg_sig)]

# Keep genes with recognized symbols
deg_annotated <- deg_sig[
  !is.na(deg_sig$GeneSymbol) &
    deg_sig$GeneSymbol != "",
]

# Separate upregulated and downregulated genes
up_genes <- deg_annotated[
  deg_annotated$logFC >= 1,
  "GeneSymbol"
]

down_genes <- deg_annotated[
  deg_annotated$logFC <= -1,
  "GeneSymbol"
]

# Check numbers
length(up_genes)
length(down_genes)

head(deg_annotated)
# ============================================================
# 13. GO Biological Process Enrichment
# ============================================================

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
# ============================================================
# 15. KEGG Pathway Enrichment - Upregulated Genes
# ============================================================

kegg_up <- bitr(
  up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ekegg_up <- enrichKEGG(
  gene = kegg_up$ENTREZID,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# View top KEGG pathways
head(as.data.frame(ekegg_up), 10)


# ============================================================
# 16. KEGG Enrichment Plot
# ============================================================

dotplot(
  ekegg_up,
  showCategory = 15,
  title = "KEGG Pathway Enrichment - Upregulated Genes"
)
# ============================================================
# 17. Prepare Ranked Gene List for GSEA
# ============================================================

gene_rank <- data.frame(
  GeneID = rownames(deg_results),
  logFC = deg_results$logFC
)

gene_rank <- gene_rank[
  !is.na(gene_rank$logFC),
]

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

# Remove duplicate Entrez IDs
gene_rank_entrez <- gene_rank_entrez[
  !duplicated(gene_rank_entrez$ENTREZID),
]

# Create ranked gene list
gene_list <- gene_rank_entrez$logFC
names(gene_list) <- gene_rank_entrez$ENTREZID

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)


# ============================================================
# 18. GSEA - GO Biological Process
# ============================================================

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


# ============================================================
# 19. GSEA Plot
# ============================================================

dotplot(
  gsea_go,
  showCategory = 15,
  title = "GSEA - GO Biological Processes"
)
# ============================================================
# 20. Prepare Top 50 Upregulated Genes for PPI Analysis
# ============================================================

top50_up <- deg_annotated[
  deg_annotated$logFC >= 1,
]

top50_up <- top50_up[
  order(top50_up$FDR),
][1:50, ]

ppi_genes <- top50_up$GeneSymbol

ppi_genes


# ============================================================
# 21. Prepare Top 50 Downregulated Genes for PPI Analysis
# ============================================================

top50_down <- deg_annotated[
  deg_annotated$logFC <= -1,
]

top50_down <- top50_down[
  order(top50_down$FDR),
][1:50, ]

ppi_genes_down <- top50_down$GeneSymbol

ppi_genes_down
# ============================================================
# 22. Create Results Directory
# ============================================================

dir.create(
  "RNAseq_results",
  showWarnings = FALSE
)


# ============================================================
# 23. Save Differential Expression Results
# ============================================================

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
# 24. Save GO Enrichment Results
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


# ============================================================
# 25. Save KEGG Enrichment Results
# ============================================================

write.csv(
  as.data.frame(ekegg_up),
  "RNAseq_results/KEGG_Upregulated.csv",
  row.names = FALSE
)


# ============================================================
# 26. Save GSEA Results
# ============================================================

write.csv(
  as.data.frame(gsea_go),
  "RNAseq_results/GSEA_GO_Biological_Process.csv",
  row.names = FALSE
)


# ============================================================
# 27. Save PPI Gene Lists
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
# ============================================================
# 28. Save Important Figures
# ============================================================

# PCA plot
pca_plot <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, label = sample)
) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  labs(
    title = "PCA of RNA-seq Samples",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()

ggsave(
  "RNAseq_results/PCA_plot.pdf",
  pca_plot,
  width = 8,
  height = 6
)


# Volcano plot
volcano_plot <- ggplot(
  deg_results,
  aes(x = logFC, y = -log10(FDR), color = significance)
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

ggsave(
  "RNAseq_results/Volcano_plot.pdf",
  volcano_plot,
  width = 8,
  height = 6
)


# GO plots
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


# KEGG plot
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


# GSEA plot
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
# ============================================================
# 29. Final Analysis Summary
# ============================================================

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
