# GSE341745-RNAseq-Analysis

**About this project**

I worked on this project to analyze the RNA-seq dataset GSE341745 and understand how gene expression changes between the 12 h and 24 h time points.

I used R for the analysis, starting with basic data checks and filtering, followed by normalization, PCA, differential expression analysis, and pathway-level analysis. I also explored the significant genes using GO, KEGG, GSEA, and STRING PPI analysis.

**Dataset**

GEO accession: GSE341745
Organism: Mus musculus
Samples: 8
4 samples at 12 h
4 samples at 24 h
Genes in the original expression matrix: 42,003

The expression matrix was obtained from the GEO dataset and used as the starting point for my analysis.

**My analysis**

**Data preparation and filtering
**
I first imported the expression matrix into R and checked the structure, dimensions, missing values, and expression distributions.

I filtered out genes with very low expression using CPM values. I kept genes with CPM ≥ 1 in at least 4 samples.

This reduced the dataset from 42,003 genes to 11,207 genes.

I then used TMM normalization from edgeR to account for differences in library size between samples.

**PCA**

I used log-CPM transformed expression values for PCA to look at the overall variation between samples.

This helped me check how the samples grouped and whether there was separation between the 12 h and 24 h conditions.

**Differential expression analysis**

For differential expression, I used the quasi-likelihood framework in edgeR to compare 24 h vs 12 h.

I considered a gene significant when:

FDR < 0.05
|log2FC| ≥ 1

I obtained:

1,899 significant DEGs
868 upregulated genes
1,031 downregulated genes

I used a volcano plot to visualize the overall differential expression and a heatmap to look at the expression patterns of the top 50 DEGs.

**Gene annotation**

The original data contained Ensembl gene IDs, so I converted these to mouse gene symbols using org.Mm.eg.db.

This made the significant genes easier to interpret and use in the downstream enrichment analyses.

**GO enrichment**

I performed GO Biological Process enrichment separately for the upregulated and downregulated genes using clusterProfiler.

The upregulated genes were strongly associated with immune and antiviral-related processes. Some of the major terms included interferon responses, response to virus, innate immune response, and antigen processing and presentation.

For the downregulated genes, some of the enriched processes were related to reactive oxygen species metabolism, oxidative stress, endocytosis, receptor-mediated endocytosis, cilium organization, and lymphocyte differentiation.

**KEGG pathway analysis
**
I also performed KEGG pathway enrichment for the upregulated genes.

Some of the pathways that came up strongly were:

Herpes simplex virus 1 infection
Influenza A
Hepatitis C
Apoptosis
Epstein-Barr virus infection
Proteasome
Lysosome-related pathways
Antigen processing and presentation

The viral infection pathways are mainly pathway annotations containing shared host immune and antiviral genes. I therefore would not interpret their enrichment as evidence that the samples were infected with those specific viruses.

**GSEA**

After the DEG analysis, I also wanted to look at pathway-level changes without relying only on the genes that passed my DEG cutoff.

For this, I ranked the genes according to their log2 fold change and performed GSEA using GO Biological Process gene sets.

**PPI analysis**

Finally, I explored the protein–protein interaction relationships among the top differentially expressed genes using STRING.

The main connected part of the network included genes such as:

Ccl2, Kdr, Axl, Vwf, P2rx7, Pycard, Bak1, and Prok2.

From the network structure, Ccl2 and Kdr appeared to have relatively high connectivity in the displayed network.

I treated these as candidate hub genes based on the network structure rather than considering them definitive hub proteins.

**Main observations**

One of the clearest patterns I observed was the enrichment of interferon, antiviral, and innate immune-related processes among the upregulated genes at 24 h.

There were also changes involving oxidative stress, endocytosis, apoptosis, lysosomal processes, and cellular organization.

The PPI analysis gave me another way to look at the genes identified by the differential expression analysis and highlighted a connected group involving Ccl2, Kdr, Axl, and Vwf.

**Results**

I have included the main output files in the RNAseq_results folder:

Differential expression results
Annotated significant DEGs
GO enrichment results
KEGG enrichment results
GSEA results
PPI gene lists
PCA plot
Volcano plot
GO enrichment plots
KEGG enrichment plot
GSEA plot

**Tools I used**

R | RStudio | edgeR | clusterProfiler | org.Mm.eg.db | AnnotationDbi | readxl | ggplot2 | pheatmap

I used STRING for the PPI network analysis.

## Project structure

* GSE341745_RNAseq_Analysis.R — Main R script containing the complete analysis
* RNAseq_results/ — Contains the results, plots, and output files generated during the analysis
* README.md — Project overview, methods, and main findings

 Dataset source

**GSE341745 — NCBI Gene Expression Omnibus (GEO)
**
I used the GEO accession GSE341745 to identify the dataset used in this project.

