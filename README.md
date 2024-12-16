# GeneStream
This Shiny application, **GeneStream**, provides an integrated workflow for analyzing and visualizing high-throughput gene expression data. It includes multiple tabs for:
1. Samples: Examine sample metadata, view summary statistics, tables, and histograms to understand the distribution of sample-related variables.
2. Counts: Load and filter normalized count matrices by variance and non-zero criteria, then visualize the effects through diagnostic plots, clustered heatmaps, and PCA.
3. Differential Expression Analysis: Upload DE results (e.g., from DESeq2), visualize them in sortable tables, and generate volcano plots to identify genes with significant expression changes.
4. Gene Set Enrichment Analysis (GSEA): Upload GSEA results, filter pathways by adjusted p-values and direction of effect, display top pathways in a barplot, and explore data through tables and scatter plots of NES vs. -log10(padj).

Through an interactive interface, users can explore and summarize their data, apply various filters, and visualize key results in an intuitive manner. The application supports uploading multiple CSV files, adjusting parameters using sliders and dropdowns, and downloading filtered results.
