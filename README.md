# GeneStream
Welcome to **GeneStream**! Here you can:
- Examine sample metadata to understand dataset variables.
- Filter and visualize normalized gene counts.
- Explore differential expression results, including volcano plots.
- Perform and interpret gene set enrichment analysis (GSEA) to identify enriched pathways.

## app.R
This Shiny application, **GeneStream**, provides an integrated workflow for analyzing and visualizing high-throughput gene expression data. It includes multiple tabs for:
1. Samples: Examine sample metadata, view summary statistics, tables, and histograms to understand the distribution of sample-related variables.
2. Counts: Load and filter normalized count matrices by variance and non-zero criteria, then visualize the effects through diagnostic plots, clustered heatmaps, and PCA.
3. Differential Expression Analysis: Upload DE results (e.g., from DESeq2), visualize them in sortable tables, and generate volcano plots to identify genes with significant expression changes.
4. Gene Set Enrichment Analysis (GSEA): Upload GSEA results, filter pathways by adjusted p-values and direction of effect, display top pathways in a barplot, and explore data through tables and scatter plots of NES vs. -log10(padj).

Through an interactive interface, users can explore and summarize their data, apply various filters, and visualize key results in an intuitive manner. The application supports uploading multiple CSV files, adjusting parameters using sliders and dropdowns, and downloading filtered results.

<img width="1470" alt="Screenshot 2024-12-16 at 13 50 07" src="https://github.com/user-attachments/assets/5f0cc6c0-5541-44ee-8010-54cd47e4fb25" />

## process_data.R
This Shiny application provides a comprehensive workflow for processing, analyzing, and visualizing gene expression data. It features the following functionalities:
1. Series Matrix to CSV: Upload GEO series matrix files, extract sample metadata, and convert them into a structured CSV format. Preview extracted sample characteristics before downloading.
2. TSV to CSV Conversion: Transform tab-delimited TSV files into CSV format, ensuring compatibility with other downstream analyses. Preview the converted table to validate the data.
3. Differential Expression Analysis (DESeq2): Perform DESeq2-based differential expression analysis. Users can upload raw count matrices and sample metadata, define group comparisons, and generate a results table with log2 fold changes, p-values, and adjusted p-values (padj). Results can be previewed and downloaded.
4. Gene Set Enrichment Analysis (FGSEA): Perform FGSEA using differential expression results and gene set files in GMT format. Users can rank genes based on DESeq2 statistics (e.g., log2 fold change) and identify enriched pathways. Results, including adjusted p-values, normalized enrichment scores (NES), and leading-edge genes, are available for preview and download.

Through an intuitive interface, ExpressionExplorer allows users to preprocess data, apply filters, and generate key insights with diagnostic plots and enrichment analyses. This app supports seamless uploads, customizable parameters, and downloadable outputs to streamline transcriptomic workflows.

<img width="1139" alt="Screenshot 2024-12-16 at 13 56 22" src="https://github.com/user-attachments/assets/c4e1ad91-f0d2-4122-ba84-ec1f7fb57516" />





