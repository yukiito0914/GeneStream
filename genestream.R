## Author: Yuki Ito
## yuki@bu.edu

## This Shiny application provides an integrated workflow for analyzing and visualizing 
## high-throughput gene expression data. It includes multiple tabs for:
## 1. Samples: Examine sample metadata, view summary statistics, tables, and histograms 
##    to understand the distribution of sample-related variables.
## 2. Counts: Load and filter normalized count matrices by variance and non-zero criteria, 
##    then visualize the effects through diagnostic plots, clustered heatmaps, and PCA.
## 3. Differential Expression Analysis: Upload DE results (e.g., from DESeq2), 
##    visualize them in sortable tables, and generate volcano plots to identify genes 
##    with significant expression changes.
## 4. Gene Set Enrichment Analysis (GSEA): Upload GSEA results, filter pathways by 
##    adjusted p-values and direction of effect, display top pathways in a barplot, 
##    and explore data through tables and scatter plots of NES vs. -log10(padj).

## Through an interactive interface, users can explore and summarize their data, 
## apply various filters, and visualize key results in an intuitive manner.
## The application supports uploading multiple CSV files, adjusting 
## parameters using sliders and dropdowns, and downloading filtered results.

options(shiny.maxRequestSize = 1024 * 1024 * 1024)

library(shiny)
library(dplyr)
library(bslib)
library(ggplot2)
library(DT)
library(ComplexHeatmap)  
library(circlize)
library(colourpicker) 
library(fgsea)

ui <- fluidPage(
  titlePanel("GeneStream"),
  "Welcome to GeneStream! This is Yuki Ito. Here you can:",
  tags$ul(
    tags$li("Examine sample metadata to understand dataset variables."),
    tags$li("Filter and visualize normalized gene counts."),
    tags$li("Explore differential expression results, including volcano plots."),
    tags$li("Perform and interpret gene set enrichment analysis (GSEA) to identify enriched pathways.")
  ),
  tabsetPanel(
    tabPanel("Samples", 
             fluidRow(
               column(3,
                      "The distinct values and distributions of sample information are important to understand before conducting analysis of corresponding sample data. You can load and examine a sample information matrix.",
                      br(),
                      fileInput("sample_info", "Input: Sample file (.csv)", accept = c(".csv")),
                      actionButton("submit_samples", "Submit", style = "width: 100%;")),
               column(9,
                      tabsetPanel(
                        tabPanel("Summary", tableOutput("summary_table")),
                        tabPanel("Table", DTOutput("sample_data_table")),
                        tabPanel("Plots", 
                                 fluidRow(
                                   column(4,
                                          "You can generate histograms of continuous variables.",
                                          selectInput("plot_column", "Choose a column to plot:", choices = NULL, selected = NULL),
                                          selectInput("group_column", "Choose a column to group by:", choices = NULL, selected = NULL)),
                                   column(8,
                                          plotOutput("sample_data_plot"))
                                 ),
                                 )
                      ))
             )),
    tabPanel("Counts",
             fluidRow(
               column(3,
                      "Exploring and visualizing counts matrices can aid in selecting gene count filtering strategies and understanding counts data structure. You can choose different gene filtering thresholds and assess their effects using diagnostic plots of the counts matrix.",
                      fileInput("normalized_count", "Input: Normalized counts matrix (.csv)", accept = c(".csv")),
                      actionButton("submit_counts", "Submit", style = "width: 100%;"),
                      # Conditional panel to show sliderInput after button is clicked
                      conditionalPanel(
                        condition = "input.submit_counts > 0",
                        sliderInput("var_percentile", "Genes with at least X percentile of variance", min=0, max=100, value=50),
                        sliderInput("non_zero_samples", "genes with at least X samples that are non-zero", min=0, max=100, value=50),
                        actionButton("filter_counts", "Filter", style = "width: 100%;")
                      )),
               column(9,
                      tabsetPanel(
                        tabPanel("Filtering", tableOutput("filtering_table")),
                        tabPanel("Plots",
                                 "You can generate diagnostic scatter plots, where genes passing filters are marked in a darker color, and genes filtered out are lighter. (Median count vs Variance & Median count vs Number of zeros)",
                                 actionButton("plot_med", "Generate plots"),
                                 plotOutput("med_var_plot"), 
                                 plotOutput("med_zeros_plot")
                        ),
                        tabPanel("Heatmap", 
                                 "You can generate a clustered heatmap of counts remaining after filtering",
                                 checkboxInput("log_transform", "Apply log transformation", value = TRUE),
                                 actionButton("plot_heatmap", "Plot"),
                                 plotOutput("heatmap")),
                        tabPanel("PCA", 
                                 actionButton("run_pca", "Run　PCA"),
                                 selectInput("x_pc", "Select X-axis PC:", choices = NULL),
                                 selectInput("y_pc", "Select Y-axis PC:", choices = NULL),
                                 actionButton("plot_pca", "Plot PCA"),
                                 plotOutput("pca"))
                      ))
             )),
    tabPanel("Different Expression Analysis",
             fluidRow(
               column(3,
                      "Differential expression identifies which genes, if any, are implicated in a specific biological comparison. You can load and explore a differential expression dataset.",
                      fileInput("dea_res", "Input: Results of a differential expression analysis (.csv)", accept = c(".csv")),
                      actionButton("submit_dea", "Submit", style = "width: 100%;")),
               column(9,
                      tabsetPanel(
                        tabPanel("Show Result", DTOutput("dea_table")),
                        tabPanel("Explore", 
                                 fluidRow(
                                   column(4,
                                          'A volcano plot can be generated with "log2 fold-change" on the x-axis and "p-adjusted" on the y-axis.',
                                          radioButtons("x_name", "Choose the column for the x-axis",
                                                       choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                                                       selected = "log2FoldChange"),
                                          radioButtons("y_name", "Choose the column for the y-axis",
                                                       choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                                                       selected = "padj"),
                                          colourInput("color1", "Base point color","blue"),
                                          colourInput("color2", "Highlight point color","red"),
                                          sliderInput("slider", "Select the magnitude of the p adjusted coloring:", min=-40, max=0, value=-20),
                                          actionButton("plotButton", "Plot", style = "width: 100%;")),
                                   column(8,
                                          tabsetPanel(
                                            tabPanel("Plot", plotOutput("volcano")),
                                            tabPanel("Table", tableOutput("de_table"))
                                          ))
                                 )
                        )
                      )
               )
             )),
    tabPanel("Gene Set Enrichment Analysis",
             fluidRow(
               column(3,
                      "Gene Set Enrichment Analysis (GSEA) provides a systematic approach to identify significantly enriched gene sets from high-throughput data, offering insights into biological pathways or processes associated with experimental conditions. You can upload differential expression results, rank genes by a chosen statistic, and test for enrichment using curated gene set databases. ",
                      fileInput("fgsea_res", "Input: A table of fgsea results from the differential expression data (.csv)", accept = c(".csv")),
                      actionButton("submit_fgsea", "Submit", style = "width: 100%;")),
               column(9,
                      tabsetPanel(
                        tabPanel("Barplot",
                                 fluidRow(column(4,
                                                 "You can generate barplots of fgsea NES for top pathways selected by slider.",
                                                 sliderInput("pathway_slider", "Number of top pathways to plot by adjusted p-value: ", min=10, max=50, value=30)),
                                          column(8,
                                                 plotOutput("fgsea_bar")))),
                        tabPanel("Table",
                                 fluidRow(column(4,
                                                 "You can generate a sortable data table displaying the results.",
                                                 sliderInput("filter_padj_tbl", "Filter table by adjusted p-value: ", min=0, max=1, value=0.5),
                                                 radioButtons("neg_or_pos", "What kind of pathway",choices = c("All", "Positive", "Negative"), selected = "All"),
                                                 downloadButton("download_fgsea_filtered", "Download", style = "width: 100%;")),
                                          column(8,
                                                 DTOutput("fgsea_table")))),
                        tabPanel("Scatter Plot",
                                 fluidRow(column(4,
                                                 "You can generate a scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis, with gene sets below threshold in grey color.",
                                                 sliderInput("filter_padj_plt", "Filter table by adjusted p-value: ", min=0, max=1, value=0.5)),
                                          column(8,
                                                 plotOutput("fgsea_plot"))))
                      ))
             ))
  )
)


server <- function(input, output, session) {
  
  # Samples tab
  
  load_sample <- reactive({
    dataf <- read.csv(req(input$sample_info$datapath))
    return(dataf)
  })
  
  # Function to create a table that includes a summary of the type and values in each column
  summary_table <- function(sample_info) {
    summary_list <- list()
    for (col in names(sample_info)) {
      col_data <- sample_info[[col]]
      col_type <- class(col_data)[1]  # Get the primary class (e.g., "character", "numeric")
      # Check if the column is numeric
      if (is.numeric(col_data)) {
        # Calculate mean and sd for numeric columns
        m <- mean(col_data, na.rm = TRUE)
        s <- sd(col_data, na.rm = TRUE)
        val_str <- sprintf("%.1f (+/- %.1f)", m, s)
      } else {
        # For non-numeric, list unique values
        unique_vals <- unique(col_data)
        val_str <- paste(unique_vals, collapse = ", ")
      }
      # Add summary to the list
      summary_list[[col]] <- list(
        col_name = col,
        type = col_type,
        mean_or_dist_val = val_str
      )
    }
    summary_df <- do.call(rbind, lapply(summary_list, as.data.frame))
    summary_df <- dplyr::rename(summary_df,
                                "Column Name" = col_name,
                                "Type" = type,
                                "Mean (sd) or Distinct Values" = mean_or_dist_val)
    return(summary_df)
  }
  
  # Function to create histograms of continuous variables
  draw_histogram <- function(sample_info, plot_col, group_col){
    sample_hist <- sample_info %>%
      ggplot(aes_string(x = plot_col, fill = group_col)) +
      geom_histogram(position = "dodge", alpha = 0.7) +
      theme_minimal() +
      labs(title = paste("Histogram of", plot_col, "grouped by", group_col),
           x = plot_col,
           fill = group_col)
    return(sample_hist)
  }
  
  output$summary_table <- renderTable({
    req(input$submit_samples)
    summary_table(load_sample())
  })
  
  output$sample_data_table <- renderDT({
    req(input$submit_samples)
    datatable(load_sample(), options = list(lengthChange = TRUE))
    load_sample()
  })
  
  # Update the choices for the grouping variable dropdown
  observe({
    req(input$submit_samples)
    continuous_cols <- sapply(load_sample(), function(col) is.numeric(col) && !is.factor(col))
    updateSelectInput(session, "plot_column", choices = names(load_sample())[continuous_cols], selected = NULL)
    grouping_cols <- sapply(load_sample(), function(col) is.character(col) || is.factor(col))
    updateSelectInput(session, "group_column", choices = names(load_sample())[grouping_cols], selected = NULL)
  })
  
  output$sample_data_plot <- renderPlot({
    req(input$plot_column, input$group_column) 
    draw_histogram(load_sample(), input$plot_column, input$group_column)
  })
  
  
  # Counts tab
  
  load_counts <- reactive({
    dataf <- read.csv(req(input$normalized_count$datapath), stringsAsFactors = FALSE)
    rownames(dataf) <- dataf$GeneID
    dataf$GeneID <- NULL
    return(dataf)
  })
  
  # Update the max value of slider input based on the number of samples
  observeEvent(input$submit_counts, {
    req(load_counts())
    col_count <- max(0, ncol(load_counts())) 
    updateSliderInput(session, "non_zero_samples", max = col_count, value = col_count)
  })
  
  # Function to filter count data
  filtering <- function(normalized_count, var_percentile, non_zero_samples){
    gene_variances <- apply(normalized_count, 1, var)
    threshold <- quantile(gene_variances, probs = var_percentile / 100)
    filtered <- normalized_count[gene_variances >= threshold, ]
    non_zero_counts <- rowSums(filtered != 0)
    filtered <- filtered[non_zero_counts >= non_zero_samples, ] 
    return(filtered)
  }
  
  filter_count <- reactive({
    req(input$filter_counts) 
    filtering(load_counts(), input$var_percentile, input$non_zero_samples)
  })
  
  # Function to create a table summarizing the effect of the filtering
  filtering_table <- function(normalized_count, filtered_count){
    filtering_tbl <- data.frame(
      num_samples = as.integer(ncol(normalized_count)),
      total_genes = nrow(normalized_count),
      filtered_genes = paste(nrow(filtered_count), "(", round((nrow(filtered_count)/nrow(normalized_count))*100, 1), "%)", sep=""),
      non_filtered_genes = paste(nrow(normalized_count)-nrow(filtered_count), "(", round((nrow(normalized_count)-nrow(filtered_count))/nrow(normalized_count)*100, 1), "%)", sep="")
    )
    filtering_tbl<- dplyr::rename(filtering_tbl,
                           "Number of samples" = num_samples,
                           "Total number of genes" = total_genes,
                           "Number of genes passing filter" = filtered_genes,
                           "Number of genes not passing filter" = non_filtered_genes)
    return(filtering_tbl)
  }
  
  output$filtering_table <- renderTable({
    req(input$filter_counts)
    filtering_table(load_counts(), filter_count())
  })
  
  # Helper function to add metadata
  add_metadata <- function(normalized_count, filtered_count) {
    normalized_count$Median <- apply(normalized_count, 1, median) + 1
    normalized_count$Variance <- apply(normalized_count, 1, var) + 1
    normalized_count$ZeroCounts <- apply(normalized_count, 1, function(x) sum(x == 0))
    normalized_count$Status <- ifelse(rownames(normalized_count) %in% rownames(filtered_count), "Passed", "Filtered")
    return(normalized_count)
  }
  
  # Function to create median count vs variance scatter plot
  med_var_plot <- function(normalized_count, filtered_count){
    normalized_count <- add_metadata(normalized_count, filtered_count)
    med_var <- ggplot(normalized_count, aes(x = Median, y = Variance, color = Status)) +
      geom_point(alpha = 0.7) +
      scale_x_log10() +  
      scale_y_log10() + 
      theme_minimal() +
      scale_color_manual(values = c("Passed" = "darkblue", "Filtered" = "lightblue")) +
      labs(title = "Median Count vs Variance",
           x = "Median Count (log scale)",
           y = "Variance (log scale)",
           color = "Filter Status")
    return(med_var)
  }
  
  # Function to create median count vs number of zeros scatter plot
  med_zeros_plot <- function(normalized_count, filtered_count){
    normalized_count <- add_metadata(normalized_count, filtered_count)
    med_zeros <- ggplot(normalized_count, aes(x = Median, y = ZeroCounts, color = Status)) +
      geom_point(alpha = 0.7) +
      scale_x_log10() +  
      theme_minimal() +
      scale_color_manual(values = c("Passed" = "darkblue", "Filtered" = "lightblue")) +
      labs(title = "Median Count vs Number of Zeros",
           x = "Median Count (log scale)",
           y = "Number of Zeros",
           color = "Filter Status")
    return(med_zeros)
  }
  
  output$med_var_plot <- renderPlot({
    req(input$plot_med) 
    med_var_plot(load_counts(), filter_count())
  })
  
  output$med_zeros_plot <- renderPlot({
    req(input$plot_med) 
    med_zeros_plot(load_counts(), filter_count())
  })
  
  # Function to create a heatmap
  plot_heatmap <- function(filtered_count, log_transform = TRUE){
    filtered_count <- as.matrix(filtered_count)
    if (log_transform) {
      data_to_plot <- log2(filtered_count + 1)  # log-transform with a pseudocount
    } else {
      data_to_plot <- filtered_count
    }
    Heatmap(data_to_plot, name = "Expression", cluster_rows = TRUE, cluster_columns = TRUE, col = colorRamp2(c(0, max(data_to_plot)/2, max(data_to_plot)), c("blue", "white", "red")), show_column_names = TRUE)
  }
  
  output$heatmap <- renderPlot({
    req(input$plot_heatmap) 
    plot_heatmap(filter_count(), log_transform = input$log_transform)
  })
  
  # Function to perform PCA
  run_pca <- function(filtered_count, scale_data = TRUE){
    # Remove non-numeric data
    numeric_data <- filtered_count[sapply(filtered_count, is.numeric)]
    pca <- prcomp(numeric_data, scale. = scale_data)
    list(
      pca = pca,
      explained_variance = summary(pca)$importance[2, ] * 100,
      pcs = as.data.frame(pca$x)
    )
  }
  
  pca_results <- reactiveVal(NULL)
  
  # Run PCA and update dropdown menu for principal components to plot
  observeEvent(input$run_pca, {
    req(input$run_pca)
    pca_result <- run_pca(filter_count())
    pca_results(pca_result)
    pcs <- colnames(pca_result$pcs)
    updateSelectInput(session, "x_pc", choices = pcs, selected = pcs[1])
    updateSelectInput(session, "y_pc", choices = pcs, selected = pcs[2])
  })
  
  # Function to generate PCA plot
  plot_pca <- function(pca_result, x_pc = "PC1", y_pc = "PC2", color_column = NULL) {
    pcs <- pca_result$pcs
    explained_variance <- pca_result$explained_variance
    x_label <- paste(x_pc, "(", round(explained_variance[which(names(explained_variance) == x_pc)], 1), "%)", sep = "")
    y_label <- paste(y_pc, "(", round(explained_variance[which(names(explained_variance) == y_pc)], 1), "%)", sep = "")
    pca_plot <- ggplot(pcs, aes_string(x = x_pc, y = y_pc, color = color_column)) +
      geom_point(size = 3, alpha = 0.7) +
      theme_minimal() +
      labs(title = "PCA Projection", x = x_label, y = y_label, color = "Groups") +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
      coord_fixed()
    return(pca_plot)
  }
  
  output$pca <- renderPlot({
    req(input$x_pc, input$y_pc, input$plot_pca)
    plot_pca(
      pca_result = pca_results(),
      x_pc = input$x_pc,
      y_pc = input$y_pc
    )
  })
  
  # DE tab
  
  load_DE <- reactive({
    dataf <- read.csv(req(input$dea_res$datapath), stringsAsFactors = FALSE)
    return(dataf)
  }) 
  
  output$dea_table <- renderDT({
    req(input$submit_dea)
    datatable(load_DE(), options = list(lengthChange = TRUE))
  })
  
  # Function to generate volcano plot
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    out <- dataf %>%
      mutate(volc_plot_status  = case_when(
        is.na(!!sym(y_name)) ~ NA_character_, 
        !!sym(y_name) < 10**slider ~"TRUE",
        !!sym(y_name) >= 10**slider ~"FALSE")) %>%
      ggplot(aes(x=!!sym(x_name), y=-log10(!!sym(y_name)), color=volc_plot_status)) +
      geom_point() +
      scale_color_manual(values=c("TRUE" = color2, "FALSE" = color1), 
                         name=paste0(y_name, " < 1 x 10^", slider)) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.box.just = "center") +
      guides(color = guide_legend(nrow = 1, byrow = TRUE)) 
    return(out)
  }
  
  # Function to generate a table including filtered genes based on a selected p-adjusted value threshold
  draw_table <- function(dataf, slider) {
    out <- dataf %>%
      filter(padj < 10**slider) %>%
      mutate(across(c(baseMean, log2FoldChange, lfcSE, stat), ~ formatC(.x, format = "f", digits = 2))) %>%
      mutate(across(c(pvalue, padj), ~ formatC(.x, format = "e", digits = 3)))
    colnames(out)[1] <- "gene"
    return(out)
  }
  
  output$volcano <- renderPlot({
    req(input$plotButton)
    volcano_plot(load_DE(), input$x_name, input$y_name, input$slider, input$color1, input$color2)
  })
  
  output$de_table <- renderTable({
    req(input$plotButton)
    draw_table(load_DE(), input$slider)
  })
  
  # Gene Set Enrichment Analysis tab
  
  fgsea_res <- reactiveVal(NULL)
  
  # Observe an event when the "submit_fgsea" button is clicked
  observeEvent(input$submit_fgsea, {
    req(input$fgsea_res)
    df <- read.csv(input$fgsea_res$datapath, stringsAsFactors = FALSE)
    
    if("leadingEdge" %in% colnames(df)){
      if(is.list(df$leadingEdge)){
        df$leadingEdge <- sapply(df$leadingEdge, function(x) paste(x, collapse=", "))
      }
    }
    fgsea_res(df) # Store the processed data frame in a reactive variable
  })
  
  # Generate a barplot of fgsea NES for top pathways selected by slider
  output$fgsea_bar <- renderPlot({
    req(fgsea_res())
    top_n <- input$pathway_slider
    df <- fgsea_res() %>%
      arrange(desc(abs(NES))) %>% 
      head(top_n) %>%
      mutate(color_group = NES > 0)  
    ggplot(df, aes(x = reorder(pathway, abs(NES)), y = NES, fill = color_group)) +
      geom_col() +
      coord_flip() +
      theme_minimal() +
      scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
      labs(
        title = paste("Top", top_n, "Pathways by NES"),
        x = "Pathway",
        y = "NES",
        fill = "Sign"
      )
  })
  
  filtered_table <- reactive({
    req(fgsea_res())
    df <- fgsea_res()
    # filter by adjusted p-value
    df <- df %>% filter(padj <= input$filter_padj_tbl)
    # filter by NES sign
    if (input$neg_or_pos == "Positive") {
      df <- df %>% filter(NES > 0)
    } else if (input$neg_or_pos == "Negative") {
      df <- df %>% filter(NES < 0)
    } # No filtering if input$neg_or_pos == "All"
    df
  })
  
  output$fgsea_table <- renderDT({
    req(filtered_table())
    datatable(filtered_table(), options = list(pageLength=10, lengthChange=TRUE))
  })
  
  # Download filtered fgsea result
  output$download_fgsea_filtered <- downloadHandler(
    filename = function(){"fgsea_filtered.csv"},
    content = function(file){
      write.csv(filtered_table(), file, row.names=FALSE)
    }
  )
  
  # Generate a scatter plot of NES
  output$fgsea_plot <- renderPlot({
    req(fgsea_res())
    df <- fgsea_res()
    # filter by padj
    df <- df %>% mutate(logpadj = -log10(padj))
    df <- df %>% mutate(status = ifelse(padj <= input$filter_padj_plt, "In", "Out"))
    ggplot(df, aes(x=NES, y=logpadj, color=status))+
      geom_point(alpha=0.7)+
      scale_color_manual(values=c("In"="red","Out"="grey"))+
      theme_minimal()+
      labs(title="NES vs -log10(padj)", x="NES", y="-log10(padj)", color="")
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
