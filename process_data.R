## Author: Yuki Ito
## yuki@bu.edu

## This Shiny application provides a comprehensive workflow for processing, analyzing, and visualizing gene expression data. It features the following functionalities:

## 1. Series Matrix to CSV: Upload GEO series matrix files, extract sample metadata, and convert them into a structured CSV format. Preview extracted sample characteristics before downloading.
## 2. TSV to CSV Conversion: Transform tab-delimited TSV files into CSV format, ensuring compatibility with other downstream analyses. Preview the converted table to validate the data.
## 3. Differential Expression Analysis (DESeq2): Perform DESeq2-based differential expression analysis. Users can upload raw count matrices and sample metadata, define group comparisons, and generate a results table with log2 fold changes, p-values, and adjusted p-values (padj). Results can be previewed and downloaded.
## 4. Gene Set Enrichment Analysis (FGSEA): Perform FGSEA using differential expression results and gene set files in GMT format. Users can rank genes based on DESeq2 statistics (e.g., log2 fold change) and identify enriched pathways. Results, including adjusted p-values, normalized enrichment scores (NES), and leading-edge genes, are available for preview and download.

## Through an intuitive interface, ExpressionExplorer allows users to preprocess data, apply filters, and generate key insights with diagnostic plots and enrichment analyses. This app supports seamless uploads, customizable parameters, and downloadable outputs to streamline transcriptomic workflows.

options(shiny.maxRequestSize = 1024 * 1024 * 1024)

library(shiny)
library(dplyr)
library(DESeq2)
library(fgsea)

ui <- fluidPage(
  titlePanel("Series Matrix, TSV Converter, DESeq2, and FGSEA"),
  tabsetPanel(
    tabPanel("Series Matrix to CSV",
             sidebarLayout(
               sidebarPanel(
                 fileInput("matrix_file", "Upload Series Matrix File (.txt)", accept = ".txt"),
                 actionButton("process", "Process File", style = "width: 100%; margin-top: 10px;"),
                 downloadButton("download_csv", "Download Sample Information")
               ),
               mainPanel(
                 h4("Preview of Sample Information"),
                 tableOutput("table_preview"),
                 verbatimTextOutput("debug_output")
               )
             )
    ),
    tabPanel("TSV to CSV",
             sidebarLayout(
               sidebarPanel(
                 fileInput("tsv_file", "Upload TSV File (.tsv)", accept = c(".tsv")),
                 actionButton("process_tsv", "Process TSV", style = "width: 100%; margin-top: 10px;"),
                 downloadButton("download_tsv_csv", "Download TSV as CSV")
               ),
               mainPanel(
                 h4("Preview of TSV Data"),
                 tableOutput("tsv_table_preview")
               )
             )
    ),
    tabPanel("DE (DESeq2)",
             sidebarLayout(
               sidebarPanel(
                 fileInput("de_count", "Upload raw count (.csv)", accept = ".csv"),
                 fileInput("de_coldata", "Upload colData (.csv)", accept = ".csv"),
                 selectInput("de_group", "Group column for DE", choices = NULL),
                 actionButton("run_de", "Run DESeq2", style="width:100%; margin-top:10px;"),
                 downloadButton("download_de_results", "Download DE results")
               ),
               mainPanel(
                 h4("DE Results"),
                 tableOutput("de_table")
               )
             )
    ),
    tabPanel("FGSEA",
             sidebarLayout(
               sidebarPanel(
                 fileInput("fgsea_gmt", "Upload GMT file", accept = ".gmt"),
                 fileInput("fgsea_de", "Upload DESeq2 result (.csv)", accept = ".csv"),
                 selectInput("fgsea_stat_col", "Select stat column for ranking:", choices = NULL),
                 actionButton("run_fgsea", "Run fgsea", style="width:100%; margin-top:10px;"),
                 downloadButton("download_fgsea", "Download FGSEA Results")
               ),
               mainPanel(
                 h4("FGSEA Results"),
                 tableOutput("fgsea_table")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  parsed_data <- reactiveVal(NULL)
  debug_messages <- reactiveVal("")
  
  # Series Matrix processing function
  process_file <- function(filepath) {
    lines <- readLines(filepath)
    characteristics_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
    
    if (length(characteristics_lines) == 0) {
      stop("No !Sample_characteristics_ch1 lines found. Check the file format.")
    }
    
    characteristics_values <- lapply(characteristics_lines, function(line) {
      vals <- unlist(strsplit(line, "\t"))
      if (length(vals) <= 1) return(NULL)
      vals <- vals[-1]
      gsub("^\"|\"$", "", vals)
    })
    
    characteristics_values <- Filter(Negate(is.null), characteristics_values)
    if (length(characteristics_values) == 0) {
      stop("No valid characteristics values found. Check the file format.")
    }
    
    lengths_check <- sapply(characteristics_values, length)
    if (length(unique(lengths_check)) != 1) {
      stop("Characteristics columns do not have consistent lengths.")
    }
    
    characteristics_matrix <- do.call(rbind, characteristics_values)
    
    extract_features <- function(values) {
      parsed <- lapply(values, function(value) {
        if (!nzchar(value) || !grepl(":", value)) return(NULL)
        parts <- unlist(strsplit(value, ":\\s*"))
        if (length(parts) == 2) parts else NULL
      })
      parsed <- Filter(Negate(is.null), parsed)
      if (length(parsed) == 0) return(NULL)
      do.call(rbind, parsed)
    }
    
    feature_data <- lapply(seq_len(ncol(characteristics_matrix)), function(i) {
      extract_features(characteristics_matrix[, i])
    })
    feature_data <- Filter(Negate(is.null), feature_data)
    if (length(feature_data) == 0) {
      stop("No valid 'Key: Value' pairs extracted.")
    }
    
    combined_features <- do.call(rbind, feature_data)
    unique_features <- unique(combined_features[, 1])
    
    result <- data.frame(matrix(NA, nrow = ncol(characteristics_matrix), ncol = length(unique_features)))
    colnames(result) <- unique_features
    
    for (i in seq_along(feature_data)) {
      parsed <- feature_data[[i]]
      for (j in seq_len(nrow(parsed))) {
        key <- parsed[j, 1]
        val <- parsed[j, 2]
        result[i, key] <- val
      }
    }
    
    return(result)
  }
  
  observeEvent(input$process, {
    req(input$matrix_file)
    tryCatch({
      parsed <- process_file(input$matrix_file$datapath)
      parsed_data(parsed)
      showNotification("File processed successfully.", type = "message")
    }, error = function(e) {
      debug_messages(paste(debug_messages(), "Error:", e$message, sep = "\n"))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$table_preview <- renderTable({
    req(parsed_data())
    head(parsed_data(), 10)
  })
  
  output$debug_output <- renderText({
    debug_messages()
  })
  
  output$download_csv <- downloadHandler(
    filename = function() { "sample_information.csv" },
    content = function(file) {
      write.csv(parsed_data(), file, row.names = FALSE)
    }
  )
  
  # TSV -> CSV
  tsv_data <- reactiveVal(NULL)
  
  observeEvent(input$process_tsv, {
    req(input$tsv_file)
    df <- read.csv(input$tsv_file$datapath, sep = "\t", stringsAsFactors = FALSE)
    tsv_data(df)
    showNotification("TSV processed successfully.", type = "message")
  })
  
  output$tsv_table_preview <- renderTable({
    req(tsv_data())
    head(tsv_data(), 10)
  })
  
  output$download_tsv_csv <- downloadHandler(
    filename = function() { "converted_from_tsv.csv" },
    content = function(file) {
      write.csv(tsv_data(), file, row.names = FALSE)
    }
  )
  
  # DE (DESeq2)
  de_counts <- reactiveVal(NULL)
  de_coldata <- reactiveVal(NULL)
  de_results <- reactiveVal(NULL)
  
  observeEvent(input$de_count, {
    req(input$de_count)
    df <- read.csv(input$de_count$datapath, row.names = 1, stringsAsFactors = FALSE)
    de_counts(df)
  })
  
  observeEvent(input$de_coldata, {
    req(input$de_coldata)
    cd <- read.csv(input$de_coldata$datapath, row.names = 1, stringsAsFactors = FALSE)
    de_coldata(cd)
    updateSelectInput(session, "de_group", choices = colnames(cd))
  })
  
  observeEvent(input$run_de, {
    req(de_counts(), de_coldata(), input$de_group)
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(de_counts()),
                                  colData = de_coldata(),
                                  design = as.formula(paste0("~", input$de_group)))
    dds <- DESeq(dds)
    res <- results(dds)
    res_df <- as.data.frame(res)
    res_df$GeneID <- rownames(res_df)
    # rearrange columns so GeneID first
    res_df <- res_df[, c("GeneID", colnames(res_df)[colnames(res_df)!="GeneID"])]
    de_results(res_df)
    showNotification("DESeq2 analysis completed.", type = "message")
    
    stat_cols <- c("log2FoldChange", "stat", "pvalue", "padj")
    stat_cols <- stat_cols[stat_cols %in% colnames(res_df)]
    # updateSelectInput(session, "fgsea_stat_col", choices = stat_cols)
  })
  
  output$de_table <- renderTable({
    req(de_results())
    head(de_results(), 10)
  })
  
  output$download_de_results <- downloadHandler(
    filename = function(){ "DE_results.csv" },
    content = function(file) {
      write.csv(de_results(), file, row.names = FALSE)
    }
  )
  
  # FGSEA tab
  fgseaRes <- reactiveVal(NULL)
  fgsea_data <- reactiveVal(NULL)
  
  observeEvent(input$fgsea_de, {
    req(input$fgsea_de)
    fg_df <- read.csv(input$fgsea_de$datapath, stringsAsFactors = FALSE)
    fgsea_data(fg_df)
    stat_cols <- c("log2FoldChange", "stat", "pvalue", "padj")
    stat_cols <- stat_cols[stat_cols %in% colnames(fg_df)]
    updateSelectInput(session, "fgsea_stat_col", choices = stat_cols)
  })
  
  observeEvent(input$run_fgsea, {
    req(fgsea_data(), input$fgsea_gmt, input$fgsea_stat_col)
    pathways <- gmtPathways(input$fgsea_gmt$datapath)
    fg_df <- fgsea_data()
    ranks <- setNames(fg_df[[input$fgsea_stat_col]], fg_df$GeneID)
    ranks <- ranks[!is.na(ranks)] 
    
    fgsea_out <- fgsea(pathways = pathways,
                       stats = ranks,
                       minSize = 15,
                       maxSize = 500,
                       nperm = 1000)
    fgsea_out <- fgsea_out %>% arrange(padj)
    fgsea_out$leadingEdge <- sapply(fgsea_out$leadingEdge, function(x) paste(x, collapse = ", "))
    fgseaRes(fgsea_out)
    showNotification("FGSEA analysis completed.", type = "message")
  })
  
  output$fgsea_table <- renderTable({
    req(fgseaRes())
    head(fgseaRes(), 10)
  })
  
  output$download_fgsea <- downloadHandler(
    filename = function(){"fgsea_results.csv"},
    content = function(file) {
      write.csv(fgseaRes(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
