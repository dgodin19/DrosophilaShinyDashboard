options(shiny.maxRequestSize = 30 * 1024^2)

library(shiny)
library(DT)
library(limma)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(fgsea)

ui <- fluidPage(
  
  sidebarLayout(
    
    sidebarPanel(
      
      fileInput("counts_file","Upload normalized counts matrix CSV"),
      fileInput("sample_metadata","Upload sample metadata CSV"),
      fileInput("raw_counts","Upload raw counts matrix"),
      fileInput("fgsea_results_file","Upload fgsea results CSV"),
      
      conditionalPanel(
        condition = "input.main_tabs == 'counts_tab'",
        
        hr(),
        h4("Counts filtering"),
        
        sliderInput(
          "var_filter",
          "Variance percentile",
          min = 0,
          max = 100,
          value = 50,
          step = 1
        ),
        
        sliderInput(
          "nonzero_filter",
          "Min non-zero samples",
          min = 0,
          max = 80,
          value = 10,
          step = 1
        )
      ),
      
      conditionalPanel(
        condition = "input.main_tabs == 'fgsea_tab'",
        
        hr(),
        h4("fgsea controls"),
        
        sliderInput(
          "num_paths",
          "Number of pathways",
          min = 5,
          max = 50,
          value = 20
        ),
        
        sliderInput(
          "padj_filter",
          "Adjusted p-value threshold",
          min = 0,
          max = 0.25,
          value = 0.05
        ),
        
        radioButtons(
          "nes_direction",
          "NES direction",
          choices = c("All","Positive","Negative")
        ),
        
        downloadButton("download_fgsea")
      )
      
    ),
    
    mainPanel(
      
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel(
          "Sample Information Exploration",
          
          tabsetPanel(
            
            tabPanel(
              "Summary",
              
              fluidRow(
                column(9, tableOutput("sample_summary"))
              ),
              
              fluidRow(
                column(3,"Number of rows:"),
                column(9,textOutput("sample_rows"))
              ),
              
              fluidRow(
                column(3,"Number of columns:"),
                column(9,textOutput("sample_cols"))
              )
            ),
            
            tabPanel(
              "Sample Data Table",
              fluidRow(
                column(12,DTOutput("sample_table"))
              )
            ),
            
            tabPanel(
              "Variable Visualization",
              
              fluidRow(
                column(
                  4,
                  selectInput(
                    "variable",
                    "Select variable",
                    choices = c("treatment","sex","timepoint","lifestage","All"),
                    selected = "treatment"
                  )
                )
              ),
              
              fluidRow(
                column(12,plotOutput("sample_barplot"))
              )
            )
          )
        ),
        
        tabPanel(
          "Counts Matrix Exploration",
          value = "counts_tab",
          
          tabsetPanel(
            
            tabPanel(
              "Counts Table",
              fluidRow(
                column(12,DTOutput("counts_table"))
              )
            ),
            
            tabPanel(
              "Filtering Summary",
              fluidRow(
                column(12,tableOutput("summary_of_counts_filtering"))
              )
            ),
            
            tabPanel(
              "Filtering Diagnostics",
              
              fluidRow(
                column(6,plotOutput("median_variance_plot")),
                column(6,plotOutput("median_zero_plot"))
              )
            ),
            
            tabPanel(
              "Clustered Heatmap",
              
              fluidRow(
                column(12,plotOutput("counts_heatmap",height="700px"))
              )
            ),
            
            tabPanel(
              "PCA",
              
              fluidRow(
                
                column(
                  3,
                  selectInput(
                    "pca_x",
                    "X-axis PC",
                    choices=paste0("PC",1:10),
                    selected="PC1"
                  )
                ),
                
                column(
                  3,
                  selectInput(
                    "pca_y",
                    "Y-axis PC",
                    choices=paste0("PC",1:10),
                    selected="PC2"
                  )
                ),
                
                column(
                  3,
                  selectInput(
                    "pca_color",
                    "Color by",
                    choices=c("treatment","sex","timepoint","lifestage"),
                    selected="treatment"
                  )
                )
              ),
              
              fluidRow(
                column(12,plotOutput("pca_scatterplot",height="700px"))
              )
            )
            
          )
        ),
        
        tabPanel(
          "Differential Expression",
          
          tabsetPanel(
            
            tabPanel(
              "Results Table",
              
              fluidRow(
                column(
                  4,
                  selectInput(
                    "design_matrix",
                    "Select design",
                    choices=c(
                      "treatment","sex","timepoint","lifestage",
                      "treatment + sex",
                      "treatment + timepoint",
                      "treatment + lifestage",
                      "treatment + sex + timepoint",
                      "treatment + sex + lifestage",
                      "treatment + timepoint + lifestage",
                      "sex + timepoint + lifestage",
                      "timepoint + lifestage",
                      "timepoint + sex",
                      "sex + lifestage",
                      "treatment + sex + timepoint + lifestage"
                    ),
                    selected="treatment"
                  )
                )
              ),
              
              fluidRow(
                column(12,DTOutput("differential_expression"))
              )
            ),
            
            tabPanel(
              "P-value Histogram",
              fluidRow(
                column(12,plotOutput("pvalue_histogram",height="700px"))
              )
            ),
            
            tabPanel(
              "Log2FC Distribution",
              fluidRow(
                column(12,plotOutput("log2fc",height="700px"))
              )
            ),
            
            tabPanel(
              "MA Plot",
              fluidRow(
                column(12,plotOutput("ma_plot",height="700px"))
              )
            ),
            
            tabPanel(
              "Volcano Plot",
              fluidRow(
                column(12,plotOutput("volcano_plot",height="700px"))
              )
            )
            
          )
        ),
        
        tabPanel(
          "Gene Set Enrichment Analysis",
          value = "fgsea_tab",
          
          tabsetPanel(
            
            tabPanel(
              "Top Pathways",
              
              fluidRow(
                column(
                  12,
                  plotOutput("fgsea_barplot",click="pathway_click")
                )
              ),
              
              fluidRow(
                column(
                  12,
                  tableOutput("selected_pathway")
                )
              )
            ),
            
            tabPanel(
              "Results Table",
              
              fluidRow(
                column(12,DTOutput("fgsea_table"))
              )
            ),
            
            tabPanel(
              "NES vs Significance",
              
              fluidRow(
                column(12,plotOutput("fgsea_scatter",height="700px"))
              )
            )
            
          )
        )
        
      )
      
    )
    
  )
  
)

server <- function(input, output, session) {
  
  # Inputs
  
  counts_data <- reactive({
    req(input$counts_file)
    
    file_ext <- tools::file_ext(input$counts_file$name)
    
    validate(
      need(
        file_ext %in% c("csv", "tsv"),
        "Counts file must be a CSV or TSV file."
      )
    )
    
    counts_df <- tryCatch({
      
      if (file_ext == "csv") {
        read.csv(input$counts_file$datapath, check.names = FALSE)
      } else {
        read.delim(input$counts_file$datapath, check.names = FALSE)
      }
      
    }, error = function(e) {
      validate(
        need(FALSE, "Counts file could not be read. Ensure it is properly formatted CSV/TSV.")
      )
    })
    
    validate(
      need(
        ncol(counts_df) == 81,
        paste(
          "Counts matrix must contain 1 gene column + 80 samples. Found",
          ncol(counts_df)
        )
      )
    )
    
    validate(
      need(
        tolower(trimws(colnames(counts_df)[1])) %in%
          c("gene", "gene_id", "geneid", "ensembl_gene_id", "x", ""),
        paste("First column must contain gene IDs. Found:", colnames(counts_df)[1])
      )
    )
    
    colnames(counts_df)[1] <- "gene_id"
    
    counts_df
  })
  
  
  sample_metadata <- reactive({
    req(input$sample_metadata)
    
    file_ext <- tools::file_ext(input$sample_metadata$name)
    
    validate(
      need(
        file_ext %in% c("csv", "tsv"),
        "Sample metadata must be a CSV or TSV file."
      )
    )
    
    sample_df <- tryCatch({
      
      if (file_ext == "csv") {
        read.csv(input$sample_metadata$datapath, check.names = FALSE)
      } else {
        read.delim(input$sample_metadata$datapath, check.names = FALSE)
      }
      
    }, error = function(e) {
      validate(
        need(FALSE, "Sample metadata file is not a properly formatted CSV/TSV.")
      )
    })
    
    validate(
      need(
        nrow(sample_df) == 80,
        paste(
          "Sample metadata must contain exactly 80 sample rows. Found",
          nrow(sample_df)
        )
      )
    )
    
    required_cols <- c(
      "sample",
      "treatment",
      "sex",
      "timepoint",
      "lifestage",
      "replicate"
    )
    
    validate(
      need(
        all(required_cols %in% colnames(sample_df)),
        paste(
          "Metadata must contain the following columns:",
          paste(required_cols, collapse = ", ")
        )
      )
    )
    
    sample_df
  })
  
  
  raw_counts <- reactive({
    req(input$raw_counts)
    
    file_ext <- tools::file_ext(input$raw_counts$name)
    
    validate(
      need(
        file_ext %in% c("csv", "tsv"),
        "Raw counts file must be a CSV or TSV file."
      )
    )
    
    counts_matrix <- tryCatch({
      
      if (file_ext == "csv") {
        read.csv(input$raw_counts$datapath, check.names = FALSE)
      } else {
        read.delim(input$raw_counts$datapath, check.names = FALSE)
      }
      
    }, error = function(e) {
      validate(
        need(FALSE, "Raw counts file is not a properly formatted CSV/TSV.")
      )
    })
    
    counts_matrix
  })
  
  
  # Sample Metadata Functions 
  
  sample_metadata_summary_table <- function(sample_metadata) {
    
    num_rows_whole <- nrow(sample_metadata)
    num_cols_whole <- ncol(sample_metadata)
    
    summary_df <- data.frame(
      Column_Name = character(),
      Type = character(),
      Distinct_Values = character(),
      stringsAsFactors = FALSE
    )
    
    for (col in colnames(sample_metadata)) {
      
      values <- unique(sample_metadata[[col]])
      
      summary_df <- rbind(
        summary_df,
        data.frame(
          Column_Name = col,
          Type = class(sample_metadata[[col]])[1],
          Distinct_Values = paste(values, collapse = ", "),
          stringsAsFactors = FALSE
        )
      )
    }
    
    result <- list(
      num_rows = num_rows_whole,
      num_cols = num_cols_whole,
      column_summary = summary_df
    )
    
    return(result)
  }
  
  
  sample_summary_data <- reactive({
    req(sample_metadata())
    sample_metadata_summary_table(sample_metadata())
  })
  
  
  output$sample_summary <- renderTable({
    sample_summary_data()$column_summary
  })
  
  
  output$sample_rows <- renderText({
    sample_summary_data()$num_rows
  })
  
  
  output$sample_cols <- renderText({
    sample_summary_data()$num_cols
  })
  
  
  output$sample_table <- DT::renderDT({
    req(sample_metadata())
    
    datatable(
      sample_metadata(),
      options = list(pageLength = 10)
    )
  })
  
  
  sample_barplot <- function(sample_metadata, variable) {
    
    if (variable == "All") {
      
      df_long <- stack(
        sample_metadata[, c("treatment", "sex", "timepoint", "lifestage")]
      )
      
      result <- ggplot(df_long, aes(x = values, fill = ind)) +
        geom_bar(position = "dodge") +
        labs(
          x = "Category",
          y = "Count",
          fill = "Variable",
          title = "Sample metadata distribution"
        )
      
    } else {
      
      result <- ggplot(sample_metadata, aes_string(x = variable, fill = variable)) +
        geom_bar() +
        labs(
          x = variable,
          y = "Count",
          title = paste("Barplot of", variable),
          fill = variable
        )
    }
    
    return(result)
  }
  
  
  output$sample_barplot <- renderPlot({
    req(sample_metadata(), input$variable)
    
    sample_barplot(
      sample_metadata(),
      input$variable
    )
  })
  
  
  # Counts Functions
  
  create_counts_table <- function(counts_df, var_filter, nonzero_filter) {
    
    gene_ids <- counts_df[, 1]
    
    counts <- counts_df[, -1]
    counts <- as.data.frame(lapply(counts, as.numeric))
    
    gene_var <- apply(counts, 1, var)
    
    var_threshold <- quantile(
      gene_var,
      probs = var_filter / 100
    )
    
    nonzero_counts <- rowSums(counts > 0)
    
    keep <- gene_var >= var_threshold &
      nonzero_counts >= nonzero_filter
    
    filtered <- counts_df[keep, ]
    
    filtered
  }
  
  
  display_filtering_effects <- function(filtered, original) {
    
    num_samples <- ncol(original) - 1
    total_genes <- nrow(original)
    
    num_filtered_genes <- nrow(filtered)
    num_unfiltered_genes <- total_genes - num_filtered_genes
    
    percent_passing_genes <- (num_filtered_genes / total_genes) * 100
    percent_failing_genes <- (num_unfiltered_genes / total_genes) * 100
    
    result <- data.frame(
      Metric = c(
        "Number of samples",
        "Total genes",
        "Genes passing filter",
        "Genes failing filter"
      ),
      Value = c(
        num_samples,
        total_genes,
        num_filtered_genes,
        num_unfiltered_genes
      ),
      Percent = c(
        NA,
        NA,
        round(percent_passing_genes, 2),
        round(percent_failing_genes, 2)
      )
    )
    
    result
  }
  
  
  compute_gene_stats <- function(counts_df, var_filter, nonzero_filter) {
    
    counts <- counts_df[, -1]
    counts <- as.data.frame(lapply(counts, as.numeric))
    
    gene_median <- apply(counts, 1, median)
    gene_var <- apply(counts, 1, var)
    
    zero_counts <- rowSums(counts == 0)
    nonzero_counts <- rowSums(counts > 0)
    
    var_threshold <- quantile(gene_var, probs = var_filter / 100)
    
    pass <- gene_var >= var_threshold &
      nonzero_counts >= nonzero_filter
    
    data.frame(
      median = gene_median,
      variance = gene_var,
      zeros = zero_counts,
      pass_filter = pass
    )
  }
  
  
  gene_stats <- reactive({
    req(counts_data(), input$var_filter, input$nonzero_filter)
    
    compute_gene_stats(
      counts_data(),
      input$var_filter,
      input$nonzero_filter
    )
  })
  
  
  observe({
    req(counts_data())
    print(head(counts_data()))
  })
  
  
  counts_filtered <- reactive({
    req(counts_data(), input$var_filter, input$nonzero_filter)
    
    create_counts_table(
      counts_data(),
      input$var_filter,
      input$nonzero_filter
    )
  })
  
  
  output$counts_table <- DT::renderDT({
    
    req(counts_filtered())
    
    DT::datatable(
      head(counts_filtered(), 20),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })
  
  
  filter_summary <- reactive({
    
    req(counts_data(), counts_filtered())
    
    display_filtering_effects(
      filtered = counts_filtered(),
      original = counts_data()
    )
  })
  
  clustered_heatmap <- function(filtered, log_transform = TRUE){
    
    gene_ids <- filtered[,1]
    
    counts <- filtered[,-1]
    counts <- as.data.frame(lapply(counts, as.numeric))
    
    rownames(counts) <- gene_ids
    
    if(log_transform){
      counts <- log2(counts + 1)
    }
    
    gene_var <- apply(counts, 1, var)
    
    top_idx <- order(gene_var, decreasing = TRUE)[1:30]
    
    top30 <- counts[top_idx, ]
    
    top30_z <- t(scale(t(top30)))
    
    pheatmap(
      top30_z,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 6,
      fontsize_col = 7,
      color = colorRampPalette(c("navy","white","firebrick3"))(50),
      main = "Top 30 Variable Genes (Z-scored)"
    )
  }
  
  pca_scatterplot <- function(filtered_counts, sample_metadata, x_pc, y_pc, color_var){
    
    counts <- filtered_counts[,-1]
    counts <- as.data.frame(lapply(counts, as.numeric))
    
    counts_t <- t(counts)
    
    pca_results <- prcomp(counts_t, scale. = TRUE)
    
    var_explained <- (pca_results$sdev^2) / sum(pca_results$sdev^2)
    
    pca_df <- as.data.frame(pca_results$x)
    pca_df$sample <- rownames(pca_df)
    
    plot_df <- merge(pca_df, sample_metadata, by = "sample")
    
    x_index <- as.numeric(sub("PC","",x_pc))
    y_index <- as.numeric(sub("PC","",y_pc))
    
    ggplot(
      plot_df,
      aes(
        x = .data[[x_pc]],
        y = .data[[y_pc]],
        color = .data[[color_var]]
      )
    ) +
      geom_point(size = 3) +
      labs(
        x = paste0(x_pc," (",round(var_explained[x_index]*100,2),"%)"),
        y = paste0(y_pc," (",round(var_explained[y_index]*100,2),"%)"),
        color = color_var,
        title = "PCA of Filtered Counts"
      )
  }
  
  
  # Differential Expression
  
  run_limma <- function(counts_matrix, sample_metadata, design_matrix) {
    
    gene_ids <- counts_matrix[, 1]
    
    counts <- counts_matrix[, -1]
    counts <- as.data.frame(lapply(counts, as.numeric))
    rownames(counts) <- gene_ids
    
    sample_metadata <- droplevels(sample_metadata)
    
    design_formula <- as.formula(paste("~", design_matrix))
    
    design <- model.matrix(design_formula, data = sample_metadata)
    
    dge <- DGEList(counts = counts)
    
    keep <- filterByExpr(dge, design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    dge <- calcNormFactors(dge)
    
    v <- voom(dge, design, plot = FALSE)
    
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    
    result <- topTable(
      fit,
      coef = ncol(design),
      number = Inf,
      sort.by = "P"
    )
    
    result$gene_id <- rownames(result)
    
    result <- result[, c(      "logFC",      "AveExpr",      "t",      "P.Value",      "adj.P.Val",      "B"    )]
    
    return(result)
  }
  
  
  # Plots 
  
  plot_pvals <- function(results) {
    
    ggplot(results, aes(x = P.Value)) +
      geom_histogram(
        binwidth = 0.02,
        fill = "steelblue",
        color = "black"
      ) +
      theme_minimal() +
      labs(
        title = "Histogram of limma p-values",
        x = "P-value",
        y = "Count"
      )
  }
  
  
  plot_log2fc <- function(results) {
    
    ggplot(results, aes(x = logFC)) +
      geom_histogram(
        binwidth = 0.25,
        fill = "steelblue",
        color = "black"
      ) +
      theme_minimal() +
      labs(
        title = "Distribution of log2 Fold Changes",
        x = "log2 Fold Change",
        y = "Count"
      )
  }
  
  
  plot_ma <- function(results) {
    
    ggplot(results, aes(x = AveExpr, y = logFC)) +
      geom_point(alpha = 0.4) +
      geom_hline(yintercept = 0, color = "red") +
      theme_minimal() +
      labs(
        title = "MA Plot",
        x = "Average Expression",
        y = "log2 Fold Change"
      )
  }
  
  
  plot_volcano <- function(results) {
    
    results$log10padj <- -log10(results$adj.P.Val)
    
    results$volc_plot_status <- "NS"
    
    results$volc_plot_status[      results$adj.P.Val < 0.05 & results$logFC > 1    ] <- "UP"
    
    results$volc_plot_status[      results$adj.P.Val < 0.05 & results$logFC < -1    ] <- "DOWN"
    
    ggplot(results, aes(
      x = logFC,
      y = log10padj,
      color = volc_plot_status
    )) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c(
        "UP" = "red",
        "DOWN" = "blue",
        "NS" = "gray"
      )) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_minimal() +
      labs(
        title = "Volcano Plot",
        x = "log2 Fold Change",
        y = "-log10 Adjusted P-value"
      )
  }
  
  
  # FGSEA
  
  fgsea_results <- reactive({
    req(input$fgsea_results_file)
    
    df <- read.csv(input$fgsea_results_file$datapath)
    
    df
  })
  
  
  top_pathways <- function(fgsea_results, num_paths) {
    
    top_paths <- fgsea_results %>%
      arrange(padj) %>%
      slice_head(n = num_paths) %>%
      mutate(
        direction = ifelse(NES > 0, "Up", "Down"),
        pathway = factor(pathway, levels = pathway[order(NES)])
      )
    
    ggplot(top_paths, aes(x = pathway, y = NES, fill = direction)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c(
        "Up" = "firebrick",
        "Down" = "steelblue"
      )) +
      theme_minimal() +
      labs(
        title = "Top enriched pathways",
        x = "Pathway",
        y = "Normalized Enrichment Score"
      )
  }
  
  
  filtered_fgsea <- reactive({
    
    df <- fgsea_results()
    
    df <- df[df$padj <= input$padj_filter, ]
    
    if (input$nes_direction == "Positive") {
      df <- df[df$NES > 0, ]
    }
    
    if (input$nes_direction == "Negative") {
      df <- df[df$NES < 0, ]
    }
    
    df
  })
  
  
  fgsea_scatter <- function(df, padj_threshold) {
    
    df$significant <- df$padj <= padj_threshold
    
    ggplot(df, aes(
      x = NES,
      y = -log10(padj),
      color = significant
    )) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("grey70", "red")) +
      theme_minimal() +
      labs(
        x = "Normalized Enrichment Score",
        y = "-log10 adjusted p-value",
        title = "FGSEA pathway significance"
      )
  }
  
  
  # Outputs 
  
  output$summary_of_counts_filtering <- renderTable({
    filter_summary()
  })
  
  
  output$median_variance_plot <- renderPlot({
    
    df <- gene_stats()
    
    plot(
      df$median,
      df$variance,
      col = ifelse(df$pass_filter, "black", "grey80"),
      pch = 16,
      log = "xy",
      xlab = "Median expression",
      ylab = "Variance"
    )
  })
  
  
  output$median_zero_plot <- renderPlot({
    
    df <- gene_stats()
    
    plot(
      df$median,
      df$zeros,
      col = ifelse(df$pass_filter, "black", "grey80"),
      pch = 16,
      log = "x",
      xlab = "Median expression",
      ylab = "Number of zero samples"
    )
  })
  
  
  output$counts_heatmap <- renderPlot({
    req(counts_filtered())
    
    clustered_heatmap(
      counts_filtered(),
      log_transform = TRUE
    )
  })
  
  
  output$pca_scatterplot <- renderPlot({
    
    req(counts_filtered(), sample_metadata())
    
    pca_scatterplot(
      counts_filtered(),
      sample_metadata(),
      input$pca_x,
      input$pca_y,
      input$pca_color
    )
  })
  
  
  de_results <- reactive({
    
    req(raw_counts(), sample_metadata(), input$design_matrix)
    
    run_limma(
      raw_counts(),
      sample_metadata(),
      input$design_matrix
    )
  })
  
  
  output$differential_expression <- DT::renderDT({
    
    req(de_results())
    
    DT::datatable(
      de_results(),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })
  
  
  output$pvalue_histogram <- renderPlot({
    req(de_results())
    plot_pvals(de_results())
  })
  
  
  output$log2fc <- renderPlot({
    req(de_results())
    plot_log2fc(de_results())
  })
  
  
  output$ma_plot <- renderPlot({
    req(de_results())
    plot_ma(de_results())
  })
  
  
  output$volcano_plot <- renderPlot({
    req(de_results())
    plot_volcano(de_results())
  })
  
  
  output$fgsea_barplot <- renderPlot({
    req(fgsea_results())
    top_pathways(fgsea_results(), input$num_paths)
  })
  
  
  output$selected_pathway <- renderTable({
    
    req(input$pathway_click)
    
    df <- fgsea_results()
    
    df %>%
      arrange(padj) %>%
      slice(input$pathway_click$y)
  })
  
  
  output$fgsea_scatter <- renderPlot({
    
    req(fgsea_results())
    
    fgsea_scatter(
      fgsea_results(),
      input$padj_filter
    )
  })
  
  
  output$fgsea_table <- renderDT({
    
    req(filtered_fgsea())
    
    datatable(filtered_fgsea())
  })
  
  
  output$download_fgsea <- downloadHandler(
    
    filename = function() {
      "filtered_fgsea_results.csv"
    },
    
    content = function(file) {
      write.csv(filtered_fgsea(), file, row.names = FALSE)
    }
  )
  
}

shinyApp(ui, server)