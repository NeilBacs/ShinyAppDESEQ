
library(shiny)
library(DESeq2)
library(tidyverse)
library(DT)

#
# ------------------------------ UI SECTION ------------------------------
#

ui <- fluidPage(
  titlePanel("DESeq2 RNA-seq Analysis Shiny App"),
  
  sidebarLayout(
    sidebarPanel(
      # User selects which count matrix to load, added dummy choice
      selectInput("count_matrix", "Choose dataset:",
                  choices = c("GSE164073", "empty count")),
      
      # User selects which metadata file to load,added dummy choice
      selectInput("metadata", "Choose dataset:",
                  choices = c("metadata.csv","empty metadata")),
      
      # Placeholder UI element for infection reference dropdown
      uiOutput("infection_ui"),
      
      # Button to run DESeq2
      actionButton("runDE", "Run DESeq2")
    ),
    
    mainPanel(
      tabsetPanel(
        # MA plot results
        tabPanel("MA Plot", plotOutput("maPlot")),
        # Top significant genes table
        tabPanel("Top Genes", DTOutput("topGenes")),
        # Text summary of DESeq2 run
        tabPanel("DESeq2 Summary", verbatimTextOutput("summary"))
      )
    )
  )
)

#
# ------------------------------ SERVER SECTION ------------------------------
#

server <- function(input, output, session) {
  
  # ---------------- Load count matrix based on dropdown -----------------
  counts_data <- reactive({
    req(input$count_matrix)
    
    # If user selected empty dataset, show text error
    validate(
      need(input$count_matrix == "GSE164073", 
           "❌ No count matrix selected. Please choose a valid dataset.")
    )
    
    read.csv("GSE164073/GSE164073_Eye_count_matrix.csv.gz",
             row.names = 1)
  })
  
  
  
  # ---------------- Load metadata based on dropdown -----------------
  meta_data <- reactive({
    req(input$metadata)
    
    validate(
      need(input$metadata == "metadata.csv", 
           "❌ No metadata selected. Please choose a valid metadata file.")
    )
    
    metadata <- read.csv("GSE164073/metadata.csv", row.names = 1)
    rownames(metadata) <- metadata$description 
    
    metadata <- metadata[, c("tissue.ch1", "infection.ch1")]
    colnames(metadata) <- c("tissue", "infection")
    
    metadata$infection <- factor(metadata$infection)
    metadata
  })
  

  
  
  # ------------- Dynamic UI: infection reference level dropdown ----------
  output$infection_ui <- renderUI({
    req(meta_data())
    
    selectInput(
      "refLevel",
      "Reference level for infection:",
      # Levels of infection factor become choices
      choices = levels(meta_data()$infection),
      selected = levels(meta_data()$infection)[1]  # default selection
    )
  })
  
  
  # ---------------------- RUN DESEQ2 WHEN BUTTON IS CLICKED -------------------
  deseq_res <- eventReactive(input$runDE, {
    
    # Make sure all required inputs are available
    req(counts_data(), meta_data(), input$refLevel)
    
    # Copy metadata and change reference level of infection factor
    meta <- meta_data()
    meta$infection <- relevel(meta$infection, ref = input$refLevel)
    
    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(
      countData = counts_data(),
      colData = meta,
      design = ~ infection
    )
    
    # Remove genes with extremely low counts
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    
    # Run DESeq2 normalization + dispersion + differential expression
    dds <- DESeq(dds)
    
    # Return results table
    results(dds, alpha = 0.05)
  })
  
  
  # -------------------------- MA PLOT OUTPUT -------------------------
  output$maPlot <- renderPlot({
    req(deseq_res())
    plotMA(deseq_res())
  })
  
  
  # -------------------------- TOP GENES TABLE -------------------------
  output$topGenes <- renderDT({
    req(deseq_res())
    
    res <- deseq_res()
    
    # Filter for statistically significant genes (padj < 0.05)
    sig_genes <- subset(res, padj < 0.05)
    
    # Convert to dataframe for table
    sig_df <- data.frame(
      gene = rownames(sig_genes),
      baseMean = sig_genes$baseMean,
      log2FC = sig_genes$log2FoldChange,
      padj = sig_genes$padj
    )
    
    # Sort by padj and take top 100
    top_genes <- sig_df[order(sig_df$padj), ] |> head(100)
    
    # Display using DataTable
    datatable(top_genes)
  })
  
  
  # -------------------------- SUMMARY TEXT -------------------------
  output$summary <- renderPrint({
    req(deseq_res())
    summary(deseq_res())
  })
}

# -------------------------- RUN SHINY APP -------------------------
shinyApp(ui = ui, server = server)
