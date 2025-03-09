library(shiny)
library(bslib)
library(reactable)
library(DESeq2)
library(plotly)
library(htmlwidgets)

extract_annotations <- function(file_path) {
  # Read the GCT file
  lines <- readLines(file_path)

  # Extract annotation rows
  num_annotation_rows <- as.numeric(strsplit(lines[2], "\\s+")[[1]][4])
  annotation_rows <- lines[3:(3 + num_annotation_rows)]
  annotations <- setNames(
    sapply(seq_along(annotation_rows), function(i) strsplit(annotation_rows[i], "\t")[[1]][1]),
    seq_along(annotation_rows)
  )
  return(annotations)
}

convert_gct <- function(file_path, annot_row) {
  # convert gct to dataframe
  lines <- readLines(file_path)
  num_data_rows <- as.numeric(strsplit(lines[2], "\\s+")[[1]][1])
  id_row_index <- which(grepl("^id", lines))
  header <- c(lines[id_row_index])
  header <- gsub(" ", "_", header)
  annotation_row <- c(lines[annot_row])
  annotation_row <- gsub(" ", "_", annotation_row)
  data_rows <- lines[(length(lines) - num_data_rows + 1):length(lines)]
  final_content <- c(header, annotation_row, data_rows)
  table <- read.table(text = paste(final_content, collapse = "\n"), header = TRUE)
  #return(table)
  return(list(
    table = table,
    header = header,
    annotation_row = annotation_row
  ))
}

normalize <- function(table, header, annotation_row) {
  # normalize gct dataframe using vst
  count_data <- table[-1,]
  row.names(count_data) <- count_data$id
  count_data <- count_data[,-1]
  row_ids <- row.names(count_data)
  count_data <- as.data.frame(lapply(count_data, function(x) as.numeric(as.character(x))))
  row.names(count_data) <- row_ids
  header <- unlist(strsplit(header, "\\t"))
  annotation_row <- unlist(strsplit(annotation_row, "\\t"))
  metadata <- data.frame(sample = header[-1], condition = annotation_row[-1])
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = metadata,
                                design = ~ condition)
  vsd <- vst(dds, blind = TRUE)
  vst_matrix <- assay(vsd)
  norm_table <- vst_matrix
  row.names(norm_table) <- row.names(count_data)
  return(list(
    table = norm_table,
    matrix = vst_matrix,
    sample_info = metadata
  ))
}

do_pca <- function(norm_table, metadata) {
  pca <- prcomp(t(norm_table), center = TRUE, scale. = FALSE)
  ## variance explained plot
  variance_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)
  fig_var <- plot_ly(
    x = seq_along(variance_explained),
    y = variance_explained,
    type = 'scatter',
    mode = 'markers+lines',
    marker = list(size = 10, color = 'orange', line = list(color = 'black', width = 1)),
    name = 'Variance Explained'
  ) %>%
    layout(
      title = "Explained Variance Ratio",
      xaxis = list(title = "Principal Components", tickvals = seq_along(variance_explained)),
      yaxis = list(title = "Explained Variance Ratio (%)"),
      hovermode = "x unified"
    )
  var_df <- data.frame(
    PC = seq_along(variance_explained),
    Variance = variance_explained
  )
  # PLOT PCA
  ## axis labels
  pc1_label <- paste0("PC1 (", round(variance_explained[1], 2), "%)")
  pc2_label <- paste0("PC2 (", round(variance_explained[2], 2), "%)")
  pc3_label <- paste0("PC3 (", round(variance_explained[3], 2), "%)")
  pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], metadata)
  ## 2D plot
  fig_2D <- plot_ly(
    data = pca_df,
    x = ~PC1, 
    y = ~PC1,
    color = ~condition,
    text = row.names(pca_df),
    marker = list(size = 7, line = list(color = 'black', width = 1)),
    type = 'scatter',
    mode = 'markers',
    textposition = 'top center'
  ) %>%
    layout(
      title = "Principal Component Analysis (PCA)",
      xaxis = list(title = pc1_label),
      yaxis = list(title = pc2_label)
    )
  ## 3D plot
  fig_3D <- plot_ly(
    data = pca_df,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~condition,
    text = row.names(pca_df),
    marker = list(size = 8, line = list(color = 'black', width = 1)),
    type = "scatter3d",
    mode = "markers"
  )
  fig_3D <- fig_3D %>%
    layout(
      title = "Principal Component Analysis (PCA)",
      scene = list(
        xaxis = list(title = pc1_label),
        yaxis = list(title = pc2_label),
        zaxis = list(title = pc3_label)
      )
    )
  return(list(
    pca = pca_df,
    var_exp = var_df,
    var_2D = fig_var,
    pca_2D = fig_2D,
    pca_3D = fig_3D
  ))
}

# Shiny app UI
ui <- fluidPage(
  titlePanel("Shiny PCA"),
  sidebarLayout(
    sidebarPanel(
      fileInput("gct_file", "Upload GCT File", accept = ".gct"),
      uiOutput("annotation_selector"),
      actionButton("submit", "Submit")
      #verbatimTextOutput("selected_key")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Instructions", 
          p("1. Upload a GCT file."),
          p("2. Select an annotation from the dropdown menu."),
          p("3. Click 'Submit' to display the converted table.")
        ),
        tabPanel("Expression Table",
          downloadButton("save_table", "Save as CSV"),
          reactableOutput("gct_table")),
        tabPanel("Normalized Table",
          p("variance-stabilizing transformation (VST) transformation is applied to gene expression count data to stabilize the variance across the range of expression levels."),
          downloadButton("save_norm_table", "Save as CSV"),
          reactableOutput("norm_table")),
        tabPanel("PCA Sceeplot", 
          p("A scree plot is a line graph that shows the eigenvalues of principal components (PCs)"),
          br(),
          plotlyOutput("var_img"),
          reactableOutput("var_table"),
          downloadButton("save_var_table", "Save as CSV")),
        tabPanel("PCA plot",
          p("Principal Component Analysis is a statistical technique introduced by mathematician Karl Pearson in 1901. It works by transforming high-dimensional data into a lower-dimensional space while maximizing the variance (or spread) of the data in the new space. This helps preserve the most important patterns and relationships in the data."),
          br(),
          br(),
          plotlyOutput("PCA_2D"),
          br(),
          p("PC1 and PC2 are the first and second principal components that capture the most and second most variation in a dataset, followed by third and so on."),
          br(),
          plotlyOutput("PCA_3D"),
          br(),
          reactableOutput("pca_scores"),
          downloadButton("save_pca_table", "Save as CSV"))
      )
    )
  )
)

# Shiny app server
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 500 * 1024^2)
  # Reactive value to store annotations
  annotations <- reactiveVal(NULL)
  observeEvent(input$gct_file, {
    req(input$gct_file)
    annotations(extract_annotations(input$gct_file$datapath))
  })

  # Generate dropdown menu for annotations
  output$annotation_selector <- renderUI({
    req(annotations())
    annotation_choices <- setNames(as.list(names(annotations())), annotations())
    selectInput("selected_annotation", "Select Annotation:", choices = annotation_choices, selected = NULL)
  })

  # Display selected annotation key
  output$selected_key <- renderText({
    req(input$selected_annotation)
  
    # Display the selected annotation key (which is already the index)
    paste("Selected Annotation Key:", input$selected_annotation)
  })
  observeEvent(input$submit, {
    req(input$gct_file, input$selected_annotation)
    
    # Get annotation index from the selection
    selected_index <- as.numeric(input$selected_annotation) + 2
    
    # Convert GCT and display as table
    gct_convert <- convert_gct(input$gct_file$datapath, selected_index)
    gct_table <- gct_convert$table
    output$gct_table <- renderReactable({
      reactable(gct_table, searchable = TRUE, filterable = TRUE, highlight = TRUE)
    })
    # Normalize and display
    norm_out <- normalize(gct_table, gct_convert$header, gct_convert$annotation_row)
    norm_table <- norm_out$table
    output$norm_table <- renderReactable({
      reactable(norm_table, searchable = TRUE, filterable = TRUE, highlight = TRUE)
    })
    # PCA and display
    pca_result <- do_pca(norm_table, norm_out$sample_info)
    output$var_img <- renderPlotly({
      pca_result$var_2D
    })
    var_table <- pca_result$var_exp
    output$var_table <- renderReactable({
      reactable(var_table, searchable = TRUE, filterable = TRUE, highlight = TRUE)
    })
    output$PCA_2D <- renderPlotly({
      pca_result$pca_2D
    })
    output$PCA_3D <- renderPlotly({
      pca_result$pca_3D
    })
    pca_scores <- pca_result$pca
    output$pca_scores <- renderReactable({
      reactable(pca_scores, searchable = TRUE, filterable = TRUE, highlight = TRUE)
    })
    # save dataframes
    output$save_table <- downloadHandler(
      filename = function() {
        paste0("table_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(gct_table, file, row.names = FALSE)
      }
    )
    output$save_norm_table <- downloadHandler(
      filename = function() {
        paste0("norm_table_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(norm_table, file, row.names = FALSE)
      }
    )
    output$save_var_table <- downloadHandler(
      filename = function() {
        paste0("pc_variance_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(var_table, file, row.names = FALSE)
      }
    )
    output$save_pca_table <- downloadHandler(
      filename = function() {
        paste0("pca_scores_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(pca_scores, file, row.names = FALSE)
      }
    )
  })
}

# Run the Shiny app
shinyApp(ui, server)