libraries <- c("shiny", "scDHA", "torch")
install.packages(libraries, dependencies = TRUE)

library(shiny)
library(scDHA)
library(torch)

ui <- fluidPage(
  tags$style(HTML("
    .section-container {
      margin: 10px; /* Add margin to each section */
    }
    .button-container {
      margin-top: 20px;
      margin-bottom: 20px;
    }
    .plot-container {
      margin-top: 20px;
      margin-bottom: 20px;
    }
    .table-container {
      margin-top: 20px;
      margin-bottom: 20px;
      max-height: 300px; /* Set max height for scrollable tables */
      overflow-y: auto; /* Enable vertical scrolling for tables */
    }
    .output-image {
      margin: 20px; /* Add margin to the image */
    }
    .grid-container {
      display: grid;
      grid-template-columns: repeat(2, 1fr); /* Update the number of columns to 2 */
      grid-template-rows: repeat(2, 1fr);
      gap: 20px;
    }
    .grid-item {
      border: 1px solid #ccc;
      padding: 10px;
    }
  ")),

  fluidRow(
    #Data Upload Section
    column(width = 6,
           div(class = "section-container",
               div(class = "grid-item",
                   fileInput(inputId = 'data',
                             label = 'Please upload your data',
                             accept = '.csv'),
                   helpText("The input is a matrix with rows as samples and columns as genes or transcripts"),
                   numericInput(inputId = 'genenumbers',
                                label = 'Choose the number of desired genes',
                                value = 5000),
                   actionButton(inputId = 'runButton', label = 'Run Analysis'),
                   div(class = "button-container",
                       downloadButton("downloadData", "Log2 Tranformed Data"),
                       downloadButton("downloadCluster", "Clustering Result"),
                       downloadButton("downloadLatent", "Bottleneck"),
                       downloadButton("downloadPlotDataPNG", "scDHA Plot")
                   )
               ),
               div(class = "grid-item",
                   h4("About"),
                   p("This shiny app preforms single-cell RNA-seq analysis using the scDHA package.",
                     "For more information please visit : ",
                     a(href = "https://github.com/duct317/scDHA", "https://github.com/duct317/scDHA")
                   )
               )
           )
    ),

    #Plot section
    column(width = 6,
           div(class = "section-container",
               div(class = "grid-item",
                   div(class = "plot-container",
                       plotOutput("scatterPlot"), class = "output-image")
               )
           )
    )
  ),

  fluidRow(
    # Clustering result section
    column(width = 6,
           div(class = "section-container",
               div(class = "grid-item",
                   div(class = "table-container",
                       tableOutput("clusterTable")
                   )
               )
           )
    ),

    #Bottleneck section
    column(width = 6,
           div(class = "section-container",
               div(class = "grid-item",
                   div(class = "table-container",
                       tableOutput("latentTable")
                   )
               )
           )
    )
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize = 10 * 1024^2)

  filtered_data <- reactiveVal(NULL)
  result <- reactiveVal(NULL)
  cluster_data <- reactiveVal(NULL)
  latent_data <- reactiveVal(NULL)

  observeEvent(input$runButton, {
    req(input$data)
    data <- read.csv(input$data$datapath, header = TRUE, row.names = 1)
    data_matrix <- as.matrix(data)
    data <- log2(data_matrix + 1)
    filtered_data(data)
    result(scDHA(data, seed = 1))

    cluster_data_df <- data.frame(cluster = result()$cluster)
    rows <- rownames(data)
    rownames(cluster_data_df) <- rows
    cluster_data(cluster_data_df)

    latent_data_df <- as.data.frame(result()$latent)
    rownames(latent_data_df) <- rows
    latent_data(latent_data_df)
  })

  output$scatterPlot <- renderPlot({
    if (!is.null(result())) {
      scDHA_result <- scDHA.vis(result(), seed = 1)
      plot(scDHA_result$pred, col = factor(label), xlab = "scDHA1", ylab = "scDHA2")
    }
  })

  output$clusterTable <- renderTable({
    if (!is.null(cluster_data())) {
      cluster_data()
    }
  }, rownames = TRUE)

  output$latentTable <- renderTable({
    if (!is.null(latent_data())) {
      latent_data()
    }
  }, rownames = TRUE)

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("log_2_transformed_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      filtered_data_df <- as.data.frame(filtered_data())
      write.csv(filtered_data_df, file, row.names = TRUE)
    }
  )

  output$downloadCluster <- downloadHandler(
    filename = function() {
      paste("cluster_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(cluster_data())) {
        write.csv(cluster_data(), file, row.names = TRUE)
      }
    }
  )

  output$downloadLatent <- downloadHandler(
    filename = function() {
      paste("scDHA_latent_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(latent_data())) {
        write.csv(latent_data(), file, row.names = TRUE)
      }
    }
  )

  output$downloadPlotDataPNG <- downloadHandler(
    filename = function() {
      paste("scatter_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      if (!is.null(result())) {
        scDHA_result <- scDHA.vis(result(), seed = 1)
        png(file, width = 800, height = 600)
        plot(scDHA_result$pred, col = factor(label), xlab = "scDHA1", ylab = "scDHA2")
        dev.off()
      }
    }
  )
}

shinyApp(ui, server)
