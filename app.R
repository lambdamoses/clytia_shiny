library(shiny)
library(tidyverse)
library(loomR)
library(plotly)
library(shinythemes)

cell_attrs <- readRDS("clytia_cell_attrs.Rds")
gene_names <- readRDS("clytia_gene_names.Rds")
clytia_loom <- connect("clytia.loom")
#show1 <- readRDS("clytia_show.Rds")
#velo <- readRDS("clytia_velocity.Rds")

# Default file name
default_fn <- function(ext) {
  paste0("plot-", Sys.time(), ".", ext)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   theme = shinytheme("flatly"),
   # Application title
   titlePanel("Clytia data explorer"),
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         radioButtons("dim_reduction", "Dimension reduction",
                      choices = c("PCA", "tSNE", "UMAP")),
         fluidRow(
           column(6, numericInput("dim_use_x", "Dimension for x", 1, min = 1, step = 1)),
           column(6, numericInput("dim_use_y", "Dimension for y", 2, min = 1, step = 1))),
         selectInput("color_by", "Color by", 
                     choices = c("none", "cluster", "gene", "nGene", "nUMI", "cell_density"),
                     selected = "none"),
         conditionalPanel(condition = "input.color_by == 'gene'",
                          selectInput("gene", "Gene symbol", choices = gene_names)),
         sliderInput("pt_size", "Point size", min = 0, max = 2, value = 1, step = 0.1),
         sliderInput("alpha", "Alpha (opacity)", 0, 1, 1, step = 0.05),
         radioButtons("theme", "Plot theme", choices = c("light", "dark")),
         checkboxInput("velo", "Plot RNA velocity"),
         checkboxInput("interactive", "Interactive plot", value = FALSE),
         actionButton("submit", "Make plot"),
         checkboxInput("save_plot", "Save plot? (only applies for non-interactive plots)",
                       value = FALSE),
         conditionalPanel(condition = "input.save_plot == true && input.interactive == false",
                          selectInput("fig_unit", "Plot size unit",
                                      choices = c("in", "cm", "mm")),
                          fluidRow(column(6, numericInput("fig_width", 
                                                          "Plot width", 
                                                          value = 6)),
                                   column(6, numericInput("fig_height",
                                                          "Plot height",
                                                          value = 4))),
                          radioButtons("plot_format", "Format",
                                       choices = c("jpeg", "png", "pdf")),
                          downloadButton("Save", label = "Save"))
      ),
      
      # Show the plot
      mainPanel(
         plotOutput("Plot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # The data frame used for plotting
  df <- eventReactive(input$submit, {
    names_get <- c(paste0(input$dim_reduction, c(input$dim_use_x, input$dim_use_y)))
    if (!input$color_by %in% c("none", "cell_density", "gene")) {
      names_get <- c(names_get, input$color_by)
    }
    df <- cell_attrs[,names_get]
    if (input$color_by == "gene") {
      ind <- which(gene_names == input$gene)
      df[[input$gene]] <- clytia_loom[["matrix"]][,ind]
    }
    df
  })
  # Generate the basic plot
  g <- eventReactive(input$submit, {
    dr_names <- paste0(input$dim_reduction, c(input$dim_use_x, input$dim_use_y))
    p <- ggplot(df(), aes_string(dr_names[1], dr_names[2]))
    if (input$color_by == "none") {
      p <- p +
        geom_point(size = input$pt_size, alpha = input$alpha)
    } else if (input$color_by == "cell_density") {
      p <- p +
        geom_hex(bins = 150) +
        scale_fill_viridis_c()
    } else if (input$color_by == "gene") {
      p <- p +
        geom_point(aes_string(color = input$gene), size = input$pt_size, alpha = input$alpha) +
        scale_color_viridis_c()
    } else {
      p <- p +
        geom_point(aes_string(color = input$color_by), 
                   size = input$pt_size, alpha = input$alpha)
      if (input$color_by != "cluster") {
        p <- p +
          scale_color_viridis_c()
      }
    }
    if (input$theme == "light") {
      p <- p +
        theme_bw()
    } else {
      p <- p + theme_dark()
    }
    p
  })
  output$Plot <- renderPlot(g())
  # Default file name
  fn <- reactive({
    default_fn(input$plot_format)
  })
  output$Save <- downloadHandler(filename = fn(), 
                                 content = function (file) {
                                   ggsave(file, plot = g(),
                                          width = input$fig_width,
                                          height = input$fig_height,
                                          units = input$fig_unit,
                                          device = input$plot_format)
                                 })
}

# Run the application 
shinyApp(ui = ui, server = server)

