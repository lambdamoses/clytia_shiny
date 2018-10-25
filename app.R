library(shiny)
library(ggplot2)
library(tibble)
library(dplyr)
library(loomR)
library(plotly)
library(shinythemes)
library(scales)
library(velocyto.R)
# For async programming
library(future)
library(promises)
# For loader for slow processes
library(shinycssloaders)

# To do:
# Use plot caching
# Use async programming
# Display table showing marker genes for each cluster
# It would be cool if I can let the users click on a gene to plot it

cell_attrs <- readRDS("clytia_cell_attrs.Rds")
gene_names <- readRDS("clytia_gene_names.Rds")
clytia_loom <- connect("clytia.loom")
theme_set(theme_bw())
# Set cell colors for velocyto plot
velo_set_colors <- function(mode = "discrete", vec, alpha) {
  switch(mode,
         "discrete" = {
           setNames(ac(cell_attrs$cluster_colors, alpha = alpha), 
                    cell_attrs$cell_names)
         },
         "continuous" = {
           color_vec <- viridis_pal()(256)[as.numeric(cut(vec, 256))]
           setNames(ac(color_vec, alpha = alpha), cell_attrs$cell_names)
         })
}

# Define UI for plot settings
ui <- fluidPage(
   theme = shinytheme("flatly"),
   # Application title
   titlePanel("Clytia data explorer"),
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        tags$head(tags$script('
                                var dimension = [0, 0];
                              $(document).on("shiny:connected", function(e) {
                              dimension[0] = window.innerWidth;
                              dimension[1] = window.innerHeight;
                              Shiny.onInputChange("dimension", dimension);
                              });
                              $(window).resize(function(e) {
                              dimension[0] = window.innerWidth;
                              dimension[1] = window.innerHeight;
                              Shiny.onInputChange("dimension", dimension);
                              });
                              ')),
         radioButtons("dim_reduction", "Dimension reduction",
                      choices = c("PCA", "tSNE", "UMAP")),
         fluidRow(
           column(6, numericInput("dim_use_x", "Dimension for x", 1, min = 1, step = 1)),
           column(6, numericInput("dim_use_y", "Dimension for y", 2, min = 1, step = 1))),
         selectInput("color_by", "Color by", 
                     choices = c("none", "cluster", "gene", "nGene", "nUMI", "cell_density"),
                     selected = "none"),
         # options depend on what to use to color
         uiOutput("cont_params", inline = TRUE),
         radioButtons("theme", "Plot theme", choices = c("light", "dark")),
         helpText("Plotting RNA velocity may take a while (about 1 minute)"),
         checkboxInput("velo", "Plot RNA velocity"),
         uiOutput("velo_params"),
         checkboxInput("interactive", "Interactive plot", value = FALSE),
         actionButton("submit", "Make plot"),
         helpText("Save plot only works for static plots"),
         checkboxInput("save_plot", "Save plot?",
                       value = FALSE),
         conditionalPanel(condition = "input.save_plot == true && input.interactive == false",
                          selectInput("fig_unit", "Plot size unit",
                                      choices = c("in", "cm", "mm"),
                                      label = c("inch", "centimeter", 
                                                "millimeter")),
                          fluidRow(column(6, numericInput("fig_width", 
                                                          "Plot width", 
                                                          value = 6)),
                                   column(6, numericInput("fig_height",
                                                          "Plot height",
                                                          value = 4))),
                          radioButtons("plot_format", "Format",
                                       choices = c("jpeg", "png", "pdf", "tiff")),
                          downloadButton("Save", label = "Save"))
      ),
      
      # Show the plot
      mainPanel(
        uiOutput("Plot")
      )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$cont_params <- renderUI({
    if (input$color_by == "cell_density") {
      sliderInput("n_bins", "Number of bins", 30, 200, 150, step = 5)
    } else if (input$color_by == "gene") {
      column(12,
        textInput("gene", "Gene symbol"),
        sliderInput("pt_size", "Point size", min = 0, max = 2, value = 1, step = 0.1),
        sliderInput("alpha", "Alpha (opacity)", 0, 1, 1, step = 0.05))
    } else {
      column(12,
        sliderInput("pt_size", "Point size", min = 0, max = 2, value = 1, step = 0.1),
        sliderInput("alpha", "Alpha (opacity)", 0, 1, 1, step = 0.05))
    }
  })
  
  # Settings for velocyto plot
  output$velo_params <- renderUI({
    if (input$velo) {
      column(12,
             sliderInput("n_grid", "Number of grid points along each axis",
                         10, 100, value = 40, step = 1),
             sliderInput("arrow_scale", "Arrow scale",1, 5, value = 3, step = 0.5))
    } else {
      return()
    }
    
  })
  # Get these plot settings: interactive, velocyto
  if_interactive <- eventReactive(input$submit, input$interactive)
  if_velo <- eventReactive(input$submit, input$velo)
  # Load data required for velocyto plot only when velo is TRUE
  observeEvent(input$submit, {
    if (if_velo() && !exists("show1") && !exists("velo") && !if_interactive()) {
      showNotification("Loading RNA velocity results")
      show1 <<- readRDS("clytia_show.Rds")
      velo <<- readRDS("clytia_velocity.Rds")
    }
  })
  
  # The data frame used for plotting
  df <- eventReactive(input$submit, {
    names_get <- c(paste0(input$dim_reduction, c(input$dim_use_x, input$dim_use_y)))
    if (input$velo) {
      df <- as.matrix(cell_attrs[,names_get])
      rownames(df) <- cell_attrs$cell_names
    } else {
      if (!input$color_by %in% c("none", "cell_density", "gene")) {
        names_get <- c(names_get, input$color_by)
      }
      df <- cell_attrs[, c(names_get, "barcode")]
      if (input$color_by == "gene") {
        ind <- which(gene_names == input$gene)
        gene_vals <- clytia_loom[["layers/scale_data"]][,ind]
        # Truncate at 10
        gene_vals[gene_vals > 10] <- 10
        df[[input$gene]] <- gene_vals
      }
    }
    df
  })
  # Generate the basic plot
  g <- eventReactive(input$submit, {
    if (input$velo) {
      # I really look forward to the next release of velocyto.R, when ggplot2 will be used
      if (input$color_by == "cluster") {
        col_mode <- "discrete"
      } else {
        col_mode <- "continuous"
      }
      if (input$color_by %in% c("cell_density", "none")) {
        if (input$color_by == "cell_density") {
          showNotification("Can't color by cell density in RNA velocity plot",
                           duration = NULL)
        }
        colors_use <- NULL
      } else if (input$color_by != "gene") {
        col_vec <- cell_attrs[,input$color_by]
        colors_use <- velo_set_colors(col_mode, col_vec, input$alpha)
      } else {
        ind <- which(gene_names == input$gene)
        # Come back here if truncation is the culprit
        col_vec <- clytia_loom[["layers/scale_data"]][,ind]
        colors_use <- velo_set_colors(col_mode, col_vec, input$alpha)
      }
      showNotification("Computing arrows", duration = NULL)
      show.velocity.on.embedding.cor(emb = df(), 
                                     vel = velo, show.grid.flow = TRUE, 
                                     arrow.scale = 3, grid.n = 40, cc = show1$cc,
                                     cell.colors = colors_use,
                                     cex = input$pt_size)
    } else {
      dr_names <- paste0(input$dim_reduction, c(input$dim_use_x, input$dim_use_y))
      p <- ggplot(df(), aes_string(dr_names[1], dr_names[2], label = "barcode"))
      if (input$color_by == "none") {
        p <- p +
          geom_point(size = input$pt_size, alpha = input$alpha)
      } else if (input$color_by == "cell_density") {
        p <- p +
          geom_hex(bins = input$n_bins) +
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
      if (input$theme == "dark") {
        p <- p + theme_dark()
      }
      p
    }
  })
  output$Plot <- renderUI({
    if (if_interactive()) {
      if (if_velo()) {
        textOutput("text_out")
      } else {
        withSpinner(plotlyOutput("Plotly_out", height = input$dimension[1] / 2))
      }
    } else {
      withSpinner(plotOutput("Plot_out", height = input$dimension[1] / 2))
    }
  })
  output$text_out <- renderText("Interactive mode is not available for RNA velocity plot yet.
                                It will be available for the next release of velocyto.R.")
  output$Plotly_out <- renderPlotly(toWebGL(ggplotly(g())))
  output$Plot_out <- renderPlot(g())
  # Default file name
  output$Save <- downloadHandler(filename = "plot", 
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
