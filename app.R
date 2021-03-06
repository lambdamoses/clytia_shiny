library(shiny)
library(ggplot2)
library(loomR)
library(plotly)
library(shinythemes)
library(scales)
library(velocyto.R)
# For async programming
#library(future)
#library(promises)
# For loader for slow processes
library(shinycssloaders)

# To do:
# Use plot caching (need to wait for the next shiny release)
# Use async programming
# Display table showing marker genes for each cluster
# It would be cool if I can let the users click on a gene to plot it

cell_attrs <- readRDS("clytia_cell_attrs.Rds")
gene_names <- readRDS("clytia_gene_names.Rds")
clytia_loom <- connect("clytia_scaled.loom")
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
         uiOutput("cont_params"),
         radioButtons("theme", "Plot theme", choices = c("light", "dark")),
         helpText("Plotting RNA velocity may take a while (about 1 minute)"),
         checkboxInput("velo", "Plot RNA velocity"),
         uiOutput("velo_params"),
         checkboxInput("interactive", "Interactive plot", value = FALSE),
         actionButton("submit", "Make plot"),
         conditionalPanel("input.velo == false && input.interactive == false",
                          checkboxInput("save_plot", "Save plot?",
                                        value = FALSE)),
         uiOutput("save_ui", inline = TRUE)
      ),
      
      # Show the plot
      mainPanel(
        uiOutput("Plot"),
        uiOutput("Plot2") # Put the violin plot in a separate place
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
  # Get other plot settings
  color_by <- eventReactive(input$submit, input$color_by)
  pt_size <- eventReactive(input$submit, input$pt_size)
  alpha <- eventReactive(input$submit, input$alpha)
  gene <- eventReactive(input$submit, input$gene)
  n_bins <- eventReactive(input$submit, input$n_bins)
  theme_plt <- eventReactive(input$submit, input$theme)
  dim_reduction <- eventReactive(input$submit, input$dim_reduction)
  dim_use_x <- eventReactive(input$submit, input$dim_use_x)
  dim_use_y <- eventReactive(input$submit, input$dim_use_y)
  n_grid <- eventReactive(input$submit, input$n_grid)
  arrow_scale <- eventReactive(input$submit, input$arrow_scale)
  
  # Load data required for velocyto plot only when velo is TRUE
  observeEvent(input$submit, {
    if (if_velo() && !exists("show1") && !exists("velo") && !if_interactive() &&
        color_by() != "cell_density") {
      showNotification("Loading RNA velocity results")
      show1 <<- readRDS("clytia_show.Rds")
      velo <<- readRDS("clytia_velocity.Rds")
    }
  })
  
  # The data frame used for plotting
  df <- eventReactive(input$submit, {
    max_dims <- if (input$dim_reduction == "PCA") 75 else 2
    validate(need(input$dim_use_x <= max_dims && input$dim_use_y <= max_dims,
                  paste(input$dim_reduction, "was only done up to dimension", max_dims)))
    names_get <- c(paste0(input$dim_reduction, c(input$dim_use_x, input$dim_use_y)))
    if (input$velo) {
      df <- as.matrix(cell_attrs[,names_get])
      rownames(df) <- cell_attrs$cell_names
    } else {
      if (!input$color_by %in% c("none", "cell_density", "gene")) {
        names_get <- c(names_get, input$color_by)
      }
      df <- cell_attrs[, c(names_get, "barcode", "cluster")]
      if (input$color_by == "gene") {
        validate(need(input$gene %in% gene_names,
                      paste("Gene", input$gene, "is absent from this dataset.")))
        ind <- which(gene_names == input$gene)
        gene_vals <- clytia_loom[["matrix"]][,ind]
        # Truncate at 10
        gene_vals[gene_vals > 10] <- 10
        df[[input$gene]] <- gene_vals
        df <- df[,c(1,2,5,3,4)]
      }
    }
    df
  })
  
  output$Plot <- renderUI({ 
    if (if_interactive()) {
        withSpinner(plotlyOutput("Plotly_out", height = input$dimension[1] * 3/8))
    } else {
      withSpinner(plotOutput("Plot_out", height = input$dimension[1] * 3/8))
    }
  })
  
  output$Plot2 <- renderUI({
    validate(need(!if_velo() && !color_by() %in% c("none", "cell_density"),
                  ""))
    if (color_by() %in% c("none", "cell_density") || if_velo()) {
      return()
    } else {
      if (if_interactive()) {
        withSpinner(plotlyOutput("Plotly_out2", height = input$dimension[1] * 3/8))
      } else {
        withSpinner(plotOutput("Plot_out2", height = input$dimension[1] * 3/8))
      }
    }
  })
  output$save_ui <- renderUI({
    if (!input$velo && input$save_plot && !input$interactive) {
      if (!color_by() %in% c("none", "cell_density")) {
        wellPanel(radioButtons("which_save", "Which plot to save?",
                            choiceNames = c("Top plot", "Bottom plot"),
                            choiceValues = c("p1", "p2")),
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
               downloadButton("Save", label = "Save")
        )
      } else {
        wellPanel(selectInput("fig_unit", "Plot size unit",
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
      }
    } else {
      return()
    }
  })
  
  output$text_out <- renderText("Interactive mode is not available for RNA velocity plot yet.
                                It will be available for the next release of velocyto.R.")
  output$Plotly_out <- renderPlotly({
    validate(need(!if_velo(), 
                  "Interactive mode is not available for RNA velocity plots at present.
                  It should be available with the next release of velocyto.R."))
    # Using just plotly is faster than ggplotly
    if (color_by() == "cell_density") {
      nms <- names(df())
      p <- plot_ly(x = df()[,1], y = df()[,2]) %>% 
        add_histogram2d(nbinsx = n_bins(), nbinsy = n_bins()) %>% 
        layout(xaxis = list(title = nms[1]), yaxis = list(title = nms[2]))
    } else {
      nms <- names(df())
      p <- plot_ly(x = df()[,1], y = df()[,2], hoverinfo = "text",
                   text = paste0(nms[1], ": ", round(df()[,1], 2), "\n",
                                 nms[2], ": ", round(df()[,2],2), "\n",
                                 "cluster: ", df()$cluster, "\n",
                                 "barcode: ", df()$barcode)) %>% 
        layout(xaxis = list(title = nms[1], zeroline = FALSE), 
               yaxis = list(title = nms[2], zeroline = FALSE))
      if (color_by() == "none") {
        p <- p %>% 
          add_markers(type = "scattergl", opacity = alpha(), color = I("black"),
                      marker = list(size = pt_size() * 3, sizemode = "diameter"))
      } else if (color_by() == "cluster") {
        p <- p %>% 
          add_markers(type = "scattergl", opacity = alpha(), 
                      color = df()[,3], colors = hue_pal()(20),
                      marker = list(size = pt_size() * 3, sizemode = "diameter"))
      } else {
        p <- p %>% 
          add_markers(type = "scattergl", opacity = alpha(), color = df()[,3],
                      marker = list(size = pt_size() * 3, sizemode = "diameter"))
      }
    }
    if (theme_plt() == "dark") {
      p <- p %>% 
        layout(plot_bgcolor = "#7F7F7F")
    }
    p
  })
  
  output$Plotly_out2 <- renderPlotly({
    validate(need(!if_velo(), 
                  "Interactive mode is not available for RNA velocity plots at present.
                  It should be available with the next release of velocyto.R."))
    if (color_by() == "cluster") {
      counts <- tapply(rep(1,12165), df()$cluster, length)
      p2 <- plot_ly(x = 0:19, y = counts,
                    color = as.factor(0:19), colors = hue_pal()(20),
                    type = "bar", hoverinfo = "text", text = paste(counts)) %>% 
        layout(xaxis = list(title = "cluster"),
               yaxis = list(title = "count"),
               showlegend = FALSE)
    } else {
      nms <- names(df())
      p2 <- plot_ly(x = df()$cluster, y = df()[,3], color = df()$cluster,
                    colors = hue_pal()(20), type = "violin",
                    box = list(visible = FALSE)) %>% 
        layout(xaxis = list(title = "cluster"), 
               yaxis = list(title = nms[3]))
    }
    p2
  })
  
  output$Plot_out <- renderPlot({
    if (if_velo()) {
      validate(need(color_by() != "cell_density",
                    "Can't color by cell density in RNA-velocity plot."))
      # I really look forward to the next release of velocyto.R, when ggplot2 will be used
      if (color_by() == "cluster") {
        col_mode <- "discrete"
      } else {
        col_mode <- "continuous"
      }
      if (color_by() %in% c("cell_density", "none")) {
        colors_use <- NULL
      } else if (color_by() != "gene") {
        col_vec <- cell_attrs[,color_by()]
        colors_use <- velo_set_colors(col_mode, col_vec, alpha())
      } else {
        ind <- which(gene_names == gene())
        # Come back here if truncation is the culprit
        col_vec <- clytia_loom[["matrix"]][,ind]
        colors_use <- velo_set_colors(col_mode, col_vec, alpha())
      }
      nms <- colnames(df())
      showNotification("Computing arrows", duration = 45)
      show.velocity.on.embedding.cor(emb = df(), 
                                     vel = velo, show.grid.flow = TRUE, 
                                     arrow.scale = arrow_scale(), grid.n = n_grid(), cc = show1$cc,
                                     cell.colors = colors_use,
                                     cex = pt_size(), xlab = nms[1], ylab = nms[2])
    } else {
      dr_names <- paste0(dim_reduction(), c(dim_use_x(), dim_use_y()))
      p <- ggplot(df(), aes_string(dr_names[1], dr_names[2], label = "barcode"))
      if (color_by() == "none") {
        p <- p +
          geom_point(size = pt_size(), alpha = alpha())
      } else if (color_by() == "cell_density") {
        p <- p +
          geom_hex(bins = n_bins()) +
          scale_fill_viridis_c()
      } else if (color_by() == "gene") {
        p <- p +
          geom_point(aes_string(color = gene()), size = pt_size(), alpha = alpha()) +
          scale_color_viridis_c()
      } else {
        p <- p +
          geom_point(aes_string(color = color_by()), 
                     size = pt_size(), alpha = alpha())
        if (color_by() != "cluster") {
          p <- p +
            scale_color_viridis_c()
        }
      }
      if (theme_plt() == "dark") {
        p <- p + theme_dark()
      }
      p1 <<- p
      p
    }
  })
  
  output$Plot_out2 <- renderPlot({
    validate(need(!color_by() %in% c("none", "cell_density") && !if_velo(), ""))
    # Add violin plot or barplot
    if (color_by() == "cluster") {
      p2 <<- ggplot(df(), aes(cluster, fill = cluster)) +
        geom_bar() +
        theme(legend.position = "none")
    } else if (color_by() == "gene") {
      p2 <<- ggplot(df(), aes_string("cluster", gene(), fill = "cluster")) +
        geom_violin() +
        theme(legend.position = "none")
    } else {
      p2 <<- ggplot(df(), aes_string("cluster", color_by(), fill = "cluster")) +
        geom_violin() +
        theme(legend.position = "none")
    }
    p2
  })
  # Default file name
  output$Save <- downloadHandler(filename = "plot", 
                                 content = function (file) {
                                   pl <- if (color_by() %in% c("none", "cell_density")) {
                                     p1
                                   } else if (input$which_save == "p2") {
                                     p2
                                   } else {
                                     p1
                                   }
                                   ggsave(file, pl,
                                          width = input$fig_width,
                                          height = input$fig_height,
                                          units = input$fig_unit,
                                          device = input$plot_format)
                                 })
}

# Run the application 
shinyApp(ui = ui, server = server)
