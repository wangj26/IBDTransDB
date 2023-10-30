# IBDTransDB Explore Application
# An application used to explore omics data

library(RSQLite)

library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyBS)
library(shinybusy)
library(shinyalert)
#library(shinyjs)

library(DT)
library(kableExtra)

library(dplyr)
library(tidyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)

library(WebGestaltR)


#==============================================================================================================
# Import functions used to access the DB, create PCA/volcano plots, perform signature analysis, etc
#==============================================================================================================
source("./app_functions/dataset_database_query.R")
source("./app_functions/pca.R")
source("./app_functions/DGE_viewer.R")
source("./app_functions/signature_viewer.R")
source("./app_functions/enrichment.R")
source("./app_functions/cell_deconvolution.R")
source("./app_functions/dataset_download.R")
#==============================================================================================================


#==============================================================================================================
# Create the connection to the ImmunoVerse DB
#==============================================================================================================
db <- dbConnect(RSQLite::SQLite(), dbname = "./IBDTransDB.db")
#==============================================================================================================

global_dataset_acc <- NULL

# Define the UI for the application
ui <- navbarPage(
    id = "dataset-navbar",
    
    theme = "style.css",
    tags$head(
      tags$style("
        body {
          overflow-x: hidden; overflow-y: hidden; overflow: visible !important;
        } 
        .targets_button_wrapper:before, .targets_button_wrapper:after {
          content: ' ';
          width: 100px;
          height: 1px;
          margin: 0 10px;
          vertical-align: super;
          background-color: #4fbfa8;
          display: inline-block;
        }
        .radio-inline {
          padding-left: 20px;
          padding-right: 50px;
        }
        .shiny-output-error { 
          visibility: hidden;
        }
        .shiny-output-error:before { 
          visibility: hidden;
        }
        .bttn-unite.bttn-md { 
          width: 450px !important;
        }
        "
      )
    ),

   #useShinyjs(),
    
    title = "IBDExplore",
    
    ### DATASET DESCRIPTION TAB ###        
    tabPanel("Dataset Description",
             id = "dataset-description",
             add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
             column(10, offset = 1,
                    column(10, offset = 1,
                           wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                                     uiOutput("dataset_description"),
                                     br(), hr(), br(),
                                     fluidRow(
                                       column(3, uiOutput("pca_plot_options")),
                                       column(5, plotOutput("pca_plot")),
                                       column(4, 
                                              uiOutput("pc_test_options"),
                                              br(),
                                              uiOutput("pc_test_output"), 
                                              plotOutput("pc_test_plot", height = "225px"))
                                     ),
                                     fluidRow(
                                       column(3, br(), uiOutput("scree_plot_button"))
                                     )
                           ))
             )
    ),
    
    ### DIFFERENTIAL GENE EXPRESSION VIEWER TAB ###
    tabPanel("DGE Viewer",
             id = "dge-viewer",
             add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
             column(8, offset = 2,
                    fluidRow(
                      wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 200px;",
                                uiOutput("volcano_and_comparison_options"), 
                                br(),
                                div(style = "height: 50px;",
                                    selectizeInput(inputId = "volcano_genes",
                                                   label = "Select genes of interest to be displayed in the DE Gene Table Viewer",
                                                   choices = NULL,
                                                   multiple = TRUE,
                                                   width = "500px"))
                      )
                    ),
                    fluidRow(
                      column(6, 
                             wellPanel(style = "background-color: #fff; border-color: #2c3e50; margin-left: -1em; height: 725px;",
                                       h3("DE Gene Table Viewer", style = "text-align: center;"),
                                       DT::dataTableOutput("comparison_data_table"),
                                       br(),
                                       uiOutput("comparison_data_table_buttons"))),
                      column(6, 
                             wellPanel(style = "background-color: #fff; border-color: #2c3e50; margin-right: -1em; height: 725px;",
                                       h3("Volcano Viewer", style = "text-align: center;"),
                                       div(style = "display: flex; justify-content: space-between;",
                                           dropdownButton(h5("volcano plot options"),
                                                          radioGroupButtons(inputId = "volcano_pval_option",
                                                                            label = NULL,
                                                                            choices = c("p-value" = "p_value", "adjusted p-value" = "p_value_adj"),
                                                                            selected = "p_value",
                                                                            width = "300px",
                                                                            checkIcon = list(
                                                                              yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                                              no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))),
                                                          numericInput(inputId = "volcano_p_value_threshold",
                                                                       label = "Threshold for the p-value",
                                                                       value = 0.05,
                                                                       min = 0,
                                                                       max = 1,
                                                                       step = 0.01,
                                                                       width = "200px"),
                                                          numericInput(inputId = "volcano_log_fc_threshold",
                                                                       label = "Threshold for the logFC",
                                                                       value = 1,
                                                                       min = 0,
                                                                       max = 10,
                                                                       step = 0.1,
                                                                       width = "200px"),
                                                          circle = TRUE, 
                                                          status = "primary",
                                                          size = "sm",
                                                          icon = icon("cog"), 
                                                          width = "350px",
                                                          right = FALSE,
                                                          tooltip = tooltipOptions(title = "Click to see plot options", placement = "right")),
                                           downloadButton(outputId = "download_volcano_plot",
                                                          label = "Download Volcano Plot",
                                                          style = "text-align center; width: 200px")),
                                       br(),
                                       plotOutput("volcano_plot"), 
                                       uiOutput("volcano_significant_genes")))
                    ),
                    fluidRow(
                      wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                                h3("Gene Expression Box Plot Viewer", style = "text-align: center;"),
                                p("Select rows in the comparison metrics table to display associated expression plots",
                                  "for the given comparison. A maximum of 9 rows can be selected at a time."),
                                div(style = "display: flex; justify-content: space-between;",
                                    # dropdownButton(inputId = "box_plot_options",
                                    #                h5("plot options"),
                                    #                radioGroupButtons(inputId = "gene_plot_display_pval",
                                    #                                  label = NULL,
                                    #                                  choices = c("p-value" = "p_value", "adjusted p-value" = "p_value_adj"),
                                    #                                  selected = "p_value",
                                    #                                  width = "300px",
                                    #                                  checkIcon = list(
                                    #                                    yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                    #                                    no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))),
                                    #                switchInput(inputId = "gene_plot_display_jitter",
                                    #                            label = "show data points",
                                    #                            value = FALSE,
                                    #                            width = "300px"),
                                    #                circle = TRUE,
                                    #                status = "primary",
                                    #                size = "sm",
                                    #                icon = icon("cog"),
                                    #                width = "350px",
                                    #                right = FALSE,
                                    #                tooltip = tooltipOptions(title = "Click to see plot options",
                                    #                                         placement = "right")),
                                    downloadButton(outputId = "download_expression_boxplots",
                                                   label = "Download Plot(s)",
                                                   style = "text-align center; width: 175px")),
                                br(),
                                uiOutput("gene_plots"))
                    )
             )
    ),
    
    ### SIGNATURE VIEWER TAB ###
    tabPanel("Signature Viewer",
             id = "signature-viewer",
             add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
             column(10, offset = 1,
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 1050px;",
                              column(5, uiOutput("signature_plot_options")),
                              column(7, 
                                     div(style = "display: flex; justify-content: space-between;",
                                         actionButton(inputId = "generate_signature_plot",
                                                      label = "Generate Plot(s)",
                                                      style = "text-align center; width: 150px"),
                                         # dropdownButton(p("plot options"),
                                         #                radioGroupButtons(inputId = "signature_display_pval",
                                         #                                  label = NULL,
                                         #                                  choices = c("p-value" = "p_value", "adjusted p-value" = "p_value_adj"),
                                         #                                  selected = "p_value",
                                         #                                  width = "300px",
                                         #                                  checkIcon = list(
                                         #                                    yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                         #                                    no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))),
                                         #                switchInput(inputId = "signature_display_jitter",
                                         #                            label = "show data points",
                                         #                            value = FALSE,
                                         #                            width = "300px"),
                                         #                circle = TRUE,
                                         #                status = "primary",
                                         #                size = "sm",
                                         #                icon = icon("cog"),
                                         #                width = "350px",
                                         #                right = FALSE,
                                         #                tooltip = tooltipOptions(title = "Click to see plot options", placement = "left"), inputId="signatureDropDown")
                                         ),
                                     br(),
                                     uiOutput("signature_plots"), 
                                     br(),
                                     DT::dataTableOutput("signature_plot_table"))
             ))
    ),
    
    ### ENRICHMENT ANALYSIS TAB ###
    tabPanel("Enrichment Analysis",
             id = "enrichment-analysis",
             add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
             column(8, offset = 2,
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                              h2("Enrichment Analysis",
                                 style = "text-align: center;"),
                              br(),
                              radioGroupButtons(inputId = "enrichment_option",
                                                label = NULL,
                                                choices = c("ORA", "GSEA"),
                                                selected = "ORA",
                                                status = "primary"),
                              br(),
                              uiOutput("enrichment_options")
                    ),
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                              uiOutput("enrichment_report")
                    )
             )
    ),
    
    ### CELL DECONVOLUTION TAB ###
    tabPanel("Cell Deconvolution",
             id = "cell-deconvolution",
             add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
             column(10, offset = 1,
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 300px;",
                              uiOutput("cell_deconvolution_options")),
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                              uiOutput("cell_deconvolution_results"),
                              #DT::dataTableOutput("cell_deconvolution_table"),
                              br(),
                              uiOutput("cell_deconvolution_plot_options"),
                              uiOutput("cell_deconvolution_plots"))
             )
    )
    
)


# Define the server functionality for the application
server <- function(input, output, session) {
  
    ### REACTIVE VALUES ###
    active_comparison_data <- reactiveVal(NULL)
    active_enrichment_data <- reactiveValues(gene_set_collection = NULL,
                                             ranked_genes = NULL,
                                             gsea_results = NULL,
                                             top_gsea_results = NULL)
    values <- reactiveValues()
    ##########
    
    observe({
      # Import the selected dataset ID from the URL parameter
      #dataset_acc <- parseQueryString(session$clientData$url_search)[["dataset_acc"]]
      #global_dataset_acc <<- dataset_acc
      dataset_acc <- "GSE16879"
      global_dataset_acc <<- "GSE16879"
      if (!is.null(dataset_acc)) {
        # Pull the necessary dataset description from the DB
        dataset_description <- get_dataset_description(dataset_acc, db)
        showNotification("Dataset description done loading", type = "default")
        
        dataset_acc_identifier <- dataset_description$id
        
        #==============================================================================================================
        ### DATASET DESCRIPTION
        #==============================================================================================================
        
        # Generate the output for the dataset description (e.g. title, data type)
        output$dataset_description <- renderUI({
          div(h3(paste0(dataset_description$title, " (",  dataset_acc,")")),
              br(),
              splitLayout(cellWidths = c("33%", "33%", "33%"),
                          list(p("Disease:", dataset_description$disease),
                               p("Organism:", dataset_description$organism),
                               p("Data Type:", dataset_description$experiment_type)),
                          list(p("Source:", dataset_description$source),
                               p("Cell Type:", dataset_description$cell_type),
                               p("Sample Number:", dataset_description$sample_number)),
                          list(p("Treatment:", dataset_description$treatment),
                               p("Timepoint:", dataset_description$timepoint),
                               p("Dose:", dataset_description$dose))),
              p("Platform:", dataset_description$platform),
              p("Normalization Method:", dataset_description$normalization_method),
              p("Summary:", dataset_description$summary),
              p("Design:", dataset_description$design))
        })
        #==============================================================================================================
    
        # Pull the mapping of genes identifiers to gene names
        gene_map <- get_gene_map(db)
        
        # Pull the necessary data (exp data, sample annotations, comparisons, etc) from the DB
        sample_ann <- get_sample_ann(dataset_acc, db)
        showNotification("Sample annotations done loading", type = "default")
        comparisons <- get_comparisons(dataset_acc, db)
        showNotification("Comparison descriptions done loading", type = "default")
        comparison_data <- get_comparison_data(comparisons$id[1], gene_map, db)
        showNotification("Comparison data done loading", type = "default")
        # Create a sample to group mapping from the comparison description
        sample_group_map <- map_samples_to_groups(comparisons)
        
        #==============================================================================================================
        ### PCA
        #==============================================================================================================
        
        # Pull the PCA data from the DB
        pca_list <- get_pca(dataset_acc, sample_ann, db)
        pca_df <- pca_list[["pca_df"]]  
        pc_variance <- pca_list[["pc_variance"]]
        
        # Generate the PCA options (i.e. selections for feature, x axis, y axis)
        output$pca_plot_options <- renderUI({
          div(
            selectInput(inputId = "pca_group",
                        label = "The PCA plot will be colored by the selected variable",
                        choices = c("---", unlist(str_split(unique(sample_ann$sample_ann_type), ";"))),
                        multiple = FALSE,
                        selected = "---"),
            selectInput(inputId = "pca_filter_feature",
                        label = "PCA plot filter",
                        choices = c("---", unlist(str_split(unique(sample_ann$sample_ann_type), ";"))),
                        multiple = FALSE,
                        selected = "---"),
            selectInput(inputId = "pca_filter_selections",
                        label = "Filter choices",
                        choices = NULL,
                        multiple = TRUE),
            div(style = "display: flex; justify-content: space-between;",
                selectInput(inputId = "pca_x_axis",
                            label = "x axis",
                            choices = list("PC1" = 1, "PC2" = 2, "PC3" = 3, "PC4" = 4, "PC5" = 5),
                            selected = 1,
                            multiple = FALSE,
                            width = "100px"),
                selectInput(inputId = "pca_y_axis",
                            label = "y axis",
                            choices = list("PC1" = 1, "PC2" = 2, "PC3" = 3, "PC4" = 4, "PC5" = 5),
                            selected = 2,
                            multiple = FALSE,
                            width = "100px")),
            div(style = "display: flex; justify-content: space-between;",
                prettyCheckbox(inputId = "plot_ellipse",
                               label = "Circle groups",
                               value = FALSE,
                               status = "primary",
                               shape = "curve",
                               inline = TRUE),
                prettyCheckbox(inputId = "hide_legend",
                               label = "Hide legend",
                               value = FALSE,
                               status = "primary",
                               shape = "curve",
                               inline = TRUE))
          )
        })
        
        # Generate the initial PCA plot (not colored by any group)
        output$pca_plot <- renderPlot({
          generate_pca_plot(pca_df, pc_variance,
                            input$pca_group, 
                            input$pca_filter_feature, input$pca_filter_selections,
                            input$pca_x_axis, input$pca_y_axis,
                            input$plot_ellipse, input$hide_legend)
        })
        
        # Generate the PC test options (i.e. selections for feature, group 1, group 2, PC)
        output$pc_test_options <- renderUI({
          div(
            div(style = "display: flex; justify-content: space-between;",
                p("Association between selected feature and PC."),
                dropdownButton(p("The currently plotted feature and the PC specified on the x axis are used for the test.",
                                 "A Wilcoxon rank sum test is performed between the two groups specified."),
                               circle = TRUE, 
                               status = "primary",
                               size = "xs",
                               icon = icon("question"), 
                               width = "300px")),
            div(style = "display: inline-block; vertical-align: top; padding-right: 10px;",
                selectInput(inputId = "pc_test_group1",
                            label = "Group 1",
                            choices = NULL,
                            multiple = FALSE,
                            width = "150px")),
            div(style = "display: inline-block; vertical-align: top; padding-left: 10px;",
                selectInput(inputId = "pc_test_group2",
                            label = "Group 2",
                            choices = NULL,
                            multiple = FALSE,
                            width = "150px")),
            actionButton(inputId = "run_pc_test",
                         label = "Run Test",
                         style = "text-align center; width: 150px")
          )
        })
        
        # Generate the button used to display the scree plot
        output$scree_plot_button <- renderUI({
          actionButton(inputId = "display_scree_plot",
                       label = "Show percent variance explained",
                       style = "text-align center; width: 250px")
        })
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### Observer for PCA plot and PC test options
        #==============================================================================================================
        
        # Update options in PCA plot filter
        observeEvent(input$pca_filter_feature, {
          if (input$pca_filter_feature != "---") {
            # Update the options for the PCA filter selections
            updateSelectInput(session,
                              inputId = "pca_filter_selections",
                              choices = unique(pca_df[[input$pca_filter_feature]]),
                              selected = NULL)
          } else {
            # Update the options for the PCA filter selections
            updateSelectInput(session,
                              inputId = "pca_filter_selections",
                              choices = NULL)
          }
        })
        
        # Update options in PC test
        observeEvent(input$pca_group, {
          if (input$pca_group != "---") {
            # Update the options for group 1
            updateSelectInput(session,
                              inputId = "pc_test_group1",
                              choices = unique(pca_df[input$pca_group]),
                              selected = unique(pca_df[input$pca_group])[1])
            # Update the options for group 2
            updateSelectInput(session,
                              inputId = "pc_test_group2",
                              choices = unique(pca_df[input$pca_group]),
                              selected = ifelse(length(unique(pca_df[input$pca_group])) > 1, 
                                                unique(pca_df[input$pca_group])[2],
                                                unique(pca_df[input$pca_group])[1]))
          } else {
            # Update the options for group 1
            updateSelectInput(session,
                              inputId = "pc_test_group1",
                              choices = NULL)
            # Update the options for group 2
            updateSelectInput(session,
                              inputId = "pc_test_group2",
                              choices = NULL)
          }
        })
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### Observer for test between selected feature and selected PC
        #==============================================================================================================
        observeEvent(input$run_pc_test, {
          
          if ((nchar(input$pc_test_group1) > 0) & (nchar(input$pc_test_group2) > 0)) {
            # Run a Wilcoxon rank sum test between group 1 and group 2 for the selected feature and PC
            p_value <- run_pc_test(pca_df, input$pca_x_axis, input$pca_group, 
                                   input$pca_filter_feature, input$pca_filter_selections,
                                   input$pc_test_group1, input$pc_test_group2)
            
            if (!is.null(p_value)) {
              # Display the p-value for the test between the selected feature and the selected PC
              output$pc_test_output <- renderUI({
                p(paste0(isolate(input$pca_group), " and PC", isolate(input$pca_x_axis), ": p = ", p_value))
              })
              
              # Display a box plot of the PC coordinates
              output$pc_test_plot <- renderPlot({
                generate_pc_boxplot(pca_df, input$pca_x_axis, input$pca_group, 
                                    input$pca_filter_feature, input$pca_filter_selections,
                                    input$pc_test_group1, input$pc_test_group2, p_value)
              })
            } else {
              # Display a warning message
              output$pc_test_output <- renderUI({
                p("Cannot run test due to insufficient groups. Please change your filters.", style = "color: red;")
              })
              
              output$pc_test_plot <- renderPlot({
                NULL
              })
            }
          }
          
        })
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### Observer for scree plot modal
        #==============================================================================================================
        # Event reactive for scree plot generation
        reactive_scree_plot <- eventReactive(input$display_scree_plot, {
          generate_scree_plot(pc_variance)
        })
        
        # Generate the scree plot
        output$scree_plot <- renderPlot({
          reactive_scree_plot()
        })
        
        # Display the modal with the scree plot
        observeEvent(input$display_scree_plot, {
          showModal(modalDialog(
            plotOutput("scree_plot"),
            footer = NULL,
            easyClose = TRUE
          ))
        })
        #==============================================================================================================
        
    
        #==============================================================================================================
        ### DGE VIEWER
        #==============================================================================================================
        
        # Generate the options for the volcano plot and comparison metrics table (comparison and gene selection)
        output$volcano_and_comparison_options <- renderUI({
            div(style = "display: inline-block; vertical-align: top; padding-left: 10px; padding-right: 10px;",
                selectInput(inputId = "volcano_comparison_input",
                            label = "Select a comparison of interest",
                            choices = setNames(comparisons$comparison, gsub(";", "-", comparisons$comparison)),
                            multiple = FALSE,
                            width = "500px"))
        })
        
        # Update the genes options
        updateSelectizeInput(session, "volcano_genes", choices = unique(comparison_data$gene), server = TRUE)
        
        # Display the table with comparison metrics for the top 50 DE genes
        output$comparison_data_table <- DT::renderDataTable(server = FALSE, {
          generate_comparison_table(comparisons, comparisons$comparison[1], comparison_data, active_comparison_data)
        })
        
        # Display the buttons under with comparison metrics table (download table, show missing genes, jump to signature viewer)
        output$comparison_data_table_buttons <- renderUI({
          div(style = "display: flex; justify-content: space-between;",
              downloadButton(outputId = "comparison_data_table_download_button",
                             label = "Download Table",
                             style = "text-align center; width: 200px"),
              # actionButton(inputId = "comparison_data_table_show_missing_genes",
              #              label = "Show Missing Genes",
              #              style = "text-align center; width: 200px; padding-right: 5px;")
              )
        })
        
        # DOWNLOAD HANDLER for comparison data table
        output$comparison_data_table_download_button <- downloadHandler(
          filename = function() {
            paste0(dataset_acc, "_", input$volcano_comparison_input, ".txt")
          },
          content = function(file) {
            # Download the zip file
            write.table(active_comparison_data(), file = file, sep = "\t", row.names = FALSE)
          },
          contentType = "txt"
        )
        
        # Generate the volcano plot
        output$volcano_plot <- renderPlot({
          values$volcano_plot <- generate_volcano_plot(comparisons, comparisons$comparison[1], comparison_data, "p_value", 1, 0.05)
          values$volcano_plot
        })
    
        # Display the number of DE genes
        output$volcano_significant_genes <- renderUI({
          determine_num_signif_genes(comparisons, comparisons$comparison[1], comparison_data, "p_value", 1, 0.05)
        })
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### Observers for volcano plot/comparison table/expression box plot viewer
        #==============================================================================================================
        
        ### OBSERVER FOR COMPARISON INPUT ###
        observeEvent(input$volcano_comparison_input, {
          
          if (all(!is.null(c(input$volcano_comparison_input,input$volcano_pval_option, input$volcano_log_fc_threshold, input$p_value_threshold)))) {
            # Pull the comparison data for the selected comparison and append to the existing comparison data
            selected_comparison_id <- comparisons %>%
              filter(comparison == input$volcano_comparison_input) %>%
              pull(id)
            if (!(selected_comparison_id %in% unique(comparison_data$comparison_id))) {
              selected_comparison_data <- get_comparison_data(selected_comparison_id, gene_map, db)
              comparison_data <<- rbind(comparison_data, selected_comparison_data)
              showNotification("Comparison data done loading", type = "default")
            }
            
            if (length(input$volcano_genes) == 0) {
              # Update the table with comparison metrics for the top 50 DE genes
              output$comparison_data_table <- DT::renderDataTable(server = FALSE, {
                generate_comparison_table(comparisons, input$volcano_comparison_input, comparison_data, active_comparison_data)
              })
              # Update the volcano plot
              output$volcano_plot <- renderPlot({
                values$volcano_plot <- generate_volcano_plot(comparisons, input$volcano_comparison_input, comparison_data,
                                      input$volcano_pval_option, input$volcano_log_fc_threshold, input$volcano_p_value_threshold)
              
                values$volcano_plot  
                #print("No update!")
              })
            } else {
              # Update the table with comparison metrics for the selected genes
              output$comparison_data_table <- DT::renderDataTable(server = FALSE, {
                update_comparison_table(input$volcano_genes, comparisons, input$volcano_comparison_input, 
                                        comparison_data, active_comparison_data)
              })
              # Update the volcano plot and label the selected genes
              output$volcano_plot <- renderPlot({
                update_volcano_plot(input$volcano_genes, comparisons, input$volcano_comparison_input, 
                                    comparison_data, input$volcano_pval_option, 
                                    input$volcano_log_fc_threshold, input$volcano_p_value_threshold)
                
                #print("Update1!!!!")
                })
            }
            
            # Update the number of DE genes
            output$volcano_significant_genes <- renderUI({
              determine_num_signif_genes(comparisons, 
                                         input$volcano_comparison_input, 
                                         comparison_data, 
                                         input$volcano_pval_option,
                                         input$volcano_log_fc_threshold, 
                                         input$volcano_p_value_threshold)
            })
          }
          
        }, ignoreInit = TRUE)
        ##########
        
        
        ## OBSERVER FOR GENES ###
        observeEvent(input$volcano_genes, {
          # Update the table with comparison metrics and volcano plot to reflect the selected genes
          if (is.null(input$volcano_genes)) {
            if(is.null(input$volcano_comparison_input)){
              output$comparison_data_table <- DT::renderDataTable(server = FALSE, {
                generate_comparison_table(comparisons, comparisons$comparison[1], comparison_data, active_comparison_data)
              })
            }else{
              output$comparison_data_table <- DT::renderDataTable(server = FALSE, {
                generate_comparison_table(comparisons, input$volcano_comparison_input, comparison_data, active_comparison_data)
              })
            }
            
          } else {
            output$comparison_data_table <- DT::renderDataTable(server = FALSE, {
              update_comparison_table(input$volcano_genes, comparisons, input$volcano_comparison_input, comparison_data, active_comparison_data)
            })
          }
        }, ignoreNULL = FALSE)
        
        
        
        #########
        
        #shinyjs::hide("box_plot_options")
        ### OBSERVER FOR COMPARISON METRICS TABLE ROW SELECTION ###
        
        
        observeEvent(input$comparison_data_table_rows_selected, {
          #print(input$comparison_data_table_rows_selected)
          # Pull the expression data given the selected genes
          selected_genes <- active_comparison_data()$gene[input$comparison_data_table_rows_selected]
          exp_data <- get_exp_data(selected_genes, gene_map, dataset_acc_identifier, sample_ann, db)
          
          # Label the selected genes in the volcano plot
          # output$volcano_plot <- renderPlot({
          #   update_volcano_plot(selected_genes,
          #                       comparisons,
          #                       input$volcano_comparison_input,
          #                       comparison_data,
          #                       input$volcano_pval_option,
          #                       input$volcano_log_fc_threshold,
          #                       input$volcano_p_value_threshold)
          # 
          #   #print("Update!")
          # })
          # 
          # Generate a list of expression plots based on selected rows
          # print(gene_plots)
          gene_plots <- generate_gene_plots(input$comparison_data_table_rows_selected,
                                            comparisons,
                                            input$volcano_comparison_input,
                                            active_comparison_data(),
                                            exp_data,
                                            sample_group_map,
                                            "p_value",
                                            FALSE)

          #shinyjs::show("box_plot_options")
          values$box_plots<-gene_plots

          # Define the plot layout (index for each plot) and generate the plotOutput objects for the expression plots
          plot_layout <- determine_plot_layout(length(gene_plots))
          output$gene_plots <- renderUI({
            splitLayout(cellWidths = c("33%", "33%", "33%"),
                        lapply(names(gene_plots)[plot_layout[[1]]], function(x) {
                          plotOutput(x)
                        }),
                        lapply(names(gene_plots)[plot_layout[[2]]], function(x) {
                          plotOutput(x)
                        }),
                        lapply(names(gene_plots)[plot_layout[[3]]], function(x) {
                          plotOutput(x)
                        }))
          })

          # Generate the expression plots for the selected genes
          sapply(names(gene_plots), function(x) {
            output[[x]] <- renderPlot({
              gene_plots[[x]]
            })
          })

          # Label the selected genes in the volcano plot
          output$volcano_plot <- renderPlot({
            update_volcano_plot(selected_genes,
                                comparisons,
                                input$volcano_comparison_input,
                                comparison_data,
                                input$volcano_pval_option,
                                input$volcano_log_fc_threshold,
                                input$volcano_p_value_threshold)
            #print("Update!")
          })

          
        })
        
        
        ##########
        
        
        ### Change between p-value/adjusted p-value and jitter/no jitter display in the expression plots ###
        observeEvent(input$gene_plot_display_pval, {
          # Pull the expression data given the selected genes
          selected_genes <- active_comparison_data()$gene[input$comparison_data_table_rows_selected]
          exp_data <- get_exp_data(selected_genes, gene_map, dataset_acc_identifier, sample_ann, db)
          
          # Generate a list of expression plots based on selected rows
          gene_plots <- generate_gene_plots(input$comparison_data_table_rows_selected, 
                                            comparisons, 
                                            input$volcano_comparison_input, 
                                            active_comparison_data(), 
                                            exp_data,
                                            sample_group_map,
                                            input$gene_plot_display_pval,
                                            input$gene_plot_display_jitter)
          
          # Generate the expression plots for the selected genes
          sapply(names(gene_plots), function(x) {
            output[[x]] <- renderPlot({
              gene_plots[[x]]
            })
          })
        }, ignoreInit = TRUE)
        
        observeEvent(input$gene_plot_display_jitter, {
          # Pull the expression data given the selected genes
          selected_genes <- active_comparison_data()$gene[input$comparison_data_table_rows_selected]
          exp_data <- get_exp_data(selected_genes, gene_map, dataset_acc_identifier, sample_ann, db)
          
          # Generate a list of expression plots based on selected rows
          gene_plots <- generate_gene_plots(input$comparison_data_table_rows_selected, 
                                            comparisons, 
                                            input$volcano_comparison_input, 
                                            active_comparison_data(), 
                                            exp_data,
                                            sample_group_map,
                                            input$gene_plot_display_pval,
                                            input$gene_plot_display_jitter)
          # Generate the expression plots for the selected genes
          sapply(names(gene_plots), function(x) {
            output[[x]] <- renderPlot({
              gene_plots[[x]]
            })
          })
        }, ignoreInit = TRUE)
        ##########
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### Observer for missing genes modal
        #==============================================================================================================
        # Event reactive for missing genes (when a gene set is selected)
        reactive_missing_genes <- eventReactive(input$comparison_data_table_show_missing_genes, {
          if (input$signature_geneset_selection) {
            selected_genes <- get_geneset_data(input$signature_geneset_input, db)
          } else {
            selected_genes <- strsplit(input$signature_text_input, "\n")
          }
          selected_genes <- unlist(selected_genes)
          
          # Pull the expression data and display the missing genes
          exp_data <- get_exp_data(selected_genes, gene_map, dataset_acc_identifier, sample_ann, db)
          setdiff(selected_genes, exp_data$gene)
        })
        
        # Generate the missing genes
        output$comparison_data_table_missing_genes <- renderUI({
          div(
            h5(paste("The following", length(reactive_missing_genes()), "genes were not found in this dataset")),
            p(paste(reactive_missing_genes(), collapse = ", "))
          )
        })
        
        # Display the modal with the missing genes
        observeEvent(input$comparison_data_table_show_missing_genes, {
          showModal(modalDialog(
            uiOutput("comparison_data_table_missing_genes"),
            footer = NULL,
            easyClose = TRUE
          ))
        })
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### SIGNATURE VIEWER
        #==============================================================================================================
    
        # Generate the options for the expression plots (group selections and gene specifications)
        output$signature_plot_options <- renderUI({
          div(
            h4("Group Selections"),
            div(style = "border-style: solid; border-width: 1px; border-radius: 5px; padding: 5px;",
                div(style = "display: flex; justify-content: space-between;",
                    selectInput(inputId = "signature_comparison_input",
                                label = HTML("Select one or multiple comparisons of interest.<br/>P-values will be displayed in the figure(s)."),
                                choices = setNames(comparisons$comparison, gsub(";", "-", comparisons$comparison)),
                                multiple = TRUE,
                                width = "450px"),
                    prettyCheckbox(inputId = "signature_comparison_selection",
                                   label = NULL, 
                                   value = TRUE,
                                   status = "danger",
                                   shape = "curve")),
                p(span("OR", style = "background: #fff; padding: 0 10px;"),
                  style = "font-size: 16; font-weight: bold; text-align: center; border-bottom: 1px solid #000; line-height: 0.1em; margin: 10px 0 20px;"),
                div(style = "display: flex; justify-content: space-between;",
                    selectInput(inputId = "signature_group_input",
                                label = HTML("Select one or multiple groups of interest.<br/>P-values will NOT be displayed in the figure."),
                                choices = sample_group_map$group,
                                multiple = TRUE,
                                width = "450px"),
                    prettyCheckbox(inputId = "signature_group_selection",
                                   label = NULL, 
                                   value = FALSE,
                                   status = "danger",
                                   shape = "curve"))
            ),
            br(),
            h4("Gene Specifications"),
            div(style = "border-style: solid; border-width: 1px; border-radius: 5px; padding: 5px;",
              textAreaInput(inputId = "signature_text_input",
                            label = "Paste a list of gene names (each gene must be on a new line)",
                            placeholder = "paste genes here",
                            width = "450px",
                            height = "150px")
            )
          )
        })
        #==============================================================================================================
    
    
        #==============================================================================================================
        ### ObserveEvents for SIGNATURE VIEWER plotting
        #==============================================================================================================
        
        # UPDATE COMPARISON/GROUP SELECTION
        observeEvent(input$signature_comparison_selection, {
          updatePrettyCheckbox(
            session,
            "signature_group_selection",
            value = ifelse(input$signature_comparison_selection, FALSE, TRUE)
          )
        })
        observeEvent(input$signature_group_selection, {
          updatePrettyCheckbox(
            session,
            "signature_comparison_selection",
            value = ifelse(input$signature_group_selection, FALSE, TRUE)
          )
        })
        
        ### GENERATE SIGNATURE PLOTS ###
        observeEvent(input$generate_signature_plot, {
          
          isError <- FALSE
          if((input$signature_comparison_selection==TRUE && is.null(input$signature_comparison_input)) || (input$signature_group_selection==TRUE && is.null(input$signature_group_input)) || (input$signature_text_input=="")){
            shinyalert(title = "Insufficient Input Supplied",
                       text = "Please select at least one comparison/condition and input at least one gene symbol!",
                       type = "error")
            isError <- TRUE
          }
         
          if(isError == FALSE){
            # Pull the appropriate genes
            selected_genes <- strsplit(input$signature_text_input, "\n")
            selected_genes <- unlist(selected_genes)
            
            # Pull the expression data
            exp_data <- get_exp_data(selected_genes, gene_map, dataset_acc_identifier, sample_ann, db)
            
            # Display a notification if at least one specified gene is missing from the dataset
            genes_not_found <- setdiff(selected_genes, exp_data$gene)
            if (length(genes_not_found) != 0) {
              showNotification(paste(length(genes_not_found), "genes were not fonund in the dataset:",
                                     str_trunc(paste(selected_genes, collapse = ", "), width = 50, ellipsis = "...")),
                               type = "warning")
            }
            
            # Generate the expression plot(s)
            if (input$signature_comparison_selection) {
              # Pull the comparison data for the selected comparison and append to the existing comparison data
              selected_comparison_ids <- comparisons %>%
                filter(comparison == input$signature_comparison_input) %>%
                pull(id)
              selected_comparison_ids <- selected_comparison_ids[!(selected_comparison_ids %in% unique(comparison_data$comparison_id))]
              if (length(selected_comparison_ids) != 0 & length(selected_genes) == 1) {
                selected_comparison_data <- get_comparison_data(selected_comparison_ids, gene_map, db)
                comparison_data <<- rbind(comparison_data, selected_comparison_data)
                showNotification("Comparison data done loading", type = "default")
              }
              
              # Generate a list of expression plots (one per selected comparison)
              exp_plots <- generate_comparison_exp_plots(selected_genes, 
                                                         input$signature_comparison_input, 
                                                         exp_data, 
                                                         sample_group_map, 
                                                         comparisons, 
                                                         comparison_data,
                                                         "p_value",
                                                         FALSE)
              # Define the plot layout (index for each plot)
              plot_layout <- determine_plot_layout(length(exp_plots))
              
              # Generate the plotOutput objects for the expression plots
              output$signature_plots <- renderUI({
                splitLayout(cellWidths = c("33%", "33%", "33%"),
                            lapply(names(exp_plots)[plot_layout[[1]]], function(x) {
                              plotOutput(x, height = "250px")
                            }),
                            lapply(names(exp_plots)[plot_layout[[2]]], function(x) {
                              plotOutput(x, height = "250px")
                            }),
                            lapply(names(exp_plots)[plot_layout[[3]]], function(x) {
                              plotOutput(x, height = "250px")
                            }))
              })
              
              # Generate the expression plots for the selected comparisons
              sapply(names(exp_plots), function(x) {
                output[[x]] <- renderPlot({
                  exp_plots[[x]]
                })
              })
              
              # Remove the table of comparisons if it exists
              output$signature_plot_table <- renderDataTable({
                NULL
              })
            } else {
              num_groups <- length(input$signature_group_input)
              
              # Generate the plotOutput object for the grouped expression plot
              output$signature_plots <- renderUI({
                plotOutput("signature_grouped_plot")
              })
              # Generate the grouped expression plot (one figure with all selected groups)
              output$signature_grouped_plot <- renderPlot({
                generate_grouped_exp_plot(selected_genes, 
                                          input$signature_group_input, 
                                          exp_data, 
                                          sample_group_map,
                                          FALSE)
              })
              
              # Generate the table of comparisons for the generated plot
              output$signature_plot_table <- renderDataTable({
                generate_signature_plot_table(selected_genes, 
                                              input$signature_group_input, 
                                              comparisons,
                                              exp_data, 
                                              sample_group_map)
              })
            }
          }
          
        })
        ##########
        
        ### Change between p-value/adjusted p-value and jitter/no jitter display in the expression plots ###
        observeEvent(input$signature_display_pval, {
          # Pull the appropriate genes based on whether the user wants to use the pasted genes, uploaded genes, or gene set
          if (input$signature_text_selection) {
            selected_genes <- strsplit(input$signature_text_input, "\n")
          } else {
            selected_genes <- get_geneset_data(input$signature_geneset_input, db)
          }
          selected_genes <- unlist(selected_genes)
          
          # Pull the expression data
          exp_data <- get_exp_data(selected_genes, gene_map, dataset_acc_identifier, sample_ann, db)
          
          if (input$signature_comparison_selection) {
            # Pull the comparison data for the selected comparison and append to the existing comparison data
            selected_comparison_ids <- comparisons %>%
              filter(comparison == input$signature_comparison_input) %>%
              pull(id)
            selected_comparison_ids <- selected_comparison_ids[!(selected_comparison_ids %in% unique(comparison_data$comparison_id))]
            if (length(selected_comparison_ids) != 0 & length(selected_genes) == 1) {
              selected_comparison_data <- get_comparison_data(selected_comparison_ids, gene_map, db)
              comparison_data <<- rbind(comparison_data, selected_comparison_data)
              showNotification("Comparison data done loading", type = "default")
            }
            
            # Generate a list of expression plots (one per selected comparison)
            exp_plots <- generate_comparison_exp_plots(selected_genes, 
                                                       input$signature_comparison_input, 
                                                       exp_data, 
                                                       sample_group_map, 
                                                       comparisons, 
                                                       comparison_data,
                                                       input$signature_display_pval,
                                                       input$signature_display_jitter)
            
            # Generate the expression plots for the selected comparisons
            sapply(names(exp_plots), function(x) {
              output[[x]] <- renderPlot({
                exp_plots[[x]]
              })
            })
            
            # Remove the table of comparisons if it exists
            output$signature_plot_table <- renderDataTable({
              NULL
            })
          }
        }, ignoreInit = TRUE, ignoreNULL = TRUE)
        
        observeEvent(input$signature_display_jitter, {
          # Pull the appropriate genes based on whether the user wants to use the pasted genes, uploaded genes, or gene set
          if (input$signature_text_selection) {
            selected_genes <- strsplit(input$signature_text_input, "\n")
          } else {
            selected_genes <- get_geneset_data(input$signature_geneset_input, db)
          }
          selected_genes <- unlist(selected_genes)
          
          # Pull the expression data
          exp_data <- get_exp_data(selected_genes, gene_map, dataset_acc_identifier, sample_ann, db)
          
          if (input$signature_comparison_selection) {
            # Pull the comparison data for the selected comparison and append to the existing comparison data
            selected_comparison_ids <- comparisons %>%
              filter(comparison == input$signature_comparison_input) %>%
              pull(id)
            selected_comparison_ids <- selected_comparison_ids[!(selected_comparison_ids %in% unique(comparison_data$comparison_id))]
            if (length(selected_comparison_ids) != 0 & length(selected_genes) == 1) {
              selected_comparison_data <- get_comparison_data(selected_comparison_ids, gene_map, db)
              comparison_data <<- rbind(comparison_data, selected_comparison_data)
              showNotification("Comparison data done loading", type = "default")
            }
            
            # Generate a list of expression plots (one per selected comparison)
            exp_plots <- generate_comparison_exp_plots(selected_genes, 
                                                       input$signature_comparison_input, 
                                                       exp_data, 
                                                       sample_group_map, 
                                                       comparisons, 
                                                       comparison_data,
                                                       input$signature_display_pval,
                                                       input$signature_display_jitter)
            
            # Generate the expression plots for the selected comparisons
            sapply(names(exp_plots), function(x) {
              output[[x]] <- renderPlot({
                exp_plots[[x]]
              })
            })
            
            # Remove the table of comparisons if it exists
            output$signature_plot_table <- renderDataTable({
              NULL
            })
          } else {
            num_groups <- length(input$signature_group_input)
            # Generate the grouped expression plot (one figure with all selected groups)
            output$signature_grouped_plot <- renderPlot({
              generate_grouped_exp_plot(selected_genes, 
                                        input$signature_group_input, 
                                        exp_data, 
                                        sample_group_map,
                                        input$signature_display_jitter)
            })
          }
        }, ignoreInit = TRUE, ignoreNULL = TRUE)
        ##########
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### ENRICHMENT ANALYSIS
        #==============================================================================================================
        
        # Set the default enrichment analysis to ORA
        output$enrichment_options <- renderUI({
          div(
            p("This enrichment analysis option uses WebGestAlt to perform Over-Representation Analysis on the genes",
              "pertaining to the selected comparison."),
            br(),
            selectInput(inputId = "enrichment_comparison_input",
                        label = "Select a comparison of interest for the enrichment analysis",
                        choices = setNames(comparisons$comparison, gsub(";", "-", comparisons$comparison)),
                        multiple = FALSE,
                        width = "575px"),
            div(
                div(style = "display: inline-block;",
                    radioGroupButtons(inputId = "enrichment_pval_option",
                                      label = NULL,
                                      choices = c("p-value" = "p_value", "adjusted p-value" = "p_value_adj"),
                                      selected = "p_value",
                                      width = "300px",
                                      direction = "horizontal",
                                      checkIcon = list(
                                        yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                        no = tags$i(class = "fa fa-square-o", style = "color: steelblue")))),
                div(style = "display: inline-block; padding-right: 15px;",
                    numericInput(inputId = "enrichment_p_value_threshold",
                                 label = "p-value cutoff",
                                 value = 0.05,
                                 min = 0,
                                 max = 1,
                                 step = 0.01,
                                 width = "150px")),
                div(style = "display: inline-block; padding-left: 15px;",
                    radioGroupButtons(inputId = "enrichment_logfc_option",
                                      label = NULL,
                                      choices = c("up", "down"),
                                      selected = "up",
                                      width = "250px",
                                      direction = "horizontal",
                                      checkIcon = list(
                                        yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                        no = tags$i(class = "fa fa-square-o", style = "color: steelblue")))),
                div(style = "display: inline-block;",
                    numericInput(inputId = "enrichment_log_fc_threshold",
                                 label = "logFC cutoff",
                                 value = 1,
                                 min = 0,
                                 max = 10,
                                 step = 0.1,
                                 width = "150px"))),
            splitLayout(cellWidths = c("60%", "40%"),
                        list(
                          selectInput(inputId = "enrichment_database_type_input",
                                      label = "Select a type of database (changing this selection will alter available options on the right)",
                                      choices = names(list_enrichment_databases()),
                                      multiple = FALSE,
                                      width = "100%"),
                          br(), br(), br(), br(), br()
                        ),
                        list(
                          selectInput(inputId = "enrichment_database_input",
                                      label = "Select a database to use",
                                      choices = list_enrichment_databases(),
                                      multiple = FALSE,
                                      width = "100%"),
                          br(), br(), br(), br(), br()
                        )
            ),
            div(style = "display: flex; justify-content: space-between;",
                p(""),
                actionBttn(inputId = "run_ora",
                           label = "Run ORA",
                           class = "btn enrichment-button",
                           width = "400px",
                           style = "unite",
                           color = "danger"),
                p(""))
          )
        })
        
        # Generate the options for ORA or GSEA if the associated option is selected
        observeEvent(input$enrichment_option, {
          if (input$enrichment_option == "ORA") {
            output$enrichment_options <- renderUI({
              div(
                p("This enrichment analysis option uses WebGestalt to perform Over-Representation Analysis on the genes",
                  "pertaining to the selected comparison."),
                br(),
                selectInput(inputId = "enrichment_comparison_input",
                            label = "Select a comparison of interest for the enrichment analysis",
                            choices = setNames(comparisons$comparison, gsub(";", "-", comparisons$comparison)),
                            multiple = FALSE,
                            width = "575px"),
                div(
                  div(style = "display: inline-block;",
                      radioGroupButtons(inputId = "enrichment_pval_option",
                                        label = NULL,
                                        choices = c("p-value" = "p_value", "adjusted p-value" = "p_value_adj"),
                                        selected = "p_value",
                                        width = "300px",
                                        direction = "horizontal",
                                        checkIcon = list(
                                          yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                          no = tags$i(class = "fa fa-square-o", style = "color: steelblue")))),
                  div(style = "display: inline-block; padding-right: 15px;",
                      numericInput(inputId = "enrichment_p_value_threshold",
                                   label = "p-value cutoff",
                                   value = 0.05,
                                   min = 0,
                                   max = 1,
                                   step = 0.01,
                                   width = "150px")),
                  div(style = "display: inline-block; padding-left: 15px;",
                      radioGroupButtons(inputId = "enrichment_logfc_option",
                                        label = NULL,
                                        choices = c("up", "down"),
                                        selected = "up",
                                        width = "250px",
                                        direction = "horizontal",
                                        checkIcon = list(
                                          yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                          no = tags$i(class = "fa fa-square-o", style = "color: steelblue")))),
                  div(style = "display: inline-block;",
                      numericInput(inputId = "enrichment_log_fc_threshold",
                                   label = "logFC cutoff",
                                   value = 1,
                                   min = 0,
                                   max = 10,
                                   step = 0.1,
                                   width = "150px"))),
                splitLayout(cellWidths = c("60%", "40%"),
                            list(
                              selectInput(inputId = "enrichment_database_type_input",
                                          label = "Select a type of database (changing this selection will alter available options on the right)",
                                          choices = names(list_enrichment_databases()),
                                          multiple = FALSE,
                                          width = "100%"),
                              br(), br(), br(), br(), br()
                            ),
                            list(
                              selectInput(inputId = "enrichment_database_input",
                                          label = "Select a database to use",
                                          choices = list_enrichment_databases(),
                                          multiple = FALSE,
                                          width = "100%"),
                              br(), br(), br(), br(), br()
                            )
                ),
                div(style = "display: flex; justify-content: space-between;",
                    p(""),
                    actionBttn(inputId = "run_ora",
                               label = "Run ORA",
                               class = "btn enrichment-button",
                               width = "400px",
                               style = "unite",
                               color = "danger"),
                    p(""))
              )
            })
          } else {
            output$enrichment_options <- renderUI({
              div(
                p("This enrichment analysis option uses WebGestalt to perform Gene Set Eenrichment Analysis on the genes",
                  "pertaining to the selected comparison."),
                p("The signed -log10(p-value) is used as the ranking metric."),
                br(),
                selectInput(inputId = "enrichment_comparison_input",
                            label = "Select a comparison of interest for the enrichment analysis",
                            choices = setNames(comparisons$comparison, gsub(";", "-", comparisons$comparison)),
                            multiple = FALSE,
                            width = "700px"),
                splitLayout(cellWidths = c("60%", "40%"),
                            list(
                              selectInput(inputId = "enrichment_database_type_input",
                                          label = "Select a type of database (changing this selection will alter available options on the right)",
                                          choices = names(list_enrichment_databases()),
                                          multiple = FALSE,
                                          width = "100%"),
                              br(), br(), br(), br(), br()
                            ),
                            list(
                              selectInput(inputId = "enrichment_database_input",
                                          label = "Select a database to use",
                                          choices = list_enrichment_databases(),
                                          multiple = FALSE,
                                          width = "100%"),
                              br(), br(), br(), br(), br()
                            )
                ),
                div(style = "display: flex; justify-content: space-between;",
                    p(""),
                    actionBttn(inputId = "run_gsea",
                               label = "Run GSEA",
                               class = "btn enrichment-button",
                               width = "400px",
                               style = "unite",
                               color = "danger"),
                    p(""))
              )
            })
          }
        })
        #==============================================================================================================
    
    
        #==============================================================================================================
        ### ObserveEvent for ENRICHMENT ANALYSIS
        #==============================================================================================================
        
        # Update the database options when the user selects a database type(s)
        observeEvent(input$enrichment_database_type_input, {
          enrichment_databases <- list_enrichment_databases()
          updateSelectInput(session, 
                            "enrichment_database_input",
                            choices = enrichment_databases[[input$enrichment_database_type_input]])
        })
        
        # RUN ORA AND PRODUCE ENRICHMENT RESULTS TABLE
        observeEvent(input$run_ora, {
          # Pull the comparison data for the selected comparison and append to the existing comparison data
          selected_comparison_id <- comparisons %>%
            filter(comparison == input$enrichment_comparison_input) %>%
            pull(id)
          if (!(selected_comparison_id %in% unique(comparison_data$comparison_id))) {
            selected_comparison_data <- get_comparison_data(selected_comparison_id, gene_map, db)
            comparison_data <<- rbind(comparison_data, selected_comparison_data)
            showNotification("Comparison data done loading", type = "default")
          }
          
          # Generate the organism map between names in the DB and names used by Webgestalt
          organism_map <- generate_organism_map()
          
          # Generate the enrichment results and return the directory with the html report
          project_name <- perform_ora(input$enrichment_comparison_input,
                                      organism_map[[dataset_description$organism]],
                                      input$enrichment_database_input,
                                      dataset_acc,
                                      input$enrichment_pval_option, 
                                      input$enrichment_p_value_threshold,
                                      input$enrichment_logfc_option,
                                      input$enrichment_log_fc_threshold, 
                                      comparisons, 
                                      comparison_data)
          
          if (is.null(project_name)) {
            # Display a message stating no genes met the thresholds
            showModal(
              modalDialog(
                h4("No genes meet the specified cutoffs or no significant gene set. Please make adjustments to run ORA."),
                easyClose = TRUE,
                footer = NULL
              )
            )
          } else {
            # Modify the html report so that it is ready for viewing
            # modify_enrichment_report(project_name)
            # Move the associated files to the www directory so they can be accessed
            file.copy(paste0("./enrichment_results/Project_", project_name), "./www/", recursive = TRUE)
            
            # Display the html report
            output$enrichment_report <- renderUI({
              tags$iframe(seamless = "seamless",
                          src = paste0("Project_", project_name, "/Report_", project_name, ".html"),
                          width = "100%",
                          style ="height: 100vh;")
            })
          }
        })
        
        # RUN GSEA AND PRODUCE ENRICHMENT RESULTS TABLE
        observeEvent(input$run_gsea, {
          # Pull the comparison data for the selected comparison and append to the existing comparison data
          selected_comparison_id <- comparisons %>%
            filter(comparison == input$enrichment_comparison_input) %>%
            pull(id)
          if (!(selected_comparison_id %in% unique(comparison_data$comparison_id))) {
            selected_comparison_data <- get_comparison_data(selected_comparison_id, gene_map, db)
            comparison_data <<- rbind(comparison_data, selected_comparison_data)
            showNotification("Comparison data done loading", type = "default")
          }
          
          # Generate the organism map between names in the DB and names used by Webgestalt
          organism_map <- generate_organism_map()
    
          # Generate the enrichment results and return the directory with the html report
          project_name <- perform_gsea(input$enrichment_comparison_input,
                                       organism_map[[dataset_description$organism]],
                                       input$enrichment_database_input,
                                       dataset_acc,
                                       comparisons,
                                       comparison_data)
          # Modify the html report so that it is ready for viewing
          # modify_enrichment_report(project_name)
          # Move the associated files to the www directory so they can be accessed
          file.copy(paste0("./enrichment_results/Project_", project_name), "./www/", recursive = TRUE)
          
          # Display the html report
          output$enrichment_report <- renderUI({
            tags$iframe(seamless = "seamless",
                        src = paste0("Project_", project_name, "/Report_", project_name, ".html"),
                        width = "100%",
                        style ="height: 100vh;")
          })
        })
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### CELL DECONVOLUTION
        #==============================================================================================================
        
        # List the diseases that have trained cell deconvolution models
        trained_models <- list.dirs("cell_deconvolution", full.names = FALSE)[-1]
        
        # Display the cell deconvolution options if cell deconvolution is available for the given dataset
        cell_deconvolution_available <- TRUE
        formatted_disease <- gsub(";Healthy|;nonIBD|;nonRA", "", dataset_description$disease)
        if (!(formatted_disease %in% c("CD", "UC", "CD;UC", "AD", "RA"))) {
          cell_deconvolution_available <- FALSE
        }
        if (dataset_description$organism != "Human") {
          cell_deconvolution_available <- FALSE
        }
        if (!(dataset_description$experiment_type %in% c("RNASeq", "Microarray"))) {
          cell_deconvolution_available <- FALSE
        }
        if (dataset_description$source == "Blood") {
          cell_deconvolution_available <- FALSE
        }
        
        if (cell_deconvolution_available) {
          output$cell_deconvolution_options <- renderUI({
            div(
              p("The cell deconvolution function uses SCADEN in the backend."),
              p("The model used was trained on the appropriate single cell and bulk transcriptomics data given this dataset's features."),
              tags$ul(
                tags$li(paste("Disease:", dataset_description$disease)),
                tags$li(paste("Organism:", dataset_description$organism)),
                tags$li(paste("Source:", dataset_description$source))
              ),
              p("Performing cell deconvolution will generate estimated cell fractions for each sample in this dataset."),
              br(), br(),
              div(style = "text-align: center;",
                  actionButton(inputId = "run_cell_deconvolution",
                               label = "Perform Cell Deconvolution",
                               style = "text-align center; width: 300px;"))
            )
          })
        } else {
          output$cell_deconvolution_options <- renderUI({
            div(
              h4("Cell deconvolution is currently not supported for this dataset"),
              p("Cell deconvolution is supported for datasets that meet the following criteria."),
              tags$ul(
                tags$li("Disease: AD, CD, UC, IBD"),
                tags$li("Organism: Human"),
                tags$li("Data Type: Microarray, RNA-Seq"),
              )
            )
          })
        }
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### Observers for CELL DECONVOLUTION
        #==============================================================================================================
        
        # Create a reactive value for the cell fractions table (this is necessary when running tests between groups)
        reactive_cell_fractions <- reactiveVal()
        
        # RUN CELL DECONVOLUTION BUTTON
        observeEvent(input$run_cell_deconvolution, {
          
          # Generate the header for the cell deconvolution results and the options for the cell fractions table
          output$cell_deconvolution_results <- renderUI({
            div(
              h3("Cell Deconvolution Results", style = "text-align: center;"),
              br(),
              div(style = "display: flex; justify-content: space-between;",
                  # materialSwitch(inputId = "display_cell_fractions",
                  #                label = "Display cell fractions table", 
                  #                status = "primary"),
                  downloadButton(outputId = "cell_fractions_table_download_button",
                                 label = "Download Table",
                                 style = "width: 200px")),
              br(),
              hr(),
              br(),
              p("Use the options below to visualize estimated cell fractions.",
                "Cell fraction box plots (one per cell type) will be generated for the selected comparison.",
                "Only plots for cell types that meet the specified thresholds for the p-value and minimum cell fraction will be generated.",
                "If you would like to visualize cell fractions across all cell types, set the thresholds to 1 and 0 respectively.")
            )
          })
          
          # Perform cell deconvolution
          #a <- get_cell_fractions(dataset_acc, sample_group_map, db)
          reactive_cell_fractions(get_cell_fractions(dataset_acc, sample_group_map, db))
          
          # Display the styled cell fractions table
          # output$cell_deconvolution_table <- DT::renderDataTable(server = FALSE, {
          #   generate_cell_fractions_table(reactive_cell_fractions())
          # })
          
          # Generate the options for the cell fraction plots
          output$cell_deconvolution_plot_options <- renderUI({
            div(style = "display: flex; justify-content: space-between;",
                div(style = "display: inline-block; vertical-align: middle;",
                    div(style = "display: inline-block; vertical-align: middle; padding-left: 5px; padding-right: 5px;",
                        pickerInput(inputId = "cell_deconvolution_comparison",
                                    label = "Select a comparison to visualize",
                                    choices = comparisons$comparison,
                                    selected = NULL,
                                    multiple = FALSE,
                                    inline = TRUE,
                                    width = "fit",
                                    choicesOpt = list(style = "font-size: 75%;"))),
                    div(style = "display: inline-block; vertical-align: middle; padding-left: 5px; padding-right: 5px;",
                        numericInput(inputId = "cell_deconvolution_p_value",
                                     label = "Threshold for the p-value",
                                     value = 0.05,
                                     min = 0,
                                     max = 1,
                                     step = 0.01,
                                     width = "200px")),
                    div(style = "display: inline-block; vertical-align: middle; padding-left: 5px; padding-right: 5px;",
                        numericInput(inputId = "cell_deconvolution_min_fraction",
                                     label = "Minimum cell fraction",
                                     value = 0.001,
                                     min = 0,
                                     max = 1,
                                     step = 0.001,
                                     width = "200px"))),
                div(style = "display: inline-block; vertical-align: middle;",
                    actionButton(inputId = "plot_cell_fractions",
                                 label = "Plot Cell Fractions",
                                 style = "text-align center; width: 150px"))
              
            )
          })
          
        })
        
        # Display or hide the cell fractions table 
        # observeEvent(input$display_cell_fractions, {
        #   if (input$display_cell_fractions) {
        #     show("cell_deconvolution_table")
        #   } else {
        #     hide("cell_deconvolution_table")
        #   }
        # })
        
        # GENERATE CELL FRACTION PLOTS
        observeEvent(input$plot_cell_fractions, {
          
          # Generate a list of cell fraction plots based on the specified comparisons and selected column
          cell_fraction_plots <- generate_cell_fraction_plots(input$cell_deconvolution_comparison,
                                                              reactive_cell_fractions(),
                                                              input$cell_deconvolution_p_value,
                                                              input$cell_deconvolution_min_fraction)
          
          if (is.null(cell_fraction_plots)) {
            output$cell_deconvolution_plots <- renderUI({
              h3("No cell types meet the specified thresholds")
            })
          } else {
            # Define the plot layout (index for each plot)
            plot_layout <- determine_cell_fraction_layout(length(cell_fraction_plots))
            # Generate the plotOutput objects for the cell fraction plots
            output$cell_deconvolution_plots <- renderUI({
              splitLayout(cellWidths = c("25%", "25%", "25%", "25%"),
                          lapply(names(cell_fraction_plots)[plot_layout[[1]]], function(x) {
                            plotOutput(x)
                          }),
                          lapply(names(cell_fraction_plots)[plot_layout[[2]]], function(x) {
                            plotOutput(x)
                          }),
                          lapply(names(cell_fraction_plots)[plot_layout[[3]]], function(x) {
                            plotOutput(x)
                          }),
                          lapply(names(cell_fraction_plots)[plot_layout[[4]]], function(x) {
                            plotOutput(x)
                          }))
            })
            # Generate the cell fraction plots for the selected cell type
            sapply(names(cell_fraction_plots), function(x) {
              output[[x]] <- renderPlot({
                cell_fraction_plots[[x]]
              })
            })
          }

        })
        
        # DOWNLOAD HANDLER for cell fractions table
        output$cell_fractions_table_download_button <- downloadHandler(
          filename = function() {
            paste0(dataset_acc, "_cell_fractions.txt")
          },
          content = function(file) {
            # Download the txt file
            cell_fractions_table <- reactive_cell_fractions()
            cell_fractions_table$sample_id <- rownames(cell_fractions_table)
            cell_fractions_table <- relocate(cell_fractions_table, sample_id)
            write.table(cell_fractions_table, file = file, sep = "\t", row.names = FALSE)
          },
          contentType = "txt"
        )
        
        # # CELL FRACTIONS TABLE
        # observeEvent(input$cell_deconvolution_comparison_input, {
        #   output$cell_deconvolution_table <- DT::renderDataTable(server = FALSE, {
        #     generate_cell_fractions_table(cell_fractions, 
        #                                   input$cell_deconvolution_comparison_input, 
        #                                   sample_group_map)
        #   })
        # })
        # 
        # # OBSERVER FOR CELL DECONVOLUTION TABLE COLUMN SELECTION
        # observeEvent(input$cell_deconvolution_table_columns_selected, {
        #   
        #   # Generate a list of cell fraction plots based on the specified comparisons and selected column
        #   cell_fraction_plots <- generate_cell_fraction_plots(input$cell_deconvolution_comparison, 
        #                                                       input$cell_deconvolution_table_columns_selected, 
        #                                                       cell_fractions)
        #   # Define the plot layout (index for each plot)
        #   plot_layout <- determine_cell_fraction_layout(length(cell_fraction_plots))
        #   
        #   # Generate the plotOutput objects for the cell fraction plots
        #   output$cell_deconvolution_plots <- renderUI({
        #     splitLayout(cellWidths = c("25%", "25%", "25%", "25%"),
        #                 lapply(names(cell_fraction_plots)[plot_layout[[1]]], function(x) {
        #                   plotOutput(x)
        #                 }),
        #                 lapply(names(cell_fraction_plots)[plot_layout[[2]]], function(x) {
        #                   plotOutput(x)
        #                 }),
        #                 lapply(names(cell_fraction_plots)[plot_layout[[3]]], function(x) {
        #                   plotOutput(x)
        #                 }),
        #                 lapply(names(cell_fraction_plots)[plot_layout[[4]]], function(x) {
        #                   plotOutput(x)
        #                 }))
        #   })
        #   
        #   # Generate the cell fraction plots for the selected cell type
        #   sapply(names(cell_fraction_plots), function(x) {
        #     output[[x]] <- renderPlot({
        #       cell_fraction_plots[[x]]
        #     })
        #   })
        #   
        # })
        #==============================================================================================================
        
        
        #==============================================================================================================
        ### DOWNLOAD DATASET
        #==============================================================================================================
        # Generate the comparison options and download button for the download dataset section
        output$download_options <- renderUI({
          div(
            pickerInput(inputId = "download_comparison_selection",
                        label = "Selected comparisons will be included in the download",
                        choices = comparisons$comparison,
                        selected = comparisons$comparison,
                        multiple = TRUE,
                        options = list(`actions-box` = TRUE,
                                       `multiple-separator` = " / ",
                                       style = "font-size: 75%;"),
                        choicesOpt = list(style = rep("font-size: 75%;", length(comparisons$comparison))),
                        width = "700px"),
            br(), br(), br(),
            downloadButton(outputId = "download_button",
                           label = "Download Dataset",
                           style = "text-align center; width: 200px")
          )
        })
        #==============================================================================================================
        
        #==============================================================================================================
        ### DownloadHandler for Volcano Plots
        #==============================================================================================================
        output$download_volcano_plot = downloadHandler(
          filename = function(){
            paste0(input$volcano_comparison_input, Sys.Date(), ".pdf")
          },
          content = function(file){
            l <- list()
            # pdf(file,width =14 )
            # device = function(..., width, height) {
            # grDevices::png(..., width = width, height = height, res = 600, units = "in")
            # }
            # ggsave(file, plot = plot_output(), device = device, height = 7, width = 14)
            
            #l<-generate_volcano_plot(comparisons, comparisons$comparison[1], comparison_data, "p_value", 1, 0.05)
            ggsave(file,values$volcano_plot)
            # arrangeGrob(l, nrow = 1)  
            # dev.off()
          }
        )
        #==============================================================================================================
        ### DownloadHandler for Box Plots
        #==============================================================================================================
        output$download_expression_boxplots = downloadHandler(
          filename = function(){
            paste0(input$volcano_comparison_input, Sys.Date(), ".pdf")
          },
          content = function(file){
            l <- list()
            # pdf(file,width =14 )
            # device = function(..., width, height) {
            # grDevices::png(..., width = width, height = height, res = 600, units = "in")
            # }
            # ggsave(file, plot = plot_output(), device = device, height = 7, width = 14)
            
            #l<-generate_volcano_plot(comparisons, comparisons$comparison[1], comparison_data, "p_value", 1, 0.05)
            ggplot2::ggsave(filename = file,
                            plot = gridExtra::marrangeGrob(values$box_plots, nrow = 1, ncol = 1), 
                            device = "pdf")
            # arrangeGrob(l, nrow = 1)  
            # dev.off()
          }
        )
        #==============================================================================================================
        ### DownloadHandler for dataset download
        #==============================================================================================================
        output$download_button <- downloadHandler(
          filename = function() {
            paste0(dataset_acc, ".zip")
          },
          content = function(file) {
            showModal(modalDialog(h4("Downloading Data"), footer = NULL))
            on.exit(removeModal())
            
            # Generate a named list with the files to be downloaded
            files_to_download <- prepare_files(dataset_acc,
                                               input$download_comparison_selection,
                                               dataset_description,
                                               comparisons,
                                               sample_ann)
            # Create a temporary directory to store the files
            dir.create(dataset_acc)
            # Download the zip file
            for (file_name in names(files_to_download)) {
              path <- paste0(dataset_acc, "/", file_name, ".csv")
              if (file_name == "data") {
                write.csv(files_to_download[[file_name]], file = path, row.names = TRUE)
                next
              }
              write.csv(files_to_download[[file_name]], file = path, row.names = FALSE)
            }
            zip(zipfile = file, files = paste0(dataset_acc, "/", names(files_to_download), ".csv"))
          },
          contentType = "application/zip"
        )
        #==============================================================================================================
      }
      
      session$onSessionEnded(function() {
        
        stopApp()
      })
    })
    
}

onStop(function() {
  # Close the connection to the DB
  dbDisconnect(db)
  
  # Remove the directory associated with the dataset download
  system(paste("rm -r", global_dataset_acc))
})

# Run the application 
shinyApp(ui = ui, server = server)
