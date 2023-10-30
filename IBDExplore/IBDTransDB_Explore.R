# IBDTransDB Explore Application
# An application used to explore omics data

library(RSQLite)

library(shiny)
library(shinyWidgets)
library(shinyBS)
library(shinybusy)
library(rintrojs)
library(shinyjs)

library(DT)

library(dplyr)
library(tidyr)
library(stringr)

library(rstudioapi)


#==============================================================================================================
# Import functions used to generate table summaries, table outputs, expression plots, etc
#==============================================================================================================
source("./app_functions/database_query.R")
source("./app_functions/dataset_table.R")
#==============================================================================================================


#==============================================================================================================
# Create the connection to the Immunoverse DB
#==============================================================================================================
db <- dbConnect(RSQLite::SQLite(), dbname = "./IBDTransDB.db")
#==============================================================================================================

                  
# Define the UI for the application
ui <- navbarPage(
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
        "
      )
    ),
    tags$script(HTML("var header = $('.navbar> .container-fluid');
                       header.append('<div style=\"float:right;color:white;padding-top:8px;\"><h4>Developed by Immunology Computational Biology - GRC</h3></div>');
                       console.log(header)")),
    
    useShinyjs(),
    extendShinyjs(text = "shinyjs.browseURL = function(url) {window.open(url,'_blank');}", 
                  functions = "browseURL"),
    
    title = "IBDExplore",
    
    ### DATASETS TAB ###
    tabPanel("Datasets",
             column(3,
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 850px;",
                              
                              h4("The options presented below can be used to filter for datasets that match the specified criteria"),
                              br(),
                              
                              selectInput(
                                inputId = "disease_input",
                                label = "Select diseases of interest",
                                choices = get_keywords("Disease", db),
                                multiple = TRUE
                              ),
                              br(),
                              
                              selectInput(
                                inputId = "tissue_input",
                                label = "Select tissues of interest",
                                choices = get_keywords("Source", db),
                                multiple = TRUE
                              ),
                              selectInput(
                                inputId = "cell_type_input",
                                label = "Select cell types of interest",
                                choices = get_keywords("CellType", db),
                                multiple = TRUE
                              ),
                              br(),
                              
                              selectInput(
                                inputId = "treatment_input",
                                label = "Select treatments of interest",
                                choices = get_keywords("Treatment", db),
                                multiple = TRUE
                              ),
                              selectInput(
                                inputId = "timepoint_input",
                                label = "Select timepoints of interest",
                                choices = get_keywords("Timepoint", db),
                                multiple = TRUE
                              )

                    )),
             
             column(9,
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 850px;",
                              column(12, 
                                     p("To explore a dataset of interest, select the corresponding row from the table below."),
                                     br(),
                                     DT::dataTableOutput("dataset_table"))
                              
                    ))
    ),
    
    ### ABOUT TAB ###        
    # tabPanel("About",
    #          column(8, offset = 2,
    #                 wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 775px;",
    #                           p("ImmunoExplore is a web application targeted for biologists and", 
    #                             "bioinformaticians that allows single dataset interrogation.",
    #                             "Dataset exploration functionality includes sinlge target, multi-target,",
    #                             "and gene set level queries, dimensionality reduction, differential",
    #                             "gene expression analysis for dataset relevant comparisons, expression",
    #                             "profile visualization, enrichment analysis on the fly, and cell",
    #                             "deconvolution on the fly."),
    #                           p("A list of functionality for bulk datasets and single cell datasets can be found below."),
    #                           br(),
    #                           h4("Bulk Datasets"),
    #                           tags$ul(
    #                             tags$li("A description of the data across relevant features including a study summary and the experimental design"),
    #                             tags$li("Dimensionality reduction including customizable PCA plot, feature-PC tests of association, and percent variance explained plotting"),
    #                             tags$li("Easy querying of differential gene expression analysis results across study-related comparisons"),
    #                             tags$li("Flexible gene-level and signature-level expression visualizations"),
    #                             tags$li("Enrichment analysis on the fly using ORA or preranked GSEA"),
    #                             tags$li("Cell deconvolution on the fly using SCADEN")
    #                           ),
    #                           br(),
    #                           h4("Single Cell Datasets"),
    #                           tags$ul(
    #                             tags$li("A description of the data across relevant features including a study summary and the experimental design"),
    #                             tags$li("Easy querying of differential gene expression analysis results across study-related comparisons"),
    #                             tags$li("Enrichment analysis on the fly using ORA or preranked GSEA")
    #                           )))
    # )
    # 
    ### UPLOAD DATASET TAB ###
    # tabPanel("Upload Dataset",
    #          column(8, offset = 2,
    #                 h3("Upload Your Dataset"),
    #                 p("You can upload your dataset using this interface.",
    #                   "Accepted data types include microarray, bulk RNA-Seq, single cell RNA-Seq, and mass spec.",
    #                   "Please refer to the dataset upload guidlines in order to properly format the necessary files.",
    #                   "A complete submission will include all files referenced in the guidelines.",
    #                   "We may reach out if we need additional information."),
    #                 wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 400px;",
    #                           column(6,
    #                                  textAreaInput(inputId = "dataset_upload_first_name",
    #                                                label = "First Name",
    #                                                width = "250px"),
    #                                  textAreaInput(inputId = "dataset_upload_last_name",
    #                                                label = "Last Name",
    #                                                width = "250px"),
    #                                  textAreaInput(inputId = "dataset_upload_email",
    #                                                label = "Email",
    #                                                width = "300px"),
    #                                  textAreaInput(inputId = "dataset_upload_data_type",
    #                                                label = "Data Type",
    #                                                width = "2500px"),
    #                                  fileInput(inputId = "dataset_upload_input", 
    #                                            label = "Upload a zip file with associated files for the given dataset", 
    #                                            multiple = FALSE,
    #                                            accept = c("text/csv", ".txt/.csv"),
    #                                            width = "100%"))
    #                           
    #          ))
    # )
    
)


# The active datasets table and selected row (defined globally)
active_datasets_table <- reactiveVal(NULL)

# Define the server functionality for the application
server <- function(input, output, session) {
    
    # Pull the datasets
    datasets <- get_datasets(db)
    
    # Display the table with all datasets
    output$dataset_table <- DT::renderDataTable(server = FALSE, {
      style_table(datasets)
    })

    
    #==============================================================================================================
    ### Observer for filtering the dataset table (updates whenever the user changes an option)
    #==============================================================================================================
    observe({
      input$disease_input
      input$tissue_input
      input$cell_type_input
      input$treatment_input
      input$dose_input
      
      # Update the table with all datasets
      output$dataset_table <- DT::renderDataTable(server = FALSE, {
        # Filter the datasets table for the selected options
        updated_table <- update_table(datasets, input)
        # Update the active datasets table in the global environment
        active_datasets_table(updated_table)
        
        style_table(updated_table)
      })
    })
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### Observer for launching a dataset specific application (when user clicks on a table row)
    #==============================================================================================================
    observeEvent(input$dataset_table_rows_selected, {
      # Pull the selected dataset ID and experiment type
      dataset_acc <- active_datasets_table()$dataset_acc[input$dataset_table_rows_selected]
      js$browseURL(stringr::str_interp("https://abbviegrc.shinyapps.io/ibdexplore_data/?dataset_acc=${dataset_acc}"))
      
      #shiny::tags$a(href = "https://www.example.com", "Visit Example Website")

    })
    #==============================================================================================================
    
}

# Close the connection to the DB
onStop(function() {
  dbDisconnect(db)
})

# Run the application 
shinyApp(ui = ui, server = server)
