# IBDTransDB Integrate Application
# An application used to integrate data from multiple data sources for target identification
# Produces a list of targets that can be ranked using several metrics


library(RSQLite)

library(DT)
library(kableExtra)
library(reactable)

library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyBS)
library(shinybusy)
library(shinyalert)
library(shinyjs)

library(tidyr)
library(dplyr)
library(stringr)

library(poolr)
library(stats)

library(parallel)


#==============================================================================================================
# Import functions used to generate table summaries, table outputs, expression plots, etc
#==============================================================================================================
source("./app_functions/database_query.R")
source("./app_functions/query_options.R")
source("./app_functions/query_summary.R")
source("./app_functions/data_selection.R")
source("./app_functions/target_ranking.R")
#==============================================================================================================


#==============================================================================================================
# Create the connection to the Immunoverse DB
#==============================================================================================================
db <- dbConnect(RSQLite::SQLite(), dbname = "./IBDTransDB.db")
#==============================================================================================================

                  
# Define the UI for the application
ui <- navbarPage(
    title = "IBDIntegrate",
    id = "navbar-tabset",
    includeCSS("www/style.css"),
    
    useShinyjs(),
    extendShinyjs(text = "shinyjs.browseURL = function(url) {window.open(url,'_blank');}", 
                  functions = "browseURL"),
    
    # Alert if the user does not specify any genes or diseases
    useShinyalert(),
    
    tags$head(
      tags$script(HTML('
        $(document).ready(function() {
          $(".section .header").click(function() {
            var section = $(this).parent();
            if (section.hasClass("expanded")) {
              section.removeClass("expanded");
              section.find(".content").slideUp();
              section.find(".expand-icon").removeClass("fa-chevron-down").addClass("fa-chevron-right");
            } else {
              $(".section.expanded").removeClass("expanded");
              $(".section.expanded .content").slideUp();
              $(".section.expanded .expand-icon").removeClass("fa-chevron-down").addClass("fa-chevron-right");
              section.addClass("expanded");
              section.find(".content").slideDown();
              section.find(".expand-icon").removeClass("fa-chevron-right").addClass("fa-chevron-down");
            }
          });
        });
      ')),
      tags$script(src = 'popover.js'),
      tags$script(HTML("
        $(document).ready(function(){
          $('body').popover({
            selector: '[data-toggle=\"popover\"]',
            html: true,
            trigger: 'hover',
            container: 'body'
          });
        });
        $(document).ready(function() {
          $('header').find('nav').append(\'<span class='development-header'> Developed by Immunology Computational Biology - GRC </span>\');
        })
      "))
    ),
    tags$script(HTML("var header = $('.navbar> .container-fluid');
                       header.append('<div style=\"float:right;color:white;padding-top:8px;\"><h4>Developed by Immunology Computational Biology - GRC</h3></div>');
                       console.log(header)")),
    
    ### DATA SELECTION TAB ###
    tabPanel("Data Selection",
             column(3,
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 1250px;",
                              
                              h4("REQUIRED", style = "text-align: center;"),
                              
                              div(
                                class = "section",
                                div(
                                  class = "header",
                                  icon("dna", class = "icon"),
                                  h4("TARGET INPUT")
                                ),
                                radioGroupButtons(
                                  inputId = "target_option",
                                  label = NULL,
                                  size = "sm",
                                  direction = "horizontal",
                                  choices = c("GENE" = "GENE-level", "PATHWAY" = "PATHWAY-level"),
                                  selected = "GENE-level",
                                  checkIcon = list(yes = tags$i(class = "fa fa-check-square",
                                                                style = "color: steelblue"),
                                                   no = tags$i(class = "fa fa-square-o", 
                                                               style = "color: steelblue"))),
                                uiOutput("target_input")
                              ),
                              br(),
                              
                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "disease_input",
                                  label = "Select diseases of interest",
                                  choices = get_keywords("Disease", db),
                                  multiple = TRUE
                                )
                              ),
                              br(),
                              
                              h4("OPTIONAL", style = "text-align: center;"),

                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "tissue_input",
                                  label = "Select tissues of interest",
                                  choices = get_keywords("Source", db),
                                  multiple = TRUE
                                )
                              ),
                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "cell_type_input",
                                  label = "Select cell types of interest",
                                  choices = get_keywords("CellType", db),
                                  multiple = TRUE
                                )
                              ),
                              br(),
                              
                              # div(
                              #   style = "margin-bottom: 5px;",
                              #   selectInput(
                              #     inputId = "platform_input",
                              #     label = "Select platforms of interest",
                              #     choices = get_keywords("Platform", db),
                              #     multiple = TRUE
                              #   )
                              # ),
                              # br(),
                              # 
                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "treatment_input",
                                  label = "Select treatments of interest",
                                  choices = get_keywords("Treatment", db),
                                  multiple = TRUE
                                )
                              ),
                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "timepoint_input",
                                  label = "Select timepoints of interest",
                                  choices = get_keywords("Timepoint", db),
                                  multiple = TRUE
                                )
                              ),
                              br(),
                              
                              div(class = "fixed-button",
                                  actionBttn(inputId = "confirm_selections",
                                             label = "confirm selections",
                                             style = "material-flat",
                                             color = "danger"))
                    )),
             
             column(9,
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 1250px; overflow-y: scroll; max-height: 1250px;",
                              add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
                              fluidRow(
                                uiOutput("summary_title"),
                                br(),
                                column(6,
                                       uiOutput("user_query_summary")),
                                column(6,
                                       uiOutput("data_summary"))
                              ),
                              br(), br(),
                              fluidRow(
                                uiOutput("selections_section"),
                                br(),
                                column(6, 
                                       div(style = "padding-left: 10px; padding-right: 10px",
                                           uiOutput("dataset_buttons"),
                                           DT::dataTableOutput("dataset_table"))),
                                column(6, 
                                       div(style = "padding-left: 10px; padding-right: 10px",
                                           uiOutput("comparison_buttons"),
                                           DT::dataTableOutput("comparison_table")))
                              ),
                              br(), br(),
                              div(class = "targets_button_wrapper", 
                                  style = "text-align: center", 
                                  uiOutput("generate_ranked_targets_button"))
                              
                    )
                    
            )
    ),
    
    ### RESULTS PANEL TAB ###
    tabPanel("Results Panel",
             column(10, offset = 1,
                    fluidRow(
                      wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                                uiOutput("ranked_targets_title"),
                                reactableOutput("ranked_targets_table"),
                                br(),
                                uiOutput("generate_download_button"))
                    ),
                    fluidRow(
                      column(6, style = "padding-left: 0px; padding-right: 10px;",
                             wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                                       uiOutput("immunoexplore_section"))),
                      column(6, style = "padding-left: 10px; padding-right: 0px;",
                             wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                                       uiOutput("immunocompare_section")))
                    )
             )
    )
    
    # ### HELP TAB ###
    # tabPanel("Help",
    #          column(10, offset = 1,
    #                 wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
    #                           add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
    #                           includeHTML("homepage.html")))
    # )
    
)

# Define the server functionality for the application
server <- function(input, output, session) {
    
    ### REACTIVE VALUES ###
    datasets <- reactiveVal(NULL)
    comparisons <- reactiveVal(NULL)
    comparison_data <- reactiveVal(NULL)
    comparison_table <- reactiveVal(NULL)
    targets <- reactiveVal(NULL)
    pathway_database <- reactiveVal(NULL)
    pathway_data <- reactiveVal(NULL)
    selected_dataset_acc_list <- reactiveVal(NULL)
    integration_table <- reactiveVal(NULL)
    gene_map <- reactiveVal(NULL)
    ##########
    
    
    #==============================================================================================================
    ### ObserveEvents for updating query options (updates whenever the user changes an option)
    #==============================================================================================================
    # observeEvent(input$disease_input, {
    #   change_query_options(session, "disease", input,
    #                        input$disease_input, input$organism_input,
    #                        input$tissue_input, input$treatment_input, 
    #                        input$experiment_type_input, input$cell_type_input, input$timepoint_input, 
    #                        db)
    # })
    # observeEvent(input$tissue_input, {
    #   change_query_options(session, "tissue", input, 
    #                        input$disease_input, input$organism_input,
    #                        input$tissue_input, input$treatment_input, 
    #                        input$experiment_type_input, input$cell_type_input, input$timepoint_input,
    #                        db)
    # })
    # observeEvent(input$cell_type_input, {
    #   change_query_options(session, "cell_type", input, 
    #                        input$disease_input, input$organism_input,
    #                        input$tissue_input, input$treatment_input, 
    #                        input$experiment_type_input, input$cell_type_input, input$timepoint_input,
    #                        db)
    # })
    # observeEvent(input$experiment_type_input, {
    #   change_query_options(session, "experiment_type", input, 
    #                        input$disease_input, input$organism_input,
    #                        input$tissue_input, input$treatment_input, 
    #                        input$experiment_type_input, input$cell_type_input, input$timepoint_input,
    #                        db)
    # })
    # observeEvent(input$treatment_input, {
    #   change_query_options(session, "treatment", input, 
    #                        input$disease_input, input$organism_input,
    #                        input$tissue_input, input$treatment_input, 
    #                        input$experiment_type_input, input$cell_type_input, input$timepoint_input, 
    #                        db)
    # })
    # 
    # observeEvent(input$timepoint_input, {
    #   change_query_options(session, "timepoint", input, 
    #                        input$disease_input, input$organism_input,
    #                        input$tissue_input, input$treatment_input, 
    #                        input$experiment_type_input, input$cell_type_input, input$timepoint_input, 
    #                        db)
    # })
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### Set initial target option and ObserveEvent for TARGET SWITCH
    #==============================================================================================================
    
    # # Set the initial target selection option (Option 1)
    # output$gene_selection <- renderUI({
    #   p("Perform integration analysis on ALL targets", align = "center")
    # })
    
    # # Update the input for targets when a different option is selected
    # observeEvent(input$gene_switch, {
    #   if (input$gene_switch == "Option 1") {
    #     output$gene_selection <- renderUI({
    #       p("Perform integration for ALL targets", align = "center", style = "color: white;")
    #     })
    #   } else if (input$gene_switch == "Option 2") {
    #     output$gene_selection <- renderUI({
    #       div(
    #         p("Perform integration for predefined list(s) of targets", align = "center", style = "color: white;"),
    #         pickerInput(inputId = "target_lists",
    #                     label = NULL, 
    #                     choices = get_predefined_groups(db),
    #                     selected = "all",
    #                     multiple = TRUE,
    #                     width = "275px"))
    #     })
    #   } else {
    #     output$gene_selection <- renderUI({
    #       div(
    #         p("Perform integration for an uploaded list of targets", align = "center", style = "color: white;"),
    #         fileInput(inputId = "target_file", 
    #                   label = NULL, 
    #                   multiple = FALSE, 
    #                   accept = c("text/csv", ".txt/.csv"),
    #                   width = "auto"))
    #     })
    #   }
    # })
    
    # Set the initial target input (GENE-level)
    output$target_input <- renderUI({
      div(
        radioGroupButtons(
          inputId = "gene_switch",
          label = NULL,
          size = "sm",
          choices = c("Option 1", "Option 2"),
          selected = "Option 1",
          checkIcon = list(yes = tags$i(class = "fa fa-check-square",
                                        style = "color: steelblue"),
                           no = tags$i(class = "fa fa-square-o", 
                                       style = "color: steelblue"))
        ),
        uiOutput("gene_selection")
      )
    })
    # Set the initial target selection option (Option 1)
    output$gene_selection <- renderUI({
      div(
        p("Perform integration for predefined target list(s)", align = "center"),
        selectInput(inputId = "target_lists",
                    label = NULL,
                    choices = get_predefined_groups(db),
                    selected = NULL,
                    multiple = TRUE)
      )
    })
    
    # Update the input for targets when switching between GENE-level, SIGNATURE-level, and PATHWAY-level
    observeEvent(input$target_option, {
      if (input$target_option == "GENE-level") {
        output$target_input <- renderUI({
          div(
            radioGroupButtons(
              inputId = "gene_switch",
              label = NULL,
              size = "sm",
              choices = c("Option 1", "Option 2"),
              selected = "Option 1",
              checkIcon = list(yes = tags$i(class = "fa fa-check-square",
                                            style = "color: steelblue"),
                               no = tags$i(class = "fa fa-square-o", 
                                           style = "color: steelblue"))
            ),
            uiOutput("gene_selection")
          )
        })
      } else if (input$target_option == "SIGNATURE-level") {
        output$target_input <- renderUI({
          div(
            radioGroupButtons(
              inputId = "signature_switch",
              label = NULL,
              size = "sm",
              choices = c("Option 1", "Option 2"),
              selected = "Option 1",
              checkIcon = list(yes = tags$i(class = "fa fa-check-square",
                                            style = "color: steelblue"),
                               no = tags$i(class = "fa fa-square-o", 
                                           style = "color: steelblue"))
            ),
            uiOutput("signature_selection")
          )
        })
      } else {
        # Define the pathway choices
        pathways <- list("pathway" = c("KEGG" = "pathway_KEGG", "Reactome" = "pathway_Reactome"))
        output$target_input <- renderUI({
          div(
            p("Perform integration for pathways", align = "center"),
            selectInput(inputId = "pathway_database_selection",
                        label = "Select a pathway database",
                        choices = pathways,
                        selected = NULL,
                        multiple = FALSE,
                        width = "auto"))
        })
      }
    })
    
    # Update the input for targets when a different option is selected
    observeEvent(input$gene_switch, {
      if (input$gene_switch == "Option 1") {
        output$gene_selection <- renderUI({
          div(
            p("Perform integration for predefined target list(s)", align = "center"),
            selectInput(inputId = "target_lists",
                        label = NULL,
                        choices = get_predefined_groups(db),
                        selected = NULL,
                        multiple = TRUE)
          )
        })
      } else {
        output$gene_selection <- renderUI({
          div(
            p("Perform integration for an uploaded list of targets", align = "center"),
            fileInput(inputId = "target_file", 
                      label = NULL, 
                      multiple = FALSE, 
                      accept = c("text/csv", ".txt/.csv"),
                      width = "auto"))
        })
      }
    })
    
    # Update the input for signatures when a different option is selected
    observeEvent(input$signature_switch, {
      if (input$signature_switch == "Option 1") {
        output$signature_selection <- renderUI({
          div(
            p("Integration for a single custom signature", align = "center"),
            p("(each gene must be entered on a separate line)", align = "center"),
            textAreaInput(inputId = "signature_paste",
                          label = NULL,
                          placeholder = "paste genes here")
          )
        })
      } else {
        output$signature_selection <- renderUI({
          div(
            p("Upload a GMT file for multiple signatures", align = "center"),
            p("(only gene names are accepted)", align = "center"),
            fileInput(inputId = "signature_file", 
                      label = NULL, 
                      multiple = FALSE, 
                      accept = c("text/csv", ".txt/.csv"),
                      width = "auto")
          )
        })
      }
    })
    #==============================================================================================================

    
    #==============================================================================================================
    ### Confirm selections and produce output when the user clicks the "Confirm Selections" button
    #==============================================================================================================
    observeEvent(input$confirm_selections, {
      
       
        # Generate an alert if the user does not specify any genes or diseases
        if (!isTruthy(input$disease_input)) {
          shinyalert(title = "Insufficient Input Supplied",
                     text = "You need to specify at least one disease",
                     type = "error")
        } else {
          
          # Define the user query
          diseases <- unlist(input$disease_input)
          organisms <- unlist(input$organism_input)
          treatments <- unlist(input$treatment_input)
          sources <- unlist(input$tissue_input)
          cell_types <- unlist(input$cell_type_input)
          experiment_types <- unlist(input$experiment_type_input)
          timepoints <- unlist(input$timepoint_input)
          
          # Pull the dataset descriptions, comparison descriptions, and gene set descriptions
          datasets <- get_datasets(diseases, organisms, sources, treatments, experiment_types, platforms, cell_types, timepoints, db)
          comparisons <- get_comparisons(datasets, db)
          # Update the reactive values
          datasets(datasets)
          comparisons(comparisons)
          
          # Pull the gene map
          gene_map <- get_gene_map(db)
          # Update the reactive value
          gene_map(gene_map)
          
          isError <- FALSE
          
          # Determine the targets to rank based on the user option specified
          if (input$target_option == "GENE-level") {
            if (input$gene_switch == "Option 1") {
              ### THIS WILL BE CHANGED - INCLUDE SETS OF GENES IN DB ###
              targets <- get_predefined_group_genes(input$target_lists, db)
              
            } else {
              # Pull the targets from the uploaded target list
              targets <- read.csv(input$target_file$datapath, header = FALSE)
              targets <- toupper(unique(targets$V1))
            }
            
            targets <- intersect(targets,gene_map()$gene)
            
            if(length(targets)==0){
              shinyalert(title = "Invalid gene symbols",
                         text = "Please input valid gene symbols!",
                         type = "error")
              isError <- TRUE
            }else{
              # Update the reactive value
              targets(targets)
              # Pull the comparison data
              comparison_data <- get_comparison_data(comparisons, targets, gene_map, db)
              # Update the reactive value
              comparison_data(comparison_data)
            }
          } else if (input$target_option == "SIGNATURE-level") {
            if (input$signature_switch == "Option 1") {
              # Pull the targets for the given signature
              targets <- toupper(str_trim(unlist(strsplit(input$signature_paste, "\n"))))
              # Update the reactive value
              targets(targets)
              # Pull the comparison data
              comparison_data <- get_comparison_data(comparisons, targets, gene_map, db)
              # Update the reactive value
              comparison_data(comparison_data)
            } else {
              
            }
          } else {
            # Pull the selected pathway databases and pathways
            # Update the reactive value
            pathway_database <- input$pathway_database_selection
            pathway_database(pathway_database)
            
            # Pull the pathway data
            pathway_data <- get_pathway_data(comparisons, pathway_database, db)
            # Update the reactive value
            pathway_data(pathway_data)
          }
          
          if(isError==FALSE){
            
            # Jump to the Comparison Summary tab
            updateTabsetPanel(session, "main_body_tabset", selected = "Summary and Selections")
            
            
            #==============================================================================================================
            ### SUMMARY SECTION
            #==============================================================================================================
            
            # Generate the title of the summary section
            output$summary_title <- renderUI({
              h2("Summary of Matching Data", align = "center")
            })
            
            # Generate the table with summaries of datasets, comparisons, and genesets matching the user query
            output$user_query_summary <- renderUI({
              tags$div(
                tags$p("The following output meets the user specified critieria listed below."),
                tags$ul(
                  tags$li(paste0("diseases - ", paste0(diseases, collapse = ", "))),
                  tags$li(paste0("sources - ", paste0(sources, collapse = ", "))),
                  tags$li(paste0("celltypes - ", paste0(experiment_types, collapse = ", "))),
                  tags$li(paste0("treatments - ", paste0(treatments, collapse = ", "))),
                  tags$li(paste0("timepoints - ", paste0(timepoints, collapse = ", ")))
                )
              )
            })
            # Generate the table with summaries of datasets, comparisons, and genesets matching the user query
            output$data_summary <- renderText({
              generate_summary_table(datasets, comparisons, diseases, experiment_types, db)
            })
            #==============================================================================================================
            
            
            #==============================================================================================================
            ### DATA SELECTIONS
            #==============================================================================================================
            
            # Generate the title and description of the data selection section
            output$selections_section <- renderUI({
              div(
                h2("Data Selection", align = "center"),
                br(), br(),
                div(style = "border: solid 1px black; border-radius: 5px; padding: 10px;",
                    p("Guidlines", style = "font-weight: bold"),
                    p("Choose the datasets that you would like to include in the integration analysis.",
                      "Selected rows in the table correspond to datasets that will be included.",
                      "By default all datasets that match the user specified options are selected.", 
                      "Deselecting datasets will filter the available options for comparisons.",
                      "A given comparison (e.g. Responder vs NonResponder) may be present across multiple datasets.",
                      "Selecting a row from the comparison table will include that type of comparison",
                      "across all selected datasets that have the given comparison."),
                    p("Note", style = "font-weight: bold"),
                    p("We recommend selecting one or a few comparisons that pertain to your specific question as opposed to selecting all."))
              )
            })
            
            # Display the select all button and deselect all button for the datasets
            output$dataset_buttons <- renderUI({
              div(
                h3("Datasets", align = "center"),
                div(style = "display: inline-block; vertical-align: top;",
                    actionButton("select_all_datasets",
                                 label = "select all",
                                 width = "100px"),
                    actionButton("deselect_all_datasets",
                                 label = "deselect all",
                                 width = "100px")
                )
              )
            })
            # Update reactive value
            selected_dataset_acc_list(datasets$dataset_acc)
            #a <- generate_dataset_table(datasets, TRUE)
            # Generate the dataset table
            output$dataset_table <- renderDataTable(server = FALSE, {
              generate_dataset_table(datasets, TRUE)
            })
            
            # Display the select all button and deselect all button for the comparisons
            output$comparison_buttons <- renderUI({
              div(
                h3("Comparisons", align = "center"),
                div(style = "display: inline-block; vertical-align: top;",
                    actionButton("select_all_comparisons",
                                 label = "select all",
                                 width = "100px"),
                    actionButton("deselect_all_comparisons",
                                 label = "deselect all",
                                 width = "100px")
                )
              )
            })
            # Generate the comparison table
            comparison_table <- generate_comparison_table(datasets$dataset_acc, comparisons)
            # Update the reactive value
            comparison_table(comparison_table)
            # Display the styled comparison table (for selecting comparisons)
            output$comparison_table <- renderDataTable(server = FALSE, {
              style_comparison_table(comparison_table, FALSE)
            })
            #==============================================================================================================
            
            
            #==============================================================================================================
            ### RUN IMMUNOINTEGRATE BUTTON
            #==============================================================================================================
            
            # Display the button used to generate the ranked targets
            output$generate_ranked_targets_button <- renderUI({
              actionButton("run_immunointegrate", 
                           label = "Run Integration",
                           style = "width: 250px; font-size: 18px") 
            })
            #==============================================================================================================
            
          }
        }
    })
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### ObserveEvent: filtering comparisons on selected datasets and selecting/deselecting all datasets/comparisons
    #==============================================================================================================
    # DATASETS TABLE SELECTION
    observeEvent(input$dataset_table_rows_selected, {
      # Pull the selected dataset IDs
      selected_dataset_acc_list <- datasets() %>%
        slice(input$dataset_table_rows_selected) %>%
        pull(dataset_acc)
      # Update reactive value
      selected_dataset_acc_list(selected_dataset_acc_list)
      
      # Generate the comparison table
      comparison_table <- generate_comparison_table(selected_dataset_acc_list, comparisons())
      # Update the reactive value
      comparison_table(comparison_table)
      # Display the styled comparison table (for selecting comparisons)
      output$comparison_table <- renderDataTable(server = FALSE, {
        style_comparison_table(comparison_table, FALSE)
      })
    })
    
    # SELECT ALL DATASETS
    observeEvent(input$select_all_datasets, {
      # Display the styled dataset table (for selecting datasets)
      output$dataset_table <- renderDataTable(server = FALSE, {
        generate_dataset_table(datasets(), TRUE)
      })
    })
    
    # DESELECT ALL DATASETS
    observeEvent(input$deselect_all_datasets, {
      # Display the styled dataset table (for selecting datasets)
      output$dataset_table <- renderDataTable(server = FALSE, {
        generate_dataset_table(datasets(), FALSE)
      })
    })
    
    # SELECT ALL COMPARISONS
    observeEvent(input$select_all_comparisons, {
      # Generate the comparison table
      comparison_table <- generate_comparison_table(selected_dataset_acc_list(), comparisons())
      # Update the reactive value
      comparison_table(comparison_table)
      # Display the styled comparison table (for selecting comparisons)
      output$comparison_table <- renderDataTable(server = FALSE, {
        style_comparison_table(comparison_table, TRUE)
      })
    })
    
    # DESELECT ALL COMPARISONS
    observeEvent(input$deselect_all_comparisons, {
      # Generate the comparison table
      comparison_table <- generate_comparison_table(selected_dataset_acc_list(), comparisons())
      # Update the reactive value
      comparison_table(comparison_table)
      # Display the styled comparison table (for selecting comparisons)
      output$comparison_table <- renderDataTable(server = FALSE, {
        style_comparison_table(comparison_table, FALSE)
      })
    })
    #==============================================================================================================

    
    #==============================================================================================================
    ### ObserveEvent for RUN INTEGRATION
    #==============================================================================================================
    observeEvent(input$run_immunointegrate, {

      # Generate an alert if the user does not select any comparisons
      if (length(input$comparison_table_rows_selected) == 0 || length(input$geneset_table_rows_selected == 0)) {
        shinyalert(title = "Insufficient Input Supplied",
                   text = "You need to select at least one comparison or at least one gene set",
                   type = "error")
      } else {
        # Jump to the Integration Results tab
        updateTabsetPanel(session, "navbar-tabset", selected = "Results Panel")
        
        # Determine the comparisons that are selected
        selected_comparisons <- comparison_table()[input$comparison_table_rows_selected,]
        selected_comparisons <- pull(selected_comparisons, full_comparison)
        selected_comparison_ids <- format_comparison(comparisons()) %>%
          filter(full_comparison %in% selected_comparisons, dataset_acc %in% selected_dataset_acc_list()) %>%
          pull(id)

        
        #==========
        if (input$target_option == "GENE-level") {
          # Pull the comparison data based on the selected comparisons
          selected_comparison_data <- filter(comparison_data(), comparison_id %in% selected_comparison_ids)
          # Dataframe used to store the ranked targets and associated metrics
          rank_df <- rank_targets(targets(), comparisons(), selected_comparison_data)
          targets(rank_df$target) ##reorder the targets based on rank_df
          extra_info <- generate_target_extra_info(targets(), comparisons(), selected_comparison_data)
          
        } else if (input$target_option == "SIGNATURE-level") {
          # Pull the comparison data based on the selected comparisons
          selected_comparison_data <- filter(comparison_data(), comparison_id %in% selected_comparison_ids)
          # Dataframe used to store the ranked signatures and associated metrics
          rank_df <- rank_signature(targets(), selected_comparison_ids, comparisons(), db)
        } else {
          # Pull the pathway data based on the selected comparisons
          selected_pathway_data <- filter(pathway_data(), comparison_id %in% selected_comparison_ids)
          # Dataframe used to store the ranked pathways and associated metrics
          rank_df <- rank_pathways(comparisons(),selected_pathway_data)
          extra_info <- generate_pathway_extra_info(rank_df,comparisons(), selected_pathway_data)
        }
        # Update the reactive value
        integration_table(rank_df)
        showNotification(ui = "Please Note: Table loading may take few seconds to fully load",type = "message")
        # Display the title for the results section
        output$ranked_targets_title <- renderUI({
          h3("Integration Results", align = "center")
        })
        
        # Display the content for the ranked targets table
        output$ranked_targets_table <- renderReactable({
          style_rank_df(rank_df, extra_info)
        })
        
        # Display the download button for the table
        output$generate_download_button <- renderUI({
          downloadButton("download_button",
                         label = "Download Table",
                         style = "color: black; background-color: #DCDCDC; border-color: #505050; border-width: 1px;")
        })
        #==========
        
        
        #==========
        # Display the "Run ImmunoExplore" section
        output$immunoexplore_section <- renderUI({
          div(style = "display: flex; justify-content: space-between;",
              h3("Connect to IBDExplore"),
              actionButton("run_immunoexplore", 
                           label = "IBDExplore",
                           style = "width: 250px; font-size: 18px; text-align; center")
            # p("Clicking the", span("Run ImmunoExplore", style = "font-weight: bold;"), "button below will bring you to",
            #   "the ImmunoExplore main page where you can further interrogate a dataset of interest.",
            #   "By default, datasets shown in the dataset selection tab will be those included in this integration analysis."),
            # br(),
          )
        })
        #==========
        
        
        #==========
        # Display the "Run ImmunoCompare" section
        output$immunocompare_section <- renderUI({
          div(style = "display: flex; justify-content: space-between;",
              h3("Connect to IBDCompare", align = "center"),
              actionButton("run_immunocompare", 
                           label = "IBDCompare",
                           style = "width: 250px; font-size: 18px; text-align; center")
            # p("You can select targets of interest from the drop down menu.",
            #   "Clicking the", span("Run ImmunoCompare", style = "font-weight: bold;"), "button will navigate you to ImmunoCompare",
            #   "where you can futher interrogate these targets across the data sources used in the integration analysis.",
            #   "Limit your selection to 8 targets."),
            # div(style = "display: flex; justify-content: space-between;",
            #     selectizeInput(inputId = "volcano_comparison_input",
            #                    label = "Select targets of interest",
            #                    choices = integration_table()$target,
            #                    multiple = TRUE,
            #                    width = "250px",
            #                    options = list(maxItems = 8)),
            #     actionButton("run_immunocompare", 
            #                  label = "ImmunoCompare",
            #                  style = "width: 250px; font-size: 18px; text-align; center"))
          )
        })
        #==========
      }
    })
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### DownloadHandler for download button
    #==============================================================================================================
    output$download_button <- downloadHandler(
      filename = function() { 
        "integration_results.csv"
      },
      content = function(file) {
        write.csv(integration_table(), file, row.names = FALSE)
      }
    )
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### ObserveEvent for RUN IMMUNOEXPLORE
    #==============================================================================================================
    observeEvent(input$run_immunoexplore, {
      js$browseURL(stringr::str_interp("https://abbviegrc1.shinyapps.io/ibdexplore/"))
    })
    #==============================================================================================================
    
    #==============================================================================================================
    ### ObserveEvent for RUN IMMUNOCOMPARE
    #==============================================================================================================
    observeEvent(input$run_immunocompare, {
      js$browseURL(stringr::str_interp("https://abbviegrc2.shinyapps.io/ibdcompare/"))
    })
    #==============================================================================================================
    
}

# Close the connection to the DB
onStop(function() {
  dbDisconnect(db)
})

# Run the application 
shinyApp(ui = ui, server = server)
