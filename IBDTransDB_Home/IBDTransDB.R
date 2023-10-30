# IBDTransDB Application
# A platform of several tools for multi-omics exploration, visualization, and integration


library(shiny)
library(shinyWidgets)
library(shinyBS)
library(shinybusy)
library(rintrojs)
library(shinyjs)


# Define the UI for the application
ui <- navbarPage(title = "IBDTransDB",
                 theme = "style.css",
                 fluid = TRUE, 
                 collapsible = TRUE,
                 tags$head(
                   tags$style(HTML(".navbar-brand {width: 200px;}
                              .navbar {background-color: white;}
                              .navbar-default .navbar-brand {color: #E0E0E0;}
                              .tab-panel {background-color: #F0F0F0; color: #0B0C42}
                              .navbar-default .navbar-nav > .active > a, 
                              .navbar-default .navbar-nav > .active > a:focus, 
                              .navbar-default .navbar-nav > .active > a:hover {
                                background-color: #E0E0E0;
                              }"))
                 ),
                 
                 ### HOME ###
                 tabPanel("Home",
                          includeHTML("homepage.html"),
                          tags$head(
                            tags$link(rel = "stylesheet",
                                      type = "text/css",
                                      href = "plugins/font-awesome-4.7.0/css/font-awesome.min.css")
                          )
                 )
                 
)


# Define the server functionality for the application
server <- function(input, output, session) {

  # Response form for users
  # formServer(form_info)
  
}

# Run the application 
shinyApp(ui = ui, server = server)
