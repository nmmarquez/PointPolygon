#.libPaths(c("~/R3.5/", .libPaths()))
library(shiny)
library(shinydashboard)
library(leaflet)
library(sp)
library(dplyr)
library(DT)

rWidthSamples <- c("3"=3, "5"=5, "10"=10)
samps <- paste0(
    rep(paste0("rwidth_", rWidthSamples), 3),
    rep(c(" poly", " mix", " ov"), each=length(rWidthSamples)))


header <- dashboardHeader(
  title = 'Point Polygon Estimates'
)

body <- dashboardBody(
    fluidRow(
      column(width=12,
             tabBox(id='tabvals', width=NULL,
                    tabPanel('Model Estimates', plotOutput('mest'), value=1),
                    tabPanel('Model Stats Summary', DTOutput('dt1'), value=2),
                    tabPanel('Model Stats', DTOutput('dt2'), value=3),
                    tabPanel('Description', htmlOutput('desc'), value=4)
             )
      )
    ),
    status="danger",
    tags$head(tags$style(HTML('
                              /* logo */
                              .skin-blue .main-header .logo {
                              background-color: #070B19;
                              }
                              
                              /* logo when hovered */
                              .skin-blue .main-header .logo:hover {
                              background-color: #070B19;
                              }
                              
                              /* navbar (rest of the header) */
                              .skin-blue .main-header .navbar {
                              background-color: #070B19;
                              }        
                              
                              /* main sidebar */
                              .skin-blue .main-sidebar {
                              background-color: #070B19;
                              }
                              
                              /* active selected tab in the sidebarmenu */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #ff0000;
                              }
                              
                              /* other links in the sidebarmenu */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: #00ff00;
                              color: #000000;
                              }
                              
                              /* other links in the sidebarmenu when hovered */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #DF0101;
                              }
                              /* toggle button when hovered  */                    
                              .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #DF0101;
                              }
                              /* Highlighted Tab Color*/
                              .nav-tabs-custom .nav-tabs li.active {
                              border-top-color: #DF0101;
                              }')))
)

paramDF <- expand.grid(
  rangeE = c(.3, .5, .7),
  covVal = c(2, .4, -.5, .2, -2),
  covType = c("random", "spatial", "cluster"),
  M = seq(50, 300, by=50),
  seed = 1:10)

sidebar <- dashboardSidebar(
    conditionalPanel(
        condition='input.tabvals==1',
        selectInput('sampling', 'Sampling', samps),
        selectInput('range', 'Range', c(.7, .5, .3)),
        selectInput('cov', 'Cov Value', c(2, .4, -.5, .2, -2)),
        selectInput('ct', 'Cov Type',  c("random", "spatial", "cluster")),
        selectInput('M', 'M', seq(50, 300, by=50)),
        selectInput('seed', 'Seed', 1:10),
        selectInput('sd', 'SD', c(FALSE, TRUE))
    )
)

dashboardPage(
  header,
  sidebar,
  body
)