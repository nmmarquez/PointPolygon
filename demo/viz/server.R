.libPaths(c("~/R3.5/", .libPaths()))
library(shiny)
library(shinydashboard)
library(ggplot2)
library(sp)
library(dplyr)
library(PointPolygon)
library(readr)
library(DT)

dt1 <- read_csv("~/Data/utaziResults/aggRes.csv")
dt2 <- read_csv("~/Data/utaziResults/results.csv")

shinyServer(function(input,output){
    output$mest <- renderPlot({
        fn <- modelname <- paste0(
            "~/Data/utaziTest/",
            "range=", input$range,
            ",cov=", input$cov,
            ",covtype=", input$ct,
            ",M=", input$M,
            ",seed=", input$seed, ".Rds"
        )
        modelRez <- readRDS(fn)
        ggFieldEst(modelRez$sim, modelRez$pred, sd=input$sd)
    })
    
    output$dt1 <- renderDT({
        dt1
    })
    
    output$dt2 <- renderDT({
        dt2
    })
  
})