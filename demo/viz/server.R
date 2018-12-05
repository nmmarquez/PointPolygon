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
fileName <- '~/Documents/PointPolygon/demo/ppExplain.html'
descHTML <- readChar(fileName, file.info(fileName)$size)

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
        predList <- list(
            riemann = modelRez$pred$riemann[[input$sampling]],
            utazi = modelRez$pred$utazi[[input$sampling]],
            resample = modelRez$pred$resample[[input$sampling]],
            point = modelRez$pred$point$point
        )
        ggFieldEst(modelRez$sim, predList, sd=input$sd) +
            labs(fill=ifelse(input$sd, "SD", "Probability"))
    })
    
    output$dt1 <- renderDT({
        dt1
    })
    
    output$dt2 <- renderDT({
        dt2
    })
    
    output$desc <- renderUI({
        withMathJax(HTML(descHTML))
    })
  
})