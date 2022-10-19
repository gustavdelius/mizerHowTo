

#' shiny app to tweak erepro
#'
#' @export
shiny_erepro <- function(input, dat = NULL) {
  params_shiny <- input
  ui=fluidPage(

    # Application title
    titlePanel("erepro calibration"),

    fluidRow(
      column(4, wellPanel(
        #sliderInput("kappa", "log10 Resource Carrying Capacity:", min = 8, max = 12, value = log10(params_shiny@resource_params$kappa),
        #            step = 0.1),
        #   sliderInput("Rmax", "log10 Maximum Recruitment:", min = 1, max = 12, value = 12,
        #              step = 0.1),
        sliderInput("erepro", "log10 Reproductive Efficiency:", min = -8, max = 1, value = log10(params_shiny@species_params$erepro[1]),
                    step = 0.1)
      )),
      column(6,
             plotOutput("plot1", width = 600, height = 600),
             if(!is.null(dat)) plotOutput("plot2", width = 600, height = 600)
      ))



  )
  server = function(input, output) {

#    sim <- observeEvent(input$erepro, {
#       print("yo")
#       params_shiny@species_params$erepro <- rep(10^input$erepro,12)
#       params_shiny <- setParams(params_shiny)
#       sim_shiny <- project(params_shiny, effort = 1, t_max = 50)
# })
#     output$plot1 <- renderPlot({
#       plot(sim())
#     })

    output$plot1 <- renderPlot({
      # set up params using values given, need check and change parameter values so units work in days units
      params_shiny@species_params$erepro <- rep(10^input$erepro,length(params_shiny@species_params$erepro))
      # params@species_params$Rmax <- rep(10^input$Rmax,12)
      params_shiny <- setParams(params_shiny)#,kappa=10^input$kappa)
      # run without fishing
      sim_shiny <- project(params_shiny, effort = 1, t_max = 50)
      plot(sim_shiny)
    })

    if(!is.null(dat))
    {
    output$plot2 <- renderPlot({
      # set up params using values given, need check and change parameter values so units work in days units
      params_shiny@species_params$erepro <- rep(10^input$erepro,length(params_shiny@species_params$erepro))
      # params@species_params$Rmax <- rep(10^input$Rmax,12)
      params_shiny <- setParams(params_shiny)#,kappa=10^input$kappa)
      # run without fishing
      sim_shiny <- project(params_shiny, effort = 1, t_max = 100)
      plotPredObsYield(sim_shiny,dat)
      # plotBiomass(sim_shiny)
    })
    }
  }
  shinyApp(ui, server)
}



#' shiny app to tweak erepro per species and see the effects on Fmsy
#'
#' @export
shiny_fmsy <- function(params,dat) {


params_shiny <- params
FmsyDat <- dat

ui=fluidPage(

  # Application title
  titlePanel("North Sea Fmsy"),

  fluidRow(
    column(4, wellPanel(
      # sliderInput("kappa", "log10 Resource Carrying Capacity:", min = 8, max = 12, value = log10(params_optim@resource_params$kappa),
      #             step = 0.1),
      #   sliderInput("Rmax", "log10 Maximum Recruitment:", min = 1, max = 12, value = 12,
      #              step = 0.1),
      sliderInput("erepro", "log10 Reproductive Efficiency:", min = -8, max = 1, value = log10(params_shiny@species_params$erepro[1]),
                  step = 0.1),
      sliderTextInput(
        inputId = "species",
        label = "Species name",
        choices = params_shiny@species_params$species,
        selected = params_shiny@species_params$species[1],
        grid = T
      ),
    )),

    column(6,
           plotOutput("distPlot", width = 600, height = 600)
    ))



)

server = function(input, output) {

  output$distPlot <- renderPlot({
    # set up params using values given, need check and change parameter values so units work in days units
    params_shiny@species_params$erepro[which(params_shiny@species_params$species == input$species)] <-10^(input$erepro)

    plotFmsy(params_shiny, speciesData = list(input$species,FmsyDat), effortRes = 10)
  })
}
# add species slider to display jsut one species at a time

shinyApp(ui = ui, server = server)

}



#' @title Gamma calibration app

#' @description Shiny app which focus on tweaking gamma per species and see the
#' effects on growth, feeding level, diet and shape of the size spectrum / biomass
#' through time
#'
#' @param params An object of class \linkS4class{MizerParams}.
#' @param dat A dataframe containing the yield value, to pass to `plotPredObsYield`.
#' Default is null (plot won't be displayed).
#'
#' @export

shiny_gamma <- function(params, dat = NULL)
{

    params_shiny <- params

    ui=fluidPage(
        titlePanel("Gamma calibration"),
        fluidRow(
            column(2,
                   wellPanel(
                       actionButton("run", "Run the simulation")
                   ),
                   wellPanel(
                       actionButton("done", "Return", icon = icon("check"),
                                    onclick = "setTimeout(function(){window.close();},500);")
                   ),
                   wellPanel(
                       uiOutput("toggle_gamma"))
            )
            ,
            column(5,
                   plotOutput("plotGrowth", height = 600),
                   plotOutput("plotFeed", height = 600),
                   if(!is.null(dat)) plotOutput("plotYield", height = 600)
            ),
            column(5,
                   plotOutput("plotDiet", height = 600),
                   plotOutput("plotCal", height = 600)
            )
        )
    )

    server <- function(input, output) {
        output$toggle_gamma <- renderUI({
            n_spec <- length(params_shiny@species_params$species)
            lapply(1:n_spec, function(i) {
                div(
                    paste(params_shiny@species_params$species[i],
                          'Initial Value:',
                          params_shiny@species_params$gamma[i]),
                    numericInput(inputId = paste0('gamma', i),
                                 "Search volume:",
                                 value = params_shiny@species_params$gamma[i])
                )
            })
        })
        # reactive expression
        sim <- eventReactive(input$run, {
            n_spec <- length(params_shiny@species_params$species)
            new_gamma <- c(sapply(1:n_spec, function(i) {
                inputName <- paste("gamma", i, sep = "")
                input[[inputName]]
            }))

            shiny_gamma_output <<- new_gamma

            params_shiny@species_params$gamma <- new_gamma
            params_shiny <- setParams(params_shiny)
            sim_shiny <- project(params_shiny, effort = 1, t_max = 100)
        })

        output$plotGrowth <- renderPlot({
            plotGrowthCurves2(sim(), species_panel = T)
        })

        output$plotFeed <- renderPlot({
            plotFeedingLevel(sim(), include_critical = T)
        })

        output$plotCal <- renderPlot({
            plotCalibration(sim())
        })

        if(!is.null(dat))
        {output$plotYield <- renderPlot({
            plotPredObsYield(sim(),dat)
        })
        }

        output$plotDiet <- renderPlot({
            plotDiet2(sim())
        })

        observeEvent(input$done, {
            stopApp(sim()@params@species_params$gamma)
        })

    }

    runGadget(ui, server, viewer = browserViewer())
}


#' shiny app to tweak gamma per species and see the effects on growth
#'
#' @export
shiny_kappa <- function(params, dat = NULL)
{
params_shiny <- params

ui=fluidPage(

  # Application title
  titlePanel("kappa calibration"),

  fluidRow(
    column(4, wellPanel(
      sliderInput("kappa", "log10 Resource Carrying Capacity:", min = 8, max = 12, value = log10(params_shiny@resource_params$kappa),
                  step = 0.1),
    )),
    column(6,
           plotOutput("plot1", width = 600, height = 600),
           if(!is.null(dat)) plotOutput("plot2", width = 600, height = 600)
    ))



)

server = function(input, output) {

  output$plot1 <- renderPlot({

    params_shiny <- setParams(params_shiny,kappa=10^input$kappa)
    sim_shiny <- project(params_shiny, effort = 1, t_max = 100)
    plotGrowthCurves2(sim_shiny, species_panel = T)
  })


  if(!is.null(dat))
  {
    output$plot2 <- renderPlot({
      params_shiny <- setParams(params_shiny,kappa=10^input$kappa)
      sim_shiny <- project(params_shiny, effort = 1, t_max = 100)
      plotPredObsYield(sim_shiny,dat)
    })
  }
}


shinyApp(ui, server)
}
