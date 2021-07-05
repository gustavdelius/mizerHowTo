

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




#' shiny app to tweak gamma per species and see the effects on growth
#'
#' @export

shiny_gamma <- function(params, dat = NULL)
{

params_shiny <- params

ui=fluidPage(
  titlePanel("Gamma calibration"),
  fluidRow(
    column(4,
           wellPanel(
             actionButton("run", "Run the simulation")
           ),
           wellPanel(
             "Sprat's intitial value:",
             textOutput("SpratGamma"),
             numericInput("gamma1", "Search volume:", min = .1*params_shiny@species_params$gamma[1], max = 10*params_shiny@species_params$gamma[1],
                          value = params_shiny@species_params$gamma[1]),
             "Sandeel's intitial value:",
             textOutput("SandeelGamma"),
             numericInput("gamma2", "Search volume:", min = .1*params_shiny@species_params$gamma[2], max = 10*params_shiny@species_params$gamma[2],
                          value = params_shiny@species_params$gamma[2], step = 0.01*params_shiny@species_params$gamma[2]
             ),
             "N.pout's intitial value:",
             textOutput("N.poutGamma"),
             numericInput("gamma3", "Search volume:", min = .1*params_shiny@species_params$gamma[3], max = 10*params_shiny@species_params$gamma[3],
                          value = params_shiny@species_params$gamma[3], step = 0.01*params_shiny@species_params$gamma[3]
             ),
             "Dab's intitial value:",
             textOutput("DabGamma"),
             numericInput("gamma4", "Search volume:", min = .1*params_shiny@species_params$gamma[4], max = 10*params_shiny@species_params$gamma[4],
                          value = params_shiny@species_params$gamma[4], step = 0.01*params_shiny@species_params$gamma[4]
             ),
             "Herring's intitial value:",
             textOutput("HerringGamma"),
             numericInput("gamma5", "Search volume:", min = .1*params_shiny@species_params$gamma[5], max = 10*params_shiny@species_params$gamma[5],
                          value = params_shiny@species_params$gamma[5], step = 0.01*params_shiny@species_params$gamma[5]
             ),
             "Gurnard's intitial value:",
             textOutput("GurnardGamma"),
             numericInput("gamma6", "Search volume:", min = .1*params_shiny@species_params$gamma[6], max = 10*params_shiny@species_params$gamma[6],
                          value = params_shiny@species_params$gamma[6], step = 0.01*params_shiny@species_params$gamma[6]
             ),
             "Sole's intitial value:",
             textOutput("SoleGamma"),
             numericInput("gamma7", "Search volume:", min = .1*params_shiny@species_params$gamma[7], max = 10*params_shiny@species_params$gamma[7],
                          value = params_shiny@species_params$gamma[7], step = 0.01*params_shiny@species_params$gamma[7]
             ),
             "Whiting's intitial value:",
             textOutput("WhitingGamma"),
             numericInput("gamma8", "Search volume:", min = .1*params_shiny@species_params$gamma[8], max = 10*params_shiny@species_params$gamma[8],
                          value = params_shiny@species_params$gamma[8], step = 0.01*params_shiny@species_params$gamma[8]
             ),
             "Plaice's intitial value:",
             textOutput("PlaiceGamma"),
             numericInput("gamma9", "Search volume:", min = .1*params_shiny@species_params$gamma[9], max = 10*params_shiny@species_params$gamma[9],
                          value = params_shiny@species_params$gamma[9], step = 0.01*params_shiny@species_params$gamma[9]
             ),
             "Haddock's intitial value:",
             textOutput("HaddockGamma"),
             numericInput("gamma10", "Search volume:", min = .1*params_shiny@species_params$gamma[10], max = 10*params_shiny@species_params$gamma[10],
                          value = params_shiny@species_params$gamma[10], step = 0.01*params_shiny@species_params$gamma[10]
             ),
             "Saithe's intitial value:",
             textOutput("SaitheGamma"),
             numericInput("gamma11", "Search volume:", min = .1*params_shiny@species_params$gamma[11], max = 10*params_shiny@species_params$gamma[11],
                          value = params_shiny@species_params$gamma[11], step = 0.01*params_shiny@species_params$gamma[11]
             ),
             "Cod's intitial value:",
             textOutput("CodGamma"),
             numericInput("gamma12", "Search volume:", min = .1*params_shiny@species_params$gamma[12], max = 10*params_shiny@species_params$gamma[12],
                          value = params_shiny@species_params$gamma[12], step = 0.01*params_shiny@species_params$gamma[12]
             ),
           )
    ),
    column(6,
           plotOutput("plot1", width = 600, height = 600),
           plotOutput("plot2", width = 600, height = 600),
           plotOutput("plot3", width = 600, height = 600),
           if(!is.null(dat)) plotOutput("plot4", width = 600, height = 600)
    )
  )

)

server = function(input, output) {
  output$SpratGamma <- renderText({
    params_shiny@species_params$gamma[1]
  })
  output$SandeelGamma <- renderText({
    params_shiny@species_params$gamma[2]
  })
  output$N.poutGamma <- renderText({
    params_shiny@species_params$gamma[3]
  })
  output$DabGamma <- renderText({
    params_shiny@species_params$gamma[4]
  })
  output$HerringGamma <- renderText({
    params_shiny@species_params$gamma[5]
  })
  output$GurnardGamma <- renderText({
    params_shiny@species_params$gamma[6]
  })
  output$SoleGamma <- renderText({
    params_shiny@species_params$gamma[7]
  })
  output$WhitingGamma <- renderText({
    params_shiny@species_params$gamma[8]
  })
  output$PlaiceGamma <- renderText({
    params_shiny@species_params$gamma[9]
  })
  output$HaddockGamma <- renderText({
    params_shiny@species_params$gamma[10]
  })
  output$SaitheGamma <- renderText({
    params_shiny@species_params$gamma[11]
  })
  output$CodGamma <- renderText({
    params_shiny@species_params$gamma[12]
  })
  # reactive expression
  sim <- eventReactive(input$run, {
    params_shiny@species_params$gamma <- c(input$gamma1,input$gamma2,input$gamma3,input$gamma4,input$gamma5,input$gamma6,
                                           input$gamma7,input$gamma8,input$gamma9,input$gamma10,input$gamma11,input$gamma12)
    params_shiny <- setParams(params_shiny)
    # print(params_shiny@species_params$gamma) # check if everything is going well
    sim_shiny <- project(params_shiny, effort = 1, t_max = 100)
  })

  output$plot1 <- renderPlot({
    plotGrowthCurves2(sim(), species_panel = T)
  })

  output$plot2 <- renderPlot({
    plotFeedingLevel2(sim(), include_critical = T)
  })

  output$plot3 <- renderPlot({
    plotCalibration(sim())
  })
  if(!is.null(dat))
    {output$plot4 <- renderPlot({
    plotPredObsYield(sim(),dat)
  })
  }
}

shinyApp(ui, server)
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
