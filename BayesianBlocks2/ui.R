
require(shiny)
require(shinydashboard)
require(leaflet)

header <- dashboardHeader(title = "Bayesian Blocks (Demo)",titleWidth = 450
)

body <- dashboardBody(
  fluidRow(
    column(width = 3,
           box(title="Block Summary",width=NULL,status = "warning",
               textOutput("constants"),
               textOutput("linears"),
               textOutput("pruned"),
               tags$head(tags$style("#constants{color: black;font-size: 18px;text-align: right;}")),
               tags$head(tags$style("#linears{color: black;font-size: 18px;text-align: right;}")),
               tags$head(tags$style("#pruned{color: black;font-size: 18px;text-align: right;}"))
               ),
           box(title="Algorithm Options",width=NULL,
               checkboxGroupInput(inputId="type", label="Desired block shapes:",
                            choices=list("Constant","Linear"),
                            selected=c("Constant")),
               radioButtons("c_choice", "Penalization type:",
                            choices = list("AIC","BIC","Likelihood Ratio Test"),
                            selected = c("Likelihood Ratio Test")),
               conditionalPanel(
                 condition = "input.c_choice == 'Likelihood Ratio Test' && (input.type !='Linear' && input.type !='Constant')",
                 sliderInput("alpha",
                             "Likelihood Ratio Test Rejection Alpha:",
                             min = 0,  max = 1, value = 0.05))),
           
           box(title="Plot Options", width = NULL,
               checkboxGroupInput("show","What to show:",
                                  choiceNames =
                                    list("Binned Data","Data Points","Bayesian Blocks"),
                                  choiceValues =
                                    list("hist","points","blocks"),
                                  selected =list("points","blocks")),
               checkboxInput("true",label="True Intensity",value=TRUE),
               conditionalPanel(
                 condition = "input.show.indexOf('hist') != -1",sliderInput("binwidth",
                           "Histogram Binwidth (Milliseconds):",
                           min = 10,  max = 1000, value = 300,step=10))
           )
    ),
    column(width = 9,
           box(width = NULL, solidHeader = TRUE,
               plotOutput("plot",width = "100%", height = "100%")
           ),
           box(width = NULL,
               sliderInput("c",
                           "Block penalization value:",
                           min = 0,  max = 12, value = 4,step=0.01)
           )
    )
  )
)

dashboardPage(
  header,
  dashboardSidebar(disable = TRUE),
  body
)