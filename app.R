library(shiny)
source("utmToPlanas.R")

# Interfaz de usuario
ui <- fluidPage(
  titlePanel("Aplicación Independiente: UTM a Planas"),
  UTMPlanasModuleUi("utmModule")
)

# Servidor
server <- function(input, output, session) {
  callModule(UTMPlanasModuleServer, "utmModule")
}

# Lanzar la aplicación
shinyApp(ui, server)