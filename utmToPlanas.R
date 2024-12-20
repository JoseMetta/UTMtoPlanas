# modules/utmToPlanas.R

library(shiny)
library(sf)
library(geosphere)
library(dplyr)
library(DT)

# Función de transformación UTM a Planas
funcion_UTM_planas <- function(p, crs_input) {
  a <- 6378137.00
  b <- 6356752.314
  e2 <- 0.00669438
  et2 <- 0.006739497
  c <- 6399593.626
  k_0 <- 0.9996
  
  # Convertir a sf usando columnas E y N
  h <- st_as_sf(p, coords = c("ESTE", "NORTE"), crs = crs_input)
  q <- st_transform(h, crs = 4326)
  m <- st_coordinates(q)
  
  ## Función que calcula azimuth entre puntos
  st_azimuth = function(x, y) {
    # Checks
    stopifnot(all(st_is(x, "POINT")))
    stopifnot(all(st_is(y, "POINT")))
    
    # Extract geometry
    x = st_geometry(x)
    y = st_geometry(y)
    
    # Recycle 'x' or 'y' if necessary
    if (length(x) < length(y)) {
      ids = rep(1:length(x), length.out = length(y))
      x = x[ids]
    }
    if (length(y) < length(x)) {
      ids = rep(1:length(y), length.out = length(x))
      y = y[ids]
    }
    
    # Get coordinate matrices
    x_coords = st_coordinates(x)
    y_coords = st_coordinates(y)
    
    # Calculate azimuths
    x1 = x_coords[, 1]
    y1 = x_coords[, 2]
    x2 = y_coords[, 1]
    y2 = y_coords[, 2]
    az = (180 / pi) * atan2(x2 - x1, y2 - y1)
    names(az) = NULL
    az[az < 0] = az[az < 0] + 360
    
    # Replace with 'NA' for identical points
    az[x1 == x2 & y1 == y2] = NA
    
    return(az)
  }
  
  # Inicializar matrices
  k_esc <- matrix(nrow = nrow(m), ncol = 1)
  k_ele <- matrix(nrow = nrow(m), ncol = 1)
  k_com <- matrix(nrow = nrow(m), ncol = 1)
  azimut_grad <- matrix(nrow = nrow(m), ncol = 1)
  azimut_rad <- matrix(nrow = nrow(m), ncol = 1)
  distance_UTM <- matrix(nrow = nrow(m), ncol = 1)
  distance_top <- matrix(nrow = nrow(m), ncol = 1)
  e_top <- matrix(nrow = nrow(m), ncol = 1)
  n_top <- matrix(nrow = nrow(m), ncol = 1)
  
  for (i in 1:nrow(m)) {
    x <- 500000 - p[i, 1]
    q <- 0.000001 * x
    
    n_k <- a / ((1 - e2 * (sin(m[i, 2] * pi / 180))^2)^(1/2))
    p_k <- (10^(12)) * ((1 + et2 * (cos(m[i, 2] * pi / 180))^2) / (2 * n_k^2 * k_0^2))
    k_esc[i, 1] <- k_0 * (1 + p_k * q^2 + 0.00003 * q^4)
    k_ele[i, 1] <- (a * (1 - e2)) / (a * (1 - e2) + p[i, 3] * (1 - e2 * (sin(m[i, 2] * pi / 180))^2)^(3/2))
    k_com[i, 1] <- k_esc[i, 1] * k_ele[i, 1]
    if (i == 1) {
      azimut_grad[i, 1] <- 90
    } else {
      azimut_grad[i, 1] <- st_azimuth(h[1, 1], h[i, 1])
    }
    azimut_rad[i, 1] <- azimut_grad[i, 1] * pi / 180
    distance_UTM[i, 1] <- st_distance(h)[i, 1]
    distance_top[i, 1] <- distance_UTM[i, 1] / mean(c(k_com[i, 1], k_com[1, 1]))
    e_top[i, 1] <- p[1, 1] + distance_top[i, 1] * sin(azimut_rad[i, 1])
    n_top[i, 1] <- p[1, 2] + distance_top[i, 1] * cos(azimut_rad[i, 1])
  }
  
  # Crear data.frame de resultados
  datos <- data.frame(
    "Lat" = m[, 2], "Long" = m[, 1], 
    "k_esc" = k_esc, "K_ele" = k_ele, 
    "K_com" = k_com, "dist_UTM" = distance_UTM, 
    "Acimut_rad" = azimut_rad, "acimut_grad" = azimut_grad, 
    "dist_top" = distance_top, "E_top" = e_top, "N_top" = n_top
  )
  return(datos)
}


# Leer el archivo CSV
epsgGeodesicUtmFile <- read.csv("www/epsgGeodesicUtm.csv", stringsAsFactors = FALSE)

# UI del módulo
UTMPlanasModuleUi <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    titlePanel("Transformación UTM a Planas"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput(ns("archivoC"), label = h3("Subir archivo con datos"), accept = ".csv"),
        selectInput(ns("crs"), "Seleccione el CRS", choices = epsgGeodesicUtmFile$SRC),  # Aquí debes cargar tus CRS
        actionButton(ns("process"), "Procesar")
      ),
      mainPanel(
        h4("Resultados"),
        uiOutput(ns("panel_utm_planas")),
        uiOutput(ns("results_ui"))
      )
    )
  )
}

# Servidor del módulo
UTMPlanasModuleServer <- function(input, output, session) {
  ns <- session$ns
  datos <- reactiveValues()
  
  ## Genera el panel de procesos al momento de cargar los datos
  # Cargar CRS al inicio
  observe({
    # Leer el archivo CSV
    epsgGeodesicUtmFile <- read.csv("www/epsgGeodesicUtm.csv", stringsAsFactors = FALSE)
    
    # Validar que tenga las columnas esperadas
    if (!all(c("SRC", "EPSG") %in% colnames(epsgGeodesicUtmFile))) {
      showNotification("El archivo epsgGeodesicUtm.csv no tiene las columnas SRC y EPSG.", type = "error")
      return()
    }
    
    # Actualizar los valores del selector
    updateSelectInput(session, ns("crs"), choices = epsgGeodesicUtmFile$SRC)
    
    # Guardar el archivo en datos reactivos
    datos$epsgGeodesicUtmFile <- epsgGeodesicUtmFile
    
    print("Archivo base de CRS guardado exitosamente")
  })
  
  observeEvent(input$process, {
    req(input$archivoC)
    print("Archivo cargado por usuario:")
    print(input$archivoC)
  })
  
  observeEvent(input$process, { #INICIO
    req(input$archivoC)
  
  output$panel_utm_planas <- renderUI(expr = if (!is.null(input$archivoC)) {
    print("Renderizando panel UTM Planas")
    fluidPage(
      sidebarPanel(
        width = 3,
        #selectInput(ns("crs"), "CRS", epsgGeodesicUtmFile$SRC),
        h4("UTM correction"),
        actionButton(ns("inicio_correccion_UTM"), label = "Start", class = "btn-info")
      ),
      mainPanel(
        class = "well",
        h4("Select a row in the table (pivot point)"),
        dataTableOutput(ns("tabla_inicio_utm"))
      )
    )   
  } else {
    NULL
  })
  
  }) #FINAL
  
  ## Genera la tabla al momento de cargar los datos
  output$tabla_inicio_utm <- renderDataTable({
    req(input$archivoC)
    datos$puntos_coordenadasUTM <- read.csv(input$archivoC$datapath, sep = ",", header = TRUE)
    datatable(datos$puntos_coordenadasUTM, options = list(
      pageLength = 5,
      scrollX = TRUE,
      searching = FALSE
    ), selection = 'single')
  })
  
  ####### Realiza el ajuste de UTM planas
  observeEvent(input$inicio_correccion_UTM, {
    print("Iniciando ajuste de correcciones")
    if (length(input$tabla_inicio_utm_rows_selected) == 0) {
      showNotification(
        h4("Select a row regarding the pivot point"), 
        action = NULL, duration = 5, type = "warning"
      )
      return()
    }
    
    datos$datos_utm_ordenados <- rbind(
      datos$puntos_coordenadasUTM[input$tabla_inicio_utm_rows_selected, ], 
      datos$puntos_coordenadasUTM[-input$tabla_inicio_utm_rows_selected, ]
    )
    datos_utm <- datos$datos_utm_ordenados[, c(1:4)] # Ajusta índices si es necesario
    epsgGeodesicUtmFile <- read.csv("www/epsgGeodesicUtm.csv")
    crs_UTM_input <- epsgGeodesicUtmFile[epsgGeodesicUtmFile$SRC %in% input$crs, "EPSG"]
    print("Entrando a función UTM planas")
    datos$correccion_utm <- funcion_UTM_planas(datos_utm, as.numeric(crs_UTM_input))
  })
  
  output$results_ui <- renderUI({
    if (!is.null(datos$correccion_utm)) {
      fluidRow(
        column(
          12, 
          h3("Resultados"),
          dataTableOutput(ns("results_table")),
          downloadLink(ns('descarga_utm_correccion'), 'Descargar Resultados')
        )
      )
    } else {
      h4("No hay resultados aún.")
    }
  })
  
  output$results_table <- renderDataTable({
    req(datos$correccion_utm)
    datatable(datos$correccion_utm, options = list(
      pageLength = 5,
      scrollX = TRUE
    ))
  })
  
  output$descarga_utm_correccion <- downloadHandler(
    filename = "resultados_proceso_UTM-planas.csv",
    content = function(file) {
      write.csv(datos$correccion_utm, file, row.names = FALSE)
    }
  )
}
