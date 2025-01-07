# modules/utmToPlanas.R

library(shiny)
library(sf)
library(geosphere)
library(dplyr)
library(DT)
library(shinycssloaders)
library(leaflet)
library(leaflet.extras)
library(leafem)
library(sodium)

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
  
  for (i in 1:nrow(m)) { # O(n)
    cat(sprintf("Procesando fila %d\n", i))
    flush.console()
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
      azimut_grad[i, 1] <- st_azimuth(h[1, 1], h[i, 1]) # O(n)
    } #hacia arriba O(1)
    azimut_rad[i, 1] <- azimut_grad[i, 1] * pi / 180
    distance_UTM[i, 1] <- st_distance(h)[i, 1] # O(n)
    distance_top[i, 1] <- distance_UTM[i, 1] / mean(c(k_com[i, 1], k_com[1, 1]))
    n_top[i, 1] <- p[1, 1] + distance_top[i, 1] * sin(azimut_rad[i, 1]) # swap entre columnas Norte y Este
    e_top[i, 1] <- p[1, 2] + distance_top[i, 1] * cos(azimut_rad[i, 1])
  }
  
  # Crear data.frame de resultados
  datos <- data.frame(
    "Lat" = m[, 2], "Long" = m[, 1], 
    "k_esc" = k_esc, "K_ele" = k_ele, 
    "K_com" = k_com, "dist_UTM" = distance_UTM, 
    "Acimut_rad" = azimut_rad, "acimut_grad" = azimut_grad, 
    "dist_top" = distance_top, "N_top" = e_top, "E_top" = n_top
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
        helpText("Nota: Asegúrese de que su archivo CSV tenga el siguiente orden y nombres de columnas: Punto, ESTE, NORTE, h"),
        selectInput(ns("crs"), "Seleccione el CRS", choices = epsgGeodesicUtmFile$SRC, "Seleccione el CRS para el mapa generado"),  # Aquí debes cargar tus CRS
        actionButton(ns("process"), "Procesar")
      ),
      mainPanel(
        h4("Resultados"),
        uiOutput(ns("panel_utm_planas")),
        #uiOutput(ns("results_ui"))
      )
    ),
    # Sección adicional condicional, al mismo nivel que el título
    uiOutput(ns("results_ui"))
  )
}

# Servidor del módulo
UTMPlanasModuleServer <- function(input, output, session) {
  ns <- session$ns
  #datos <- reactiveValues()
  correccion_iniciada <- reactiveVal(FALSE)
  datos <- reactiveValues(mapa_datos_input = NULL, mapa_datos_resultados = NULL, puntos_coordenadasUTM = NULL, correccion_utm = NULL, datos_utm_ordenados = NULL)
  
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
  
  observeEvent(input$process, { #INICIO
    req(input$archivoC)
    print("Archivo cargado por usuario:")
    print(input$archivoC)
  
  output$panel_utm_planas <- renderUI(expr = if (!is.null(input$archivoC)) {
    print("Renderizando panel UTM Planas")
    fluidPage(
      sidebarPanel(
        width = 3,
        #selectInput(ns("crs"), "CRS", epsgGeodesicUtmFile$SRC),
        h4("UTM correction"),
        helpText("Previamente escoja un punto pivote en la tabla"),
        actionButton(ns("inicio_mensaje_UTM"), label = "Start", class = "btn-info")
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
  
  observeEvent(input$inicio_mensaje_UTM, {
    # Mostrar el mensaje modal inicial
    showModal(modalDialog(
      title = "Importante",
      fluidRow(h4("Advertencia: Considere un tiempo elevado de conversión para archivos con más de 1500 filas de entrada")),
      footer = actionButton(ns("inicio_correccion_UTM"), "Aceptar"),
      easyClose = FALSE # Obliga al usuario a aceptar para continuar
    ))
  })
  
  ####### Realiza el ajuste de UTM planas
  observeEvent(input$inicio_correccion_UTM, {
    removeModal()
    correccion_iniciada(TRUE)
    if (correccion_iniciada()) {
      print("correccion_iniciada cambiada a TRUE")
    }
    
    selected_rows <- isolate(input$tabla_inicio_utm_rows_selected)
    crs_input <- isolate(input$crs)
    puntos_coordenadasUTM <- isolate(datos$puntos_coordenadasUTM)
    datos_correccion_utm<- isolate(datos$correccion_utm)
    datos_utm_ordenados <- isolate(datos$datos_utm_ordenados)
    
    invalidateLater(100, session)
    later::later(function() {
      print("Iniciando ajuste de correcciones")
      if (length(selected_rows) == 0) {
        showNotification(
          h4("Select a row regarding the pivot point"), 
          action = NULL, duration = 5, type = "warning"
        )
        correccion_iniciada(FALSE)  # Detener el spinner si hay error
        return()
      }
      
      # Procesar datos con los valores capturados
      datos_utm_ordenados <- rbind(
        puntos_coordenadasUTM[selected_rows, ], 
        puntos_coordenadasUTM[-selected_rows, ]
      )
      datos_utm <- datos_utm_ordenados[, c(1:4)] # Ajusta índices si es necesario
      
      epsgGeodesicUtmFile <- read.csv("www/epsgGeodesicUtm.csv")
      crs_UTM_input <- epsgGeodesicUtmFile[epsgGeodesicUtmFile$SRC %in% crs_input, "EPSG"]
      print("Entrando a función UTM planas")
      datos_correccion_utm <- funcion_UTM_planas(datos_utm[,-1], as.numeric(crs_UTM_input))
      # Añadir nombres de puntos
      datos_correccion_utm <- cbind(Punto = datos_utm[,1], datos_correccion_utm) # Unir columnas de datos$correccion_utm con la primera columna de datos_utm que se llama "Punto"
      
      correccion_iniciada(FALSE)  # Marca que la corrección ha finalizado
      print("correccion_iniciada cambiada a FALSE")
      datos$correccion_utm <- datos_correccion_utm
      datos$datos_utm_ordenados <- datos_utm_ordenados
    }, delay = 0.1)  # Ejecutar después de un breve delay
    print(paste("Terminando el 'observeEvent(input$inicio_correccion_UTM' ¿correccion_utm es NULL?", is.null(datos$correccion_utm)))
  })
  
  ## Genera la ventana del mapa cuando los datos fueron creados correctamente
  ## MAPA
  
  output$empty_output <- renderUI({
    if (correccion_iniciada()) {
      print("Generando spinner")
      withSpinner(uiOutput("no_output"), color = "#0dc5c1")
    } else {
      print("Sin spinner generado")
      h4("No hay resultados aún.")  # Mensaje predeterminado
    }
  })
  
  output$results_ui <- renderUI({
    print(paste("Evaluando results_ui. ¿correccion_utm es NULL?", is.null(datos$correccion_utm)))
    if (!is.null(datos$correccion_utm)) {
      fluidRow(
        column(
          12, 
          h3("Resultados"),
          dataTableOutput(ns("results_table")),
          downloadLink(ns('descarga_utm_correccion'), 'Descargar Resultados'),
          h1("Generación de mapa interactivo"),
          column(10,
                 selectInput(ns("crs_mapa_UTM"), "Seleccione el CRS", epsgGeodesicUtmFile$SRC, "Seleccione el sistema de coordenadas CRS para la conversión"),
                 actionButton(ns("crear_mapa_utm"), "Start"),
                 uiOutput(ns("panel_mapa")) 
          )
        )
      )
    } else {
      fluidRow(
        column(
          12,
          #h4("No hay resultados aún."),
          uiOutput(ns("empty_output"))  # Este es dinámico
        )
      )
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
  
  # TODO: GENERACION DE MAPA INTERACTIVO
  
  output$panel_mapa <- renderUI({
    if (!is.null(datos$mapa_datos_input) && !is.null(datos$mapa_datos_resultados)) {
      print("Datos no nulos de entrada. Renderizando mapa")
      tags$div(
        style = "width: 100%; height: calc(100vh - 100px);", #ajusta el alto relativo a la ventana
        leafletOutput(ns("mapa"), width = "100%", height = "100%") #asegura que ocupe todo el contenedor
      )
    } else {
      tags$div(
        h4("No hay datos disponibles para generar el mapa."),
        style = "width: 100%; height: calc(100vh - 100px); background-color: #f5f5f5; border: 1px dashed #ccc; display: flex; justify-content: center; align-items: center;"  # Estilo para reservar espacio
      )
    }
  })
  
  ## Crea el mapa
  output$mapa<- renderLeaflet({
    req(datos$mapa_datos_input, datos$mapa_datos_resultados)
    
    print("Datos de entrada para el mapa:")
    print(st_geometry(datos$mapa_datos_input))
    print("Datos de resultados para el mapa:")
    print(st_geometry(datos$mapa_datos_resultados))
    
    mapa<-leaflet() %>%
      addProviderTiles(providers$OpenStreetMap.Mapnik, group = "OpenStreetMap.Mapnik") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri.WorldImagery") %>%
      addLayersControl(
        baseGroups = c("OpenStreetMap.Mapnik","Esri.WorldImagery"),
        overlayGroups = c("Control Points", "Results")
      )%>%
      addCircleMarkers(data = datos$mapa_datos_input, color = "red", group = "Control Points",
                       label = ~as.character(datos$mapa_datos_input$Punto))%>%
      addCircleMarkers(data=datos$mapa_datos_resultados, color="blue", group = "Results",
                       label = ~as.character(datos$mapa_datos_input$Punto))
    
    mapa
  })
  
  ## Genera los datos para crea el mapa
  observeEvent(input$crear_mapa_utm,{
    print("Generando mapa")
    print(datos$datos_utm_ordenados)
    datos_utm<-datos$datos_utm_ordenados[, c(1:4)]
    print(paste("Filas en datos_utm:", nrow(datos_utm)))
    print(paste("Filas en datos$correccion_utm:", nrow(datos$correccion_utm)))
    ##crs cargados
    epsgGeodesicUtmFile <- read.csv("www/epsgGeodesicUtm.csv")
    crs_mapa_UTM_input <- epsgGeodesicUtmFile[epsgGeodesicUtmFile$SRC %in% input$crs_mapa_UTM,c("EPSG")]
    print("Crs cargados")
    datos_mapa_UTM<-cbind(datos_utm,datos$correccion_utm) #referencia del error
    print("Datos_mapa_UTM")
    print(datos_mapa_UTM)
    
    sf_datos_locales<-st_as_sf(datos_mapa_UTM, coords=c("Long","Lat"), crs=4326)
    
    names(sf_datos_locales)[1]<-c("Punto")## No cambiar este nombre, de el depende el etiquetado del mapa
    datos$mapa_datos_input<-sf_datos_locales##Datos de entrada
    print("Mapa_datos_input: ")
    print(datos$mapa_datos_input)
    
    # TODO: Agregar manejo de NaN's
    
    datos_mapa_UTM <- datos_mapa_UTM[complete.cases(datos_mapa_UTM[, c("E_top", "N_top")]), ]
    
    if (any(is.na(datos_mapa_UTM$E_top) | is.na(datos_mapa_UTM$N_top))) {
      showNotification(
        "Advertencia: Se encontraron coordenadas faltantes y se han eliminado.",
        type = "warning",  # Tipo de notificación: mensaje de advertencia
        duration = 3       # Duración en segundos
      )
    }
    
    sf_datos_resultados<-st_as_sf(datos_mapa_UTM, coords=c("E_top","N_top"), crs=as.numeric(crs_mapa_UTM_input)) # orden invertido respecto al original
    
    sf_datos_resultados<-st_transform(sf_datos_resultados, 4326)
    
    names(sf_datos_resultados)[1]<-c("Punto") ## No cambiar este nombre, de el depende el etiquetado del mapa
    datos$mapa_datos_resultados<-sf_datos_resultados ## Datos de resultados
    print("Mapa_datos_resultados: ")
    print(datos$mapa_datos_resultados)
  })
  
}
