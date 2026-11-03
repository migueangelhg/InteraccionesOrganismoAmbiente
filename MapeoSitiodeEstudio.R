# ------------------------------------------------------------------------------
# 17. Mapa estático del área de estudio con OpenStreetMap
# ------------------------------------------------------------------------------
# Requiere instalar: install.packages(c("ggspatial", "maptiles", "sf"))

library(ggspatial)
library(maptiles)
library(sf)
library(ggplot2)

# Convertir los datos de ocurrencia a un objeto espacial (sf)
occ_sf <- st_as_sf(occ_its, coords = c("longitude", "latitude"), crs = 4326)

# Calcular el centro y el zoom aproximado para el mapa
bbox_occ <- st_bbox(occ_sf)
centro_lon <- mean(c(bbox_occ["xmin"], bbox_occ["xmax"]))
centro_lat <- mean(c(bbox_occ["ymin"], bbox_occ["ymax"]))

cat("Generando mapa estático con OpenStreetMap...\n")

p_mapa_estatico <- ggplot() +
  # Descargar y añadir el mapa base de OSM (el zoom 9 es bueno para escala regional)
  annotation_map_tile(type = "osm", zoom = 9) +
  
  # Añadir los puntos de muestreo
  geom_sf(data = occ_sf, aes(color = Estado, shape = Estado), 
          size = 3, alpha = 0.9) +
  
  # Definir colores y formas (mismos que en tus gráficos anteriores)
  scale_color_manual(values = c("Sano" = "forestgreen", "Enfermo" = "firebrick")) +
  scale_shape_manual(values = c("Sano" = 21, "Enfermo" = 24)) +
  
  # Añadir barra de escala y flecha del norte
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
  
  # Ajustar el encuadre (hacer un zoom al bbox de los puntos)
  coord_sf(xlim = c(bbox_occ["xmin"] - 0.5, bbox_occ["xmax"] + 0.5),
           ylim = c(bbox_occ["ymin"] - 0.5, bbox_occ["ymax"] + 0.5)) +
  
  # Títulos y temas
  labs(title = "Distribución de muestras de Quercus robur",
       subtitle = "Suroeste de Francia (Sano vs Enfermo)",
       x = "Longitud", y = "Latitud",
       color = "Estado", shape = "Estado") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12))

# Guardar el mapa en alta resolución
ruta_mapa_estatico <- file.path(out_dir, "mapa_estudio_osm_ggplot.png")
ggsave(ruta_mapa_estatico, p_mapa_estatico, width = 8, height = 7, dpi = 300)

cat("Mapa guardado en:", ruta_mapa_estatico, "\n")











###################################


# Instalar tidyterra si no lo tienes
if (!requireNamespace("tidyterra", quietly = TRUE)) install.packages("tidyterra")
library(tidyterra)
library(maptiles)

# 1. Descargar las teselas de OpenStreetMap para el area de los puntos
osm_tiles <- get_tiles(occ_sf, provider = "OpenStreetMap", zoom = 9)

# 2. Crear el mapa
p_mapa_estatico <- ggplot() +
  # Añadir el mapa base desde el raster descargado
  tidyterra::geom_spatraster_rgb(data = osm_tiles, alpha = 0.8) +
  
  # Añadir los puntos de muestreo
  geom_sf(data = occ_sf, aes(color = Estado, shape = Estado), 
          size = 3, alpha = 0.9) +
  
  # Definir colores y formas (mismos que en tus gráficos anteriores)
  scale_color_manual(values = c("Sano" = "forestgreen", "Enfermo" = "firebrick")) +
  scale_shape_manual(values = c("Sano" = 21, "Enfermo" = 24)) +
  
  # Añadir barra de escala y flecha del norte (estas si funcionan sin rosm)
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
  
  # Ajustar el encuadre
  coord_sf(xlim = c(bbox_occ["xmin"] - 0.5, bbox_occ["xmax"] + 0.5),
           ylim = c(bbox_occ["ymin"] - 0.5, bbox_occ["ymax"] + 0.5)) +
  
  # Títulos y temas
  labs(title = "Distribución de muestras de Quercus robur",
       subtitle = "Suroeste de Francia (Sano vs Enfermo)",
       x = "Longitud", y = "Latitud",
       color = "Estado", shape = "Estado") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12))

# 3. Guardar el mapa
ruta_mapa_estatico <- file.path(out_dir, "mapa_estudio_osm_ggplot.png")
ggsave(ruta_mapa_estatico, p_mapa_estatico, width = 8, height = 7, dpi = 300)

cat("Mapa guardado en:", ruta_mapa_estatico, "\n")











# ------------------------------------------------------------------------------
# 17. Mapa estático mejorado (OpenStreetMap + ggplot2)
# ------------------------------------------------------------------------------
library(ggspatial)
library(maptiles)
library(sf)
library(ggplot2)

# Convertir a objeto espacial
occ_sf <- st_as_sf(occ_its, coords = c("longitude", "latitude"), crs = 4326)
bbox_occ <- st_bbox(occ_sf)

# Mapa con estilo más limpio y puntos más pequeños
p_mapa_estatico <- ggplot() +
  # Usamos un estilo de mapa más sobrio y un zoom ligeramente más amplio
  annotation_map_tile(type = "cartolight", zoom = 8, progress = "none") +
  
  # Puntos más pequeños y con transparencia para no saturar
  geom_sf(data = occ_sf, aes(color = Estado, shape = Estado), 
          size = 1.8, alpha = 0.8, stroke = 0.5) +
  
  scale_color_manual(values = c("Sano" = "forestgreen", "Enfermo" = "firebrick")) +
  scale_shape_manual(values = c("Sano" = 21, "Enfermo" = 24)) +
  
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
  
  # Encuadre un poco más amplio para dar contexto
  coord_sf(xlim = c(bbox_occ["xmin"] - 0.8, bbox_occ["xmax"] + 0.8),
           ylim = c(bbox_occ["ymin"] - 0.8, bbox_occ["ymax"] + 0.8)) +
  
  labs(title = "Distribución de muestras de Quercus robur",
       subtitle = "Suroeste de Francia (Sano vs Enfermo)",
       x = "Longitud", y = "Latitud",
       color = "Estado", shape = "Estado") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12))

# Guardar con una resolución adecuada
ruta_mapa_estatico <- file.path(out_dir, "mapa_estudio_osm_ggplot_mejorado.png")
ggsave(ruta_mapa_estatico, p_mapa_estatico, width = 8, height = 7, dpi = 300)


library(leaflet)
library(htmlwidgets)

# 1. Paleta de colores y límites
pal_leaflet <- colorFactor(palette = c("forestgreen", "firebrick"), 
                           domain = occ_its$Estado)

lon_min <- as.numeric(bbox_occ["xmin"]) - 0.1
lon_max <- as.numeric(bbox_occ["xmax"]) + 0.1
lat_min <- as.numeric(bbox_occ["ymin"]) - 0.1
lat_max <- as.numeric(bbox_occ["ymax"]) + 0.1

# 2. Construir el mapa interactivo
mapa_leaflet <- leaflet(data = occ_its) %>%
  
  # --- Capas Base ---
  addProviderTiles(providers$Esri.WorldStreetMap, group = "Calles y Parques") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satélite") %>%
  
  # --- Puntos de muestreo ---
  addCircleMarkers(
    lng = ~longitude, lat = ~latitude,
    color = "black", weight = 1.5,
    fillColor = ~pal_leaflet(Estado), radius = 7, fillOpacity = 0.9,
    label = ~paste("Muestra:", run_accession),
    labelOptions = labelOptions(direction = "auto", textsize = "14px", style = list("font-weight" = "bold")),
    popup = ~paste("<b>Muestra:</b>", run_accession, "<br>",
                   "<b>Estado:</b>", Estado, "<br>",
                   "<b>Marcador:</b>", marcador)
  ) %>%
  
  # --- Escala y Leyenda ---
  addScaleBar(position = "bottomleft", 
              options = scaleBarOptions(imperial = FALSE, metric = TRUE, maxWidth = 150)) %>%
  addLegend(
    position = "bottomright", pal = pal_leaflet, values = ~Estado,
    title = "<b>Estado de salud</b>", opacity = 1
  ) %>%
  
  # --- Control de Capas ---
  addLayersControl(
    baseGroups = c("Calles y Parques", "Satélite"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  
  # --- INCRUSTAR EL INSET (MINI MAPA) ---
  addMiniMap(
    position = "topleft",         # Posición (arriba a la izquierda)
    width = 150, height = 150,    # Tamaño en píxeles
    zoomLevelOffset = -8,         # Qué tanto se aleja (-8 muestra Europa/África)
    toggleDisplay = TRUE,         # Permite ocultarlo con una flechita
    minimized = FALSE,            # Inicia abierto
    mapOptions = list(
      zoomControl = FALSE,        # Sin botones de zoom en el inset
      attributionControl = FALSE  # Sin créditos para que se vea limpio
    )
  ) %>%
  
  # Encuadrar al área de estudio
  fitBounds(lng1 = lon_min, lat1 = lat_min, lng2 = lon_max, lat2 = lat_max)

# 3. Guardar y abrir en el navegador
ruta_html <- file.path(out_dir, "mapa_interactivo_con_inset.html")
saveWidget(mapa_leaflet, ruta_html, selfcontained = TRUE)

# Abre el mapa en tu navegador
browseURL(ruta_html)
