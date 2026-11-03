# ============================================================
# FASE 1: Verificación + QC + filtrado DADA2 (ITS)
# ============================================================

# --- 0. Paquetes ---------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#if (!requireNamespace("dada2",      quietly = TRUE)) BiocManager::install("dada2")
#if (!requireNamespace("ShortRead",  quietly = TRUE)) BiocManager::install("ShortRead")
library(dada2)
library(ShortRead)

# --- 1. Inventario de archivos -------------------------------
dir_its   <- file.path("fastq_qrobur", "ITS")
fwd_files <- sort(list.files(dir_its, pattern = "_1.fastq\\.gz$", full.names = TRUE))
rev_files <- sort(list.files(dir_its, pattern = "_2.fastq\\.gz$", full.names = TRUE))

cat("Archivos R1 encontrados:", length(fwd_files), "\n")
cat("Archivos R2 encontrados:", length(rev_files), "\n")

# Comprobar emparejamiento R1 ↔ R2
sample_names <- sub("_1\\.fastq\\.gz$", "", basename(fwd_files))
sample_names_rev <- sub("_2\\.fastq\\.gz$", "", basename(rev_files))
stopifnot(identical(sample_names, sample_names_rev))
cat("Emparejamiento R1/R2 correcto. Muestras únicas:", length(sample_names), "\n\n")

# --- 2. Verificar integridad de 2 archivos de prueba ---------
cat("=== Verificación de 2 muestras de prueba ===\n")
for (i in c(1, 2)) {
  f1 <- fwd_files[i]; f2 <- rev_files[i]
  # Leer solo las primeras 1000 secuencias
  r1 <- readFastq(f1)
  r2 <- readFastq(f2)
  cat(basename(f1), "->", length(r1), "reads, longitud media:", 
      round(mean(width(r1)), 1), "pb\n")
  cat(basename(f2), "->", length(r2), "reads, longitud media:", 
      round(mean(width(r2)), 1), "pb\n\n")
}
rm(r1, r2); gc()

# --- 3. QC de calidad por posición (2 muestras) --------------
cat("=== Calidad por posición (muestra 1) ===\n")
qc_fwd <- dada2::plotQualityProfile(fwd_files[1:2])
print(qc_fwd)

# --- 4. Definir salida filtrada ------------------------------
dir_filt <- file.path("fastq_qrobur", "ITS_filtered")
dir.create(dir_filt, showWarnings = FALSE)
filt_fwd <- file.path(dir_filt, paste0(sample_names, "_F_filt.fastq.gz"))
filt_rev <- file.path(dir_filt, paste0(sample_names, "_R_filt.fastq.gz"))

# --- 5. Filtrado DADA2 (ajusta truncLen según QC visual) -----
# Valores iniciales conservadores para ITS:
#   truncLen = c(240, 200) deja margen para solapamiento R1-R2
#   maxN = 0, maxEE = c(2, 5), truncQ = 2
cat("\n=== Filtrando", length(fwd_files), "muestras (puede tardar 5-15 min) ===\n")
out <- filterAndTrim(
  fwd      = fwd_files,
  filt     = filt_fwd,
  rev      = rev_files,
  filt.rev = filt_rev,
  truncLen = c(240, 200),
  maxN     = 0,
  maxEE    = c(2, 5),
  truncQ   = 2,
  rm.phix  = TRUE,
  compress = TRUE,
  multithread = FALSE,   # Windows: dejar FALSE para evitar problemas
  verbose  = TRUE
)
cat("\n=== Filtrado completado ===\n")
print(head(out, 10))

# --- 6. Resumen de retención --------------------------------
out_df <- data.frame(
  muestra  = rownames(out),
  reads_in = out[, "reads.in"],
  reads_out = out[, "reads.out"],
  pct      = round(100 * out[, "reads.out"] / out[, "reads.in"], 1)
)
cat("\nRetención media:", round(mean(out_df$pct), 1), "%\n")
cat("Muestras con <50% retención:", sum(out_df$pct < 50), "\n")
cat("Muestras con <1000 reads finales:", sum(out_df$reads_out < 1000), "\n")

write.csv(out_df, "ITS_filter_summary.csv", row.names = FALSE)
cat("\nResumen guardado en ITS_filter_summary.csv\n")









##################### Identificación de la geografia de esas muestras ########

# ============================================================
# VERIFICACIÓN GEOGRÁFICA USANDO DATOS LOCALES (Híbrido)
# ============================================================
library(httr)
library(readr)
library(dplyr)

# 1. Asumimos que el objeto 'seleccion' de tus scripts anteriores 
# aún está en memoria. Si no, lo recreamos rápido:
if (!exists("seleccion")) {
  stop("El objeto 'seleccion' no está en memoria. Por favor, corre el script de descarga primero o dime para recrearlo.")
}

# 2. Extraer los IDs de muestra únicos que ya descargamos
muestras_ids <- unique(seleccion$sample_accession)
cat("Consultando geografía para", length(muestras_ids), "muestras específicas...\n")

# 3. Consultar la API del ENA SOLO para esas muestras (usando la opción 'accession')
# Nota: ENA permite consultar por 'sample_accession' directamente
geo_list <- list()
for (i in seq_along(muestras_ids)) {
  id <- muestras_ids[i]
  url <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", id,
                "&result=read_run&format=tsv&download=true&fields=lat,lon,country,location")
  
  res <- tryCatch(read_tsv(I(content(GET(url), as = "text", encoding = "UTF-8")), 
                           show_col_types = FALSE), 
                  error = function(e) NULL)
  
  if (!is.null(res) && nrow(res) > 0) {
    res$sample_accession <- id
    geo_list[[i]] <- res
  }
}

# 4. Unir los resultados
if (length(geo_list) > 0) {
  meta_geo <- bind_rows(geo_list)
  
  # 5. Cruzar con tu tabla original para mantener el alias de la muestra
  tabla_final <- seleccion %>%
    left_join(meta_geo, by = "sample_accession") %>%
    select(sample_alias, host, lat, lon, country, location) %>%
    distinct()
  
  # 6. Verificar resultados
  tiene_coords <- !is.na(tabla_final$lat) & !is.na(tabla_final$lon)
  cat("\n--- RESULTADO ---\n")
  cat("Muestras de Q. robur CON coordenadas:", sum(tiene_coords), "\n")
  cat("Muestras de Q. robur SIN coordenadas:", sum(!tiene_coords), "\n\n")
  
  if (sum(tiene_coords) > 0) {
    print(head(tabla_final[tiene_coords, ], 10))
  }
} else {
  cat("\nNo se pudo obtener información geográfica para estas muestras.\n")
}



##################################
#Extraccion de informacion geografica obtenida para Q. robur em cada uno de los 
#muestreos
####################################

library(leaflet)
library(readxl)
library(dplyr)
library(htmltools)
library(sf)

# 1. Cargar datos
ruta_excel <- "MaterialSuplementario/pbiomes-01-20-0004-r.sm1.xlsx"
meta <- read_excel(ruta_excel, sheet = "metadata")
mapas <- read_excel(ruta_excel, sheet = "maps")

# 2. Detectar columna de necrosis y limpiar metadata
col_necrosis <- grep("^Internal necrosis", names(meta), value = TRUE)
meta_limpia <- meta %>%
  select(`Samples ID`:`Bot number`) %>%
  rename(Internal_necrosis = all_of(col_necrosis))

# 3. Limpiar mapas
mapas <- mapas %>%
  filter(!is.na(Site)) %>%
  select(Site, longitude, latitude, Location) %>%
  distinct()

# 4. Filtrar SOLO Quercus robur (SANOS Y ENFERMOS)
q_robur <- meta_limpia %>%
  filter(Species == "robur") %>%
  left_join(mapas, by = "Site")

# 5. Jitter + colores + conteo
set.seed(123)
q_robur <- q_robur %>%
  mutate(
    lon_jit = longitude + runif(n(), -0.003, 0.003),
    lat_jit = latitude + runif(n(), -0.003, 0.003),
    color = case_when(
      Internal_necrosis == "A" ~ "darkgreen",
      Internal_necrosis == "S" ~ "red",
      TRUE                     ~ "#555555"   # gris oscuro visible
    ),
    estado = case_when(
      Internal_necrosis == "A" ~ "Sano (A)",
      Internal_necrosis == "S" ~ "Enfermo (S)",
      TRUE                     ~ "No determinado (ND)"
    )
  )

# 6. Conteo dinámico para verificar y para la leyenda
tabla_conteo <- table(q_robur$estado)
cat("=== Conteo de árboles por estado ===\n")
print(tabla_conteo)

# 7. Popups
q_robur$popup <- paste0(
  "<b>Árbol ID: ", q_robur$`Samples ID`, "</b><br>",
  "<i>Sitio: </i>", q_robur$Site, "<br>",
  "<i>Localidad: </i>", q_robur$Location, "<br>",
  "<i>Estado: </i>", q_robur$estado, "<br><hr>",
  "<b>Marcadores disponibles:</b> ITS y LSU"
)

# 8. sf y coords
q_robur_sf <- st_as_sf(q_robur, coords = c("lon_jit", "lat_jit"), crs = 4326)
coords <- st_coordinates(q_robur_sf)

# 9. Mapa corregido
leaflet(q_robur_sf, options = leafletOptions(minZoom = 5, maxZoom = 18)) %>%
  
  addProviderTiles(providers$OpenStreetMap,    group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satélite") %>%
  addProviderTiles(providers$OpenTopoMap,      group = "OpenTopoMap (relieve)") %>%
  
  addCircleMarkers(
    radius       = 4,           # más pequeños
    color        = ~color,
    fillColor    = ~color,
    fillOpacity  = 1,           # opacidad total
    weight       = 1,           # borde fino
    stroke       = TRUE,
    popup        = ~popup,
    label        = ~paste(`Samples ID`, "-", estado),
    labelOptions = labelOptions(noHide = FALSE, direction = "auto"),
    group        = "Árboles Q. robur"
  ) %>%
  
  addLegend(
    position = "bottomright",
    colors   = c("darkgreen", "red", "#555555"),
    labels   = c(
      paste0("Sano (A) – n = ", sum(q_robur$estado == "Sano (A)")),
      paste0("Enfermo (S) – n = ", sum(q_robur$estado == "Enfermo (S)")),
      paste0("No determinado (ND) – n = ", sum(q_robur$estado == "No determinado (ND)"))
    ),
    title    = "Estado de salud<br>de Q. robur",
    opacity  = 1
  ) %>%
  
  addLayersControl(
    baseGroups    = c("OpenStreetMap", "Esri Satélite", "OpenTopoMap (relieve)"),
    overlayGroups = c("Árboles Q. robur"),
    options       = layersControlOptions(collapsed = FALSE)
  ) %>%
  
  fitBounds(
    lng1 = min(coords[, 1]), lat1 = min(coords[, 2]),
    lng2 = max(coords[, 1]), lat2 = max(coords[, 2])
  )

x11()

#### tABLA DE DATOS MAESTRA
# ============================================================
# CREACIÓN DE LA TABLA MAESTRA DEFINITIVA (Sin duplicados)
# ============================================================
library(dplyr)
library(readxl)
library(stringr)

# --- 1. Cargar y limpiar la tabla de mapas (Solución al S23) ---
ruta_excel <- "MaterialSuplementario/pbiomes-01-20-0004-r.sm1.xlsx"
mapas <- read_excel(ruta_excel, sheet = "maps") %>%
  filter(!is.na(Site)) %>%
  # Agrupamos por sitio y promediamos las coordenadas para S23 (y cualquier otro duplicado)
  group_by(Site) %>%
  summarise(
    longitude = mean(longitude, na.rm = TRUE),
    latitude  = mean(latitude, na.rm = TRUE),
    .groups = "drop"
  )

# --- 2. Extraer la información del ENA (seleccion) ---
# Limpiamos los sufijos _ITS y _LSU para obtener el alias base (ej. S2.O1.A)
seleccion$alias_base <- gsub("_(ITS|LSU)$", "", seleccion$sample_alias)

# Extraemos el Sitio (S1, S2...) y el Estado (A o S)
seleccion$Site_ena   <- str_extract(seleccion$alias_base, "^S\\d+")
seleccion$Estado_ena <- str_extract(seleccion$alias_base, "[AS]$")

# --- 3. Cruzar la información (¡El paso clave!) ---
master_metadata <- seleccion %>%
  # Unimos por el sitio
  left_join(mapas, by = c("Site_ena" = "Site")) %>%
  # Traducimos el estado
  mutate(
    Estado = case_when(
      Estado_ena == "A" ~ "Sano",
      Estado_ena == "S" ~ "Enfermo",
      TRUE              ~ "Desconocido"
    )
  ) %>%
  # ESTA ES LA LÍNEA MÁGICA: Elimina cualquier duplicado por número de acceso al run
  distinct(run_accession, .keep_all = TRUE)

# --- 4. Verificación de los números ---
cat("=== TABLA MAESTRA CREADA EXITOSAMENTE ===\n")
cat("Total de archivos FASTQ únicos:", nrow(master_metadata), "\n\n")

cat("--- Distribución por marcador ---\n")
print(table(master_metadata$marcador))

cat("\n--- Distribución por estado de salud ---\n")
print(table(master_metadata$Estado))

# --- 5. Guardar la tabla para el futuro ---
write.csv(master_metadata, "Master_Metadata_FASTQ.csv", row.names = FALSE)
cat("\n Tabla guardada como 'Master_Metadata_FASTQ.csv'\n")


# ============================================================
# INVESTIGACIÓN DE DISCREPANCIAS
# ============================================================
library(dplyr)
library(readxl)

# 1. ¿Qué árbol del ENA no está como "robur" en el Excel?
meta_check <- read_excel(ruta_excel, sheet = "metadata")

# Lista de alias_base del ENA (63 árboles)
alias_ena <- unique(master_metadata$alias_base)

# Lista de Samples ID del Excel para robur
id_robur_excel <- meta_check %>%
  filter(Species == "robur") %>%
  pull(`Samples ID`)

# ¿Cuál alias del ENA NO está en la lista de robur del Excel?
faltantes <- setdiff(alias_ena, id_robur_excel)
cat("Árboles del ENA que NO son 'robur' en el Excel:", length(faltantes), "\n")
print(faltantes)

# ¿Cómo está etiquetado ese árbol en el Excel?
if (length(faltantes) > 0) {
  cat("\n--- Metadata del árbol faltante en el Excel ---\n")
  meta_check %>%
    filter(`Samples ID` %in% faltantes) %>%
    select(`Samples ID`, Species, Site, Location) %>%
    print()
}

# 2. ¿Por qué S24 no está en el ENA?
cat("\n--- Árboles de S24 en el Excel ---\n")
s24_excel <- meta_check %>%
  filter(Site == "S24") %>%
  select(`Samples ID`, Species, Site, Location)
print(s24_excel)

cat("\nÁrboles de S24 que SÍ están en el ENA:", 
    sum(s24_excel$`Samples ID` %in% alias_ena), "\n")
