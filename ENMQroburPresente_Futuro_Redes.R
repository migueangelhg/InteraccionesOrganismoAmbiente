# ==============================================================================
# ENM con elipsoides - Tres modelos: especie (GBIF), Sano (ITS), Enfermo (ITS)
# Q. robur
# ------------------------------------------------------------------------------
# Entradas:
#   - 0023551-261023074918325.csv (GBIF, ya filtrado en la descarga)
#   - Master_Metadata_FASTQ.csv (muestras ITS)
#   - WorldClim 2.1 a 30 arc-sec (descarga automatica de multiples tiles)
# Salidas en resultados/Elipsoides_Qrobur/:
#   - gbif_limpio.csv
#   - bio_crop.tif
#   - objetos_elipsoides.rds
#   - solapamiento_3_modelos.txt
#   - plot_solapamiento_3_modelos.png
# ==============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(ggplot2)
  library(geodata); library(CoordinateCleaner); library(ntbox)
})

root <- "C:/Users/31628/Documents/InteraccionPlantaMicroPracticas"
setwd(root)
out_dir <- file.path(root, "resultados", "Elipsoides_Qrobur")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

worldclim_res <- 0.5
margen_deg    <- 1.0


# ------------------------------------------------------------------------------
# 1. Localizacion y carga del archivo GBIF
# ------------------------------------------------------------------------------
nombre_gbif <- "0023551-261023074918325.csv"
coincidencias <- list.files(root, pattern = paste0("^", nombre_gbif, "$"),
                            recursive = TRUE, full.names = TRUE)
if (!length(coincidencias)) {
  coincidencias <- list.files(root, pattern = "^\\d+-\\d+\\.csv$",
                              recursive = TRUE, full.names = TRUE)
}
if (!length(coincidencias)) stop("No se encontro el archivo GBIF en ", root)
ruta_gbif <- coincidencias[1]
cat("Archivo GBIF:", ruta_gbif, "\n")
cat("Tamano:", round(file.info(ruta_gbif)$size / 1e6, 1), "MB\n")

gbif_raw <- tryCatch(
  rgbif::occ_download_import(ruta_gbif),
  error = function(e) {
    message("occ_download_import fallo, lectura manual: ", e$message)
    lineas <- readLines(ruta_gbif, n = 50, warn = FALSE)
    idx <- grep("gbifID|occurrenceID|scientificName", lineas)[1]
    if (is.na(idx)) idx <- 1
    read.delim(ruta_gbif, skip = idx - 1, sep = "\t",
               quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  })
cat("Registros GBIF cargados:", nrow(gbif_raw), "\n")


# ------------------------------------------------------------------------------
# 2. Depuracion de GBIF
# ------------------------------------------------------------------------------
gbif_clean <- gbif_raw %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude),
         decimalLongitude != 0, decimalLatitude != 0)

if ("coordinateUncertaintyInMeters" %in% names(gbif_clean)) {
  gbif_clean <- gbif_clean %>%
    filter(is.na(coordinateUncertaintyInMeters) |
             coordinateUncertaintyInMeters <= 10000)
  cat("Tras filtro de incertidumbre <= 10 km:", nrow(gbif_clean), "\n")
}

species_col <- if ("species" %in% names(gbif_clean)) "species" else "scientificName"
gbif_cc <- CoordinateCleaner::clean_coordinates(
  gbif_clean,
  lon = "decimalLongitude", lat = "decimalLatitude",
  species = species_col,
  tests = c("capitals", "centroids", "equal", "gbif", "institutions",
            "zeros", "seas", "duplicates"),
  seas_scale = 10, verbose = FALSE
)
cols_test <- grep("^\\.", names(gbif_cc), value = TRUE)
if (length(cols_test)) gbif_cc <- gbif_cc[apply(gbif_cc[, cols_test, drop = FALSE], 1, all), ]
cat("Tras CoordinateCleaner:", nrow(gbif_cc), "\n")

gbif_cc$grid_id <- paste(round(gbif_cc$decimalLongitude, 1),
                         round(gbif_cc$decimalLatitude, 1), sep = "_")
set.seed(42)
gbif_thin <- gbif_cc %>%
  group_by(grid_id) %>% slice_sample(n = 1) %>% ungroup() %>%
  select(-grid_id)
cat("Tras rarificacion a 0.1 grados:", nrow(gbif_thin), "\n")

if (nrow(gbif_thin) > 10000) {
  set.seed(42)
  gbif_thin <- gbif_thin[sample(nrow(gbif_thin), 10000), ]
  cat("Submuestreado a 10,000 registros\n")
}

write.csv(gbif_thin, file.path(out_dir, "gbif_limpio.csv"), row.names = FALSE)


# ------------------------------------------------------------------------------
# 3. Carga de muestras ITS
# ------------------------------------------------------------------------------
master <- read.csv("Master_Metadata_FASTQ.csv", stringsAsFactors = FALSE)
occ_its <- master %>%
  filter(marcador == "ITS", !is.na(longitude), !is.na(latitude),
         Estado %in% c("Sano", "Enfermo"))
cat("\nMuestras ITS: Sano =", sum(occ_its$Estado == "Sano"),
    "| Enfermo =", sum(occ_its$Estado == "Enfermo"), "\n")

if (any(table(occ_its$Estado) < 5)) {
  stop("Un grupo tiene menos de 5 muestras. No se puede ajustar el elipsoide.")
}


# ------------------------------------------------------------------------------
# 4. WorldClim 30 arc-sec: mosaico de tiles que cubren todo el bbox
# ------------------------------------------------------------------------------
# El bbox cubre ~67 grados de longitud por ~28 de latitud. Un solo tile de
# WorldClim abarca ~30 x 30 grados, por lo que se necesita una grilla de
# varios tiles. Se descargan, se mosaicizan y se recorta al bbox.

lon_all <- c(gbif_thin$decimalLongitude, occ_its$longitude)
lat_all <- c(gbif_thin$decimalLatitude,  occ_its$latitude)
bbox <- c(xmin = min(lon_all) - margen_deg, xmax = max(lon_all) + margen_deg,
          ymin = min(lat_all) - margen_deg, ymax = max(lat_all) + margen_deg)
cat("BBox de estudio: lon [", round(bbox["xmin"],2), ",", round(bbox["xmax"],2),
    "] lat [", round(bbox["ymin"],2), ",", round(bbox["ymax"],2), "]\n")

wc_dir <- file.path(root, "data", "worldclim_30s")
dir.create(wc_dir, recursive = TRUE, showWarnings = FALSE)

# Grilla de centros de tiles (3 x 3 = 9 cubre con solapamiento)
lon_breaks <- seq(bbox["xmin"], bbox["xmax"], length.out = 3)
lat_breaks <- seq(bbox["ymin"], bbox["ymax"], length.out = 3)
centros <- expand.grid(lon = lon_breaks, lat = lat_breaks)

cat("Descargando", nrow(centros), "tiles WorldClim 30s. Puede tardar unos minutos.\n")
tiles <- vector("list", nrow(centros))
for (i in seq_len(nrow(centros))) {
  cat("  Tile", i, "de", nrow(centros), "en (",
      round(centros$lon[i], 1), ",", round(centros$lat[i], 1), ")\n")
  tiles[[i]] <- tryCatch(
    geodata::worldclim_tile(lon = centros$lon[i], lat = centros$lat[i],
                            var = "bio", res = worldclim_res, path = wc_dir),
    error = function(e) { message("  Fallo: ", e$message); NULL })
}
tiles <- Filter(Negate(is.null), tiles)
if (!length(tiles)) stop("No se pudo descargar ningun tile.")

# Mosaicar y recortar
cat("Mosaicando", length(tiles), "tiles...\n")
bio_global <- do.call(terra::mosaic, tiles)
bio_crop <- terra::crop(bio_global,
                        terra::ext(bbox["xmin"], bbox["xmax"],
                                   bbox["ymin"], bbox["ymax"]),
                        snap = "out")
names(bio_crop) <- paste0("bio", seq_len(nlyr(bio_crop)))
cat("Raster mosaico recortado:", ncol(bio_crop), "x", nrow(bio_crop),
    "| variables:", nlyr(bio_crop), "\n")
writeRaster(bio_crop, file.path(out_dir, "bio_crop.tif"), overwrite = TRUE)
# Cargar bio_crop desde disco (por si acaso no está en memoria)
bio_crop <- terra::rast(file.path(out_dir, "bio_crop.tif"))

# ------------------------------------------------------------------------------
# 5. Extraccion y depuracion de predictores
# ------------------------------------------------------------------------------
pts_gbif <- terra::vect(as.data.frame(gbif_thin[, c("decimalLongitude","decimalLatitude")]),
                        geom = c("decimalLongitude","decimalLatitude"),
                        crs = "EPSG:4326")
env_gbif <- terra::extract(bio_crop, pts_gbif, ID = FALSE)

na_por_capa <- sapply(env_gbif, function(x) sum(is.na(x)))
cat("NA por capa:", paste(names(na_por_capa), na_por_capa, sep = "=",
                          collapse = ", "), "\n")

validas <- sapply(env_gbif, function(x)
  sum(!is.na(x)) >= 0.9 * nrow(env_gbif) && var(x, na.rm = TRUE) > 0)
cat("Variables validas:", sum(validas), "de", length(validas), "\n")
if (!any(validas)) stop("Ninguna variable paso el filtro.")

env_gbif <- env_gbif[, validas, drop = FALSE]
bio_crop <- bio_crop[[which(validas)]]

cm <- cor(env_gbif, use = "pairwise.complete.obs")
repeat {
  hi <- which(abs(cm) > 0.7 & upper.tri(cm), arr.ind = TRUE)
  if (!nrow(hi)) break
  col_drop <- hi[1, 2]
  cat("  Eliminando", colnames(env_gbif)[col_drop], "por |r| > 0.7\n")
  env_gbif <- env_gbif[, -col_drop, drop = FALSE]
  bio_crop <- bio_crop[[-col_drop]]
  cm <- cor(env_gbif, use = "pairwise.complete.obs")
}
vars_env <- names(env_gbif)
cat("Variables finales (", length(vars_env), "D): ",
    paste(vars_env, collapse = ", "), "\n", sep = "")
if (length(vars_env) < 3) stop("Menos de 3 variables. No se puede ajustar.")

ok_gbif <- complete.cases(env_gbif)
env_gbif <- env_gbif[ok_gbif, , drop = FALSE]
cat("GBIF sin NA en env:", nrow(env_gbif), "\n")

pts_its <- terra::vect(as.data.frame(occ_its[, c("longitude","latitude")]),
                       geom = c("longitude","latitude"), crs = "EPSG:4326")
env_its <- terra::extract(bio_crop, pts_its, ID = FALSE)

ok_its <- complete.cases(env_its)
env_its <- env_its[ok_its, , drop = FALSE]
occ_its <- occ_its[ok_its, ]
cat("ITS sin NA en env:", nrow(env_its), "\n")

env_sano    <- env_its[occ_its$Estado == "Sano", , drop = FALSE]
env_enfermo <- env_its[occ_its$Estado == "Enfermo", , drop = FALSE]
# ------------------------------------------------------------------------------
# 6. Ajuste de los tres elipsoides (ntbox::cov_center)
# ------------------------------------------------------------------------------
cat("\n=== Seccion 6: ajuste de elipsoides ===\n")

fit_elip <- function(env_data, etiqueta, level = 0.95) {
  df <- as.data.frame(env_data)
  tryCatch(
    ntbox::cov_center(data = df, mve = TRUE, level = level, vars = vars_env),
    error = function(e) {
      message("cov_center fallo para ", etiqueta, ": ", e$message,
              ". Reintentando con level = 0.90.")
      ntbox::cov_center(data = df, mve = TRUE, level = 0.90, vars = vars_env)
    })
}

elip_especie <- fit_elip(env_gbif,   "Especie")
elip_sano    <- fit_elip(env_sano,    "Sano")
elip_enfermo <- fit_elip(env_enfermo, "Enfermo")

cat("\nCentroide especie:\n"); print(round(elip_especie$centroid, 3))
cat("\nCentroide Sano:\n");    print(round(elip_sano$centroid, 3))
cat("\nCentroide Enfermo:\n"); print(round(elip_enfermo$centroid, 3))


# ------------------------------------------------------------------------------
# 7. Solapamiento: Jaccard sobre la grilla ambiental (muestreo 200k celdas)
# ------------------------------------------------------------------------------
cat("\n=== Seccion 7: solapamiento Jaccard ===\n")

jaccard_muestreo <- function(e1, e2, bio_crop, vars_env,
                             n_celdas = 200000, level = 0.95) {
  idx <- match(vars_env, names(bio_crop))
  if (any(is.na(idx))) stop("Alguna variable no esta en bio_crop")
  m <- do.call(cbind, lapply(idx, function(i) values(bio_crop[[i]])))
  colnames(m) <- vars_env
  ok <- complete.cases(m)
  m <- m[ok, , drop = FALSE]
  set.seed(42)
  if (nrow(m) > n_celdas) m <- m[sample(nrow(m), n_celdas), ]
  d1 <- mahalanobis(m, e1$centroid, e1$cov)
  d2 <- mahalanobis(m, e2$centroid, e2$cov)
  thr <- qchisq(level, df = ncol(m))
  in1 <- d1 <= thr
  in2 <- d2 <= thr
  inter <- sum(in1 & in2); union <- sum(in1 | in2)
  list(jaccard = ifelse(union > 0, inter / union, NA_real_),
       inter = inter, union = union, n1 = sum(in1), n2 = sum(in2),
       n_celdas = nrow(m))
}

ov_esp_sano    <- jaccard_muestreo(elip_especie, elip_sano,    bio_crop, vars_env)
ov_esp_enfermo <- jaccard_muestreo(elip_especie, elip_enfermo, bio_crop, vars_env)
ov_sano_enf    <- jaccard_muestreo(elip_sano,    elip_enfermo, bio_crop, vars_env)

cat("\n=== Solapamiento (Jaccard, 200k celdas, 95%) ===\n")
cat("Especie vs Sano    :", round(ov_esp_sano$jaccard, 4),
    " (n1 =", ov_esp_sano$n1, ", n2 =", ov_esp_sano$n2, ")\n")
cat("Especie vs Enfermo :", round(ov_esp_enfermo$jaccard, 4),
    " (n1 =", ov_esp_enfermo$n1, ", n2 =", ov_esp_enfermo$n2, ")\n")
cat("Sano    vs Enfermo :", round(ov_sano_enf$jaccard, 4),
    " (n1 =", ov_sano_enf$n1, ", n2 =", ov_sano_enf$n2, ")\n")

solapamiento_df <- data.frame(
  par       = c("Especie_vs_Sano", "Especie_vs_Enfermo", "Sano_vs_Enfermo"),
  jaccard   = c(ov_esp_sano$jaccard, ov_esp_enfermo$jaccard, ov_sano_enf$jaccard),
  inter     = c(ov_esp_sano$inter,   ov_esp_enfermo$inter,   ov_sano_enf$inter),
  union     = c(ov_esp_sano$union,   ov_esp_enfermo$union,   ov_sano_enf$union),
  n_dentro1 = c(ov_esp_sano$n1,      ov_esp_enfermo$n1,      ov_sano_enf$n1),
  n_dentro2 = c(ov_esp_sano$n2,      ov_esp_enfermo$n2,      ov_sano_enf$n2)
)
write.csv(solapamiento_df,
          file.path(out_dir, "solapamiento_3_modelos.txt"),
          row.names = FALSE)
print(solapamiento_df)

saveRDS(list(
  elip_especie = elip_especie, elip_sano = elip_sano, elip_enfermo = elip_enfermo,
  overlap = solapamiento_df, vars_env = vars_env, bbox = bbox
), file.path(out_dir, "objetos_elipsoides.rds"))

# ------------------------------------------------------------------------------
# 8. Plot con los tres elipsoides en PC1-PC2
# ------------------------------------------------------------------------------
all_env <- rbind(env_gbif, env_sano, env_enfermo)
grupo <- factor(c(rep("Especie", nrow(env_gbif)),
                  rep("Sano",    nrow(env_sano)),
                  rep("Enfermo", nrow(env_enfermo))),
                levels = c("Especie", "Sano", "Enfermo"))

pca <- prcomp(all_env, center = TRUE, scale. = TRUE)
pcs <- predict(pca, all_env)[, 1:2]

max_esp <- 3000
idx_esp <- which(grupo == "Especie")
if (length(idx_esp) > max_esp) {
  set.seed(42); idx_esp <- sample(idx_esp, max_esp)
}
idx_plot <- c(idx_esp, which(grupo != "Especie"))

df_plot <- data.frame(PC1 = pcs[idx_plot, 1], PC2 = pcs[idx_plot, 2],
                      Grupo = grupo[idx_plot])

# Paleta unica: el color de cada registro coincide con el de su elipse
colores <- c(Especie = "grey40", Sano = "forestgreen", Enfermo = "firebrick")

titulo_expr <- expression(
  paste("Elipsoides de nicho climático de ", italic("Quercus robur"))
)

p <- ggplot(df_plot, aes(PC1, PC2, color = Grupo, fill = Grupo)) +
  geom_point(data = subset(df_plot, Grupo == "Especie"),
             alpha = 0.30, size = 0.5) +
  geom_point(data = subset(df_plot, Grupo != "Especie"),
             alpha = 0.90, size = 1.8) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.15,
               show.legend = FALSE) +
  stat_ellipse(level = 0.95, geom = "path", linewidth = 1.1) +
  scale_color_manual(values = colores) +
  scale_fill_manual(values = colores) +
  labs(title = titulo_expr,
       x = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2, 1] * 100),
       y = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2, 2] * 100)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14))

ggsave(file.path(out_dir, "plot_solapamiento_3_modelos.png"),
       p, width = 9, height = 7, dpi = 150)
# ------------------------------------------------------------------------------
# 9. Proyeccion geografica de los tres elipsoides
# ------------------------------------------------------------------------------
# Proyecta cada elipsoide sobre la grilla ambiental. Para cada celda calcula
# la distancia de Mahalanobis y la idoneidad normalizada en [0,1], donde 1
# indica el centro del elipsoide y 0 indica el borde del elipsoide al 95%.
# Se guardan dos rasters por modelo: distancia e idoneidad.

cat("\n=== Seccion 9: proyeccion geografica ===\n")

proyectar_nicho <- function(env_raster, elipsoide, vars_env, nombre, level = 0.95) {
  idx <- match(vars_env, names(env_raster))
  if (any(is.na(idx))) stop("Falta alguna variable en el raster para ", nombre)
  m <- do.call(cbind, lapply(idx, function(i) values(env_raster[[i]])))
  colnames(m) <- vars_env
  
  ok <- complete.cases(m)
  d  <- rep(NA_real_, nrow(m))
  d[ok] <- mahalanobis(m[ok, , drop = FALSE], elipsoide$centroid, elipsoide$cov)
  
  thr <- qchisq(level, df = ncol(m))
  
  dist_ras <- env_raster[[1]]
  values(dist_ras) <- d
  names(dist_ras) <- paste0("dist_", nombre)
  
  idon_ras <- env_raster[[1]]
  values(idon_ras) <- 1 - pmin(d / thr, 1)
  names(idon_ras) <- paste0("idoneidad_", nombre)
  
  list(dist = dist_ras, idoneidad = idon_ras,
       n_validas = sum(ok), n_dentro = sum(!is.na(d) & d <= thr))
}

proj_esp <- proyectar_nicho(bio_crop, elip_especie, vars_env, "especie")
proj_san <- proyectar_nicho(bio_crop, elip_sano,    vars_env, "sano")
proj_enf <- proyectar_nicho(bio_crop, elip_enfermo, vars_env, "enfermo")

writeRaster(proj_esp$dist,      file.path(out_dir, "dist_especie.tif"), overwrite = TRUE)
writeRaster(proj_esp$idoneidad, file.path(out_dir, "idoneidad_especie.tif"), overwrite = TRUE)
writeRaster(proj_san$dist,      file.path(out_dir, "dist_sano.tif"),    overwrite = TRUE)
writeRaster(proj_san$idoneidad, file.path(out_dir, "idoneidad_sano.tif"), overwrite = TRUE)
writeRaster(proj_enf$dist,      file.path(out_dir, "dist_enfermo.tif"), overwrite = TRUE)
writeRaster(proj_enf$idoneidad, file.path(out_dir, "idoneidad_enfermo.tif"), overwrite = TRUE)

cat("Celdas con valores validos:", proj_esp$n_validas, "\n")
cat("Celdas dentro del elipsoide especie:", proj_esp$n_dentro, "\n")
cat("Celdas dentro del elipsoide Sano:   ", proj_san$n_dentro, "\n")
cat("Celdas dentro del elipsoide Enfermo:", proj_enf$n_dentro, "\n")


# ------------------------------------------------------------------------------
# 10. Panel visual de las tres idoneidades
# ------------------------------------------------------------------------------
# Tres mapas lado a lado con la idoneidad de cada modelo, superponiendo las
# ocurrencias correspondientes para verificacion visual.

png(file.path(out_dir, "panel_idoneidad_3_modelos.png"),
    width = 1800, height = 900, res = 110)
par(mfrow = c(1, 3), mar = c(3, 3, 3, 5))

# Especie
plot(proj_esp$idoneidad, main = "Idoneidad - Especie",
     col = hcl.colors(50, "YlGnBu", rev = TRUE),
     xlab = "lon", ylab = "lat")
points(gbif_thin$decimalLongitude, gbif_thin$decimalLatitude,
       pch = 20, cex = 0.1, col = "black")

# Sano
plot(proj_san$idoneidad, main = "Idoneidad - Sano",
     col = hcl.colors(50, "YlGnBu", rev = TRUE),
     xlab = "lon", ylab = "lat")
points(occ_its$longitude[occ_its$Estado == "Sano"],
       occ_its$latitude[occ_its$Estado == "Sano"],
       pch = 20, cex = 0.8, col = "darkgreen")

# Enfermo
plot(proj_enf$idoneidad, main = "Idoneidad - Enfermo",
     col = hcl.colors(50, "YlGnBu", rev = TRUE),
     xlab = "lon", ylab = "lat")
points(occ_its$longitude[occ_its$Estado == "Enfermo"],
       occ_its$latitude[occ_its$Estado == "Enfermo"],
       pch = 20, cex = 0.8, col = "darkred")

dev.off()

cat("\nArchivos generados en", out_dir, ":\n")
print(list.files(out_dir, pattern = "\\.(tif|png|txt|rds)$"))

# Para el panel de Sano
plot(proj_san$idoneidad,
     xlim = c(-3, 2), ylim = c(42, 47),
     main = "Idoneidad - Sano",
     col = hcl.colors(50, "YlGnBu", rev = TRUE))
points(occ_its$longitude[occ_its$Estado == "Sano"],
       occ_its$latitude[occ_its$Estado == "Sano"],
       pch = 20, cex = 0.9, col = "darkgreen")


ext_local <- terra::ext(-3, 2, 42, 47)
proj_san_local <- terra::crop(proj_san$idoneidad, ext_local)
plot(proj_san_local, main = "Idoneidad - Sano (zoom)")


rango <- range(values(proj_san$idoneidad), na.rm = TRUE)
rango <- range(rango[rango > 0])  # ignora ceros

plot(proj_san$idoneidad,
     range = rango,
     col = hcl.colors(50, "YlGnBu", rev = TRUE))

# ------------------------------------------------------------------------------
# 13. Proyeccion a futuro con CMIP6 (SSP245, 2050)
# ------------------------------------------------------------------------------
# Descarga las 19 variables bioclimaticas del modelo ACCESS-CM2 para el
# escenario SSP245 en el periodo 2041-2060. La resolucion es 2.5 arc-min
# porque los GCMs tienen resolucion nativa de ~1 grado y downscalear mas
# fino no añade informacion real. El raster futuro se recorta al mismo bbox
# de estudio y se proyectan los tres elipsoides sobre el.

cat("\n=== Seccion 13: descarga CMIP6 y proyeccion futura ===\n")

cmip6_dir <- file.path(root, "data", "cmip6")
dir.create(cmip6_dir, recursive = TRUE, showWarnings = FALSE)

message("Descargando CMIP6 (ACCESS-CM2, SSP245, 2041-2060). Puede tardar ~10 min.")
bio_2050 <- geodata::cmip6_world(
  model = "ACCESS-CM2",
  ssp   = "245",
  time  = "2041-2060",
  var   = "bio",
  res   = 2.5,
  path  = cmip6_dir
)
if (!inherits(bio_2050, "SpatRaster")) bio_2050 <- terra::rast(bio_2050)
names(bio_2050) <- paste0("bio", seq_len(nlyr(bio_2050)))

bio_2050_crop <- terra::crop(bio_2050,
                             terra::ext(bbox["xmin"], bbox["xmax"],
                                        bbox["ymin"], bbox["ymax"]),
                             snap = "out")
cat("Raster futuro recortado:", ncol(bio_2050_crop), "x", nrow(bio_2050_crop),
    "| variables:", nlyr(bio_2050_crop), "\n")

proj_esp_2050 <- proyectar_nicho(bio_2050_crop, elip_especie, vars_env, "especie_2050")
proj_san_2050 <- proyectar_nicho(bio_2050_crop, elip_sano,    vars_env, "sano_2050")
proj_enf_2050 <- proyectar_nicho(bio_2050_crop, elip_enfermo, vars_env, "enfermo_2050")

writeRaster(proj_esp_2050$idoneidad,
            file.path(out_dir, "idoneidad_especie_2050.tif"), overwrite = TRUE)
writeRaster(proj_san_2050$idoneidad,
            file.path(out_dir, "idoneidad_sano_2050.tif"), overwrite = TRUE)
writeRaster(proj_enf_2050$idoneidad,
            file.path(out_dir, "idoneidad_enfermo_2050.tif"), overwrite = TRUE)

cat("Celdas dentro del elipsoide especie 2050:", proj_esp_2050$n_dentro, "\n")
cat("Celdas dentro del elipsoide Sano 2050:   ", proj_san_2050$n_dentro, "\n")
cat("Celdas dentro del elipsoide Enfermo 2050:", proj_enf_2050$n_dentro, "\n")


# ------------------------------------------------------------------------------
# 14. Remuestreo del presente a la grilla futura
# ------------------------------------------------------------------------------
# Para que los paneles presente y futuro sean comparables celda a celda, los
# rasters presentes (30 arc-sec) se remuestrean a la grilla del futuro
# (2.5 arc-min). Se usa vecino mas cercano porque la idoneidad es una
# probabilidad y la interpolacion bilineal sesgaria los valores.

cat("\n=== Seccion 14: remuestreo del presente a grilla futura ===\n")

template <- bio_2050_crop[[1]]

proj_san_actual <- terra::resample(proj_san$idoneidad, template, method = "near")
proj_enf_actual <- terra::resample(proj_enf$idoneidad, template, method = "near")
proj_esp_actual <- terra::resample(proj_esp$idoneidad, template, method = "near")

cat("Grilla presente remuestreada:", ncol(proj_san_actual), "x",
    nrow(proj_san_actual), "\n")


# ------------------------------------------------------------------------------
# 15. Panel 2x2: redes Sano y Enfermo, presente y futuro
# ------------------------------------------------------------------------------
# A. Red Sano - presente
# B. Red Sano - futuro SSP245 2050
# C. Red Enfermo - presente
# D. Red Enfermo - futuro SSP245 2050
# Los cuatro paneles usan la misma escala de color y el mismo zoom local
# (suroeste de Francia) para que sean comparables directamente.

cat("\n=== Seccion 15: panel 2x2 Sano/Enfermo presente/futuro ===\n")

ext_local <- terra::ext(-3, 2, 42, 47)

crop_na <- function(r, ext) {
  rc <- terra::crop(r, ext)
  rc[rc == 0] <- NA
  rc
}

p_san_act <- crop_na(proj_san_actual, ext_local)
p_san_fut <- crop_na(proj_san_2050$idoneidad, ext_local)
p_enf_act <- crop_na(proj_enf_actual, ext_local)
p_enf_fut <- crop_na(proj_enf_2050$idoneidad, ext_local)

paleta <- hcl.colors(50, "YlGnBu", rev = TRUE)

ruta_panel <- file.path(out_dir, "panel_ENM_presente_futuro.png")
png(ruta_panel, width = 2200, height = 1400, res = 130)
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 4.5), oma = c(0, 0, 2, 0))

plot(p_san_act, main = "A. Red Sano - presente",
     col = paleta, xlab = "lon", ylab = "lat", cex.main = 1.0)
points(occ_its$longitude[occ_its$Estado == "Sano"],
       occ_its$latitude[occ_its$Estado == "Sano"],
       pch = 21, cex = 1.4, bg = "forestgreen", col = "black", lwd = 0.6)

plot(p_san_fut, main = "B. Red Sano - futuro SSP245 2050",
     col = paleta, xlab = "lon", ylab = "lat", cex.main = 1.0)
points(occ_its$longitude[occ_its$Estado == "Sano"],
       occ_its$latitude[occ_its$Estado == "Sano"],
       pch = 21, cex = 1.4, bg = "forestgreen", col = "black", lwd = 0.6)

plot(p_enf_act, main = "C. Red Enfermo - presente",
     col = paleta, xlab = "lon", ylab = "lat", cex.main = 1.0)
points(occ_its$longitude[occ_its$Estado == "Enfermo"],
       occ_its$latitude[occ_its$Estado == "Enfermo"],
       pch = 24, cex = 1.4, bg = "firebrick", col = "black", lwd = 0.6)

plot(p_enf_fut, main = "D. Red Enfermo - futuro SSP245 2050",
     col = paleta, xlab = "lon", ylab = "lat", cex.main = 1.0)
points(occ_its$longitude[occ_its$Estado == "Enfermo"],
       occ_its$latitude[occ_its$Estado == "Enfermo"],
       pch = 24, cex = 1.4, bg = "firebrick", col = "black", lwd = 0.6)

mtext("Nicho climatico de las redes fungicas de Quercus robur: presente y futuro",
      outer = TRUE, cex = 1.1, font = 2)
dev.off()

cat("Panel guardado en:", ruta_panel, "\n")

saveRDS(list(
  proj_esp_2050 = proj_esp_2050,
  proj_san_2050 = proj_san_2050,
  proj_enf_2050 = proj_enf_2050,
  proj_san_actual_remuestreado = proj_san_actual,
  proj_enf_actual_remuestreado = proj_enf_actual,
  escenario = list(model = "ACCESS-CM2", ssp = "245", time = "2041-2060",
                   res = 2.5)
), file.path(out_dir, "objetos_2050.rds"))

cat("\nArchivos generados en", out_dir, ":\n")
print(list.files(out_dir, pattern = "\\.(tif|png|rds)$"))

# ------------------------------------------------------------------------------
# 16. Panel de especie: presente, futuro y cambio proyectado
# ------------------------------------------------------------------------------
cat("\n=== Seccion 16: panel especie presente/futuro/cambio ===\n")

# Remuestrear el presente de especie a la grilla futura (2.5 arc-min)
template <- bio_2050_crop[[1]]
esp_actual_res <- terra::resample(proj_esp$idoneidad, template, method = "bilinear")
esp_futuro     <- proj_esp_2050$idoneidad

# Diferencia: futuro - presente
esp_cambio <- esp_futuro - esp_actual_res

# Crop al bbox completo
ext_full <- terra::ext(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])
esp_actual_c <- terra::crop(esp_actual_res, ext_full)
esp_futuro_c <- terra::crop(esp_futuro,     ext_full)
esp_cambio_c <- terra::crop(esp_cambio,     ext_full)

paleta_idon   <- hcl.colors(50, "YlGnBu", rev = TRUE)
paleta_cambio <- hcl.colors(50, "RdYlBu", rev = FALSE)

lim_cambio <- max(abs(range(values(esp_cambio_c), na.rm = TRUE)), na.rm = TRUE)

ruta_esp <- file.path(out_dir, "panel_especie_presente_futuro_cambio.png")

# Dimensiones grandes para visualizacion: 15 x 5.5 pulgadas a 300 dpi
png(ruta_esp, width = 15, height = 5.5, units = "in", res = 300)
par(mfrow = c(1, 3), mar = c(3.5, 3.5, 1.5, 4.5))

plot(esp_actual_c, main = "Presente WorldClim",
     col = paleta_idon, zlim = c(0, 1),
     xlab = "lon", ylab = "lat", cex.main = 1.1)

plot(esp_futuro_c, main = "Futuro SSP245 2050",
     col = paleta_idon, zlim = c(0, 1),
     xlab = "lon", ylab = "lat", cex.main = 1.1)

plot(esp_cambio_c, main = "Cambio (Futuro - Presente)",
     col = paleta_cambio, zlim = c(-lim_cambio, lim_cambio),
     xlab = "lon", ylab = "lat", cex.main = 1.1)

dev.off()

cat("Panel guardado en:", ruta_esp, "\n")


####
# Limpiar dispositivos graficos colgados
while (dev.cur() > 1) dev.off()

ruta_esp <- file.path(out_dir, "panel_especie_presente_futuro_cambio.png")
if (file.exists(ruta_esp)) file.remove(ruta_esp)

# Reducimos la altura a 4.2 para eliminar el espacio en blanco
ragg::agg_png(ruta_esp, width = 15, height = 4.2, units = "in", res = 300)

# Ajustamos margenes: un poco menos de espacio arriba y abajo
par(mfrow = c(1, 3), mar = c(3, 3, 1, 4))

plot(esp_actual_c, main = "Presente WorldClim",
     col = paleta_idon, zlim = c(0, 1),
     xlab = "lon", ylab = "lat", cex.main = 1.1)

plot(esp_futuro_c, main = "Futuro SSP245 2050",
     col = paleta_idon, zlim = c(0, 1),
     xlab = "lon", ylab = "lat", cex.main = 1.1)

plot(esp_cambio_c, main = "Cambio (Futuro - Presente)",
     col = paleta_cambio, zlim = c(-lim_cambio, lim_cambio),
     xlab = "lon", ylab = "lat", cex.main = 1.1)

dev.off()

cat("PNG generado con altura ajustada:", file.exists(ruta_esp), "\n")