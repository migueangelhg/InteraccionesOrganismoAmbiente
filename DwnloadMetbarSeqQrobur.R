# ============================================================
# DESCARGA UNIFICADA: ITS + LSU de Quercus robur (PRJNA508496)
# Estrategia: curl::multi_download (paralelo interno, sin 'parallel')
# Reanudable + reintentos automáticos
# ============================================================

# --- 0. Paquetes ---------------------------------------------
for (pkg in c("curl", "httr", "readr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(curl); library(httr); library(readr)

# --- 1. Parámetros -------------------------------------------
proyecto     <- "PRJNA508496"
dir_base     <- "fastq_qrobur"
solo_R1      <- FALSE       # TRUE = solo R1 (~1 GB); FALSE = R1+R2 (~2 GB)
max_reintentos <- 3         # reintentos para archivos que fallen

# --- 2. Metadatos del ENA ------------------------------------
campos <- paste("run_accession", "sample_accession", "sample_title",
                "sample_alias", "host", "library_strategy", "read_count",
                "fastq_bytes", "fastq_ftp", sep = ",")
url <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", proyecto,
              "&result=read_run&format=tsv&download=true&fields=", campos)
meta <- read_tsv(I(content(GET(url), as = "text", encoding = "UTF-8")),
                 show_col_types = FALSE)
cat("Runs totales en el proyecto:", nrow(meta), "\n")

# --- 3. Filtrar ITS + LSU de Q. robur ------------------------
its <- meta[grepl("^Quercus robur$", trimws(meta$host), ignore.case = TRUE) &
              grepl("ITS", meta$sample_title, ignore.case = TRUE), ]
its$marcador   <- "ITS"
its$alias_core <- sub("_ITS$", "", its$sample_alias)

otro <- meta[is.na(meta$host), ]
otro$alias_core <- sub("_(ITS|LSU)$", "", otro$sample_alias)
lsu <- otro[otro$alias_core %in% its$alias_core, ]
lsu$marcador <- "LSU"

cols_comunes <- intersect(colnames(its), colnames(lsu))
seleccion <- rbind(its[, cols_comunes], lsu[, cols_comunes])

cat("ITS:", sum(seleccion$marcador == "ITS"),
    "| LSU:", sum(seleccion$marcador == "LSU"),
    "| Total:", nrow(seleccion), "runs\n")
cat("Duplicados:", sum(duplicated(seleccion$run_accession)),
    "| Alias 1:1:", setequal(unique(its$alias_core), unique(lsu$alias_core)), "\n")

# --- 4. Construir tabla plana de descargas -------------------
descargas <- do.call(rbind, lapply(seq_len(nrow(seleccion)), function(i) {
  rutas <- strsplit(seleccion$fastq_ftp[i], ";")[[1]]
  if (solo_R1) rutas <- rutas[1]
  data.frame(
    url     = paste0("https://", rutas),
    destino = file.path(dir_base, seleccion$marcador[i], basename(rutas)),
    stringsAsFactors = FALSE
  )
}))

for (m in c("ITS", "LSU")) {
  dir.create(file.path(dir_base, m), showWarnings = FALSE, recursive = TRUE)
}

ya_ok      <- file.exists(descargas$destino) & file.size(descargas$destino) > 0
pendientes <- descargas[!ya_ok, ]
cat("\nArchivos a descargar:", nrow(pendientes), "de", nrow(descargas),
    "(saltando", sum(ya_ok), "ya existentes)\n")
cat("Tamaño total en disco:", round(sum(file.info(descargas$destino[ya_ok])$size,
                                        na.rm = TRUE)/1024^3, 2), "GB\n\n")

# --- 5. Test: descargar el primer archivo --------------------
ok_test <- FALSE
if (nrow(pendientes) > 0) {
  cat("=== TEST: descargando 1 archivo ===\n")
  test <- multi_download(
    urls     = pendientes$url[1],
    destfile = pendientes$destino[1],
    resume   = TRUE,
    progress = TRUE
  )
  cat("Resultado del test:\n")
  print(test[, c("url", "success", "status_code")])
  ok_test <- isTRUE(test$success[1]) && file.size(pendientes$destino[1]) > 0
} else {
  cat("=== Nada pendiente: todos los archivos ya están descargados ===\n")
  ok_test <- TRUE
}

# --- 6. Descarga completa con reintentos ---------------------
if (nrow(pendientes) > 0 && ok_test) {
  cat("\n=== TEST OK. Iniciando descarga completa ===\n")
  res <- multi_download(
    urls     = pendientes$url,
    destfile = pendientes$destino,
    resume   = TRUE,
    progress = TRUE
  )
  cat("Ronda 1: exitosas", sum(res$success, na.rm = TRUE),
      "de", nrow(res), "\n")
  
  # Reintentos para los que fallaron
  intento <- 1
  while (any(!res$success, na.rm = TRUE) && intento < max_reintentos) {
    intento <- intento + 1
    fallidos <- pendientes[!res$success, , drop = FALSE]
    cat("\n=== Reintento", intento, "para", nrow(fallidos), "archivos ===\n")
    res2 <- multi_download(
      urls     = fallidos$url,
      destfile = fallidos$destino,
      resume   = TRUE,
      progress = TRUE
    )
    # Actualizar estados
    res$success[!res$success][res2$success] <- TRUE
    res$error[!res$success] <- res2$error[res2$success]
    cat("Reintento", intento, ": exitosas",
        sum(res2$success, na.rm = TRUE), "de", nrow(res2), "\n")
  }
  
  cat("\n=== RESUMEN FINAL ===\n")
  cat("Descargas exitosas:", sum(res$success, na.rm = TRUE),
      "de", nrow(res), "\n")
  if (any(!res$success, na.rm = TRUE)) {
    cat("Fallos persistentes:\n")
    print(res[!res$success, c("url", "error")])
  }
} else if (!ok_test) {
  cat("\n=== TEST FALLÓ. No se lanza descarga completa. ===\n")
}

# --- 7. Verificación final ----------------------------------
files <- list.files(dir_base, recursive = TRUE, full.names = TRUE)
cat("\n=== VERIFICACIÓN ===\n")
cat("Archivos en disco:", length(files), "de 252 esperados\n")
cat("Tamaño total:", round(sum(file.info(files)$size, na.rm = TRUE)/1024^3, 2), "GB\n")
cat("Directorio:", normalizePath(dir_base), "\n")

# Resumen por marcador
cat("\n--- Archivos por marcador ---\n")
print(table(sub(paste0("^", dir_base, "/"), "",
                dirname(files))))

# ============================================================
# CORRECCIÓN: verificar, limpiar y re-descargar FASTQ incompletos
# ============================================================

library(curl); library(httr); library(readr)

dir_base <- "fastq_qrobur"

# --- 1. Reconstruir metadatos (igual que antes) --------------
campos <- paste("run_accession", "sample_accession", "sample_title",
                "sample_alias", "host", "library_strategy", "read_count",
                "fastq_bytes", "fastq_ftp", sep = ",")
url <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA508496",
              "&result=read_run&format=tsv&download=true&fields=", campos)
meta <- read_tsv(I(content(GET(url), as = "text", encoding = "UTF-8")),
                 show_col_types = FALSE)

its <- meta[grepl("^Quercus robur$", trimws(meta$host), ignore.case = TRUE) &
              grepl("ITS", meta$sample_title, ignore.case = TRUE), ]
its$marcador   <- "ITS"
its$alias_core <- sub("_ITS$", "", its$sample_alias)

otro <- meta[is.na(meta$host), ]
otro$alias_core <- sub("_(ITS|LSU)$", "", otro$sample_alias)
lsu <- otro[otro$alias_core %in% its$alias_core, ]
lsu$marcador <- "LSU"

seleccion <- rbind(its[, intersect(colnames(its), colnames(lsu))],
                   lsu[, intersect(colnames(its), colnames(lsu))])

# --- 2. Construir tabla con tamaño esperado por archivo ------
descargas <- do.call(rbind, lapply(seq_len(nrow(seleccion)), function(i) {
  rutas   <- strsplit(seleccion$fastq_ftp[i], ";")[[1]]
  tamanos <- as.numeric(strsplit(seleccion$fastq_bytes[i], ";")[[1]])
  data.frame(
    url      = paste0("https://", rutas),
    destino  = file.path(dir_base, seleccion$marcador[i], basename(rutas)),
    esperado = tamanos,
    stringsAsFactors = FALSE
  )
}))

cat("Archivos esperados en total:", nrow(descargas), "\n\n")

# --- 3. Verificar qué archivos están completos ---------------
existe   <- file.exists(descargas$destino)
tam_real <- ifelse(existe, file.size(descargas$destino), 0L)
completo <- existe & !is.na(tam_real) & tam_real >= descargas$esperado * 0.95

cat("=== ESTADO ACTUAL ===\n")
cat("Completos  :", sum(completo), "\n")
cat("Incompletos:", sum(existe & !completo), "\n")
cat("Faltantes  :", sum(!existe), "\n")
cat("Total      :", nrow(descargas), "\n\n")

# --- 4. Borrar los incompletos (0 bytes o truncados) ---------
incompletos <- descargas$destino[existe & !completo]
if (length(incompletos) > 0) {
  cat("Borrando", length(incompletos), "archivos incompletos/corruptos...\n")
  file.remove(incompletos)
  # limpiar variables de estado
  existe   <- file.exists(descargas$destino)
  tam_real <- ifelse(existe, file.size(descargas$destino), 0L)
  completo <- existe & !is.na(tam_real) & tam_real >= descargas$esperado * 0.95
}

# --- 5. Re-descargar en lotes pequeños, sin resume -----------
pendientes <- descargas[!completo, ]
cat("\nArchivos pendientes de descarga:", nrow(pendientes), "\n")

if (nrow(pendientes) > 0) {
  lote <- 20
  total_ok <- 0
  fallidos <- list()
  
  for (inicio in seq(1, nrow(pendientes), by = lote)) {
    fin <- min(inicio + lote - 1, nrow(pendientes))
    idx <- inicio:fin
    cat("\n--- Lote", ceiling(inicio/lote), "de",
        ceiling(nrow(pendientes)/lote), "(", length(idx), "archivos) ---\n")
    
    # Aseguramos que el archivo no exista antes de descargar
    file.remove(pendientes$destino[idx][file.exists(pendientes$destino[idx])])
    
    res <- multi_download(
      urls     = pendientes$url[idx],
      destfile = pendientes$destino[idx],
      resume   = FALSE,   # <--- CLAVE: descarga completa siempre
      progress = TRUE
    )
    
    ok <- res$success & file.exists(pendientes$destino[idx]) &
      file.size(pendientes$destino[idx]) > 0
    total_ok <- total_ok + sum(ok, na.rm = TRUE)
    
    if (any(!ok, na.rm = TRUE)) {
      fallidos[[length(fallidos) + 1]] <- pendientes[idx, ][!ok, , drop = FALSE]
    }
    cat("  OK en este lote:", sum(ok, na.rm = TRUE), "/", length(idx), "\n")
  }
  
  cat("\n=== RESUMEN DE DESCARGA ===\n")
  cat("Descargados correctamente:", total_ok, "de", nrow(pendientes), "\n")
  
  # --- 6. Reintento final de los fallidos --------------------
  if (length(fallidos) > 0) {
    fallidos <- do.call(rbind, fallidos)
    cat("\nReintentando", nrow(fallidos), "archivos fallidos en modo secuencial...\n")
    Sys.sleep(5)  # pausa para no saturar
    for (i in seq_len(nrow(fallidos))) {
      file.remove(fallidos$destino[i][file.exists(fallidos$destino[i])])
      r <- tryCatch(
        multi_download(fallidos$url[i], fallidos$destino[i],
                       resume = FALSE, progress = FALSE),
        error = function(e) NULL
      )
      if (!is.null(r) && isTRUE(r$success[1]) &&
          file.size(fallidos$destino[i]) > 0) {
        cat("  OK:", basename(fallidos$destino[i]), "\n")
      } else {
        cat("  FALLO:", basename(fallidos$destino[i]), "\n")
      }
      Sys.sleep(1)  # 1 segundo entre archivos
    }
  }
}

# --- 7. Verificación final ----------------------------------
cat("\n=== VERIFICACIÓN FINAL ===\n")
files     <- list.files(dir_base, recursive = TRUE, full.names = TRUE)
tam_total <- sum(file.info(files)$size, na.rm = TRUE)
cat("Archivos en disco:", length(files), "de 252 esperados\n")
cat("Tamaño total:", round(tam_total/1024^3, 2), "GB (esperado ~1.24 GB)\n")

# Comprobación archivo por archivo
existe_f  <- file.exists(descargas$destino)
tam_f     <- ifelse(existe_f, file.size(descargas$destino), 0L)
completo_f <- existe_f & tam_f >= descargas$esperado * 0.95
cat("Completos verificados:", sum(completo_f), "de", nrow(descargas), "\n")

if (any(!completo_f)) {
  cat("\nArchivos aún incompletos/faltantes:\n")
  falt <- descargas[!completo_f, ]
  falt$tam_real <- tam_f[!completo_f]
  print(falt[, c("destino", "esperado", "tam_real")], n = 20)
}


