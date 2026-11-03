# ============================================================
# ASIGNACIÓN DE GRUPOS ECOLÓGICOS/FUNCIONALES CON FUNGuild
# Nguyen et al. 2016 — Fungal Ecology 20: 241-248
# ============================================================
# Este script:
#   1. Instala dependencias (phyloseq → Bioconductor; MiscMetabar → CRAN)
#   2. Descarga la base FUNGuild
#   3. Verifica el mapeo de niveles taxonómicos (13 = Genus)
#   4. Extrae los géneros únicos de ITS_taxonomy.csv
#   5. Asigna modo trófico por género
#   6. Genera tabla CSV y gráfico PDF/PNG
# ============================================================
# Requisito: ITS_taxonomy.csv en el directorio de trabajo
# ============================================================

# --- 0. Limpiar objetos viejos -----------------------------------------
rm(list = t = intersect(
  c("generos_eco", "fg_genero", "funguild_db", "fuente_asignacion",
    "generos_plot", "plot_long", "p", "generos_df", "tax",
    "n_con", "presentes", "testigo"),
  ls(envir = .GlobalEnv)
), envir = .GlobalEnv)

# --- 1. Instalación de dependencias ------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  message("Instalando phyloseq desde Bioconductor (binarios, sin Rtools)...")
  BiocManager::install("phyloseq", update = FALSE, ask = FALSE)
}

if (!requireNamespace("MiscMetabar", quietly = TRUE)) {
  message("Instalando MiscMetabar desde CRAN...")
  install.packages("MiscMetabar", repos = "https://cloud.r-project.org")
}

stopifnot(
  "phyloseq no disponible"    = requireNamespace("phyloseq",    quietly = TRUE),
  "MiscMetabar no disponible" = requireNamespace("MiscMetabar", quietly = TRUE)
)

suppressPackageStartupMessages({
  library(MiscMetabar)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# --- 2. Descargar FUNGuild ---------------------------------------------
message("\nDescargando FUNGuild DB desde el servidor oficial...")
funguild_db <- MiscMetabar::get_funguild_db()

stopifnot(
  "FUNGuild vacía o incompleta" = nrow(funguild_db) > 1000,
  "Falta columna taxon"         = "taxon"          %in% names(funguild_db),
  "Falta columna taxonomicLevel"= "taxonomicLevel" %in% names(funguild_db),
  "Falta columna trophicMode"   = "trophicMode"    %in% names(funguild_db)
)

cat("Registros en FUNGuild DB:", nrow(funguild_db), "\n")
cat("Columnas:", paste(names(funguild_db), collapse = ", "), "\n")
saveRDS(funguild_db, "funguild_db.rds")

# --- 3. Verificar que 13 = Genus ---------------------------------------
cat("\n=== DIAGNÓSTICO DE FUNGuild ===\n")
cat("\nCódigos presentes en taxonomicLevel:\n")
print(table(funguild_db$taxonomicLevel, useNA = "ifany"))

cat("\nEjemplos en nivel 13 (deberían ser géneros conocidos):\n")
print(head(funguild_db$taxon[funguild_db$taxonomicLevel == 13], 10))

testigo   <- c("Aspergillus", "Penicillium", "Cladosporium", "Fusarium")
presentes <- sum(testigo %in% funguild_db$taxon[funguild_db$taxonomicLevel == 13])
cat("\nGéneros testigo encontrados en nivel 13:", presentes, "de", length(testigo), "\n")

if (presentes < 2) {
  stop("El nivel 13 no corresponde a Genus en esta versión de FUNGuild.")
}

# --- 4. Géneros únicos detectados en el estudio ------------------------
if (!file.exists("ITS_taxonomy.csv")) {
  stop("Falta ITS_taxonomy.csv en ", getwd(),
       ". Ejecuta primero el pipeline DADA2 + UNITE.")
}

tax <- readr::read_csv("ITS_taxonomy.csv", show_col_types = FALSE)

generos_df <- tax |>
  dplyr::filter(!is.na(Genus), Genus != "") |>
  dplyr::mutate(Genus = sub("^g__", "", Genus)) |>
  dplyr::count(Genus, name = "n_ASVs", sort = TRUE)

cat("\nGéneros únicos ITS:", nrow(generos_df), "\n")

# --- 5. Filtrar FUNGuild a nivel de Género -----------------------------
fg_genero <- funguild_db |>
  dplyr::filter(taxonomicLevel == 13) |>
  dplyr::filter(!is.na(taxon), taxon != "") |>
  dplyr::select(
    Genus       = taxon,
    TrophicMode = trophicMode,
    Guild       = guild,
    Confidence  = confidenceRanking
  ) |>
  dplyr::distinct(Genus, .keep_all = TRUE)

cat("Géneros únicos en FUNGuild:", nrow(fg_genero), "\n")

if (nrow(fg_genero) < 100) {
  stop("FUNGuild con muy pocos géneros (", nrow(fg_genero), ").")
}

# --- 6. Asignación ecológica por género --------------------------------
generos_eco <- generos_df |>
  dplyr::left_join(fg_genero, by = "Genus") |>
  dplyr::mutate(
    Grupo_ecologico = dplyr::case_when(
      is.na(TrophicMode)                               ~ "Sin asignar",
      grepl("Pathotroph",  TrophicMode) &
        !grepl("Saprotroph|Symbiotroph", TrophicMode)  ~ "Patótrofo",
      grepl("Saprotroph",  TrophicMode) &
        !grepl("Pathotroph|Symbiotroph", TrophicMode)  ~ "Saprótrofo",
      grepl("Simbiotroph", TrophicMode) &
        !grepl("Pathotroph|Saprotroph",  TrophicMode)  ~ "Simbiotrofo",
      TRUE                                             ~ "Múltiple"
    )
  )

n_con <- sum(generos_eco$Grupo_ecologico != "Sin asignar")
cat("\nGéneros con asignación FUNGuild:", n_con, "de", nrow(generos_eco),
    "(", round(100 * n_con / nrow(generos_eco), 1), "%)\n\n")
print(table(generos_eco$Grupo_ecologico))

readr::write_csv(generos_eco, "ITS_generos_FUNGuild.csv")

generos_eco <-readr::read_csv("ITS_generos_FUNGuild.csv")



# --- 7. Resumen por grupo ----------------------------------------------
generos_plot <- generos_eco |>
  dplyr::count(Grupo_ecologico, wt = n_ASVs, name = "n_ASVs_total") |>
  dplyr::full_join(
    generos_eco |> dplyr::count(Grupo_ecologico, name = "n_generos"),
    by = "Grupo_ecologico"
  ) |>
  tidyr::replace_na(list(n_ASVs_total = 0, n_generos = 0)) |>
  dplyr::mutate(
    pct_generos = round(100 * n_generos     / sum(n_generos), 1),
    pct_ASVs    = round(100 * n_ASVs_total / sum(n_ASVs_total), 1)
  ) |>
  dplyr::arrange(dplyr::desc(n_generos))

print(generos_plot)

# --- 8. Gráfico de proporciones ----------------------------------------
plot_long <- generos_plot |>
  dplyr::select(Grupo_ecologico, n_generos, n_ASVs_total) |>
  tidyr::pivot_longer(
    cols = c(n_generos, n_ASVs_total),
    names_to = "Unidad", values_to = "N"
  ) |>
  dplyr::mutate(Unidad = dplyr::recode(Unidad,
                                       n_generos    = "Géneros",
                                       n_ASVs_total = "ASVs"
  )) |>
  dplyr::group_by(Unidad) |>
  dplyr::mutate(pct = N / sum(N) * 100) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    generos_plot |> dplyr::select(Grupo_ecologico, n_generos),
    by = "Grupo_ecologico"
  ) |>
  dplyr::mutate(
    Grupo_ecologico = paste0(Grupo_ecologico, " (n = ", n_generos, ")")
  )

p <- ggplot2::ggplot(
  plot_long,
  ggplot2::aes(x = Unidad, y = N, fill = Grupo_ecologico)
) +
  ggplot2::geom_col(width = 0.55, color = "white", linewidth = 0.3) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(pct >= 3, paste0(round(pct, 1), "%"), "")),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 3.3, color = "white", fontface = "bold"
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::labs(
    title    = "Asignación ecológica de los géneros fúngicos detectados",
    subtitle = sprintf(
      "Riqueza (nº de géneros, total = %d) vs. abundancia (nº de ASVs)",
      nrow(generos_eco)
    ),
    x        = NULL, y = "Número", fill = "Grupo ecológico",
    caption  = sprintf(
      "Fuente: FUNGuild (Nguyen et al. 2016). %d de %d géneros del estudio asignados (%.1f%%)",
      n_con, nrow(generos_eco), 100 * n_con / nrow(generos_eco)
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position    = "right",
    panel.grid.major.y = ggplot2::element_blank(),
    plot.title         = ggplot2::element_text(face = "bold")
  )

print(p)
ggplot2::ggsave("ITS_grupos_ecologicos.pdf", p, width = 10, height = 5)
ggplot2::ggsave("ITS_grupos_ecologicos.png", p, width = 10, height = 5, dpi = 200)

# --- 9. Resumen final --------------------------------------------------
cat("\n======================================================\n")
cat("LISTO. Archivos generados:\n",
    "  - funguild_db.rds           (", nrow(funguild_db), " registros)\n",
    "  - ITS_generos_FUNGuild.csv  (", n_con, " de ", nrow(generos_eco),
    " géneros asignados, ", round(100 * n_con / nrow(generos_eco), 1), "%)\n",
    "  - ITS_grupos_ecologicos.pdf\n",
    "  - ITS_grupos_ecologicos.png\n",
    "======================================================\n", sep = "")


# ============================================================
# EXTENSIÓN: Asignación ecológica por GREMIO (columna 'guild')
# Usa funguild_db, generos_df, tax y n_con ya existentes en memoria
# Produce:
#   - ITS_generos_FUNGuild_guild.csv  (gremios finos + categorías amplias)
#   - ITS_grupos_ecologicos_guild.pdf/.png  (gráfico con ~10 categorías)
# ============================================================

# --- 1. Verificar que los objetos base están disponibles -------------
stopifnot(
  "funguild_db no está en memoria" = exists("funguild_db", envir = .GlobalEnv),
  "generos_df no está en memoria"  = exists("generos_df",  envir = .GlobalEnv),
  "tax no está en memoria"         = exists("tax",         envir = .GlobalEnv)
)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# --- 2. Extraer gremios desde FUNGuild -------------------------------
# Filtramos a nivel Genus y desdoblamos 'guild' cuando viene compuesto
# con guiones (ej. "Plant Pathogen-Wood Saprotroph" → 2 filas)
fg_gremio <- funguild_db |>
  dplyr::filter(taxonomicLevel == 13) |>
  dplyr::filter(!is.na(taxon), taxon != "") |>
  dplyr::filter(!is.na(guild), guild != "") |>
  dplyr::transmute(
    Genus      = taxon,
    Guild_raw  = guild,
    Confidence = confidenceRanking
  ) |>
  tidyr::separate_rows(Guild_raw, sep = "[-/]") |>
  dplyr::mutate(Guild_raw = trimws(Guild_raw)) |>
  dplyr::filter(Guild_raw != "")

cat("Filas por gremio (desdobladas):", nrow(fg_gremio), "\n")
cat("Gremios únicos en FUNGuild:", length(unique(fg_gremio$Guild_raw)), "\n")

# --- 3. Mapear ~100 gremios específicos → ~10 categorías amplias ----
mapear_gremio <- function(g) {
  g <- tolower(g)
  dplyr::case_when(
    grepl("endophyt",                g) ~ "Endófito",
    grepl("ectomycorrhiz|ericoid|arbuscular|mycorrhiz",
          g) ~ "Micorrícico",
    grepl("plant pathogen",          g) ~ "Patógeno de planta",
    grepl("animal pathogen|fungal parasit|lichen parasit",
          g) ~ "Patógeno (animal/fúngico/liquen)",
    grepl("wood saprotroph",         g) ~ "Saprótrofo de madera",
    grepl("soil saprotroph|litter saprotroph|dung saprotroph|undefined saprotroph",
          g) ~ "Saprótrofo de suelo/hojarasca",
    grepl("lichenized|lichen",       g) ~ "Liquen",
    grepl("epiphyt",                 g) ~ "Epífito",
    grepl("pollen saprotroph|nectar|root", g) ~ "Otro saprótrofo",
    TRUE                                 ~ "Otro / No clasificado"
  )
}

fg_gremio <- fg_gremio |>
  dplyr::mutate(Categoria = mapear_gremio(Guild_raw))

cat("\nCategorías amplias asignadas en FUNGuild:\n")
print(table(fg_gremio$Categoria))

# --- 4. Resolver múltiples gremios por género -----------------------
# Un género puede tener varios gremios. Aplicamos una jerarquía:
# Endófito > Patógeno de planta > Patógeno animal > Micorrícico >
# Saprótrofo de madera > Saprótrofo de suelo > Liquen > Epífito > Otro.
# Los gremios completos se guardan aparte para la tabla suplementaria.
jerarquia <- c(
  "Endófito",
  "Patógeno de planta",
  "Patógeno (animal/fúngico/liquen)",
  "Micorrícico",
  "Saprótrofo de madera",
  "Saprótrofo de suelo/hojarasca",
  "Liquen",
  "Epífito",
  "Otro saprótrofo",
  "Otro / No clasificado"
)

fg_colapsado <- fg_gremio |>
  dplyr::group_by(Genus) |>
  dplyr::summarise(
    Categoria  = dplyr::first(jerarquia[jerarquia %in% Categoria]),
    Gremios    = paste(sort(unique(Guild_raw)), collapse = "; "),
    Confidence = paste(sort(unique(Confidence)), collapse = "; "),
    .groups    = "drop"
  )

cat("\nGéneros únicos en FUNGuild tras colapsar:", nrow(fg_colapsado), "\n")

# --- 5. Unir con los géneros del estudio ----------------------------
generos_eco_guild <- generos_df |>
  dplyr::left_join(fg_colapsado, by = "Genus") |>
  dplyr::mutate(
    Categoria = dplyr::coalesce(Categoria, "Sin asignar")
  )

n_con_guild <- sum(generos_eco_guild$Categoria != "Sin asignar")
cat("\nGéneros con asignación FUNGuild (guild):",
    n_con_guild, "de", nrow(generos_eco_guild),
    "(", round(100 * n_con_guild / nrow(generos_eco_guild), 1), "%)\n\n")

print(table(generos_eco_guild$Categoria))

readr::write_csv(generos_eco_guild, "ITS_generos_FUNGuild_guild.csv")

# --- 6. Resumen por categoría (guild) -------------------------------
generos_plot_guild <- generos_eco_guild |>
  dplyr::count(Categoria, wt = n_ASVs, name = "n_ASVs_total") |>
  dplyr::full_join(
    generos_eco_guild |> dplyr::count(Categoria, name = "n_generos"),
    by = "Categoria"
  ) |>
  tidyr::replace_na(list(n_ASVs_total = 0, n_generos = 0)) |>
  dplyr::mutate(
    pct_generos = round(100 * n_generos     / sum(n_generos), 1),
    pct_ASVs    = round(100 * n_ASVs_total / sum(n_ASVs_total), 1)
  ) |>
  dplyr::arrange(dplyr::desc(n_generos))

print(generos_plot_guild)

# --- 7. Gráfico por gremio (categoría amplia) -----------------------
plot_long_guild <- generos_plot_guild |>
  dplyr::select(Categoria, n_generos, n_ASVs_total) |>
  tidyr::pivot_longer(
    cols = c(n_generos, n_ASVs_total),
    names_to = "Unidad", values_to = "N"
  ) |>
  dplyr::mutate(Unidad = dplyr::recode(Unidad,
                                       n_generos    = "Géneros",
                                       n_ASVs_total = "ASVs"
  )) |>
  dplyr::group_by(Unidad) |>
  dplyr::mutate(pct = N / sum(N) * 100) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    generos_plot_guild |> dplyr::select(Categoria, n_generos),
    by = "Categoria"
  ) |>
  dplyr::mutate(
    Categoria = paste0(Categoria, " (n = ", n_generos, ")")
  )

# Mantener el orden del resumen en la leyenda
plot_long_guild$Categoria <- factor(
  plot_long_guild$Categoria,
  levels = paste0(generos_plot_guild$Categoria,
                  " (n = ", generos_plot_guild$n_generos, ")")
)

p_guild <- ggplot2::ggplot(
  plot_long_guild,
  ggplot2::aes(x = Unidad, y = N, fill = Categoria)
) +
  ggplot2::geom_col(width = 0.55, color = "white", linewidth = 0.3) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(pct >= 3, paste0(round(pct, 1), "%"), "")),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 3.3, color = "white", fontface = "bold"
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_brewer(palette = "Set3") +
  ggplot2::labs(
    title    = "Asignación ecológica por gremio (FUNGuild)",
    subtitle = sprintf(
      "Riqueza (nº de géneros, total = %d) vs. abundancia (nº de ASVs)",
      nrow(generos_eco_guild)
    ),
    x        = NULL, y = "Número", fill = "Categoría ecológica",
    caption  = sprintf(
      "Fuente: FUNGuild (Nguyen et al. 2016), columna 'guild'. %d de %d géneros asignados (%.1f%%)",
      n_con_guild, nrow(generos_eco_guild),
      100 * n_con_guild / nrow(generos_eco_guild)
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position    = "right",
    panel.grid.major.y = ggplot2::element_blank(),
    plot.title         = ggplot2::element_text(face = "bold")
  )

print(p_guild)
ggplot2::ggsave("ITS_grupos_ecologicos_guild.pdf",
                p_guild, width = 11, height = 5.5)
ggplot2::ggsave("ITS_grupos_ecologicos_guild.png",
                p_guild, width = 11, height = 5.5, dpi = 200)

# --- 8. Resumen final -----------------------------------------------
cat("\n======================================================\n")
cat("EXTENSIÓN 'guild' COMPLETADA.\n",
    "  - ITS_generos_FUNGuild_guild.csv (",
    n_con_guild, " de ", nrow(generos_eco_guild), " géneros asignados, ",
    round(100 * n_con_guild / nrow(generos_eco_guild), 1), "%)\n",
    "  - ITS_grupos_ecologicos_guild.pdf\n",
    "  - ITS_grupos_ecologicos_guild.png\n",
    "\nObjetos nuevos en memoria:\n",
    "  - fg_gremio          (FUNGuild desdoblada por gremio)\n",
    "  - fg_colapsado       (un género = una categoría + gremios finos)\n",
    "  - generos_eco_guild  (tabla final del estudio)\n",
    "  - generos_plot_guild (resumen por categoría)\n",
    "  - p_guild            (objeto ggplot)\n",
    "======================================================\n", sep = "")

# --- 7b. Gráficos de pastel por gremio (etiquetas externas sin cruces) ---
for (pkg in c("ggforce", "RColorBrewer")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(ggforce)
  library(RColorBrewer)
})

# --- 7b.1 Datos largos ------------------------------------------------
pie_data <- generos_plot_guild |>
  dplyr::select(Categoria, n_generos, n_ASVs_total) |>
  tidyr::pivot_longer(
    cols = c(n_generos, n_ASVs_total),
    names_to = "Unidad", values_to = "N"
  ) |>
  dplyr::mutate(Unidad = dplyr::recode(Unidad,
                                       n_generos    = "Riqueza de géneros",
                                       n_ASVs_total = "Abundancia de ASVs"
  )) |>
  dplyr::group_by(Unidad) |>
  dplyr::mutate(
    pct      = round(100 * N / sum(N), 1),
    prop     = N / sum(N),
    end      = cumsum(prop),
    start    = c(0, head(end, -1)),
    mid      = (start + end) / 2,
    mid_rad  = mid * 2 * pi,
    x_inner  = 0.62 * sin(mid_rad),
    y_inner  = 0.62 * cos(mid_rad),
    x_anchor = 1.02 * sin(mid_rad),
    y_anchor = 1.02 * cos(mid_rad)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Categoria = factor(Categoria, levels = generos_plot_guild$Categoria)
  )

# --- 7b.2 Paleta ------------------------------------------------------
n_cat   <- length(unique(pie_data$Categoria))
colores <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_cat)
names(colores) <- levels(pie_data$Categoria)

umbral <- 5

# --- 7b.3 Etiquetas externas con posiciones precalculadas ------------
# Todos los labels externos llevan SOLO el porcentaje entre paréntesis.
small_data <- pie_data |>
  dplyr::filter(pct < umbral) |>
  dplyr::group_by(Unidad) |>
  dplyr::arrange(mid_rad) |>
  dplyr::mutate(
    n_in_panel = dplyr::n(),
    idx        = dplyr::row_number(),
    x_label = dplyr::case_when(
      n_in_panel == 1 ~ 0,
      n_in_panel == 2 ~ c(-0.7, 0.7)[idx],
      n_in_panel == 3 ~ c(-1.1, 0, 1.1)[idx],
      n_in_panel == 4 ~ c(-1.3, -0.45, 0.45, 1.3)[idx],
      TRUE ~ seq(-1.3, 1.3, length.out = n_in_panel)[idx]
    ),
    y_label   = 1.45,
    label_txt = paste0("(", pct, "%)")
  ) |>
  dplyr::ungroup()
# --- 7b.4 Gráfico -----------------------------------------------------
p_guild <- ggplot2::ggplot(pie_data) +
  ggforce::geom_arc_bar(
    ggplot2::aes(
      x0 = 0, y0 = 0, r0 = 0, r = 1,
      start = start * 2 * pi, end = end * 2 * pi,
      fill = Categoria
    ),
    color = "white", linewidth = 0.4
  ) +
  ggplot2::geom_text(
    data = pie_data |> dplyr::filter(pct >= umbral),
    ggplot2::aes(x = x_inner, y = y_inner, label = paste0(pct, "%")),
    size = 3.3, fontface = "bold", color = "black"
  ) +
  ggplot2::geom_segment(
    data = small_data,
    ggplot2::aes(x = x_anchor, y = y_anchor,
                 xend = x_label, yend = y_label),
    color = "gray40", linewidth = 0.4
  ) +
  ggplot2::geom_text(
    data = small_data,
    ggplot2::aes(x = x_label, y = y_label, label = label_txt),
    size = 3, color = "gray20", vjust = -0.3
  ) +
  ggplot2::facet_wrap(~ Unidad, ncol = 2) +
  ggplot2::coord_fixed(clip = "off",
                       xlim = c(-1.8, 1.8),
                       ylim = c(-1.2, 1.9)) +
  ggplot2::scale_fill_manual(values = colores) +
  ggplot2::labs(
    title = "Asignación ecológica por gremio (FUNGuild)",
    fill  = "Categoría ecológica"
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    legend.position = "right",
    legend.text     = ggplot2::element_text(size = 9),
    legend.title    = ggplot2::element_text(face = "bold", size = 10),
    # --- Título con margen inferior amplio ---------------------------
    plot.title      = ggplot2::element_text(
      face  = "bold",
      hjust = 0.5,
      size  = 13,
      margin = ggplot2::margin(t = 5, b = 25)
    ),
    # --- Etiquetas de faceta con margen arriba y abajo ---------------
    strip.text      = ggplot2::element_text(
      face   = "bold",
      size   = 12,
      margin = ggplot2::margin(t = 8, b = 10)
    ),
    plot.margin     = ggplot2::margin(15, 15, 10, 15)
  )# --- 7b.5 Guardar -----------------------------------------------------
ggplot2::ggsave("ITS_grupos_ecologicos_guildPie.pdf",
                p_guild, width = 13, height = 6.5, device = "pdf")
ggplot2::ggsave("ITS_grupos_ecologicos_guildPie.png",
                p_guild, width = 13, height = 6.5, dpi = 200, device = "png")
