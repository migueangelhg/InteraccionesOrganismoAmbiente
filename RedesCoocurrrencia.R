# =============================================================================
# ANÁLISIS COMPLETO DE COMUNIDAD FÚNGICA: SANO vs ENFERMO
# Script único autocontenido
# Redes: pres=6, rho=0.40, dz=0.60, CLR pc=0.5
# =============================================================================

rm(list = ls()); gc(); set.seed(2026)

# --- 1. PAQUETES ------------------------------------------------------------
req <- c("tidyverse","vegan","betapart","igraph","ggraph","tidygraph","patchwork")
faltan <- req[!sapply(req, requireNamespace, quietly = TRUE)]
if (length(faltan)) install.packages(faltan)
suppressPackageStartupMessages({
  for (p in req) library(p, character.only = TRUE)
})

# --- 2. CARGA ---------------------------------------------------------------
asv_raw  <- readr::read_csv("ITS_ASV_table.csv",       show_col_types = FALSE)
metadata <- readr::read_csv("ITS_sample_metadata.csv", show_col_types = FALSE)
taxo     <- readr::read_csv("ITS_taxonomy.csv",        show_col_types = FALSE) %>%
  dplyr::transmute(ASV_full = ASV, Genus_tax = sub("^g__", "", Genus))
guild_map <- readr::read_csv("ITS_generos_FUNGuild_guild.csv", show_col_types = FALSE) %>%
  dplyr::select(Genus, Categoria)

cat("=== ESTRUCTURA ===\n")
cat("ASVs:", nrow(asv_raw), "| Muestras:", nrow(metadata), "\n")
print(table(metadata$Estado))

asv_raw <- asv_raw %>%
  dplyr::mutate(ASV_short = sprintf("ASV_%04d", dplyr::row_number()),
                ASV_full  = Sample)

# --- 3. MATRICES ------------------------------------------------------------
mat <- asv_raw %>%
  dplyr::select(-Sample, -ASV_full) %>%
  tibble::column_to_rownames("ASV_short") %>%
  as.matrix() %>% t()   # muestra × ASV

metadata <- metadata[match(rownames(mat), metadata$Sample), ]
stopifnot(!any(is.na(metadata$Sample)),
          all(metadata$Estado %in% c("Sano","Enfermo")))

meta_all <- metadata %>%
  dplyr::transmute(SampleID  = Sample,
                   Condicion = factor(Estado, levels = c("Sano","Enfermo")),
                   Site      = Site_ena)

idx_s <- which(meta_all$Condicion == "Sano")
idx_e <- which(meta_all$Condicion == "Enfermo")
cat(sprintf("\nMuestras: %d Sano | %d Enfermo\n", length(idx_s), length(idx_e)))

rel <- vegan::decostand(mat, "total")
pa  <- (mat > 0) * 1

# --- 4. DIVERSIDAD ALFA ----------------------------------------------------
alpha_all <- tibble::tibble(
  SampleID  = rownames(mat),
  Condicion = meta_all$Condicion,
  Shannon   = vegan::diversity(mat, "shannon"),
  Simpson   = vegan::diversity(mat, "simpson"),
  Riqueza   = vegan::specnumber(mat),
  Chao1     = vegan::estimateR(mat)["S.chao1", ]
)

cat("\n=== ALFA (media ± sd) ===\n")
print(alpha_all %>% dplyr::group_by(Condicion) %>%
        dplyr::summarise(dplyr::across(where(is.numeric),
                                       list(m = mean, s = sd)),
                         .groups = "drop"))

tests_alpha <- purrr::map_dfr(c("Shannon","Simpson","Riqueza","Chao1"),
                              function(m) {
                                f <- wilcox.test(alpha_all[[m]] ~ alpha_all$Condicion)
                                tibble::tibble(Indice = m, W = unname(f$statistic), p = f$p.value,
                                               sig = symnum(f$p.value,
                                                            cutpoints = c(0, .001, .01, .05, 1),
                                                            symbols = c("***","**","*","ns")))
                              })
cat("\n=== Wilcoxon alfa ===\n"); print(tests_alpha)

# --- 5. BETA DIVERSIDAD (Baselga) ------------------------------------------
pa_s <- pa[idx_s, , drop = FALSE]; pa_e <- pa[idx_e, , drop = FALSE]
bs   <- betapart::beta.pair(pa_s); be <- betapart::beta.pair(pa_e)
ba   <- betapart::beta.pair(rbind(pa_s, pa_e))
n_s  <- nrow(pa_s); n_e <- nrow(pa_e)
inter <- function(B) mean(as.matrix(B)[1:n_s, (n_s+1):(n_s+n_e)])

beta_resumen <- tibble::tibble(
  Condicion  = c("Sano","Enfermo","Sano↔Enfermo"),
  Beta_Sor   = c(mean(bs$beta.sor), mean(be$beta.sor), inter(ba$beta.sor)),
  Turnover   = c(mean(bs$beta.sim), mean(be$beta.sim), inter(ba$beta.sim)),
  Nestedness = c(mean(bs$beta.sne), mean(be$beta.sne), inter(ba$beta.sne))
) %>% dplyr::mutate(pct_turnover   = Turnover   / Beta_Sor,
                    pct_nestedness = Nestedness / Beta_Sor)
cat("\n=== β diversidad (Sørensen) ===\n"); print(beta_resumen)

# --- 6. PERMANOVA + BETADISPER ---------------------------------------------
dist_bray <- vegan::vegdist(mat, "bray")
perm <- vegan::adonis2(dist_bray ~ Condicion + Site, data = meta_all,
                       permutations = 999)
cat("\n=== PERMANOVA ===\n"); print(perm)

bd <- vegan::betadisper(dist_bray, meta_all$Condicion)
cat("\n=== Betadisper ===\n"); print(anova(bd))

# --- 7. NODF CON MODELO NULO -----------------------------------------------
nodf_stat <- function(x) unname(vegan::nestednodf(x)$statistic["NODF"])
null_s <- vegan::oecosimu(pa_s, nodf_stat, "swap", nsimul = 99)
null_e <- vegan::oecosimu(pa_e, nodf_stat, "swap", nsimul = 99)
cat(sprintf("\n=== NODF ===\nSano    : obs=%.3f SES=%.2f p=%.3f\n",
            null_s$oecosimu$statistic, null_s$oecosimu$z, null_s$oecosimu$pval))
cat(sprintf("Enfermo : obs=%.3f SES=%.2f p=%.3f\n",
            null_e$oecosimu$statistic, null_e$oecosimu$z, null_e$oecosimu$pval))

# =============================================================================
# 8. REDES DE CO-OCURRENCIA
# =============================================================================

# --- 8.1 CLR con pseudocount -----------------------------------------------
clr_transform <- function(m, pc = 0.5) {
  x <- m + pc
  lx <- log(x)
  sweep(lx, 1, rowMeans(lx), "-")
}

# --- 8.2 Agregar ASVs a género en proporciones crudas ----------------------
agg_genero <- function(rel_mat, asv_raw, taxo) {
  as.data.frame(rel_mat) %>%
    tibble::rownames_to_column("SampleID") %>%
    tidyr::pivot_longer(-SampleID, names_to = "ASV_short", values_to = "abund") %>%
    dplyr::left_join(asv_raw %>% dplyr::select(ASV_short, ASV_full),
                     by = "ASV_short") %>%
    dplyr::left_join(taxo, by = "ASV_full") %>%
    dplyr::filter(!is.na(Genus_tax), Genus_tax != "") %>%
    dplyr::group_by(SampleID, Genus_tax) %>%
    dplyr::summarise(abund = sum(abund), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Genus_tax, values_from = abund,
                       values_fill = 0) %>%
    tibble::column_to_rownames("SampleID") %>% as.matrix()
}

gen_all <- agg_genero(rel, asv_raw, taxo)
cat(sprintf("\nMatriz de género: %d muestras × %d géneros\n",
            nrow(gen_all), ncol(gen_all)))

# --- 8.3 Distribución de presencia ----------------------------------------
pres_all <- colSums(gen_all > 0)
print(table(cut(pres_all,
                breaks = c(0, 2, 5, 10, 20, 27, 36, 63),
                labels = c("1-2","3-5","6-10","11-20","21-27","28-36","37-63"))))

# --- 8.4 Filtro (presencia ≥ 6, como decidimos) ---------------------------
keep  <- colSums(gen_all > 0) >= 6
sub   <- gen_all[, keep, drop = FALSE]
var_s <- apply(sub[idx_s, , drop = FALSE], 2, var)
var_e <- apply(sub[idx_e, , drop = FALSE], 2, var)
sub   <- sub[, var_s > 0 & var_e > 0, drop = FALSE]

gen_s <- sub[idx_s, , drop = FALSE]
gen_e <- sub[idx_e, , drop = FALSE]
cat(sprintf("Géneros retenidos: %d (varianza > 0 en ambas)\n", ncol(sub)))

# --- 8.5 Constructor de red con CLR + Pearson + dz-avoidance --------------
construir_red_genero <- function(gen_mat, rho_min = 0.40, dz_max = 0.60,
                                 clr_pc = 0.5, aplicar_dz = TRUE) {
  clr <- clr_transform(gen_mat, pc = clr_pc)
  r   <- cor(clr, method = "pearson")
  
  mask_dz <- if (aplicar_dz) {
    ceros <- (gen_mat == 0) * 1
    frac  <- (t(ceros) %*% ceros) / nrow(gen_mat)
    m <- frac > dz_max; diag(m) <- FALSE; m
  } else matrix(FALSE, nrow(r), ncol(r))
  
  adj <- (abs(r) >= rho_min) & !mask_dz
  diag(adj) <- FALSE
  
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  tri <- upper.tri(adj) & adj
  igraph::E(g)$weight <- abs(r)[tri]
  igraph::E(g)$signo  <- sign(r)[tri]
  g
}

red_s <- construir_red_genero(gen_s, rho_min = 0.40, dz_max = 0.60)
red_e <- construir_red_genero(gen_e, rho_min = 0.40, dz_max = 0.60)

cat(sprintf("\n>>> REDES (pres=6, rho=0.40, dz=0.60) <<<\n"))
cat(sprintf("Sano    : %d nodos, %d aristas (%.1f%% pos, %.1f%% neg)\n",
            igraph::vcount(red_s), igraph::ecount(red_s),
            100*sum(igraph::E(red_s)$signo > 0)/max(igraph::ecount(red_s),1),
            100*sum(igraph::E(red_s)$signo < 0)/max(igraph::ecount(red_s),1)))
cat(sprintf("Enfermo : %d nodos, %d aristas (%.1f%% pos, %.1f%% neg)\n",
            igraph::vcount(red_e), igraph::ecount(red_e),
            100*sum(igraph::E(red_e)$signo > 0)/max(igraph::ecount(red_e),1),
            100*sum(igraph::E(red_e)$signo < 0)/max(igraph::ecount(red_e),1)))

# =============================================================================
# 9. MÉTRICAS TOPOLÓGICAS
# =============================================================================
metricas_red <- function(g, nombre) {
  if (is.null(g) || igraph::vcount(g) < 2) return(NULL)
  n    <- igraph::vcount(g); m <- igraph::ecount(g)
  comp <- igraph::components(g); deg <- igraph::degree(g)
  cl   <- igraph::cluster_louvain(g)
  safe <- function(f, d = NA_real_)
    tryCatch(suppressWarnings(f()), error = function(e) d)
  
  tibble::tibble(
    Red           = nombre,
    Nodos         = n,
    Aristas       = m,
    Conectancia   = round(m/(n*(n-1)/2), 4),
    Grado_medio   = round(mean(deg), 2),
    Grado_max     = max(deg),
    N_componentes = comp$no,
    Tam_comp_max  = max(comp$csize),
    Pct_aislados  = round(100*sum(deg == 0)/n, 1),
    Pct_positivas = round(100*sum(igraph::E(g)$signo > 0)/max(m,1), 1),
    Pct_negativas = round(100*sum(igraph::E(g)$signo < 0)/max(m,1), 1),
    Modularidad   = safe(function() igraph::modularity(cl)),
    N_modulos     = length(unique(cl$membership)),
    Tam_mod_max   = max(table(cl$membership))
  )
}

res_metricas <- dplyr::bind_rows(
  metricas_red(red_s, "Sano"),
  metricas_red(red_e, "Enfermo")
)
cat("\n=== MÉTRICAS DE RED ===\n"); print(res_metricas)

# =============================================================================
# 10. HUBS (join directo nodo ↔ Genus)
# =============================================================================
hubs_final <- function(g, nombre, top_n = 15) {
  if (is.null(g) || igraph::vcount(g) < 2) return(NULL)
  tibble::tibble(
    nodo        = igraph::V(g)$name,
    grado       = igraph::degree(g),
    betweenness = round(igraph::betweenness(g, normalized = TRUE), 3)
  ) %>%
    dplyr::arrange(dplyr::desc(grado), dplyr::desc(betweenness)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::left_join(guild_map, by = c("nodo" = "Genus")) %>%
    dplyr::mutate(Red = nombre, .before = 1)
}

hubs_s <- hubs_final(red_s, "Sano")
hubs_e <- hubs_final(red_e, "Enfermo")
cat("\n=== HUBS SANO ===\n");    print(hubs_s)
cat("\n=== HUBS ENFERMO ===\n"); print(hubs_e)

# =============================================================================
# 11. MÓDULOS (composición top 5 por módulo)
# =============================================================================
modulos_final <- function(g, nombre, top_n = 5) {
  if (is.null(g)) return(NULL)
  cl <- igraph::cluster_louvain(g)
  tibble::tibble(
    nodo   = igraph::V(g)$name,
    modulo = cl$membership,
    grado  = igraph::degree(g)
  ) %>%
    dplyr::arrange(modulo, dplyr::desc(grado)) %>%
    dplyr::group_by(modulo) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(guild_map, by = c("nodo" = "Genus")) %>%
    dplyr::mutate(Red = nombre, .before = 1)
}

mod_s <- modulos_final(red_s, "Sano")
mod_e <- modulos_final(red_e, "Enfermo")
cat("\n=== MÓDULOS SANO ===\n");    print(mod_s)
cat("\n=== MÓDULOS ENFERMO ===\n"); print(mod_e)
# =============================================================================
# 12. FIGURAS (corregido: signo_f en la capa de edges)
# =============================================================================
plot_red_final <- function(g, titulo, max_nodos = 200) {
  if (is.null(g) || igraph::vcount(g) < 2)
    return(ggplot2::ggplot() + ggplot2::annotate("text", 0, 0, label = "Red vacía") +
             ggplot2::theme_void())
  
  if (igraph::vcount(g) > max_nodos) {
    keep <- igraph::V(g)[order(igraph::degree(g), decreasing = TRUE)[seq_len(max_nodos)]]
    g <- igraph::induced_subgraph(g, keep)
  }
  
  # --- FORZAR signo_f en el grafo ANTES de convertir ---------------------
  igraph::E(g)$signo_f <- factor(igraph::E(g)$signo, levels = c(-1, 1))
  
  cl  <- igraph::cluster_louvain(g)
  com <- as.factor(cl$membership)
  grd <- igraph::degree(g)
  
  tg <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(comunidad = com, grado = grd)
  # signo_f viaja como atributo de edges (ya está en g)
  
  n_com <- length(unique(com))
  pal   <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_com)
  
  ggraph::ggraph(tg, layout = "fr") +
    ggraph::geom_edge_link(ggplot2::aes(colour = signo_f), alpha = .55) +
    ggraph::scale_edge_colour_manual(
      values = c("-1" = "#B2182B", "1" = "#2166AC"),
      labels = c("-1" = "Negativa", "1" = "Positiva"),
      name   = "Signo", drop = FALSE) +
    ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado)) +
    ggplot2::scale_size_continuous(range = c(1.5, 6), guide = "none") +
    ggplot2::scale_colour_manual(values = pal, name = "Módulo") +
    ggplot2::labs(title = titulo) +
    ggraph::theme_graph(base_family = "sans", base_size = 11) +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
}

# Regenerar + mostrar + guardar
p_s <- plot_red_final(red_s, sprintf("Red Sano — %d nodos, %d aristas",
                                     igraph::vcount(red_s), igraph::ecount(red_s)))
p_e <- plot_red_final(red_e, sprintf("Red Enfermo — %d nodos, %d aristas",
                                     igraph::vcount(red_e), igraph::ecount(red_e)))

print(p_s)
print(p_e)

ggplot2::ggsave("fig_redes_final.png",
                patchwork::wrap_plots(p_s, p_e, ncol = 2),
                width = 13, height = 7, dpi = 220)
cat("\n✓ fig_redes_final.png guardado\n")

# --- Panel principal (alfa + NMDS + beta) ---------------------------------
p_alpha <- alpha_all %>%
  tidyr::pivot_longer(c(Shannon, Simpson, Riqueza, Chao1)) %>%
  ggplot2::ggplot(ggplot2::aes(Condicion, value, fill = Condicion)) +
  ggplot2::geom_boxplot(width = .55, alpha = .85, outlier.size = 1) +
  ggplot2::geom_jitter(width = .1, size = 1, alpha = .5) +
  ggplot2::facet_wrap(~ name, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::labs(title = "Diversidad alfa", x = NULL, y = NULL) +
  ggplot2::theme_minimal(11) +
  ggplot2::theme(legend.position = "none",
                 strip.text = ggplot2::element_text(face = "bold"))

nmds <- vegan::metaMDS(mat, "bray", k = 2, trymax = 100, trace = FALSE)
sc   <- as.data.frame(vegan::scores(nmds, "sites")) %>%
  dplyr::mutate(Condicion = meta_all$Condicion)

p_nmds <- ggplot2::ggplot(sc, ggplot2::aes(NMDS1, NMDS2, colour = Condicion)) +
  ggplot2::geom_point(size = 3, alpha = .85) +
  ggplot2::stat_ellipse(level = .68) +
  ggplot2::scale_colour_brewer(palette = "Set1") +
  ggplot2::labs(title = "NMDS (Bray-Curtis)",
                subtitle = sprintf("Stress=%.3f | PERMANOVA R²=%.2f, p=%.3f",
                                   nmds$stress, perm$R2[1], perm$`Pr(>F)`[1])) +
  ggplot2::theme_minimal(11) +
  ggplot2::theme(legend.position = "bottom")

p_beta <- beta_resumen %>%
  dplyr::select(Condicion, pct_turnover, pct_nestedness) %>%
  tidyr::pivot_longer(-Condicion) %>%
  ggplot2::ggplot(ggplot2::aes(Condicion, value, fill = name)) +
  ggplot2::geom_col(width = .6) +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(value > .05,
                                                 sprintf("%.0f%%", 100*value), "")),
                     position = ggplot2::position_stack(vjust = .5),
                     colour = "white", fontface = "bold", size = 3.3) +
  ggplot2::scale_fill_manual(values = c(pct_turnover = "#4DAF4A",
                                        pct_nestedness = "#984EA3"),
                             labels = c("Recambio","Anidamiento")) +
  ggplot2::labs(title = "β diversidad (Sørensen)", y = "Proporción",
                x = NULL, fill = NULL) +
  ggplot2::theme_minimal(11) +
  ggplot2::theme(legend.position = "bottom")

panel <- (p_alpha / (p_nmds | p_beta)) +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Comunidad fúngica: Sano vs Enfermo")

# --- Guardar figuras ------------------------------------------------------
ggplot2::ggsave("fig_alpha_nmds_beta.png", panel,
                width = 12, height = 8, dpi = 220)
ggplot2::ggsave("fig_nmds.png", p_nmds, width = 7, height = 5, dpi = 220)
ggplot2::ggsave("fig_redes_final.png",
                patchwork::wrap_plots(p_s, p_e, ncol = 2),
                width = 13, height = 7, dpi = 220)

cat("\n✓ Figuras guardadas: fig_alpha_nmds_beta.png, fig_nmds.png, fig_redes_final.png\n")

# =============================================================================
# 13. EXPORTAR TABLAS
# =============================================================================
readr::write_csv(alpha_all,     "res_alpha_por_muestra.csv")
readr::write_csv(tests_alpha,   "res_alpha_wilcoxon.csv")
readr::write_csv(beta_resumen,  "res_beta_particion.csv")
readr::write_csv(res_metricas,  "res_metricas_final.csv")
readr::write_csv(hubs_s,        "res_hubs_sano.csv")
readr::write_csv(hubs_e,        "res_hubs_enfermo.csv")
readr::write_csv(mod_s,         "res_modulos_sano.csv")
readr::write_csv(mod_e,         "res_modulos_enfermo.csv")

capture.output(perm,     file = "res_permanova.txt")
capture.output(anova(bd), file = "res_betadisper_anova.txt")
capture.output(null_s,   file = "res_nodf_sano.txt")
capture.output(null_e,   file = "res_nodf_enfermo.txt")

# =============================================================================
# 14. GUARDAR RDS PARA EL RMD
# =============================================================================
saveRDS(list(
  alpha_all    = alpha_all,
  tests_alpha  = tests_alpha,
  beta_resumen = beta_resumen,
  perm         = perm,
  bd           = bd,
  nodf_s       = null_s,
  nodf_e       = null_e,
  red_s        = red_s,
  red_e        = red_e,
  metricas     = res_metricas,
  hubs_s       = hubs_s,
  hubs_e       = hubs_e,
  modulos_s    = mod_s,
  modulos_e    = mod_e,
  gen_all      = gen_all,
  params       = list(presencia = 6, rho = 0.40, dz = 0.60, clr_pc = 0.5)
), file = "resultados_completos.rds")

cat("\n=== ANÁLISIS COMPLETADO ===\n")
cat("Tablas CSV : res_alpha_*.csv, res_beta_particion.csv,\n")
cat("             res_metricas_final.csv, res_hubs_*.csv, res_modulos_*.csv\n")
cat("Texto       : res_permanova.txt, res_betadisper_anova.txt, res_nodf_*.txt\n")
cat("Figuras     : fig_alpha_nmds_beta.png, fig_nmds.png, fig_redes_final.png\n")
cat("RDS         : resultados_completos.rds\n")



# =============================================================================
# COMPARACIÓN DE REDES: pres ≥ 6 vs pres ≥ 12
# (asume que ya corriste el script completo anterior)
# =============================================================================

# --- 1. Función auxiliar para construir ambos grafos -----------------------
build_network <- function(gen_all, idx_s, idx_e, pres, rho = 0.40, dz = 0.60) {
  keep <- colSums(gen_all > 0) >= pres
  sub <- gen_all[, keep, drop = FALSE]
  var_s <- apply(sub[idx_s, , drop = FALSE], 2, var)
  var_e <- apply(sub[idx_e, , drop = FALSE], 2, var)
  sub <- sub[, var_s > 0 & var_e > 0, drop = FALSE]
  list(
    s = construir_red_genero(sub[idx_s, , drop = FALSE], rho_min = rho, dz_max = dz),
    e = construir_red_genero(sub[idx_e, , drop = FALSE], rho_min = rho, dz_max = dz),
    n = ncol(sub)
  )
}

nets_6  <- build_network(gen_all, idx_s, idx_e, 6)
nets_12 <- build_network(gen_all, idx_s, idx_e, 12)

# --- 2. Comparación de métricas --------------------------------------------
comparar_nets <- function(nets, label) {
  tibble::tibble(
    Umbral  = label,
    Generos = nets$n,
    S_Nodos = igraph::vcount(nets$s),
    S_Aristas = igraph::ecount(nets$s),
    S_Pos = sum(igraph::E(nets$s)$signo > 0),
    S_Neg = sum(igraph::E(nets$s)$signo < 0),
    E_Nodos = igraph::vcount(nets$e),
    E_Aristas = igraph::ecount(nets$e),
    E_Pos = sum(igraph::E(nets$e)$signo > 0),
    E_Neg = sum(igraph::E(nets$e)$signo < 0)
  )
}
tabla_comparativa <- dplyr::bind_rows(
  comparar_nets(nets_6,  "pres ≥ 6"),
  comparar_nets(nets_12, "pres ≥ 12")
)
cat("\n=== COMPARACIÓN pres=6 vs pres=12 ===\n")
print(tabla_comparativa)

# --- 3. Función de plot (con bug fix: signo_f en edges) -------------------
plot_red_final <- function(g, titulo, max_nodos = 200) {
  if (is.null(g) || igraph::vcount(g) < 2)
    return(ggplot2::ggplot() + ggplot2::annotate("text", 0, 0, label = "Red vacía") +
             ggplot2::theme_void())
  
  if (igraph::vcount(g) > max_nodos) {
    keep <- igraph::V(g)[order(igraph::degree(g), decreasing = TRUE)[seq_len(max_nodos)]]
    g <- igraph::induced_subgraph(g, keep)
  }
  
  igraph::E(g)$signo_f <- factor(igraph::E(g)$signo, levels = c(-1, 1))
  
  cl  <- igraph::cluster_louvain(g)
  com <- as.factor(cl$membership)
  grd <- igraph::degree(g)
  
  tg <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(comunidad = com, grado = grd)
  
  n_com <- length(unique(com))
  pal   <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_com)
  
  ggraph::ggraph(tg, layout = "fr") +
    ggraph::geom_edge_link(ggplot2::aes(colour = signo_f), alpha = .55) +
    ggraph::scale_edge_colour_manual(
      values = c("-1" = "#B2182B", "1" = "#2166AC"),
      labels = c("-1" = "Negativa", "1" = "Positiva"),
      name = "Signo", drop = FALSE) +
    ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado)) +
    ggplot2::scale_size_continuous(range = c(1.5, 6), guide = "none") +
    ggplot2::scale_colour_manual(values = pal, name = "Módulo") +
    ggplot2::labs(title = titulo) +
    ggraph::theme_graph(base_family = "sans", base_size = 11) +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
}

# --- 4. Generar las 4 figuras ---------------------------------------------
p_p6_s  <- plot_red_final(nets_6$s,  sprintf("Sano — pres ≥ 6 (%d nodos, %d aristas)",
                                             igraph::vcount(nets_6$s), igraph::ecount(nets_6$s)))
p_p6_e  <- plot_red_final(nets_6$e,  sprintf("Enfermo — pres ≥ 6 (%d nodos, %d aristas)",
                                             igraph::vcount(nets_6$e), igraph::ecount(nets_6$e)))
p_p12_s <- plot_red_final(nets_12$s, sprintf("Sano — pres ≥ 12 (%d nodos, %d aristas)",
                                             igraph::vcount(nets_12$s), igraph::ecount(nets_12$s)))
p_p12_e <- plot_red_final(nets_12$e, sprintf("Enfermo — pres ≥ 12 (%d nodos, %d aristas)",
                                             igraph::vcount(nets_12$e), igraph::ecount(nets_12$e)))

# Mostrar cada uno
print(p_p6_s);  print(p_p6_e)
print(p_p12_s); print(p_p12_e)

# --- 5. Guardar cada figura por separado ----------------------------------
ggplot2::ggsave("fig_red_pres6_sano.png",    p_p6_s,  width = 9, height = 7, dpi = 220)
ggplot2::ggsave("fig_red_pres6_enfermo.png", p_p6_e,  width = 9, height = 7, dpi = 220)
ggplot2::ggsave("fig_red_pres12_sano.png",   p_p12_s, width = 9, height = 7, dpi = 220)
ggplot2::ggsave("fig_red_pres12_enfermo.png",p_p12_e, width = 9, height = 7, dpi = 220)

# --- 6. Panel combinado 2×2 -----------------------------------------------
panel <- patchwork::wrap_plots(
  p_p6_s, p_p6_e, p_p12_s, p_p12_e,
  ncol = 2, nrow = 2, byrow = TRUE
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Redes de co-ocurrencia fúngica: Sano vs Enfermo bajo dos umbrales de presencia",
    subtitle = "Fila superior: pres ≥ 6 (exploratorio, 69 géneros) · Fila inferior: pres ≥ 12 (estándar, 37 géneros)",
    caption = "CLR + Pearson · |r| ≥ 0.40 · dz-avoidance ≤ 0.60 · Azul = cooperación, rojo = competencia"
  )

ggplot2::ggsave("fig_redes_comparativa.png", panel, width = 16, height = 13, dpi = 220)

cat("\n✓ Figuras guardadas:\n")
cat("  · fig_red_pres6_sano.png\n")
cat("  · fig_red_pres6_enfermo.png\n")
cat("  · fig_red_pres12_sano.png\n")
cat("  · fig_red_pres12_enfermo.png\n")
cat("  · fig_redes_comparativa.png (panel 2×2)\n")

# --- 7. Exportar tabla comparativa ----------------------------------------
readr::write_csv(tabla_comparativa, "res_comparativa_pres6_vs_pres12.csv")
cat("✓ res_comparativa_pres6_vs_pres12.csv exportado\n")
# =============================================================================
# v9 — Circular + save robusto + título limpio
# =============================================================================

rm(list = ls(pattern = "^plot_red"))
rm(list = c("p_p6_s","p_p6_e","p_p12_s","p_p12_e"), envir = .GlobalEnv)

plot_red_v9 <- function(g, condicion, pres, top_modulos = 10, seed = 42) {
  cat(sprintf("[v9] %s pres≥%d\n", condicion, pres))
  
  if (is.null(g) || igraph::vcount(g) < 2)
    return(ggplot2::ggplot() + ggplot2::theme_void())
  
  n_pos   <- sum(igraph::E(g)$signo > 0)
  n_neg   <- sum(igraph2::E(g)$signo < 0) |> try(silent = TRUE)
  n_neg   <- sum(igraph::E(g)$signo < 0)
  n_nodes <- igraph::vcount(g)
  n_edges <- igraph::ecount(g)
  
  igraph::E(g)$signo_f <- factor(igraph::E(g)$signo, levels = c("-1", "1"))
  
  # --- Módulos -----------------------------------------------------------
  cl <- igraph::cluster_louvain(g); mem <- cl$membership
  mod_sizes <- sort(table(mem), decreasing = TRUE)
  n_keep    <- min(top_modulos, length(mod_sizes))
  top_ids   <- as.integer(names(mod_sizes))[seq_len(n_keep)]
  fuera_top <- !mem %in% top_ids
  com_label <- ifelse(fuera_top, "Otros", paste0("M", mem))
  levels_order <- paste0("M", top_ids)
  if (any(fuera_top)) levels_order <- c(levels_order, "Otros")
  com_f <- factor(com_label, levels = levels_order)
  mod_counts <- as.integer(table(com_f))
  mod_labels <- setNames(paste0(levels(com_f), " (n=", mod_counts, ")"),
                         levels(com_f))
  
  # --- Layout CIRCULAR, ordenado por módulo -------------------------------
  grd <- igraph::degree(g)
  tg  <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(comunidad = com_f, grado = grd, nombre = name)
  
  lay <- ggraph::create_layout(tg, layout = "circle", sort.by = comunidad)
  
  # Radianes + rotación de labels para que sean radiales
  lay <- lay %>%
    dplyr::mutate(
      angle_rad = atan2(y, x),
      angle_deg = angle_rad * 180 / pi,
      # Voltear labels del lado izquierdo para que no se lean al revés
      angle_deg = ifelse(angle_deg > 90 | angle_deg < -90,
                         angle_deg + 180, angle_deg),
      hjust     = ifelse(x >= 0, 0, 1)
    )
  
  # --- Paleta -------------------------------------------------------------
  pal_base <- scales::alpha(grDevices::hcl.colors(n_keep, "Dark 3"), 0.78)
  pal      <- setNames(pal_base, paste0("M", top_ids))
  if (any(fuera_top))
    pal <- c(pal, "Otros" = scales::alpha("grey78", 0.78))
  
  # --- Info box -----------------------------------------------------------
  info_txt <- sprintf("pres \u2265 %d  \u00b7  %d nodos  \u00b7  %d interacciones",
                      pres, n_nodes, n_edges)
  
  # --- Plot ---------------------------------------------------------------
  ggraph::ggraph(lay) +
    ggraph::geom_edge_arc(
      ggplot2::aes(colour = signo_f),
      alpha = 0.40, strength = 0.6, width = 0.5
    ) +
    ggraph::scale_edge_colour_manual(
      values = c("-1" = "#B2182B", "1" = "#2166AC"),
      labels = c("-1" = sprintf("Negativa (%d)", n_neg),
                 "1"  = sprintf("Positiva (%d)", n_pos)),
      name = "Signo", limits = c("-1", "1"), drop = FALSE,
      guide = ggplot2::guide_legend(override.aes = list(
        edge_colour = c("#B2182B", "#2166AC"),
        edge_linewidth = 1.5
      ))
    ) +
    ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado),
                            alpha = 0.20, stroke = 0) +
    ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado),
                            alpha = 0.95, stroke = 0) +
    ggplot2::scale_size_continuous(range = c(3, 12), guide = "none") +
    ggplot2::scale_colour_manual(
      values = pal, labels = mod_labels, name = "Módulo",
      limits = levels(com_f), drop = FALSE
    ) +
    # Labels radiales (rotados) — todos, sin filtro
    ggrepel::geom_text_repel(
      data = lay,
      ggplot2::aes(x = x, y = y, label = nombre, angle = angle_deg,
                   hjust = hjust),
      size = 2.0,
      colour = "grey5",
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      segment.color = "grey50",
      segment.size = 0.2,
      segment.alpha = 0.6,
      force = 2, force_pull = 0.3,
      seed = seed
    ) +
    # Info box debajo del círculo
    ggplot2::annotate("label", x = 0, y = -1.6,
                      label = info_txt,
                      size = 3.2,
                      label.size = 0,
                      fill = scales::alpha("white", 0.9),
                      colour = "grey5", fontface = "bold") +
    ggplot2::coord_fixed(clip = "off",
                         xlim = c(-1.5, 1.5), ylim = c(-1.75, 1.3)) +
    ggplot2::labs(title = condicion) +   # solo "Sano" / "Enfermo"
    ggraph::theme_graph(base_family = "sans", base_size = 10) +
    ggplot2::theme(
      legend.position = "none",
      plot.title      = ggplot2::element_text(
        size = 16, face = "bold", hjust = 0.5, colour = "black",
        margin = ggplot2::margin(t = 4, b = 12)),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin     = ggplot2::margin(10, 15, 10, 15)
    )
}

# =============================================================================
# CONSTRUIR LOS 4 PANELES
# =============================================================================
cat("\n>>> Construyendo paneles con v9:\n")
p_p6_s  <- plot_red_v9(nets_6$s,  "Sano",    6)
p_p6_e  <- plot_red_v9(nets_6$e,  "Enfermo", 6)
p_p12_s <- plot_red_v9(nets_12$s, "Sano",   12)
p_p12_e <- plot_red_v9(nets_12$e, "Enfermo", 12)

# Verificación de títulos (debe ser solo "Sano" o "Enfermo")
titles <- sapply(list(p_p6_s, p_p6_e, p_p12_s, p_p12_e),
                 function(p) p$labels$title)
stopifnot(all(titles %in% c("Sano", "Enfermo")))
cat("✓ Títulos:", paste(titles, collapse = " | "), "\n")

# Mostrar el primero
print(p_p6_s)

# =============================================================================
# GUARDADO A PRUEBA DE FALLOS (ragg → cairo → PDF)
# =============================================================================
save_fig <- function(p, filename, w = 14, h = 11, dpi = 200) {
  # Intento 1: ggsave estándar
  ok <- tryCatch({
    ggplot2::ggsave(filename, plot = p, width = w, height = h,
                    dpi = dpi, bg = "white", limitsize = FALSE)
    TRUE
  }, error = function(e) FALSE)
  
  if (!ok) {
    # Intento 2: device cairo_png explícito
    ok <- tryCatch({
      grDevices::png(filename, width = w*dpi, height = h*dpi,
                     res = dpi, bg = "white", type = "cairo")
      print(p)
      grDevices::dev.off()
      TRUE
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      FALSE
    })
  }
  
  if (!ok) {
    # Intento 3: PDF
    pdf_name <- sub("\\.png$", ".pdf", filename)
    tryCatch({
      ggplot2::ggsave(pdf_name, plot = p, width = w, height = h,
                      bg = "white", limitsize = FALSE)
      cat("⚠ PDF fallback:", pdf_name, "\n")
      return(invisible(pdf_name))
    }, error = function(e) {
      cat("❌ No se pudo guardar:", filename, "-", conditionMessage(e), "\n")
    })
    return(invisible(NULL))
  }
  
  cat("✓", filename, "\n")
  invisible(filename)
}

save_fig(p_p6_s,  "fig_red_pres6_sano.png",     w = 14, h = 12)
save_fig(p_p6_e,  "fig_red_pres6_enfermo.png",  w = 14, h = 12)
save_fig(p_p12_s, "fig_red_pres12_sano.png",    w = 14, h = 12)
save_fig(p_p12_e, "fig_red_pres12_enfermo.png", w = 14, h = 12)

# =============================================================================
# PANEL 2×2 — PDF directo (PNG de 4 paneles siempre es problemático)
# =============================================================================
panel <- patchwork::wrap_plots(
  p_p6_s, p_p6_e,
  p_p12_s, p_p12_e,
  ncol = 2, nrow = 2, byrow = TRUE
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title    = "Redes de co-ocurrencia fúngica: Sano vs Enfermo",
    subtitle = "Fila superior: pres \u2265 6 \u00b7 Fila inferior: pres \u2265 12",
    caption  = "CLR + Pearson \u00b7 |r| \u2265 0.40 \u00b7 dz \u2264 0.60 \u00b7 Azul = cooperación, rojo = competencia",
    theme    = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 15, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey25"),
      plot.caption  = ggplot2::element_text(size = 8.5, colour = "grey40", hjust = 0),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    ))

# Intento como PNG primero, si falla → PDF
saved_panel <- tryCatch({
  ggplot2::ggsave("fig_redes_comparativa.png", panel,
                  width = 22, height = 18, dpi = 180, bg = "white",
                  limitsize = FALSE)
  "fig_redes_comparativa.png"
}, error = function(e) {
  ggplot2::ggsave("fig_redes_comparativa.pdf", panel,
                  width = 22, height = 18, bg = "white", limitsize = FALSE)
  "fig_redes_comparativa.pdf"
})
cat("✓ Panel:", saved_panel, "\n")















#################################################

# =============================================================================
# v12 — Dos estilos de red, ambos con etiquetas legibles
#   · Red      : layout stress + aristas CURVAS + etiquetas legibles
#   · Circular : nodos en círculo + etiquetas radiales a DOS radios alternos
# =============================================================================

rm(list = ls(pattern = "^plot_"))
rm(list = ls(pattern = "^p_"))
rm(list = ls(pattern = "^panel"))

plot_red_v12 <- function(g, cond, pres, style = c("network", "circular"),
                         top_modulos = 10, seed = 42) {
  style <- match.arg(style)
  if (is.null(g) || igraph::vcount(g) < 2)
    return(ggplot2::ggplot() + ggplot2::theme_void())
  
  n_pos   <- sum(igraph::E(g)$signo > 0)
  n_neg   <- sum(igraph::E(g)$signo < 0)
  n_nodes <- igraph::vcount(g)
  n_edges <- igraph::ecount(g)
  
  igraph::E(g)$signo_f <- factor(igraph::E(g)$signo, levels = c("-1", "1"))
  
  cl <- igraph::cluster_louvain(g); mem <- cl$membership
  mod_sizes <- sort(table(mem), decreasing = TRUE)
  n_keep <- min(top_modulos, length(mod_sizes))
  top_ids <- as.integer(names(mod_sizes))[seq_len(n_keep)]
  fuera_top <- !mem %in% top_ids
  com_label <- ifelse(fuera_top, "Otros", paste0("M", mem))
  levels_order <- paste0("M", top_ids)
  if (any(fuera_top)) levels_order <- c(levels_order, "Otros")
  com_f <- factor(com_label, levels = levels_order)
  mod_labels <- setNames(paste0(levels(com_f), " (n=",
                                as.integer(table(com_f)), ")"),
                         levels(com_f))
  
  pal_base <- scales::alpha(grDevices::hcl.colors(n_keep, "Dark 3"), 0.88)
  pal <- setNames(pal_base, paste0("M", top_ids))
  if (any(fuera_top)) pal <- c(pal, "Otros" = scales::alpha("grey75", 0.88))
  
  grd <- igraph::degree(g)
  tg <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(comunidad = com_f, grado = grd, nombre = name)
  
  info_txt <- sprintf("pres \u2265 %d  \u00b7  %d nodos  \u00b7  %d interacciones",
                      pres, n_nodes, n_edges)
  
  # ---- Paleta y escalas comunes ----
  scale_edge <- ggraph::scale_edge_colour_manual(
    values = c("-1" = "#B2182B", "1" = "#2166AC"),
    labels = c("-1" = sprintf("Negativa (%d)", n_neg),
               "1"  = sprintf("Positiva (%d)", n_pos)),
    name = "Signo", limits = c("-1", "1"), drop = FALSE,
    guide = ggplot2::guide_legend(override.aes = list(
      edge_colour = c("#B2182B", "#2166AC"),
      edge_linewidth = 1.2
    ))
  )
  scale_node <- ggplot2::scale_colour_manual(
    values = pal, name = "Módulo",
    limits = levels(com_f), drop = FALSE, labels = mod_labels
  )
  
  # ================= ESTILO 1: RED (stress + curvas) =================
  if (style == "network") {
    lay <- ggraph::create_layout(tg, layout = "stress", bbox = 25)
    
    p <- ggraph::ggraph(lay) +
      ggraph::geom_edge_arc(ggplot2::aes(colour = signo_f),
                            strength = 0.30, alpha = 0.50, width = 0.6) +
      scale_edge +
      ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado),
                              alpha = 0.25, stroke = 0) +
      ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado),
                              alpha = 0.95, stroke = 0) +
      ggplot2::scale_size_continuous(range = c(3, 12), guide = "none") +
      scale_node +
      ggrepel::geom_text_repel(
        data = lay,
        ggplot2::aes(x, y, label = nombre),
        size = 2.8, colour = "grey5", fontface = "bold",
        max.overlaps       = Inf,
        box.padding        = 0.45,
        point.padding      = 0.25,
        min.segment.length = 0,
        segment.color      = "grey40",
        segment.size       = 0.3,
        segment.alpha      = 0.6,
        force              = 3,
        force_pull         = 0.7,
        max.iter           = 100000,
        seed               = seed
      ) +
      ggplot2::labs(title = cond) +
      ggraph::theme_graph(base_family = "sans", base_size = 10) +
      ggplot2::theme(
        legend.position = "none",
        plot.title      = ggplot2::element_text(
          size = 18, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(t = 4, b = 25)),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin     = ggplot2::margin(20, 20, 20, 20)
      ) +
      ggplot2::annotate("label",
                        x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
                        label = info_txt, size = 3.5,
                        fill = scales::alpha("white", 0.9),
                        colour = "grey5", fontface = "bold",
                        label.size = NA) +
      ggplot2::coord_cartesian(clip = "off")
  }
  
  # ================= ESTILO 2: CIRCULAR (etiquetas a dos radios) =========
  if (style == "circular") {
    lay <- ggraph::create_layout(tg, layout = "circle")
    ord <- order(com_f)
    n   <- nrow(lay)
    ang <- seq(0, 2*pi, length.out = n + 1)[-1]
    lay <- lay[ord, ]
    lay$x <- cos(ang)
    lay$y <- sin(ang)
    
    # Etiquetas a DOS radios alternos → no se pisan
    lay <- lay %>%
      dplyr::mutate(
        angle_deg     = ang * 180 / pi,
        angle_deg_adj = ifelse(angle_deg > 90, angle_deg - 180,
                               ifelse(angle_deg < -90, angle_deg + 180, angle_deg)),
        hjust_lab     = ifelse(cos(ang) >= 0, 0, 1),
        label_r       = ifelse(dplyr::row_number() %% 2 == 0, 1.55, 2.15),
        x_lab         = label_r * cos(ang),
        y_lab         = label_r * sin(ang)
      )
    
    p <- ggraph::ggraph(lay) +
      ggraph::geom_edge_arc(ggplot2::aes(colour = signo_f),
                            strength = -0.12, alpha = 0.35, width = 0.5) +
      scale_edge +
      ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado),
                              alpha = 0.25, stroke = 0) +
      ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado),
                              alpha = 0.95, stroke = 0) +
      ggplot2::scale_size_continuous(range = c(3, 12), guide = "none") +
      scale_node +
      ggrepel::geom_text_repel(
        data = lay,
        ggplot2::aes(x = x_lab, y = y_lab, label = nombre,
                     angle = angle_deg_adj, hjust = hjust_lab),
        size = 2.5, colour = "grey5", fontface = "bold",
        max.overlaps       = Inf,
        box.padding        = 0.10,
        point.padding      = 0.05,
        min.segment.length = 0,
        segment.color      = "grey50",
        segment.size       = 0.2,
        segment.alpha      = 0.4,
        force              = 0.5,
        max.iter           = 100000,
        seed               = seed
      ) +
      ggplot2::coord_fixed(clip = "off",
                           xlim = c(-2.8, 2.8), ylim = c(-2.8, 2.8)) +
      ggplot2::labs(title = cond) +
      ggraph::theme_graph(base_family = "sans", base_size = 10) +
      ggplot2::theme(
        legend.position = "none",
        plot.title      = ggplot2::element_text(
          size = 18, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(t = 4, b = 25)),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin     = ggplot2::margin(20, 20, 20, 20)
      ) +
      ggplot2::annotate("label",
                        x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
                        label = info_txt, size = 3.5,
                        fill = scales::alpha("white", 0.9),
                        colour = "grey5", fontface = "bold",
                        label.size = NA)
  }
  
  p
}

# =============================================================================
# CONSTRUIR LOS 4 PANELES
# =============================================================================
p_net_s  <- plot_red_v12(nets_6$s, "Sano",    6, style = "network")
p_net_e  <- plot_red_v12(nets_6$e, "Enfermo", 6, style = "network")
p_circ_s <- plot_red_v12(nets_6$s, "Sano",    6, style = "circular")
p_circ_e <- plot_red_v12(nets_6$e, "Enfermo", 6, style = "circular")

# Mostrar el primero para verificar
print(p_net_s)

# =============================================================================
# GUARDAR CADA UNO POR SEPARADO
# =============================================================================
save_plot <- function(p, fn, w, h) {
  # PDF primero (más seguro para tamaños grandes)
  pdf_fn <- sub("\\.png$", ".pdf", fn)
  ggplot2::ggsave(pdf_fn, plot = p, width = w, height = h,
                  bg = "white", limitsize = FALSE, device = "pdf")
  cat("✓", pdf_fn, "\n")
  
  # Intento PNG a dpi moderado
  ok <- tryCatch({
    ggplot2::ggsave(fn, plot = p, width = w, height = h, dpi = 150,
                    bg = "white", limitsize = FALSE)
    TRUE
  }, error = function(e) FALSE)
  
  if (ok) cat("✓", fn, "\n") else cat("⚠ PNG falló, usar PDF:", pdf_fn, "\n")
}

save_plot(p_net_s,  "fig_net_Sano.png",     14, 12)
save_plot(p_net_e,  "fig_net_Enfermo.png",  14, 12)
save_plot(p_circ_s, "fig_circ_Sano.png",    18, 18)
save_plot(p_circ_e, "fig_circ_Enfermo.png", 18, 18)

# =============================================================================
# PANEL 2×2 FINAL — 2 estilos × 2 condiciones
# =============================================================================
panel <- patchwork::wrap_plots(
  p_net_s, p_net_e,
  p_circ_s, p_circ_e,
  ncol = 2, nrow = 2, byrow = TRUE
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Redes de co-ocurrencia fúngica: Sano vs Enfermo",
    subtitle = paste("Fila superior: red con aristas curvas (stress layout)",
                     "\u00b7 Fila inferior: topolog\u00eda circular"),
    caption = paste("CLR + Pearson \u00b7 |r| \u2265 0.40 \u00b7 dz \u2264 0.60",
                    "\u00b7 Azul = cooperaci\u00f3n, rojo = competencia",
                    "\u00b7 Etiquetas en dos radios alternos para evitar solapamiento (fila inferior)"),
    theme = ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 20, face = "bold"),
      plot.subtitle   = ggplot2::element_text(size = 13, colour = "grey25"),
      plot.caption    = ggplot2::element_text(size = 10, colour = "grey40",
                                              hjust = 0),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  )

# Guardar panel como PDF (tamaño grande seguro)
ggplot2::ggsave("fig_redes_comparativa.pdf", panel,
                width = 30, height = 26, device = "pdf",
                bg = "white", limitsize = FALSE)
cat("✓ fig_redes_comparativa.pdf\n")

# PNG a baja dpi como vista previa
tryCatch({
  ggplot2::ggsave("fig_redes_comparativa.png", panel,
                  width = 30, height = 26, dpi = 100,
                  bg = "white", limitsize = FALSE)
  cat("✓ fig_redes_comparativa.png (preview)\n")
}, error = function(e) cat("⚠ PNG panel falló, usar el PDF\n"))



# =============================================================================
# REDES CORREGIDAS — reemplaza v9 y v12
#   · Filtra nodos aislados (grado = 0) → sin "collar de cuentas"
#   · Etiqueta SOLO hubs (grado ≥ 2) → sin amontonamiento
#   · Subtítulo con conteos: N conectados · M aristas · K aislados
#   · Circular solo con subgrafo conectado
#   · annotate sin label.size (evita warnings)
#   · Guardado a prueba de Windows: cairo_pdf como fallback
# =============================================================================

plot_net_fixed <- function(g, cond, pres,
                           style           = c("stress", "circular"),
                           top_modulos     = 8,
                           min_deg_label   = 2,
                           seed            = 42) {
  style <- match.arg(style)
  
  if (is.null(g) || igraph::vcount(g) < 2)
    return(ggplot2::ggplot() + ggplot2::theme_void())
  
  # --- 1. Separar conectados / aislados ------------------------------------
  deg     <- igraph::degree(g)
  con_idx <- which(deg >= 1)
  n_iso   <- sum(deg == 0)
  n_tot   <- igraph::vcount(g)
  n_edg   <- igraph::ecount(g)
  
  if (length(con_idx) < 2)
    return(ggplot2::ggplot() + ggplot2::theme_void())
  
  g_con <- igraph::induced_subgraph(g, con_idx)
  n_con <- igraph::vcount(g_con)
  
  # --- 2. Módulos Louvain SOLO sobre conectados ----------------------------
  cl      <- igraph::cluster_louvain(g_con)
  mem     <- cl$membership
  sizes   <- sort(table(mem), decreasing = TRUE)
  n_keep  <- min(top_modulos, length(sizes))
  tops    <- as.integer(names(sizes))[seq_len(n_keep)]
  com     <- factor(ifelse(mem %in% tops, paste0("M", mem), "Otros"),
                    levels = c(paste0("M", tops), "Otros"))
  
  pal <- setNames(
    c(scales::alpha(grDevices::hcl.colors(n_keep, "Dark 3"), 0.88),
      scales::alpha("grey78", 0.85)),
    c(paste0("M", tops), "Otros")
  )
  
  # --- 3. Layout -----------------------------------------------------------
  grd <- igraph::degree(g_con)
  tg  <- tidygraph::as_tbl_graph(g_con) |>
    tidygraph::activate(nodes) |>
    dplyr::mutate(comunidad = com, grado = grd, nombre = name)
  
  if (style == "stress") {
    lay      <- ggraph::create_layout(tg, layout = "stress", bbox = 10)
    label_df <- lay |> dplyr::filter(grado >= min_deg_label)
  } else {
    n   <- igraph::vcount(g_con)
    ang <- seq(0, 2 * pi, length.out = n + 1)[-1]
    lay <- ggraph::create_layout(tg, layout = "manual",
                                 x = cos(ang), y = sin(ang))
    # ángulo radial para etiquetas (volteadas en el lado izquierdo)
    label_df <- lay |>
      dplyr::filter(grado >= min_deg_label) |>
      dplyr::mutate(
        ang_rad   = atan2(y, x),
        ang_deg   = ang_rad * 180 / pi,
        ang_deg   = ifelse(ang_deg >  90, ang_deg - 180,
                           ifelse(ang_deg < -90, ang_deg + 180, ang_deg)),
        hjust_lab = ifelse(cos(ang_rad) >= 0, 0, 1)
      )
  }
  
  # --- 4. Escalas y plot base ---------------------------------------------
  p <- ggraph::ggraph(lay) +
    ggraph::geom_edge_arc(
      ggplot2::aes(colour = factor(signo)),
      strength = if (style == "stress") 0.30 else -0.20,
      alpha    = 0.50, width = 0.5
    ) +
    ggraph::scale_edge_colour_manual(
      values = c("-1" = "#B2182B", "1" = "#2166AC"),
      labels = c("Negativa", "Positiva"),
      name   = "Signo", drop = FALSE
    ) +
    ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado),
                            alpha = 0.95, stroke = 0) +
    ggplot2::scale_size_continuous(range = c(2, 10), guide = "none") +
    ggplot2::scale_colour_manual(values = pal, name = "Módulo", drop = FALSE) +
    ggplot2::labs(
      title    = sprintf("%s · %s", cond, style),
      subtitle = sprintf(
        "pres ≥ %d  ·  %d conectados  ·  %d aristas  ·  %d aislados no mostrados",
        pres, n_con, n_edg, n_iso
      )
    ) +
    ggraph::theme_graph(base_family = "sans", base_size = 10) +
    ggplot2::theme(
      legend.position = "right",
      plot.title      = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(size = 8.5, colour = "grey30",
                                              hjust = 0.5,
                                              margin = ggplot2::margin(b = 6))
    )
  
  # --- 5. Etiquetas (solo hubs) -------------------------------------------
  if (style == "circular") {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(x = x * 1.15, y = y * 1.15, label = nombre,
                   angle = ang_deg, hjust = hjust_lab),
      size = 2.3, colour = "grey10", fontface = "bold",
      max.overlaps = Inf, segment.size = 0.2, segment.colour = "grey55",
      box.padding = 0.25, point.padding = 0.10, seed = seed
    )
  } else {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(x, y, label = nombre),
      size = 2.3, colour = "grey10", fontface = "bold",
      max.overlaps = Inf, segment.size = 0.2, segment.colour = "grey55",
      box.padding = 0.35, point.padding = 0.15, seed = seed
    )
  }
  p
}

# =============================================================================
# CONSTRUIR LOS 4 PANELES
# =============================================================================
p_ss <- plot_net_fixed(nets_6$s,  "Sano",    6, style = "stress")
p_es <- plot_net_fixed(nets_6$e,  "Enfermo", 6, style = "stress")
p_sc <- plot_net_fixed(nets_6$s,  "Sano",    6, style = "circular")
p_ec <- plot_net_fixed(nets_6$e,  "Enfermo", 6, style = "circular")

# =============================================================================
# GUARDADO INDIVIDUAL (cairo_pdf como fallback)
# =============================================================================
save_net <- function(p, base_fn, w, h, dpi = 200) {
  try(graphics.off(), silent = TRUE)
  while (!is.null(grDevices::dev.list()))
    try(grDevices::dev.off(), silent = TRUE)
  
  pdf_file <- paste0(base_fn, ".pdf")
  png_file <- paste0(base_fn, ".png")
  if (file.exists(pdf_file)) try(file.remove(pdf_file), silent = TRUE)
  if (file.exists(png_file)) try(file.remove(png_file), silent = TRUE)
  Sys.sleep(0.2)
  
  ok_pdf <- tryCatch({
    grDevices::cairo_pdf(pdf_file, width = w, height = h, bg = "white")
    print(p); grDevices::dev.off(); TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
  
  if (!ok_pdf) {
    ok_pdf <- tryCatch({
      grDevices::pdf(pdf_file, width = w, height = h, bg = "white")
      print(p); grDevices::dev.off(); TRUE
    }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
  }
  
  ok_png <- tryCatch({
    ggplot2::ggsave(png_file, plot = p, width = w, height = h,
                    dpi = dpi, bg = "white", limitsize = FALSE,
                    device = ragg::agg_png)
    TRUE
  }, error = function(e) FALSE)
  
  if (ok_pdf) cat("  ✓", pdf_file, "\n") else cat("  ✗ PDF:", pdf_file, "\n")
  if (ok_png) cat("  ✓", png_file, "\n") else cat("  ⚠ PNG:", png_file, "\n")
}

save_net(p_ss, "fig_net_Sano",      14, 11)
save_net(p_es, "fig_net_Enfermo",   14, 11)
save_net(p_sc, "fig_circ_Sano",     16, 14)
save_net(p_ec, "fig_circ_Enfermo",  16, 14)

# =============================================================================
# PANEL 2×2 — PDF grande + PNG preview a baja dpi
# =============================================================================
panel <- patchwork::wrap_plots(
  p_ss, p_es,
  p_sc, p_ec,
  ncol = 2, nrow = 2, byrow = TRUE
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title    = "Redes de co-ocurrencia fúngica: Sano vs Enfermo (pres ≥ 6)",
    subtitle = "Fila superior: stress layout  ·  Fila inferior: circular (solo subgrafo conectado)",
    caption  = paste(
      "CLR + Pearson · |r| ≥ 0.40 · dz ≤ 0.60 · Azul = cooperación · Rojo = competencia",
      "· Etiquetas solo para hubs (grado ≥ 2) · Nodos aislados no mostrados (conteo en cada subtítulo)",
      sep = "\n"
    ),
    theme = ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 18, face = "bold"),
      plot.subtitle   = ggplot2::element_text(size = 11, colour = "grey25"),
      plot.caption    = ggplot2::element_text(size = 8.5, colour = "grey40",
                                              hjust = 0),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  )

# PDF grande (seguro)
try(grDevices::cairo_pdf("fig_redes_comparativa.pdf",
                         width = 22, height = 18, bg = "white"))
print(panel); try(grDevices::dev.off(), silent = TRUE)
cat("✓ fig_redes_comparativa.pdf\n")

# PNG preview a dpi bajo
tryCatch({
  ggplot2::ggsave("fig_redes_comparativa.png", panel,
                  width = 22, height = 18, dpi = 100,
                  bg = "white", limitsize = FALSE,
                  device = ragg::agg_png)
  cat("✓ fig_redes_comparativa.png (preview)\n")
}, error = function(e) cat("⚠ PNG preview omitido:", conditionMessage(e), "\n"))






###################################################3
####################################################


# =============================================================================
# FUNCIÓN DEFINITIVA v11: Aristas visibles, cero espacio horizontal, sin líneas
# =============================================================================

plot_net_v11 <- function(g, cond, pres, style = c("stress", "circular"), 
                         top_modulos = 8, seed = 42) {
  
  style <- match.arg(style)
  if (is.null(g) || igraph::vcount(g) < 2) return(ggplot2::theme_void())
  
  degrees <- igraph::degree(g)
  n_total <- igraph::vcount(g)
  con_idx <- which(degrees >= 1)
  iso_idx <- which(degrees == 0)
  n_con   <- length(con_idx)
  n_iso   <- length(iso_idx)
  
  g_con <- igraph::induced_subgraph(g, con_idx)
  
  # --- 1. Módulos Louvain ---
  cl    <- igraph::cluster_louvain(g_con)
  mem   <- cl$membership
  sizes <- sort(table(mem), decreasing = TRUE)
  n_keep <- min(top_modulos, length(sizes))
  tops   <- as.integer(names(sizes))[seq_len(n_keep)]
  
  com_all <- character(n_total)
  com_all[con_idx] <- ifelse(mem %in% tops, paste0("M", mem), "Otros")
  com_all[iso_idx] <- "Aislado"
  
  # --- 2. Paleta de colores ---
  pal <- setNames(
    c(scales::alpha(grDevices::hcl.colors(n_keep, "Dark 3"), 0.88),
      scales::alpha("#6C3483", 0.9),  # Otros = Morado
      scales::alpha("#A6ACAF", 0.9)), # Aislado = Gris claro visible
    c(paste0("M", tops), "Otros", "Aislado")
  )
  
  # --- 3. Layout ---
  x_all <- numeric(n_total); y_all <- numeric(n_total)
  
  if (style == "circular") {
    ord <- order(factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado")))
    ang <- seq(0, 2*pi, length.out = n_total + 1)[-1]
    
    x_all <- cos(ang)[ord]
    y_all <- sin(ang)[ord]
    
    igraph::V(g)$name <- igraph::V(g)$name[ord]
    igraph::V(g)$comunidad <- factor(com_all[ord], levels = c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$grado <- degrees[ord]
    
  } else {
    lay_con <- ggraph::create_layout(tidygraph::as_tbl_graph(g_con), layout = "stress", bbox = 8)
    x_all[con_idx] <- lay_con$x
    y_all[con_idx] <- lay_con$y
    
    if (n_iso > 0) {
      x_rng <- diff(range(lay_con$x)); y_rng <- diff(range(lay_con$y))
      x_center <- mean(range(lay_con$x)); y_center <- mean(range(lay_con$y))
      ang_iso <- seq(0, 2*pi, length.out = n_iso + 1)[-1]
      x_all[iso_idx] <- x_center + (x_rng * 0.55) * cos(ang_iso)
      y_all[iso_idx] <- y_center + (y_rng * 0.55) * sin(ang_iso)
    }
    
    igraph::V(g)$comunidad <- factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$grado <- degrees
  }
  
  igraph::E(g)$signo_f <- factor(igraph::E(g)$signo, levels = c("-1", "1"))
  lay <- ggraph::create_layout(g, layout = "manual", x = x_all, y = y_all)
  
  # --- 4. Etiquetas ---
  lab_df <- data.frame(x = x_all, y = y_all, name = igraph::V(g)$name, 
                       grado = igraph::V(g)$grado, com = igraph::V(g)$comunidad, 
                       stringsAsFactors = FALSE)
  
  if (style == "circular") {
    ang_rad <- atan2(lab_df$y, lab_df$x)
    ang_deg <- ang_rad * 180 / pi
    
    lab_df$angle <- ifelse(ang_deg > 90, ang_deg - 180, 
                           ifelse(ang_deg < -90, ang_deg + 180, ang_deg))
    lab_df$hjust <- ifelse(cos(ang_rad) >= 0, 0, 1)
    
    # Aumento de radio a 1.08 para evitar empalmes con los nodos
    lab_df$x_lab <- 1.08 * cos(ang_rad)
    lab_df$y_lab <- 1.08 * sin(ang_rad)
  } else {
    lab_df$angle <- 0
    lab_df$hjust <- 0.5
    lab_df$x_lab <- lab_df$x
    lab_df$y_lab <- lab_df$y
  }
  
  hubs_df <- lab_df[lab_df$grado >= 1, , drop = FALSE] 
  aisl_df <- lab_df[lab_df$com == "Aislado", , drop = FALSE]
  
  # --- 5. Gráfico ---
  p <- ggraph::ggraph(lay) +
    # Aumento de alpha a 0.4 para que las aristas no se pierdan en el circulo
    ggraph::geom_edge_arc(
      ggplot2::aes(colour = signo_f),
      strength = if (style == "stress") 0.15 else -0.05, # Curva más sutil
      alpha    = 0.4, width = 0.5
    ) +
    ggraph::scale_edge_colour_manual(values = c("-1" = "#B2182B", "1" = "#2166AC"), guide = "none") +
    ggraph::geom_node_point(
      ggplot2::aes(colour = comunidad, size = grado), alpha = 0.95, stroke = 0
    ) +
    ggplot2::scale_size_continuous(range = c(2, 10), guide = "none") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    
    # Etiquetas Hubs (sin líneas de guía)
    ggrepel::geom_text_repel(
      data = hubs_df, 
      ggplot2::aes(x = x_lab, y = y_lab, label = name, angle = angle, hjust = hjust),
      size = 2.5, colour = "black", fontface = "bold",
      max.overlaps = Inf, 
      box.padding = if(style == "stress") 0.6 else 0.2, # Más padding en círculo
      point.padding = if(style == "stress") 0.3 else 0.1,
      force = if(style == "stress") 8 else 2, # Más fuerza en círculo
      min.segment.length = Inf, 
      segment.size = 0, 
      seed = seed
    ) +
    
    # Etiquetas Aislados (sin líneas de guía)
    ggrepel::geom_text_repel(
      data = aisl_df, 
      ggplot2::aes(x = x_lab, y = y_lab, label = name, angle = angle, hjust = hjust),
      size = 2.0, colour = "grey30", 
      max.overlaps = Inf, 
      box.padding = if(style == "stress") 0.4 else 0.2,
      point.padding = if(style == "stress") 0.2 else 0.1,
      force = if(style == "stress") 4 else 1,
      min.segment.length = Inf, 
      segment.size = 0,
      seed = seed
    ) +
    
    ggplot2::labs(
      title    = cond, 
      subtitle = sprintf("pres ≥ %d  ·  %d nodos (%d conectados, %d aislados)  ·  %d aristas",
                         pres, n_total, n_con, n_iso, igraph::ecount(g)) # Subtítulo limpio
    ) +
    ggraph::theme_graph(base_family = "sans", base_size = 12) +
    ggplot2::theme(
      legend.position = "none", 
      plot.title      = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5, margin = ggplot2::margin(b = 5)),
      plot.subtitle   = ggplot2::element_text(size = 9, colour = "grey30", hjust = 0.5, margin = ggplot2::margin(b = 10)),
      plot.margin     = ggplot2::margin(0, 0, 0, 0) 
    )
  
  # --- 6. Uso del espacio (Sin expansión) ---
  if (style == "circular") {
    p <- p + ggplot2::coord_fixed(clip = "off", 
                                  xlim = c(-1.3, 1.3), 
                                  ylim = c(-1.3, 1.3),
                                  expand = FALSE)
  } else {
    p <- p + ggplot2::coord_cartesian(clip = "off", expand = FALSE)
  }
  
  p
}

# =============================================================================
# GENERACIÓN Y GUARDADO
# =============================================================================
cat("\n>>> Generando figuras definitivas v11...\n")
p_ss <- plot_net_v11(nets_6$s, "Sano",    6, style = "stress")
p_es <- plot_net_v11(nets_6$e, "Enfermo", 6, style = "stress")
p_sc <- plot_net_v11(nets_6$s, "Sano",    6, style = "circular")
p_ec <- plot_net_v11(nets_6$e, "Enfermo", 6, style = "circular")

save_high_res <- function(p, filename, w, h) {
  try(grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, bg = "white"), silent = TRUE)
  print(p)
  try(grDevices::dev.off(), silent = TRUE)
  
  ggplot2::ggsave(paste0(filename, ".png"), plot = p, width = w, height = h, 
                  dpi = 300, bg = "white", limitsize = FALSE, device = ragg::agg_png)
  cat("  ✓", filename, "\n")
}

# Guardamos los paneles individuales en formato cuadrado para evitar estiramientos
save_high_res(p_ss, "fig_net_Sano_v11",      10, 10)
save_high_res(p_es, "fig_net_Enfermo_v11",   10, 10)
save_high_res(p_sc, "fig_circ_Sano_v11",     10, 10) 
save_high_res(p_ec, "fig_circ_Enfermo_v11",  10, 10) 

# =============================================================================
# PANEL 2x2 FINAL (Formato Cuadrado 16x16 para eliminar espacio a lo ancho)
# =============================================================================
panel <- patchwork::wrap_plots(
  p_ss, p_es, p_sc, p_ec,
  ncol = 2, nrow = 2, byrow = TRUE
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    tag_prefix = "", 
    tag_suffix = ")",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 18),
      plot.margin = ggplot2::margin(2, 2, 2, 2) 
    )
  ) & 
  ggplot2::theme(
    plot.margin = ggplot2::margin(1, 1, 1, 1)
  )

# Guardamos en 16x16 (Cuadrado perfecto)
save_high_res(panel, "fig_redes_panel_v11", 16, 16)


######################################
######################################
# ---- 1. Función de ploteo v13 (con ajustes de fuente y repulsión) -----------
plot_net_v13 <- function(g, cond, pres, style = c("stress", "circular"),
                         top_modulos = 8, seed = 42, diagnosticar = FALSE) {
  style <- match.arg(style)
  if (is.null(g) || igraph::vcount(g) < 2) return(ggplot2::theme_void())
  
  # Fix tipográfico
  igraph::V(g)$name <- gsub("^Pallilotrema$", "Papilliotrema", igraph::V(g)$name)
  
  degrees <- igraph::degree(g)
  n_total <- igraph::vcount(g)
  con_idx <- which(degrees >= 1)
  iso_idx <- which(degrees == 0)
  n_con   <- length(con_idx)
  n_iso   <- length(iso_idx)
  
  g_con <- igraph::induced_subgraph(g, con_idx)
  cl    <- igraph::cluster_louvain(g_con)
  mem   <- cl$membership
  sizes <- sort(table(mem), decreasing = TRUE)
  n_keep <- min(top_modulos, length(sizes))
  tops   <- as.integer(names(sizes))[seq_len(n_keep)]
  
  com_all <- character(n_total)
  com_all[con_idx] <- ifelse(mem %in% tops, paste0("M", mem), "Otros")
  com_all[iso_idx] <- "Aislado"
  
  pal <- setNames(
    c(scales::alpha(grDevices::hcl.colors(n_keep, "Dark 3"), 0.88),
      scales::alpha("#6C3483", 0.9),
      scales::alpha("#A6ACAF", 0.9)),
    c(paste0("M", tops), "Otros", "Aislado")
  )
  
  x_all <- numeric(n_total); y_all <- numeric(n_total)
  
  if (style == "circular") {
    ord <- order(factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado")))
    ang <- seq(0, 2*pi, length.out = n_total + 1)[-1]
    
    x_all[ord] <- cos(ang)
    y_all[ord] <- sin(ang)
    
    if (diagnosticar) {
      diagnostico_redes(g, com_all, x_all, y_all,
                        c(paste0("M", tops), "Otros", "Aislado"))
    }
    
    igraph::V(g)$comunidad <- factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$grado <- degrees
  } else {
    lay_con <- ggraph::create_layout(tidygraph::as_tbl_graph(g_con),
                                     layout = "stress", bbox = 8)
    x_all[con_idx] <- lay_con$x
    y_all[con_idx] <- lay_con$y
    if (n_iso > 0) {
      x_rng <- diff(range(lay_con$x)); y_rng <- diff(range(lay_con$y))
      x_center <- mean(range(lay_con$x)); y_center <- mean(range(lay_con$y))
      ang_iso <- seq(0, 2*pi, length.out = n_iso + 1)[-1]
      x_all[iso_idx] <- x_center + (x_rng * 0.55) * cos(ang_iso)
      y_all[iso_idx] <- y_center + (y_rng * 0.55) * sin(ang_iso)
    }
    igraph::V(g)$comunidad <- factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$grado <- degrees
  }
  
  igraph::E(g)$signo_f <- factor(igraph::E(g)$signo, levels = c("-1", "1"))
  lay <- ggraph::create_layout(g, layout = "manual", x = x_all, y = y_all)
  
  lab_df <- data.frame(x = x_all, y = y_all, name = igraph::V(g)$name,
                       grado = igraph::V(g)$grado, com = igraph::V(g)$comunidad,
                       stringsAsFactors = FALSE) |>
    dplyr::distinct(name, .keep_all = TRUE)
  
  if (style == "circular") {
    ang_rad <- atan2(lab_df$y, lab_df$x)
    ang_deg <- ang_rad * 180 / pi
    lab_df$angle <- ifelse(ang_deg > 90, ang_deg - 180,
                           ifelse(ang_deg < -90, ang_deg + 180, ang_deg))
    lab_df$hjust <- ifelse(cos(ang_rad) >= 0, 0, 1)
    lab_df$x_lab <- 1.15 * cos(ang_rad)
    lab_df$y_lab <- 1.15 * sin(ang_rad)
  } else {
    lab_df$angle <- 0; lab_df$hjust <- 0.5
    lab_df$x_lab <- lab_df$x; lab_df$y_lab <- lab_df$y
  }
  
  hubs_df <- lab_df[lab_df$grado >= 1, , drop = FALSE]
  aisl_df <- lab_df[lab_df$com == "Aislado", , drop = FALSE]
  
  edge_alpha    <- if (style == "stress") 0.3  else 0.15
  edge_strength <- if (style == "stress") 0.15 else -0.05
  
  p <- ggraph::ggraph(lay) +
    ggraph::geom_edge_arc(
      ggplot2::aes(colour = signo_f),
      strength = edge_strength, alpha = edge_alpha, width = 0.4
    ) +
    ggraph::scale_edge_colour_manual(
      values = c("-1" = "#B2182B", "1" = "#2166AC"), guide = "none"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(colour = comunidad, size = grado), alpha = 0.95, stroke = 0
    ) +
    ggplot2::scale_size_continuous(range = c(2, 10), guide = "none") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    
    # --- ETIQUETAS HUBS (Aumento de fuente y repulsión) ---
    ggrepel::geom_text_repel(
      data = hubs_df,
      ggplot2::aes(x = x_lab, y = y_lab, label = name, angle = angle, hjust = hjust),
      size = if (style == "stress") 3.2 else 2.5,       # <-- Fuente más grande en stress
      colour = "black", fontface = "bold",
      max.overlaps = Inf,
      box.padding = if (style == "stress") 0.8 else 0.4, # <-- Más separación
      point.padding = if (style == "stress") 0.5 else 0.2,
      force = if (style == "stress") 15 else 12,         # <-- Más fuerza para evitar empalmes
      min.segment.length = Inf, segment.size = 0,
      seed = seed
    ) +
    
    # --- ETIQUETAS AISLADOS (Aumento de fuente y repulsión) ---
    ggrepel::geom_text_repel(
      data = aisl_df,
      ggplot2::aes(x = x_lab, y = y_lab, label = name, angle = angle, hjust = hjust),
      size = if (style == "stress") 2.6 else 2.0,       # <-- Fuente más grande en stress
      colour = "grey30",
      max.overlaps = Inf,
      box.padding = if (style == "stress") 0.6 else 0.3, # <-- Más separación
      point.padding = if (style == "stress") 0.3 else 0.1,
      force = if (style == "stress") 10 else 5,          # <-- Más fuerza
      min.segment.length = Inf, segment.size = 0,
      seed = seed
    ) +
    ggplot2::labs(
      title = cond,
      subtitle = if (style == "stress") {
        sprintf("pres ≥ %d · %d nodos (%d conectados, %d aislados) · %d aristas",
                pres, n_total, n_con, n_iso, igraph::ecount(g))
      } else { NULL }
    ) +
    ggraph::theme_graph(base_family = "sans", base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      plot.title      = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5,
                                              margin = ggplot2::margin(b = 5)),
      plot.subtitle   = ggplot2::element_text(size = 9, colour = "grey30", hjust = 0.5,
                                              margin = ggplot2::margin(b = 10)),
      plot.margin     = ggplot2::margin(5, 5, 5, 5)
    )
  
  if (style == "circular") {
    p <- p + ggplot2::coord_fixed(clip = "off",
                                  xlim = c(-1.4, 1.4),
                                  ylim = c(-1.4, 1.4),
                                  expand = FALSE)
  } else {
    # En stress quitamos 'expand = FALSE' para darle un pequeño margen natural (5%) 
    # que permite a ggrepel tener más espacio para colocar las etiquetas grandes.
    p <- p + ggplot2::coord_cartesian(clip = "off") 
  }
  p
}

#############################################
# =============================================================================
# SCRIPT COMPLETO v13.1 — Fuentes grandes, repulsión y guardado garantizado
# =============================================================================

# 0. Paquetes necesarios
suppressPackageStartupMessages({
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

# 1. Función de ploteo optimizada
plot_net_v13 <- function(g, cond, pres, style = c("stress", "circular"),
                         top_modulos = 8, seed = 42, diagnosticar = FALSE) {
  style <- match.arg(style)
  if (is.null(g) || igraph::vcount(g) < 2) return(ggplot2::theme_void())
  
  igraph::V(g)$name <- gsub("^Pallilotrema$", "Papilliotrema", igraph::V(g)$name)
  
  degrees <- igraph::degree(g)
  n_total <- igraph::vcount(g)
  con_idx <- which(degrees >= 1)
  iso_idx <- which(degrees == 0)
  n_con   <- length(con_idx)
  n_iso   <- length(iso_idx)
  
  g_con <- igraph::induced_subgraph(g, con_idx)
  cl    <- igraph::cluster_louvain(g_con)
  mem   <- cl$membership
  sizes <- sort(table(mem), decreasing = TRUE)
  n_keep <- min(top_modulos, length(sizes))
  tops   <- as.integer(names(sizes))[seq_len(n_keep)]
  
  com_all <- character(n_total)
  com_all[con_idx] <- ifelse(mem %in% tops, paste0("M", mem), "Otros")
  com_all[iso_idx] <- "Aislado"
  
  pal <- setNames(
    c(scales::alpha(grDevices::hcl.colors(n_keep, "Dark 3"), 0.88),
      scales::alpha("#6C3483", 0.9),
      scales::alpha("#A6ACAF", 0.9)),
    c(paste0("M", tops), "Otros", "Aislado")
  )
  
  x_all <- numeric(n_total); y_all <- numeric(n_total)
  
  if (style == "circular") {
    ord <- order(factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado")))
    ang <- seq(0, 2*pi, length.out = n_total + 1)[-1]
    x_all[ord] <- cos(ang)
    y_all[ord] <- sin(ang)
    if (diagnosticar) diagnostico_redes(g, com_all, x_all, y_all, c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$comunidad <- factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$grado <- degrees
  } else {
    lay_con <- ggraph::create_layout(tidygraph::as_tbl_graph(g_con), layout = "stress", bbox = 8)
    x_all[con_idx] <- lay_con$x
    y_all[con_idx] <- lay_con$y
    if (n_iso > 0) {
      x_rng <- diff(range(lay_con$x)); y_rng <- diff(range(lay_con$y))
      x_center <- mean(range(lay_con$x)); y_center <- mean(range(lay_con$y))
      ang_iso <- seq(0, 2*pi, length.out = n_iso + 1)[-1]
      x_all[iso_idx] <- x_center + (x_rng * 0.55) * cos(ang_iso)
      y_all[iso_idx] <- y_center + (y_rng * 0.55) * sin(ang_iso)
    }
    igraph::V(g)$comunidad <- factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$grado <- degrees
  }
  
  igraph::E(g)$signo_f <- factor(igraph::E(g)$signo, levels = c("-1", "1"))
  lay <- ggraph::create_layout(g, layout = "manual", x = x_all, y = y_all)
  
  lab_df <- data.frame(x = x_all, y = y_all, name = igraph::V(g)$name,
                       grado = igraph::V(g)$grado, com = igraph::V(g)$comunidad,
                       stringsAsFactors = FALSE) |>
    dplyr::distinct(name, .keep_all = TRUE)
  
  if (style == "circular") {
    ang_rad <- atan2(lab_df$y, lab_df$x)
    ang_deg <- ang_rad * 180 / pi
    lab_df$angle <- ifelse(ang_deg > 90, ang_deg - 180, ifelse(ang_deg < -90, ang_deg + 180, ang_deg))
    lab_df$hjust <- ifelse(cos(ang_rad) >= 0, 0, 1)
    lab_df$x_lab <- 1.15 * cos(ang_rad)
    lab_df$y_lab <- 1.15 * sin(ang_rad)
  } else {
    lab_df$angle <- 0; lab_df$hjust <- 0.5
    lab_df$x_lab <- lab_df$x; lab_df$y_lab <- lab_df$y
  }
  
  hubs_df <- lab_df[lab_df$grado >= 1, , drop = FALSE]
  aisl_df <- lab_df[lab_df$com == "Aislado", , drop = FALSE]
  
  edge_alpha    <- if (style == "stress") 0.3  else 0.15
  edge_strength <- if (style == "stress") 0.15 else -0.05
  
  p <- ggraph::ggraph(lay) +
    ggraph::geom_edge_arc(
      ggplot2::aes(colour = signo_f),
      strength = edge_strength, alpha = edge_alpha, width = 0.4
    ) +
    ggraph::scale_edge_colour_manual(values = c("-1" = "#B2182B", "1" = "#2166AC"), guide = "none") +
    ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado), alpha = 0.95, stroke = 0) +
    ggplot2::scale_size_continuous(range = c(2, 10), guide = "none") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    
    # Etiquetas Hubs (Fuente grande, mucha fuerza)
    ggrepel::geom_text_repel(
      data = hubs_df,
      ggplot2::aes(x = x_lab, y = y_lab, label = name, angle = angle, hjust = hjust),
      size = if (style == "stress") 3.2 else 2.5, colour = "black", fontface = "bold",
      max.overlaps = Inf,
      box.padding = if (style == "stress") 0.8 else 0.4,
      point.padding = if (style == "stress") 0.5 else 0.2,
      force = if (style == "stress") 15 else 12,
      min.segment.length = Inf, segment.size = 0, seed = seed
    ) +
    # Etiquetas Aislados
    ggrepel::geom_text_repel(
      data = aisl_df,
      ggplot2::aes(x = x_lab, y = y_lab, label = name, angle = angle, hjust = hjust),
      size = if (style == "stress") 2.6 else 2.0, colour = "grey30",
      max.overlaps = Inf,
      box.padding = if (style == "stress") 0.6 else 0.3,
      point.padding = if (style == "stress") 0.3 else 0.1,
      force = if (style == "stress") 10 else 5,
      min.segment.length = Inf, segment.size = 0, seed = seed
    ) +
    ggplot2::labs(
      title = cond,
      subtitle = if (style == "stress") {
        sprintf("pres ≥ %d · %d nodos (%d conectados, %d aislados) · %d aristas",
                pres, n_total, n_con, n_iso, igraph::ecount(g))
      } else { NULL }
    ) +
    ggraph::theme_graph(base_family = "sans", base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      plot.title      = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5, margin = ggplot2::margin(b = 5)),
      plot.subtitle   = ggplot2::element_text(size = 9, colour = "grey30", hjust = 0.5, margin = ggplot2::margin(b = 10)),
      plot.margin     = ggplot2::margin(5, 5, 5, 5)
    )
  
  if (style == "circular") {
    p <- p + ggplot2::coord_fixed(clip = "off", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), expand = FALSE)
  } else {
    p <- p + ggplot2::coord_cartesian(clip = "off") 
  }
  p
}

# 2. Función de guardado
save_high_res <- function(p, filename, w, h) {
  try(grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, bg = "white"), silent = TRUE)
  print(p)
  try(grDevices::dev.off(), silent = TRUE)
  
  ggplot2::ggsave(paste0(filename, ".png"), plot = p, width = w, height = h,
                  dpi = 300, bg = "white", limitsize = FALSE, device = ragg::agg_png)
  cat("  ✓ Guardado:", filename, "\n")
}

# =============================================================================
# 3. AQUÍ SE GENERAN Y GUARDAN LOS ARCHIVOS
# =============================================================================
cat("\n>>> Generando gráficos...\n")

# Generar los objetos en memoria
p_ss <- plot_net_v13(nets_6$s, "Sano",    6, style = "stress")
p_es <- plot_net_v13(nets_6$e, "Enfermo", 6, style = "stress")
p_sc <- plot_net_v13(nets_6$s, "Sano",    6, style = "circular", diagnosticar = TRUE)
p_ec <- plot_net_v13(nets_6$e, "Enfermo", 6, style = "circular", diagnosticar = TRUE)

# Guardar individuales
save_high_res(p_ss, "fig_net_Sano_v13",     10, 10)
save_high_res(p_es, "fig_net_Enfermo_v13",  10, 10)
save_high_res(p_sc, "fig_circ_Sano_v13",    10, 10)
save_high_res(p_ec, "fig_circ_Enfermo_v13", 10, 10)

# Panel A) + B) con la separación extra
p_ss_ab <- p_ss + ggplot2::theme(plot.margin = ggplot2::margin(5, 18, 5, 5))
p_es_ab <- p_es + ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 18))

panel_AB <- patchwork::wrap_plots(p_ss_ab, p_es_ab, ncol = 2) +
  patchwork::plot_annotation(
    tag_levels = list(c("A)", "B)")),
    theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 18))
  )
save_high_res(panel_AB, "fig_panel_AB_v13", 16, 9)

# Panel C) + D)
panel_CD <- patchwork::wrap_plots(p_sc, p_ec, ncol = 2) +
  patchwork::plot_annotation(
    tag_levels = list(c("C)", "D)")),
    theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 18))
  )
save_high_res(panel_CD, "fig_panel_CD_v13", 16, 9)

# Panel ABCD completo
panel_ABCD <- patchwork::wrap_plots(p_ss, p_es, p_sc, p_ec, ncol = 2, nrow = 2, byrow = TRUE) +
  patchwork::plot_annotation(
    tag_levels = list(c("A)", "B)", "C)", "D)")),
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 18),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
  ) & ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 1, 1))

save_high_res(panel_ABCD, "fig_panel_ABCD_v13", 16, 16)

cat("\n>>> ¡Listo! Revisa tu carpeta de trabajo para ver los archivos .png y .pdf\n")



####################################

# 1. Función de ploteo optimizada
plot_net_v13 <- function(g, cond, pres, style = c("stress", "circular"),
                         top_modulos = 8, seed = 42, diagnosticar = FALSE) {
  style <- match.arg(style)
  if (is.null(g) || igraph::vcount(g) < 2) return(ggplot2::theme_void())
  
  igraph::V(g)$name <- gsub("^Pallilotrema$", "Papilliotrema", igraph::V(g)$name)
  
  degrees <- igraph::degree(g)
  n_total <- igraph::vcount(g)
  con_idx <- which(degrees >= 1)
  iso_idx <- which(degrees == 0)
  n_con   <- length(con_idx)
  n_iso   <- length(iso_idx)
  
  g_con <- igraph::induced_subgraph(g, con_idx)
  cl    <- igraph::cluster_louvain(g_con)
  mem   <- cl$membership
  sizes <- sort(table(mem), decreasing = TRUE)
  n_keep <- min(top_modulos, length(sizes))
  tops   <- as.integer(names(sizes))[seq_len(n_keep)]
  
  com_all <- character(n_total)
  com_all[con_idx] <- ifelse(mem %in% tops, paste0("M", mem), "Otros")
  com_all[iso_idx] <- "Aislado"
  
  pal <- setNames(
    c(scales::alpha(grDevices::hcl.colors(n_keep, "Dark 3"), 0.88),
      scales::alpha("#6C3483", 0.9),
      scales::alpha("#A6ACAF", 0.9)),
    c(paste0("M", tops), "Otros", "Aislado")
  )
  
  x_all <- numeric(n_total); y_all <- numeric(n_total)
  
  if (style == "circular") {
    ord <- order(factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado")))
    ang <- seq(0, 2*pi, length.out = n_total + 1)[-1]
    x_all[ord] <- cos(ang)
    y_all[ord] <- sin(ang)
    if (diagnosticar) diagnostico_redes(g, com_all, x_all, y_all, c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$comunidad <- factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$grado <- degrees
  } else {
    lay_con <- ggraph::create_layout(tidygraph::as_tbl_graph(g_con), layout = "stress", bbox = 8)
    x_all[con_idx] <- lay_con$x
    y_all[con_idx] <- lay_con$y
    if (n_iso > 0) {
      x_rng <- diff(range(lay_con$x)); y_rng <- diff(range(lay_con$y))
      x_center <- mean(range(lay_con$x)); y_center <- mean(range(lay_con$y))
      ang_iso <- seq(0, 2*pi, length.out = n_iso + 1)[-1]
      x_all[iso_idx] <- x_center + (x_rng * 0.55) * cos(ang_iso)
      y_all[iso_idx] <- y_center + (y_rng * 0.55) * sin(ang_iso)
    }
    igraph::V(g)$comunidad <- factor(com_all, levels = c(paste0("M", tops), "Otros", "Aislado"))
    igraph::V(g)$grado <- degrees
  }
  
  igraph::E(g)$signo_f <- factor(igraph::E(g)$signo, levels = c("-1", "1"))
  lay <- ggraph::create_layout(g, layout = "manual", x = x_all, y = y_all)
  
  lab_df <- data.frame(x = x_all, y = y_all, name = igraph::V(g)$name,
                       grado = igraph::V(g)$grado, com = igraph::V(g)$comunidad,
                       stringsAsFactors = FALSE) |>
    dplyr::distinct(name, .keep_all = TRUE)
  
  if (style == "circular") {
    ang_rad <- atan2(lab_df$y, lab_df$x)
    ang_deg <- ang_rad * 180 / pi
    lab_df$angle <- ifelse(ang_deg > 90, ang_deg - 180, ifelse(ang_deg < -90, ang_deg + 180, ang_deg))
    lab_df$hjust <- ifelse(cos(ang_rad) >= 0, 0, 1)
    lab_df$x_lab <- 1.15 * cos(ang_rad)
    lab_df$y_lab <- 1.15 * sin(ang_rad)
  } else {
    lab_df$angle <- 0; lab_df$hjust <- 0.5
    lab_df$x_lab <- lab_df$x; lab_df$y_lab <- lab_df$y
  }
  
  hubs_df <- lab_df[lab_df$grado >= 1, , drop = FALSE]
  aisl_df <- lab_df[lab_df$com == "Aislado", , drop = FALSE]
  
  edge_alpha    <- if (style == "stress") 0.3  else 0.15
  edge_strength <- if (style == "stress") 0.15 else -0.05
  
  p <- ggraph::ggraph(lay) +
    ggraph::geom_edge_arc(
      ggplot2::aes(colour = signo_f),
      strength = edge_strength, alpha = edge_alpha, width = 0.4
    ) +
    ggraph::scale_edge_colour_manual(values = c("-1" = "#B2182B", "1" = "#2166AC"), guide = "none") +
    ggraph::geom_node_point(ggplot2::aes(colour = comunidad, size = grado), alpha = 0.95, stroke = 0) +
    ggplot2::scale_size_continuous(range = c(2, 10), guide = "none") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    
    # Etiquetas Hubs (Fuente grande, mucha fuerza)
    ggrepel::geom_text_repel(
      data = hubs_df,
      ggplot2::aes(x = x_lab, y = y_lab, label = name, angle = angle, hjust = hjust),
      size = if (style == "stress") 3.2 else 2.5, colour = "black", fontface = "bold",
      max.overlaps = Inf,
      box.padding = if (style == "stress") 0.8 else 0.4,
      point.padding = if (style == "stress") 0.5 else 0.2,
      force = if (style == "stress") 15 else 12,
      min.segment.length = Inf, segment.size = 0, seed = seed
    ) +
    # Etiquetas Aislados
    ggrepel::geom_text_repel(
      data = aisl_df,
      ggplot2::aes(x = x_lab, y = y_lab, label = name, angle = angle, hjust = hjust),
      size = if (style == "stress") 2.6 else 2.0, colour = "grey30",
      max.overlaps = Inf,
      box.padding = if (style == "stress") 0.6 else 0.3,
      point.padding = if (style == "stress") 0.3 else 0.1,
      force = if (style == "stress") 10 else 5,
      min.segment.length = Inf, segment.size = 0, seed = seed
    ) +
    ggplot2::labs(
      title = cond,
      subtitle = if (style == "stress") {
        sprintf("pres ≥ %d · %d nodos (%d conectados, %d aislados) · %d aristas",
                pres, n_total, n_con, n_iso, igraph::ecount(g))
      } else { NULL }
    ) +
    ggraph::theme_graph(base_family = "sans", base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      plot.title      = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5, margin = ggplot2::margin(b = 5)),
      plot.subtitle   = ggplot2::element_text(size = 9, colour = "grey30", hjust = 0.5, margin = ggplot2::margin(b = 10)),
      plot.margin     = ggplot2::margin(5, 5, 5, 5)
    )
  
  if (style == "circular") {
    p <- p + ggplot2::coord_fixed(clip = "off", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), expand = FALSE)
  } else {
    p <- p + ggplot2::coord_cartesian(clip = "off") 
  }
  p
}

# 2. Función de guardado
save_high_res <- function(p, filename, w, h) {
  try(grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, bg = "white"), silent = TRUE)
  print(p)
  try(grDevices::dev.off(), silent = TRUE)
  
  ggplot2::ggsave(paste0(filename, ".png"), plot = p, width = w, height = h,
                  dpi = 300, bg = "white", limitsize = FALSE, device = ragg::agg_png)
  cat("  ✓ Guardado:", filename, "\n")
}

# =============================================================================
# 3. AQUÍ SE GENERAN Y GUARDAN LOS ARCHIVOS
# =============================================================================
cat("\n>>> Generando gráficos...\n")

# Generar los objetos en memoria
p_ss <- plot_net_v13(nets_6$s, "Sano",    6, style = "stress")
p_es <- plot_net_v13(nets_6$e, "Enfermo", 6, style = "stress")
p_sc <- plot_net_v13(nets_6$s, "Sano",    6, style = "circular", diagnosticar = TRUE)
p_ec <- plot_net_v13(nets_6$e, "Enfermo", 6, style = "circular", diagnosticar = TRUE)

# Guardar individuales
save_high_res(p_ss, "fig_net_Sano_v13",     10, 10)
save_high_res(p_es, "fig_net_Enfermo_v13",  10, 10)
save_high_res(p_sc, "fig_circ_Sano_v13",    10, 10)
save_high_res(p_ec, "fig_circ_Enfermo_v13", 10, 10)

# Panel A) + B) con la separación extra
p_ss_ab <- p_ss + ggplot2::theme(plot.margin = ggplot2::margin(5, 18, 5, 5))
p_es_ab <- p_es + ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 18))

panel_AB <- patchwork::wrap_plots(p_ss_ab, p_es_ab, ncol = 2) +
  patchwork::plot_annotation(
    tag_levels = list(c("A)", "B)")),
    theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 18))
  )
save_high_res(panel_AB, "fig_panel_AB_v13", 16, 9)

# Panel C) + D)
panel_CD <- patchwork::wrap_plots(p_sc, p_ec, ncol = 2) +
  patchwork::plot_annotation(
    tag_levels = list(c("C)", "D)")),
    theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 18))
  )
save_high_res(panel_CD, "fig_panel_CD_v13", 16, 9)

# Panel ABCD completo
panel_ABCD <- patchwork::wrap_plots(p_ss, p_es, p_sc, p_ec, ncol = 2, nrow = 2, byrow = TRUE) +
  patchwork::plot_annotation(
    tag_levels = list(c("A)", "B)", "C)", "D)")),
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 18),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
  ) & ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 1, 1))

save_high_res(panel_ABCD, "fig_panel_ABCD_v13", 16, 16)

cat("\n>>> ¡Listo! Revisa tu carpeta de trabajo para ver los archivos .png y .pdf\n")




