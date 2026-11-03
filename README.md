# Análisis del Microbioma Fúngico y Modelado de Nicho de *Quercus robur* (Sano vs. Enfermo)

Este repositorio contiene el código y los resultados del análisis de la comunidad fúngica asociada a *Quercus robur* en el suroeste de Francia, comparando árboles sanos y enfermos. El flujo integra secuenciación de alto rendimiento (ITS/LSU), análisis de redes de co-ocurrencia, asignación ecológica (FUNGuild) y modelado de nicho ecológico (ENM) presente y futuro.

## 1. Sitio de Estudio y Fuente de Datos

El área de muestreo se ubica en la región de Nouvelle-Aquitaine (Suroeste de Francia), un paisaje en mosaico compuesto por viñedos y bosques templados. Los datos originales de secuenciación provienen del BioProject **PRJNA508496** (Laurent et al., 2020), del cual se extrajeron exclusivamente las muestras correspondientes a *Quercus robur* para los marcadores ITS y LSU.

![Mapa de Estudio](https://raw.githubusercontent.com/migueangelhg/InteraccionesOrganismoAmbiente/main/resultados/Figuras/SitiodeEstudioQrobur.png)

**Mapa Interactivo (Leaflet):** 
El repositorio incluye un mapa interactivo generado con la librería Leaflet, que permite explorar las muestras con capas base de OpenStreetMap y satelitales, además de un minimapa de contexto. 

*Nota: GitHub no renderiza archivos HTML interactivos de forma nativa por seguridad. Para visualizarlo, haz clic en el siguiente enlace, descarga el archivo y ábrelo localmente:*

**[Descargar Mapa Interactivo de Leaflet](https://github.com/migueangelhg/InteraccionesOrganismoAmbiente/raw/main/resultados/Elipsoides_Qrobur/mapa_interactivo_con_inset.html)**

## 2. Metodología

El análisis se dividió en cinco fases principales, ejecutadas a través de scripts modulares en R:

### 2.1. Obtención y Control de Calidad de Datos
*   **Descarga:** Los archivos FASTQ crudos se descargaron desde la API del ENA mediante `curl::multi_download()`, con verificación de integridad y re-descarga automatizada de archivos truncados.
*   **Filtrado y Generación de ASVs:** Se utilizó la suite **DADA2** (Callahan et al., 2016). Las lecturas fueron filtradas con `filterAndTrim` utilizando parámetros conservadores para ITS (`truncLen = c(240, 200)`, `maxEE = c(2, 5)`, `truncQ = 2`). Se eliminaron quimeras y se generaron las tablas de Amplicon Sequence Variants (ASVs). La taxonomía se asignó utilizando la base de datos **UNITE**.

### 2.2. Asignación Ecológica y Funcional
*   **FUNGuild:** Se utilizó la base de datos de FUNGuild (Nguyen et al., 2016) mediante el paquete `MiscMetabar`. Se filtró a nivel taxonómico de Género (Nivel 13). Los gremios ecológicos complejos (ej. "Plant Pathogen-Wood Saprotroph") se desdoblaron y se mapearon a categorías funcionales amplias (Patógeno de planta, Saprótrofo de madera, Endófito, etc.) aplicando una jerarquía de prioridad para evitar redundancias.

### 2.3. Análisis de Diversidad y Estructura Comunitaria
*   **Diversidad Alfa:** Se calcularon los índices de Shannon, Simpson, Riqueza observada y Chao1. Las diferencias entre condiciones (Sano vs. Enfermo) se evaluaron mediante la prueba no paramétrica de **Wilcoxon**.
*   **Diversidad Beta:** Se calculó la distancia de Bray-Curtis. La estructura comunitaria se evaluó con un **NMDS** (Escalamiento Multidimensional No Métrico). Las diferencias significativas se testearon mediante **PERMANOVA** (`vegan::adonis2`, 999 permutaciones) y se verificó la homogeneidad de dispersión con **Betadisper** (Oksanen et al., 2020).
*   **Partición de Beta:** Se utilizó la partición de Baselga (Sørensen) para descomponer la beta diversidad en dos componentes: **Recambio** (turnover) y **Anidamiento** (nestedness) (Baselga, 2010).

### 2.4. Redes de Co-ocurrencia
*   **Transformación y Correlación:** Las abundancias de ASVs se agregaron a nivel de Género y se transformaron mediante *Centered Log-Ratio* (CLR) con un pseudo-conteo de 0.5 para manejar la composicionalidad de los datos. Se calculó la matriz de correlación de Pearson.
*   **Filtros de Red:** Se retuvieron las correlaciones con un valor absoluto `|r| >= 0.40`. Para evitar falsos positivos por ausencia de datos (sparsity), se aplicó un filtro de *dz-avoidance* (`dz <= 0.60`), eliminando conexiones entre taxones que co-ocurren en muy pocas muestras.
*   **Topología:** Las redes resultantes se analizaron con `igraph` (Csárdi & Nepusz, 2006) calculando modularidad (algoritmo de Louvain), grado medio, centralidad y detectando taxones *hub*.

### 2.5. Modelado de Nicho Ecológico (ENM)
*   **Elipsoides de Mahalanobis:** Se ajustaron elipsoides al 95% de confianza utilizando el paquete `ntbox::cov_center()` (Osorio-Olvera et al., 2020). Se generaron tres modelos: el nicho de la **especie** (usando registros globales de GBIF filtrados espacialmente y con `CoordinateCleaner`), el nicho de la **comunidad sana** y el de la **comunidad enferma**.
*   **Predictores Climáticos:** Se descargaron las 19 variables bioclimáticas de **WorldClim 2.1** (resolución de 30 arc-sec) y se aplicó un filtro de colinealidad (eliminando variables con `|r| > 0.7`) (Fick & Hijmans, 2017).
*   **Solapamiento y Proyección:** Se calculó el solapamiento geométrico mediante el índice de **Jaccard** (muestreando 200,000 celdas ambientales). Para el futuro, se proyectaron los nichos sobre el modelo **CMIP6 ACCESS-CM2** bajo el escenario **SSP245** (2041-2060), calculando el cambio en la idoneidad climática (Futuro - Presente).

## 3. Resultados

### 3.1. Diversidad Alfa y Beta
Los árboles sanos presentan valores significativamente más altos de diversidad (índices de Shannon y Simpson) en comparación con los enfermos. Sin embargo, la riqueza de especies (Chao1 y Riqueza observada) se mantiene estadísticamente similar.

El NMDS muestra una clara separación entre ambas condiciones (**PERMANOVA R² = 0.32, p = 0.001**). La partición de Sørensen revela que esta diferencia se debe casi en su totalidad al **anidamiento** (91-94%), indicando que las comunidades enfermas son un subconjunto empobrecido de las comunidades sanas, más que un recambio total de especies.

![Diversidad Alfa y Beta](https://raw.githubusercontent.com/migueangelhg/InteraccionesOrganismoAmbiente/main/resultados/Figuras/fig_alpha_nmds_beta.png)

### 3.2. Asignación Ecológica (FUNGuild)
Los gremios ecológicos dominantes, tanto en abundancia de ASVs como en riqueza de géneros, son:
1. **Patógenos de plantas** (26.1% de abundancia).
2. **Saprótrofos de madera** (19.8%).
3. **Saprótrofos de suelo/hojarasca** (17.8%).
4. **Endófitos** (15.7%).

![Gremios Ecológicos](https://raw.githubusercontent.com/migueangelhg/InteraccionesOrganismoAmbiente/main/resultados/Figuras/ITS_grupos_ecologicos_guildPie.png)

### 3.3. Redes de Co-ocurrencia
*   **Red Sana:** Presenta una topología menos conectada (53 aristas), con una estructura modular clara y presencia de nodos aislados, lo que sugiere una comunidad compartimentada y estable.
*   **Red Enferma:** Se observa una red altamente conectada (116 aristas), con un núcleo denso que sugiere una respuesta ecológica desregulada, pérdida de compartimentación o un incremento en la competencia por recursos.

![Redes de Co-ocurrencia Sano vs Enfermo](https://raw.githubusercontent.com/migueangelhg/InteraccionesOrganismoAmbiente/main/resultados/Figuras/fig_panel_AB_v13.png)

### 3.4. Modelado de Nicho Ecológico (ENM)
El análisis de solapamiento (Jaccard) indica que, aunque las comunidades sana y enferma tienen un alto solapamiento entre sí, ambas ocupan una fracción muy específica y restringida del nicho climático global de la especie hospedadora (*Q. robur*).

<img src="https://raw.githubusercontent.com/migueangelhg/InteraccionesOrganismoAmbiente/main/resultados/Figuras/plot_solapamiento_3_modelos.png" alt="Elipsoides de Nicho" width="60%">

**Proyección Futura (SSP245 - 2050):** El modelo predice una drástica reducción y desplazamiento del área climáticamente idónea. Las zonas de alta idoneidad se contraen y migran hacia el norte de Europa, lo que podría comprometer la estabilidad del holobionte planta-hongo en su distribución actual.

![Cambio Climático Futuro](https://raw.githubusercontent.com/migueangelhg/InteraccionesOrganismoAmbiente/main/resultados/Figuras/panel_especie_presente_futuro_cambio.png)

## 4. Referencias

- Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. *Global Ecology and Biogeography*, 19(1), 134–143.
- Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13(7), 581–583.
- Csárdi, G., & Nepusz, T. (2006). The igraph software package for complex network research. *InterJournal, Complex Systems*, 1695(5), 1–9.
- Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: New 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, 37(12), 4302–4315.
- Laurent, B., Marchand, M., Chancerel, E., Saint-Jean, G., Capdevielle, X., Poeydebat, C., … & Desprez-Loustau, M. L. (2020). A richer community of Botryosphaeriaceae within a less diverse community of fungal endophytes in grapevines than in adjacent forest trees revealed by a mixed metabarcoding strategy. *Phytobiomes Journal*, 4(3), 252–267.
- Nguyen, N. H., Song, Z., Bates, S. T., Branco, S., Tedersoo, L., Menke, J., … & Kennedy, P. G. (2016). FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild. *Fungal Ecology*, 20, 241–248.
- Oksanen, J., Blanchet, F. G., Friendly, M., Kindt, R., Legendre, P., McGlinn, D., … & Wagner, H. (2020). *vegan: Community Ecology Package*. R package version 2.5-7.
- Osorio-Olvera, L., Lira-Noriega, A., Soberón, J., Peterson, A. T., Falconi, M., Contreras-Díaz, R. G., … & Barve, N. (2020). ntbox: An R package with graphical user interface for modelling and evaluating multidimensional ecological niches. *Methods in Ecology and Evolution*, 11(10), 1199–1206.
- Eyring, V., Bony, S., Meehl, G. A., Senior, C. A., Stevens, B., Stouffer, R. J., & Taylor, K. E. (2016). Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. *Geoscientific Model Development*, 9(5), 1937–1958.
