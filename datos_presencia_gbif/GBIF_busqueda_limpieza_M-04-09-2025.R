#setwd("F:/Doctorado/asesorias")
setwd('C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif')


# librerias

library(readxl)
library(rgbif)
library(dplyr)
library(ggplot2)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(rnaturalearthdata)
library(elevatr)
library(sf)
library(gtools)

################# 1. Obtención de registros biológicos ############################

# carga de datos y obtención de la lista de especies (Giraldo et al. 2018)
lista_sp <- read_excel("./ganaderos_lista.xlsx", sheet="species")
lista_sp <- subset(lista_sp, lista_sp$taxonRank =="species")
lista_sp <- sort(unique(lista_sp$scientificName, decreasing = FALSE))

# consulta de posibles sinónimos de las especies usando GBIF Backbone Taxonomy
# función para buscar sinonimos de cada especie en GBIF
search_synonyms <- function(sp) {
  result <- name_lookup(query = sp, rank = c("species", "subspecies"), status = "SYNONYM")$data
  if (!is.null(result) && nrow(result) > 0) {
    result$scientificName <- sp 
    return(result)
  } else {
    return(NULL)
  }
}

# aplicación de la función a las especies (Giraldo et al. 2018)
synonym_search <- lapply(lista_sp, search_synonyms)
all_synonyms <- bind_rows(synonym_search)
#saveRDS(all_synonyms, "all_synonyms.rds") # guardar los sinÃ³nimos en RDS
all_synonyms <- readRDS("all_synonyms.rds") # cargar los datos para el siguiente procesamientos
# union de listas de especies (Giraldo et al. 2018) y sinÃ³nimos (GBIF)
synonym_list <- unique(sort(all_synonyms$canonicalName))
all_taxa <- sort(unique(c(lista_sp, synonym_list)))
#records <- bind_rows(data_list, .id = "scientificName")

# busqueda de registros geogrÃ¡ficos de cada especie (hasta 5000 registros)
# con coordenadas geogrÃ¡ficas en GBIF

result <- occ_data(scientificName = all_taxa, limit = 5000, hasCoordinate = TRUE)

# procesamiento del objeto gbif_data para obtener un data frame
data_list <- lapply(result, function(x) if (!is.null(x$data)) x$data else NULL)
data_list <- Filter(Negate(is.null), data_list) 
records <- bind_rows(data_list, .id = "scientificName")
# saveRDS(records, "GBIF_2025-03-23.rds")# activar para guardar el archivo en RDS
records <- readRDS("GBIF_2025-03-23.rds") # cargar el archivo para el siguiente proceso

# ActualizaciÃ³n del nombre de las especies segÃºn GBIF Backbone Taxonomy
canonical_names <- all_synonyms$canonicalName[match(records$scientificName, all_synonyms$canonicalName)]
accepted_names <- all_synonyms$scientificName[match(records$scientificName, all_synonyms$canonicalName)]
# Actualizar la columna scientificName en records
records <- records %>%
  mutate(
    scientificName = if_else(!is.na(canonical_names), accepted_names, scientificName)
  )
#
# Revicion de registros #
#
# subset records de Dichotomius agenor 2127 records
sp <- data.frame(subset(records, scientificName == "Dichotomius agenor"))
#
# mapas base
ecoreg <- st_read("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/Ecoregions2017/Ecoregions2017.shp") # Capa de ecoregiones WWF 2017 
ecoreg <- st_make_valid(ecoreg)
col <- ne_countries(country = "colombia", returnclass = "sf")
# col_dem
colombia <- st_as_sf(st_sfc(st_polygon(list(rbind(
  c(-80, -5), c(-80, 15), c(-60, 15), c(-60, -4), c(-80, -5)
))), crs = 4326))
col_dem <- get_elev_raster(locations = colombia, z = 6, clip = "locations")
plot(col_dem, main = "Modelo de Elevación Digital - Colombia")
col_dem <- rast(col_dem)

# Visualización de los registros  
ggplot() +
  geom_spatraster(data = col_dem) +
  scale_fill_viridis_c(name = "Elevación") +
  geom_sf(data = col, fill = "lightgrey", color = "black", alpha = 0.1) +
  geom_point(data = sp, aes(x = decimalLongitude, y = decimalLatitude), color = "darkblue")

## limpieza de puntos atípicos basados en literatura (Montoya et al. 2021)
sp <- sp[!(sp$scientificName=="Dichotomius agenor" & sp$decimalLongitude <= -81),] # remove records greater than  81 Long D agenor
sp <- sp[!(sp$scientificName=="Dichotomius agenor" & # Remove D. agenor departments based on Montoya et al. 2021
             sp$stateProvince %in% c("Boyacá",  "Meta","Casanare", "Vichada", "Boyacá")), ] # 1824 records
sp <- sp[!(sp$decimalLongitude == -76.09989 & sp$decimalLatitude == 8.039917), ] # remove UrabÃ¡ recrord, 1823 records

# limpieza de duplicados 745 records
sp$coordinates <- paste(sp$decimalLatitude, sp$decimalLongitude, sep="_")
sp <- sp[!is.na(sp$coordinates) & !duplicated(sp[c("scientificName", "coordinates")]), ]

####### ESTA PARTE SE PUEDE OMITIR Y CONTINUAR LINEA 141 ################### 
library(rgee)
Sys.setenv(RETICULATE_MINICONDA_PATH = "C:/Users/usuario/miniconda3")
library(reticulate)
use_condaenv("r-reticulate", required = TRUE)
py_config()

# Inicializar GEE usando el entorno correcto de Python
reticulate::use_condaenv("r-reticulate", required = TRUE)

library(rgee)

# Inicializar sesiÃ³n (activar Google Drive para exportaciones)
ee_Initialize(
  project = "beetle-habitat-prediction", # opcional, si tienes proyecto de GCP
  drive = TRUE
)

# Registros a formato GEE
sp_clean <- sp[, c("decimalLongitude", "decimalLatitude", "scientificName")]
sp_sf <- st_as_sf(sp_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
sp_ee <- sf_as_ee(sp_sf)

# ElevaciÃ³n  ALOS PALSAR 30m de GEE
DEM <- ee$Image("JAXA/ALOS/AW3D30/V2_2")$select("AVE_DSM")

# Extraer elevaciÃ³n
elev_ee <- DEM$reduceRegions(
  collection = sp_ee,
  reducer = ee$Reducer$first(),
  scale = 30
)

#  Tranformacion de datos de ee a R
elev_sf <- ee_as_sf(elev_ee)
#saveRDS(elev_sf, "./elev_ALOS.rds")

# Leer el archivo de elevación
elev_sf <- readRDS("./elev_ALOS.rds")

# Revisar nombres de columnas
print(names(elev_sf))  # Debe incluir "elev_ALOS"

# Convertir sp a sf
library(sf)
library(dplyr)

sp_sf <- st_as_sf(sp, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Asegurar que elev_sf está en el mismo CRS
elev_sf <- st_transform(elev_sf, crs = 4326)

# Unir elevaciones al dataset original por proximidad espacial
sp_joined <- st_join(sp_sf, elev_sf["elev_ALOS"], join = st_nearest_feature)

# Convertir a AEA, objeto con coordenadas en metros
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
sp_aea <- st_transform(sp_joined, crs = AEAstring)

# Filtrar puntos por distancia mínima de 1 km
keep <- rep(TRUE, nrow(sp_aea))
for (i in seq_len(nrow(sp_aea))) {
  if (!keep[i]) next
  dists <- st_distance(sp_aea[i, ], sp_aea[(i+1):nrow(sp_aea), ])
  close_points <- which(as.numeric(dists) < 1000)
  if (length(close_points) > 0) keep[(i + close_points)] <- FALSE
}

# Filtrar y volver a data.frame
sp_aea_filtered <- sp_aea[keep, ]
sp_aea_filtered <- st_transform(sp_aea_filtered, crs = 4326)
coords <- st_coordinates(sp_aea_filtered)
sp_aea_filtered$decimalLongitude <- coords[, "X"]
sp_aea_filtered$decimalLatitude <- coords[, "Y"]

sp_df <- as.data.frame(sp_aea_filtered)

# Confirmar que elev_ALOS sigue presente
if (!"elev_ALOS" %in% names(sp_df)) {
  stop("elev_ALOS se perdió en la conversión a data.frame. Revisa antes de filtrar.")
}

# Filtrar por elevaciÃ³n
sp_df <- sp_df %>%
  filter(elev_ALOS >= 0, elev_ALOS <= 1600)
sum(is.na(sp$elev_ALOS))
nrow(sp)  # original
nrow(sp_df) 
geom_spatraster(data = col_dem)   # <-- REVISAR

class(col_dem)
print(col_dem)

ggplot() +
    geom_spatraster(data = col_dem) +
    scale_fill_viridis_c(name = "Elevación") +
  geom_sf(data = col, fill = "gray", color = "black", alpha = 0.1) +
  geom_point(data = sp_df, aes(x = decimalLongitude, y = decimalLatitude), 
             color = "white", size = 3)


###
# extracción del "area de movilidad (M)" basado en regiones biogeográficas (WWWF ecoregions)
###
sp_sf <- st_as_sf(sp, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(ecoreg)) 
eco_col <- st_intersects(ecoreg, sp_sf)
eco_col <- ecoreg[lengths(eco_col) > 0, ]

ggplot() +
  geom_spatraster(data = col_dem) +
  scale_fill_viridis_c(name = "Elevación") +
  geom_sf(data = col, fill = NA, color = "black") +
  geom_sf(data = eco_col, fill = "darkred", color = "black", alpha=0.1) +
  geom_point(data = sp_df, aes(x = decimalLongitude, y = decimalLatitude), 
             color = "white", size = 3)

###
# Area de movilidad basado en buffer 50km 
###
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
sp_df_AEA <- st_as_sf(sp_df, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(ecoreg)) 
sp_df_AEA <- st_transform(sp_df_AEA, crs = AEAstring)

M_buffer50k <- st_buffer(sp_df_AEA, dist = 50000)

ggplot() +
  geom_spatraster(data = col_dem) +
  scale_fill_viridis_c(name = "Elevación") +
  geom_sf(data = col, fill = NA, color = "black") +
  geom_sf(data = M_buffer50k, fill = "darkred", color = "black", alpha=0.7) +
  geom_sf(data = sp_df_AEA, color = "white", size = 3)

M_buffer50k_union <- st_union(M_buffer50k)

ggplot() +
  geom_spatraster(data = col_dem) +
  scale_fill_viridis_c(name = "Elevación") +
  geom_sf(data = col, fill = NA, color = "black") +
  geom_sf(data = col, fill = NA, color = "black") +
  geom_sf(data = M_buffer50k_union, fill = "darkred", color = "black", alpha=0.7) +
  geom_point(data = sp, aes(x = decimalLongitude, y = decimalLatitude), color = "white", size=3)

#################################################################################################
######### 6. Enmascarar capas ambientales según área de calibración #############################
#################################################################################################

# --- . Guardar las áreas de calibración (M) ---
st_write(
  eco_col,
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_ecoregiones.shp",
  delete_layer = TRUE
)

st_write(
  M_buffer50k_union,
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_buffer50km.shp",
  delete_layer = TRUE
)

# --- . Cargar variables bioclimáticas de WorldClim 2.1 ---
archivos <- list.files(
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/WORLDCLIM",
  pattern = "tif$",
  full.names = TRUE
)

# Ordenar archivos de manera numérica (bio_1, bio_2, ..., bio_19)
archivos <- mixedsort(archivos)

variables <- terra::rast(archivos)

################# CARGAR VARIABLES DE SUELO 14/OCT/2025 #################################################
# =========================
# A) RECREAR bbox_vect desde tus presencias
# =========================
library(sf); library(terra)

# Asegúrate de tener sp_df_AEA en memoria:
stopifnot(exists("sp_df_AEA"))

# Pasa tus puntos a WGS84 (lon/lat)
sp_lonlat <- st_transform(sp_df_AEA, 4326)
xy <- st_coordinates(sp_lonlat)

# Min/max + acolchado
pad_deg <- 0.5
lon_min <- min(xy[,1], na.rm = TRUE) - pad_deg
lon_max <- max(xy[,1], na.rm = TRUE) + pad_deg
lat_min <- min(xy[,2], na.rm = TRUE) - pad_deg
lat_max <- max(xy[,2], na.rm = TRUE) + pad_deg

ext_bbox  <- ext(lon_min, lon_max, lat_min, lat_max)
bbox_vect <- as.polygons(ext_bbox, crs = "EPSG:4326")

# Verificación rápida
# print(bbox_vect)


# =========================
# B) RUTAS (LAS TUYAS)
# =========================
ROOT_SWC_DIR   <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/swc_fr"
ROOT_ALPHA_DIR <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/ALPHA"

OUT_SWC_TIF    <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/swc_fr/swc_fr_bbox_30s.tif"
OUT_ALPHA_TIF  <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/ALPHA/alpha_bbox_30s.tif"


# =========================
# C) FUNCIÓN ROBUSTA
# =========================

process_esri_grid_to_tif <- function(root_dir, clip_geom, out_tif) {
  message("\n== Procesando: ", root_dir, " ==")
  all_dirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
  
  # filtra tiles válidos (ignora 'info' y exige archivos de GRID)
  is_grid_dir <- function(d) {
    if (grepl("[/\\\\]info$", d, ignore.case = TRUE)) return(FALSE)
    files <- try(list.files(d, full.names = TRUE), silent = TRUE)
    if (inherits(files, "try-error")) return(FALSE)
    any(grepl("hdr\\.adf$", files, ignore.case = TRUE)) ||
      any(grepl("w001001\\.adf$", files, ignore.case = TRUE))
  }
  tiles <- Filter(is_grid_dir, all_dirs)
  
  # si el root_dir en sí mismo es el GRID (caso alpha), añádelo
  if (is_grid_dir(root_dir)) tiles <- c(root_dir, tiles)
  
  if (length(tiles) == 0) stop("No se encontraron tiles ESRI GRID válidos en: ", root_dir)
  
  cropped_list <- list()
  for (td in tiles) {
    message("  . Tile: ", td)
    r <- try(rast(td), silent = TRUE)
    if (inherits(r, "try-error")) { message("    (saltado: no se pudo leer)"); next }
    cg <- if (!identical(crs(clip_geom), crs(r))) project(clip_geom, crs(r)) else clip_geom
    
    rc <- try(crop(r, cg), silent = TRUE)
    if (inherits(rc, "try-error") || is.null(rc) || ncell(rc) == 0) {
      message("    (no cruza el área)"); next
    }
    message("    ??? cruza y se recortó")
    cropped_list[[length(cropped_list) + 1]] <- rc
  }
  
  if (length(cropped_list) == 0) stop("Ningún tile cruzó el área de recorte para: ", root_dir)
  
  message("  . Mosaicar/escribir...")
  dir.create(dirname(out_tif), showWarnings = FALSE, recursive = TRUE)
  
  # --- manejar 1 vs >1 raster ---
  if (length(cropped_list) == 1) {
    mos <- cropped_list[[1]]
  } else {
    mos <- do.call(merge, cropped_list)
  }
  
  writeRaster(mos, out_tif, overwrite = TRUE)
  message("== OK ??? ", out_tif, " ==")
  invisible(out_tif)
}

# =========================
# D) EJECUCIÓN
# =========================
# Recomendación: corre primero swc_fr
process_esri_grid_to_tif(ROOT_SWC_DIR,   bbox_vect, OUT_SWC_TIF)
# Luego alpha
process_esri_grid_to_tif(ROOT_ALPHA_DIR, bbox_vect, OUT_ALPHA_TIF)


# =========================
# E) QC RÁPIDO
# =========================
swc <- rast(OUT_SWC_TIF);  alp <- rast(OUT_ALPHA_TIF)
print(res(swc))   # ~ 0.008333333 (30")
print(res(alp))
print(global(swc, fun=c("min","max"), na.rm=TRUE))
print(global(alp, fun=c("min","max"), na.rm=TRUE))

############################## CARGANDO .TIF DE SUELOS #####################################
# Rutas a tus TIF recién creados
swc_tif <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/swc_fr/swc_fr_bbox_30s.tif"
alp_tif <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/ALPHA/alpha_bbox_30s.tif"

# 0) Plantilla: usa la malla/extent/CRS de tu stack 'variables' (WorldClim 30")
tmpl <- variables[[1]]

# 1) Leer
swc <- rast(swc_tif)  # humedad del suelo (0-1)
alp <- rast(alp_tif)  # alpha = AET/PET (0-1)

# 2) Alinear exactamente al grid de 'variables' (mismo CRS, extent, resolución)
to_tmpl <- function(r){
  if (!identical(crs(r), crs(tmpl))) r <- project(r, crs(tmpl), method = "bilinear")
  resample(r, tmpl, method = "bilinear")
}
swc30 <- to_tmpl(swc)
alp30 <- to_tmpl(alp)

# 3) Derivar estrés hídrico (más interpretable para el modelo)
stress30 <- 1 - alp30     # 0 = sin estrés; 1 = máximo estrés

# 4) Nombres y apilado
names(swc30)    <- "swc_fr"
names(stress30) <- "stress"
variables <- c(variables, swc30, stress30)  # <-- ahora 'variables' incluye 2 capas nuevas

####################################### FIN ###################################################

##############################################################################################
### BLOQUE 6.1 PROCESAMIENTO USANDO ECOREGIONES ###############################################
##############################################################################################
cat("\n=== BLOQUE 6.1: Procesamiento por ECOREGIONES ===\n")

# --- 6.1.1. Cargar M de ecoregiones ---
M_ecoreg <- st_read(
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_ecoregiones.shp"
)
M_ecoreg_vect <- vect(M_ecoreg)
M_ecoreg_vect <- terra::project(M_ecoreg_vect, crs(variables))  # Asegurar CRS

# ---6.1.2. Recortar y enmascarar ---
variables_ecoreg_crop <- crop(variables, M_ecoreg_vect)
variables_ecoreg_mask <- mask(variables_ecoreg_crop, M_ecoreg_vect)

# ---6.1.3. Guardar capas recortadas ---
dir.create("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/variables_M_ecoregiones", showWarnings = FALSE, recursive = TRUE)

terra::writeRaster(
  variables_ecoreg_mask,
  filename = paste0(
    "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/variables_M_ecoregiones/",
    names(variables_ecoreg_mask), ".tif"
  ),
  overwrite = TRUE
)

# ---6.1.4. Extraer valores de presencia ---
sp_vect <- vect(sp_df, geom = c("decimalLongitude", "decimalLatitude"), crs = crs(variables))
puntos_ecoreg <- terra::extract(variables_ecoreg_mask, sp_vect)
puntos_ecoreg <- puntos_ecoreg %>% dplyr::select(-ID) #Eliminar la columna ID que genera el desplazamiento
presencias_ecoreg <- cbind(sp_df, puntos_ecoreg)
# Guardar CSV
# --- Convertir a data.frame eliminando columnas problemáticas ---
#presencias_ecoreg_clean <- presencias_ecoreg_clean %>%
 # st_drop_geometry() %>%  # elimina columna geometry
  #dplyr::select(-dnaSequenceID, -networkKeys)  # elimina listas

### Extraemos columnas que necesitamos y guardamos en csv
# Vector con columnas a exportar
# Seleccionar solo las columnas que necesitas
df_export <- presencias_ecoreg %>%
  st_drop_geometry() %>%
  dplyr::select(
    scientificName,
    decimalLongitude,
    decimalLatitude,
    dplyr::starts_with("wc2.1_30s_bio_"),
    dplyr::any_of(c("swc_fr", "stress"))
  )

write.csv2(df_export,
           "presencias_ecoreg_limpio.csv",
           row.names = FALSE,
           fileEncoding = "UTF-8")


cat("??? Proceso de ECOREGIONES completado. Filas exportadas:",
    nrow(presencias_ecoreg), "\n")

# ---6.1.5. Verificación visual rápida ---
plot(variables_ecoreg_mask[[1]], main = "BIO1 - Recortado por ecoregiones")
plot(M_ecoreg_vect, add = TRUE)

################################################################################
### BLOQUE 6.2 PROCESAMIENTO USANDO BUFFER 50 KM ###############################
################################################################################
cat("\n=== BLOQUE 6.2: Procesamiento por BUFFER 50 km ===\n")

# ---6.2.1. Cargar M de buffer ---
M_buffer <- st_read(
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_buffer50km.shp"
)
M_buffer_vect <- vect(M_buffer)
M_buffer_vect <- terra::project(M_buffer_vect, crs(variables))  # Asegurar CRS

# ---6.2.2. Recortar y enmascarar ---
variables_buffer_crop <- terra::crop(variables, M_buffer_vect)
variables_buffer_mask <- terra::mask(variables_buffer_crop, M_buffer_vect)

# ---6.2.3. Guardar capas recortadas ---
dir.create(
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/variables_M_buffer",
  showWarnings = FALSE, recursive = TRUE
)

terra::writeRaster(
  variables_buffer_mask,
  filename = paste0(
    "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/variables_M_buffer/",
    names(variables_buffer_mask), ".tif"
  ),
  overwrite = TRUE
)

# ---6.2.4. Extraer valores de presencia ---
puntos_buffer <- terra::extract(variables_buffer_mask, sp_vect)

#  Eliminar columna ID si existe (evita desplazamiento de columnas)
if ("ID" %in% colnames(puntos_buffer)) {
  puntos_buffer <- dplyr::select(puntos_buffer, -ID)
}

presencias_buffer <- cbind(sp_df, puntos_buffer)


# --- Guardar CSV ---
# Seleccionar solo las columnas que necesitas
df_export_buffer <- presencias_buffer %>%
  st_drop_geometry() %>%
  dplyr::select(
    scientificName,
    decimalLongitude,
    decimalLatitude,
    dplyr::starts_with("wc2.1_30s_bio_"),
    dplyr::any_of(c("swc_fr", "stress"))   # <-- añade humedad y estrés
  )

# Exportar limpio
write.csv2(df_export_buffer,   # usa ; como separador (Excel en español lo abre directo)
           "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/presencias_buffer_limpio.csv",
           row.names = FALSE,
           fileEncoding = "UTF-8")


cat("??? CSV de buffer (50 km) exportado correctamente.\n")

# ---6.2.5. Verificación visual rápida ---
plot(variables_buffer_mask[[1]], main = "BIO1 - Recortado por buffer")
plot(M_buffer_vect, add = TRUE)

#VALIDACION
names(variables)  # debe listar ... "swc_fr", "stress"
global(variables[["swc_fr"]],  fun=c("min","max"), na.rm=TRUE)
global(variables[["stress"]],  fun=c("min","max"), na.rm=TRUE)


cat("CSV generados en la carpeta 'datos_presencia_gbif'\n")

############################ PROCESO DE PCA Y VIF #####################################
library(usdm)   # para VIF, si quieres probar ambos enfoques

# =========================================================
# REDUCCIÓN DE VARIABLES: CORRELACIÓN, VIF Y PCA (MEMORY-SAFE)
# =========================================================
# Paquetes -------------------------------------------------
# install.packages(c("terra","usdm","corrplot"), dependencies = TRUE)

suppressPackageStartupMessages({
  if (!requireNamespace("usdm", quietly = TRUE)) message(">> (Opcional) Instala 'usdm' para VIF: install.packages('usdm')")
  if (!requireNamespace("corrplot", quietly = TRUE)) message(">> (Opcional) Instala 'corrplot' para el gráfico de correlaciones: install.packages('corrplot')")
})

# Parámetros ----------------------------------------------
# Carpeta de salida (ajusta si prefieres otra)
OUT_DIR <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/OUT"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Tamaño de muestra para análisis tabulares (cor/VIF)
SAMPLE_SIZE <- 10000   # 10k suele ser suficiente; sube/baja según tu RAM

# Nº de componentes principales a exportar como mapas
NCOMP <- 6             # ajusta según scree/varianza explicada

# Opciones de terra para manejar memoria
terraOptions(memfrac = 0.6)   # usa ~60% de la RAM
terraOptions(todisk  = TRUE)  # operaciones por bloques con disco

# =========================
# 0) Sanity check del stack
# =========================
stopifnot(inherits(variables, "SpatRaster"))
message("Capas en 'variables': ", nlyr(variables))
message("Nombres: ", paste(names(variables), collapse = ", "))

# =========================
# 1) Muestra aleatoria (evita std::bad_alloc)
# =========================
set.seed(123)
message("Tomando muestra aleatoria de ", SAMPLE_SIZE, " píxeles...")
samp <- spatSample(variables, size = SAMPLE_SIZE, method = "random", na.rm = TRUE, as.df = TRUE)
# 'samp' es un data.frame con columnas = capas

####### PARA GUARDAR EL SAMPLE Y NO TENER QUE CORRER EL spatSample() #######
# Guardar la muestra en CSV  MAYRAAAAAAAAAAAAAAAAAAA CORRER DESDE AQUI
library(readr)
#La siguiente linea solo se corre una vez, en las proximas ejecuciones solo se hace el read()
write_csv(samp, "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/OUT/variables_sample_10k.csv")
cat("??? Muestra guardada: variables_sample_10k.csv\n")

vars_sample <- read_csv("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/OUT/variables_sample_10k.csv")
##########################################################################

# =========================
# 2) Correlaciones
# =========================
message("Calculando matriz de correlaciones (Pearson) sobre la muestra...")
cor_mat <- cor(samp, method = "pearson", use = "complete.obs")

# Guardar CSV de correlación
cor_csv <- file.path(OUT_DIR, "correlation_matrix_sample.csv")
write.csv(cor_mat, cor_csv, row.names = TRUE, fileEncoding = "UTF-8")
message(">> Matriz de correlaciones guardada: ", cor_csv)

# (Opcional) Gráfico de correlación si 'corrplot' está disponible
if (requireNamespace("corrplot", quietly = TRUE)) {
  png(file.path(OUT_DIR, "correlation_plot_sample.png"), width=1600, height=1200, res=180)
  corrplot::corrplot(cor_mat, type = "upper", tl.cex = 0.6)
  dev.off()
  message(">> Gráfico de correlación guardado: ", file.path(OUT_DIR, "correlation_plot_sample.png"))
}

# =========================
# 3) VIF (opcional)
# =========================
# Si quieres reducir manteniendo variables originales, usa VIF sobre la muestra:
vars_vif_names <- names(variables)  # por defecto, conserva todas

if (requireNamespace("usdm", quietly = TRUE)) {
  message("Aplicando VIF (umbral 10) sobre la muestra...")
  # usdm::vifstep requiere un data.frame sin NAs y solo numérico
  df_vif <- as.data.frame(stats::na.omit(samp))
  vif_obj <- usdm::vifstep(df_vif, th = 10)
  message(">> Variables seleccionadas por VIF: ", paste(vif_obj@results$Variables, collapse = ", "))
  vars_vif_names <- vif_obj@results$Variables
} else {
  message("Saltando VIF (paquete 'usdm' no disponible). Se mantienen todas las capas.")
}

# Subconjunto de raster según VIF (o todas si no se corrió VIF)
variables_vif <- variables[[vars_vif_names]]
writeRaster(variables_vif, file.path(OUT_DIR, "variables_reducidas_VIF.tif"), overwrite = TRUE)
message(">> Stack reducido por VIF (o todas si no hubo VIF) guardado en: ",
        file.path(OUT_DIR, "variables_reducidas_VIF.tif"))

# =========================
# 4) PCA con terra::pca (block-wise)
# =========================
# Nota: 'terra::pca' es eficiente en memoria y produce:
#  - $map : SpatRaster con los componentes (PC1..PCn)
#  - $model : lista con cargas (loadings) y sdev (desv. est.)
message("Calculando PCA (terra::pca) con ", NCOMP, " componentes sobre el stack reducido...")
pca_res <- terra::pca(variables_vif, ncomp = NCOMP, center = TRUE, scale = TRUE)

# Guardar mapas de componentes
pca_tif <- file.path(OUT_DIR, sprintf("PCA_components_%d.tif", NCOMP))
writeRaster(pca_res$map, pca_tif, overwrite = TRUE)
message(">> Mapas de componentes guardados: ", pca_tif)

# Varianza explicada por componente (a partir de sdev)
sdev <- pca_res$model$sdev
var_exp <- (sdev^2) / sum(sdev^2)
cum_var <- cumsum(var_exp)
pca_var_df <- data.frame(
  PC   = paste0("PC", seq_along(sdev)),
  SDev = sdev,
  VarExplained = var_exp,
  CumVarExplained = cum_var
)

pca_var_csv <- file.path(OUT_DIR, "PCA_variance_explained.csv")
write.csv(pca_var_df, pca_var_csv, row.names = FALSE, fileEncoding = "UTF-8")
message(">> Varianza explicada guardada: ", pca_var_csv)

# (Opcional) Guardar loadings (cargas) por variable para interpretar PCs
loadings <- pca_res$model$co
load_csv <- file.path(OUT_DIR, "PCA_loadings.csv")
write.csv(loadings, load_csv, row.names = TRUE, fileEncoding = "UTF-8")
message(">> Loadings guardados: ", load_csv)

# =========================
# 5) (Opcional) Elegir nº de PCs por umbral de varianza
# =========================
# Si prefieres fijar PCs por varianza acumulada (p.ej., 90%):
UMBRAL_VAR <- 0.90
k_auto <- which(cum_var >= UMBRAL_VAR)[1]
message(sprintf("Con %.0f%% de varianza acumulada, bastan %d componentes.", UMBRAL_VAR*100, k_auto))

# Puedes escribir solo las primeras k_auto capas:
pca_auto_tif <- file.path(OUT_DIR, sprintf("PCA_components_auto_%dPCs.tif", k_auto))
writeRaster(pca_res$map[[1:k_auto]], pca_auto_tif, overwrite = TRUE)
message(">> Mapas de PCs automáticos guardados: ", pca_auto_tif)

# =========================
# 6) Sugerencias para MaxEnt
# =========================
message("\nSugerencias:")
message("- Si quieres interpretar variables originales: usa 'variables_reducidas_VIF.tif'.")
message("- Si prefieres componentes (decorrelados): usa '", basename(pca_auto_tif), "' o '", basename(pca_tif), "'.")
message("- Revisa 'PCA_variance_explained.csv' y 'PCA_loadings.csv' para elegir PCs y entender sus ejes.")

