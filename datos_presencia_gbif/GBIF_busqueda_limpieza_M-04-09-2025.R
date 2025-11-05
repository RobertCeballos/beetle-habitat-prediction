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
library(readr)
library(usdm)
library(stats)

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

###################### CARGAR ARCHIVO DE ELEVACIONES ####################### 
############################################################################
elev_sf <- readRDS("./elev_ALOS.rds")

# Revisar nombres de columnas
print(names(elev_sf))  # Debe incluir "elev_ALOS"

# Convertir sp a sf
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
# B) RUTAS PARA CARGA DE ARCHIVOS Y SALIDA
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

# ===== ELIMINAR REGISTROS CON NA EN VARIABLES CLIMÁTICAS =====
# Detecta automáticamente columnas de bioclimáticas (por nombre)
cols_clima <- grep("bio", names(presencias_ecoreg), value = TRUE, ignore.case = TRUE)

if (length(cols_clima) == 0) {
  stop("No se encontraron columnas de variables bioclimáticas en presencias_ecoreg.")
}

n_bioclim <- length(cols_clima)

# Eliminar registros donde TODAS las bioclimáticas son NA
presencias_ecoreg <- presencias_ecoreg %>%
  filter(rowSums(is.na(across(all_of(cols_clima)))) < n_bioclim)

message("??? Registros válidos conservados: ", nrow(presencias_ecoreg))


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

# ===== ELIMINAR REGISTROS CON TODAS LAS BIOCLIMÁTICAS EN NA (BUFFER) =====
# Detecta columnas bioclimáticas de forma flexible
cols_clima <- grep("(^wc2\\.1_30s_)?bio", names(presencias_buffer),
                   value = TRUE, ignore.case = TRUE)

if (length(cols_clima) == 0) {
  stop("No se encontraron columnas bioclimáticas en 'presencias_buffer'.")
}

n_bioclim <- length(cols_clima)

# Identificar filas donde TODAS las bioclimáticas son NA
idx_na_all <- which(rowSums(is.na(presencias_buffer[, cols_clima])) == n_bioclim)

# Filtrar (mantiene las filas con al menos una bioclimática válida)
presencias_buffer <- presencias_buffer[rowSums(is.na(presencias_buffer[, cols_clima])) < n_bioclim, ]

message("??? BUFFER: registros conservados = ", nrow(presencias_buffer),
        if (length(idx_na_all) > 0) paste0(" | eliminados (sin clima) = ", length(idx_na_all)) else "")


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
# para VIF, SI SE QUIERE PROBAR VARIOS ENFOQUES

################################################################################
### PARTE A. REDUCCIÓN DE VARIABLES (VIF + PCA) SEGÚN EL ÁREA DE CALIBRACIÓN ###
################################################################################
cat("\n=== PARTE A: Reducción de variables ambientales por área de calibración ===\n")

# =============================
# 1. Detectar stack ambiental activo
# =============================
if (exists("variables_ecoreg_mask")) {
  RSTACK <- variables_ecoreg_mask
  M_NAME <- "ecoreg"
  message(">> Se usará el stack recortado por ECOREGIONES.")
} else if (exists("variables_buffer_mask")) {
  RSTACK <- variables_buffer_mask
  M_NAME <- "buffer"
  message(">> Se usará el stack recortado por BUFFER 50 km.")
} else {
  stop("?????? No se encontró ninguna variable recortada (ni ecoregiones ni buffer). Verifica antes de continuar.")
}

# =============================
# 2. Parámetros y carpeta de salida
# =============================
OUT_DIR <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/OUT"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SAMPLE_SIZE <- 10000   # número de píxeles a muestrear
NCOMP <- 6             # nº de componentes principales
SAMPLE_PATH <- file.path(OUT_DIR, paste0("sample_", M_NAME, "_", SAMPLE_SIZE, ".csv"))

# =============================
# 3. Funciones auxiliares
# =============================

# --- 3.1. Muestreo con caché (guarda y reutiliza) ---
sample_with_cache <- function(rstack, n, out_csv) {
  if (file.exists(out_csv)) {
    message(">> Cargando muestra cacheada: ", out_csv)
    samp <- read_csv(out_csv, show_col_types = FALSE)
  } else {
    set.seed(123)
    message(">> Generando muestra aleatoria de ", n, " píxeles...")
    samp <- spatSample(rstack, size = n, method = "random", na.rm = TRUE, as.df = TRUE)
    write_csv(samp, out_csv)
    message("??? Muestra guardada: ", out_csv)
  }
  return(samp)
}

# --- 3.2. PCA por bloques (terra::pca) ---
do_pca_blockwise <- function(rstack, ncomp, out_dir, prefix) {
  message(">> Ejecutando PCA (", prefix, ") con ", ncomp, " componentes...")
  pca_res <- terra::pca(rstack, ncomp = ncomp, center = TRUE, scale = TRUE)
  
  # Guardar mapas de componentes
  pca_tif <- file.path(out_dir, paste0(prefix, "_PCA_components_", ncomp, ".tif"))
  writeRaster(pca_res$map, pca_tif, overwrite = TRUE)
  
  # Calcular varianza explicada y guardar
  sdev <- pca_res$model$sdev
  var_exp <- (sdev^2) / sum(sdev^2)
  cum_var <- cumsum(var_exp)
  var_df <- data.frame(
    PC = paste0("PC", seq_along(sdev)),
    SDev = sdev,
    VarExplained = var_exp,
    CumVarExplained = cum_var
  )
  write.csv(var_df, file.path(out_dir, paste0(prefix, "_PCA_variance_explained.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  
  # Guardar loadings
  write.csv(pca_res$model$co, 
            file.path(out_dir, paste0(prefix, "_PCA_loadings.csv")),
            row.names = TRUE, fileEncoding = "UTF-8")
  
  message("??? PCA completado: ", pca_tif)
  return(list(model = pca_res$model, map = pca_res$map, cum_var = cum_var))
}

# =============================
# 4. Generar o cargar muestra
# =============================
samp <- sample_with_cache(RSTACK, n = SAMPLE_SIZE, out_csv = SAMPLE_PATH)

# =============================
# 5. Matriz de correlación
# =============================
cor_mat <- cor(samp, method = "pearson", use = "complete.obs")
cor_path <- file.path(OUT_DIR, paste0("correlation_", M_NAME, "_sample.csv"))
write.csv(cor_mat, cor_path, row.names = TRUE, fileEncoding = "UTF-8")
message("??? Matriz de correlación guardada: ", cor_path)

# =============================
# 6. Selección por VIF (opcional)
# =============================
vars_vif_names <- names(RSTACK)
if (requireNamespace("usdm", quietly = TRUE)) {
  message(">> Calculando VIF (umbral 10)...")
  df_vif <- as.data.frame(stats::na.omit(samp))
  vif_obj <- usdm::vifstep(df_vif, th = 10)
  vars_vif_names <- vif_obj@results$Variables
  message("??? Variables retenidas tras VIF: ", paste(vars_vif_names, collapse = ", "))
} else {
  message("?????? Paquete 'usdm' no disponible, se mantienen todas las variables.")
}

RSTACK_VIF <- RSTACK[[vars_vif_names]]
vif_path <- file.path(OUT_DIR, paste0("variables_", M_NAME, "_reducidas_VIF.tif"))
writeRaster(RSTACK_VIF, vif_path, overwrite = TRUE)
message("??? Stack reducido por VIF guardado: ", vif_path)

# =============================
# 7) PCA COMPATIBLE (sin terra::pca) REUSANDO LA MUESTRA CACHEADA
# =============================

message(">> Calculando PCA (compatibilidad moderna con terra >= 1.7)...")

# --- 7.1 Usar la muestra ya creada para correlación/VIF si existe ---
SAMPLE_SIZE_PCA <- 50000  # usa 10000 si quieres la misma muestra que VIF
if (exists("samp")) {
  samp_pca <- samp
  message(">> Usando la muestra 'samp' ya cargada (n = ", nrow(samp_pca), ").")
} else if (exists("SAMPLE_PATH") && file.exists(SAMPLE_PATH)) {
  samp_pca <- read_csv(SAMPLE_PATH, show_col_types = FALSE)
  message(">> Cargando muestra cacheada desde: ", SAMPLE_PATH, " (n = ", nrow(samp_pca), ").")
} else {
  set.seed(123)
  samp_pca <- spatSample(RSTACK_VIF, size = SAMPLE_SIZE_PCA, method = "random", na.rm = TRUE, as.df = TRUE)
  message(">> Muestra creada para PCA (n = ", nrow(samp_pca), ").")
}

# --- 7.2 Alinear nombres/orden de variables entre muestra y raster ---
#    (Evita desalineaciones por orden de capas)
if (!all(names(RSTACK_VIF) %in% colnames(samp_pca))) {
  stop("??? Hay capas en RSTACK_VIF que no aparecen como columnas en la muestra. Revisa nombres.")
}
# Reordenar columnas del sample al orden exacto del raster:
samp_pca <- samp_pca[, names(RSTACK_VIF), drop = FALSE]

# --- 7.3 PCA con prcomp() sobre la muestra (centrado y escalado) ---
pca_model <- prcomp(na.omit(samp_pca), center = TRUE, scale. = TRUE)

# --- 7.4 Extraer loadings y preparar proyección raster por bloques ---
loadings     <- pca_model$rotation[, 1:NCOMP, drop = FALSE]  # variables x PCs
center_vals  <- pca_model$center[names(RSTACK_VIF)]
scale_vals   <- pca_model$scale[names(RSTACK_VIF)]

# Chequeos de dimensiones
message(">> Verificando dimensiones: nlyr(RSTACK_VIF)=", nlyr(RSTACK_VIF),
        " | nrow(loadings)=", nrow(loadings))
if (nlyr(RSTACK_VIF) != nrow(loadings)) {
  stop("??? Dimensiones no coinciden entre stack y loadings. Revisa orden/nombres.")
}

# Función segura para aplicar PCA por bloques
pca_map_fun <- function(v) {
  # v: matriz [celdas_del_bloque x n_vars] en el mismo orden que names(RSTACK_VIF)
  v_scaled <- sweep(v, 2, center_vals, "-")
  v_scaled <- sweep(v_scaled, 2, scale_vals,  "/")
  res <- as.matrix(v_scaled) %*% loadings    # [celdas x NCOMP]
  res
}

# --- 7.5 Ejecutar terra::app() en modo seguro y exportar mapas PC ---
pc_file  <- file.path(OUT_DIR, paste0(M_NAME, "_PCA_components_", NCOMP, ".tif"))
pc_stack <- try(
  terra::app(RSTACK_VIF, fun = pca_map_fun, cores = 1, filename = pc_file, overwrite = TRUE),
  silent = TRUE
)

if (inherits(pc_stack, "try-error")) {
  message("??? Error en terra::app():"); print(pc_stack)
  stop("No se pudo generar el raster de componentes. Revisa mensajes previos.")
} else {
  names(pc_stack) <- paste0("PC", 1:NCOMP)
  message("??? Mapas de componentes guardados: ", pc_file)
}

# --- 7.6 Guardar varianza explicada y cargas (para interpretar PCs) ---
sdev    <- pca_model$sdev
var_exp <- (sdev^2) / sum(sdev^2)
cum_var <- cumsum(var_exp)

var_df <- data.frame(
  PC = paste0("PC", seq_along(sdev)),
  SDev = sdev,
  VarExplained = var_exp,
  CumVarExplained = cum_var
)

write.csv(var_df, file.path(OUT_DIR, paste0(M_NAME, "_PCA_variance_explained.csv")),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(loadings, file.path(OUT_DIR, paste0(M_NAME, "_PCA_loadings.csv")),
          row.names = TRUE, fileEncoding = "UTF-8")

message(sprintf("??? %s: %.1f%% de varianza acumulada hasta PC%d.",
                toupper(M_NAME), 100 * cum_var[NCOMP], NCOMP))


# =============================
# 8) (Opcional) Guardar PCs hasta un umbral de varianza (p. ej., 90 %)
# =============================
UMBRAL_VAR <- 0.90
k_auto <- which(cum_var >= UMBRAL_VAR)[1]

if (!is.na(k_auto) && k_auto >= 1) {
  auto_path <- file.path(OUT_DIR, sprintf("%s_PCA_components_auto_%dPCs.tif", M_NAME, k_auto))
  writeRaster(pc_stack[[1:k_auto]], auto_path, overwrite = TRUE)
  message(sprintf("??? %s: PCs hasta %.0f%% (k = %d) guardadas en: %s",
                  toupper(M_NAME), 100 * UMBRAL_VAR, k_auto, auto_path))
} else {
  message("?????? No se alcanzó el umbral de varianza especificado; no se exportaron PCs 'auto'.")
}

# =========================
# 9) Sugerencias para MaxEnt
# =========================
message("\nSugerencias:")
message("- Si quieres interpretar variables originales: usa 'variables_reducidas_VIF.tif'.")
message("- Si prefieres componentes (decorrelados): usa '", basename(pca_auto_tif), "' o '", basename(pca_tif), "'.")
message("- Revisa 'PCA_variance_explained.csv' y 'PCA_loadings.csv' para elegir PCs y entender sus ejes.")

################################################################################
### BLOQUE 8bis. ANÁLISIS VISUAL Y ESTADÍSTICO DEL PCA #########################
################################################################################
library(ggplot2)
library(reshape2)
library(viridis)

cat("\n=== BLOQUE 8bis: Visualización del PCA y soporte estadístico ===\n")

# --- 1. Gráfico Scree Plot (Varianza explicada por componente) ---
pca_var_df <- data.frame(
  PC = paste0("PC", seq_along(sdev)),
  VarExplained = var_exp,
  CumVar = cum_var
)

# Crear el gráfico
p_scree <- ggplot(pca_var_df, aes(x = PC, y = VarExplained)) +
  geom_bar(stat = "identity", fill = "#3182bd") +
  geom_line(aes(y = CumVar), group = 1, color = "#e6550d", linewidth = 1) +
  geom_point(aes(y = CumVar), color = "#e6550d", size = 2) +
  theme_minimal(base_size = 12) +
  labs(title = paste("Varianza explicada - PCA (", toupper(M_NAME), ")", sep = ""),
       y = "Proporción de varianza", x = "Componente principal")

# Mostrar en pantalla (opcional)
print(p_scree)

# Guardar el gráfico
ggsave(
  filename = file.path(OUT_DIR, paste0(M_NAME, "_PCA_screeplot.png")),
  plot = p_scree,
  width = 8, height = 5, dpi = 300
)

# --- 2. Heatmap de cargas (qué variables dominan cada PC) ---
load_df <- as.data.frame(loadings[, 1:NCOMP, drop = FALSE])
load_df$Variable <- rownames(load_df)
load_melt <- reshape2::melt(load_df, id.vars = "Variable", variable.name = "PC", value.name = "Loading")

# Opcional: usar magnitud (abs) si prefieres resaltar contribución
# load_melt$Loading <- abs(load_melt$Loading)

# Ordenar ejes para que salgan en el orden correcto
load_melt$PC <- factor(load_melt$PC, levels = colnames(loadings))
load_melt$Variable <- factor(load_melt$Variable, levels = rev(sort(unique(load_melt$Variable))))

# 2) Crear el gráfico
p_heat <- ggplot(load_melt, aes(x = PC, y = Variable, fill = Loading)) +
  geom_tile() +
  scale_fill_viridis(option = "A", name = "Carga") +
  theme_minimal(base_size = 10) +
  labs(
    title = paste("Cargas de variables en PCs (", toupper(M_NAME), ")", sep = ""),
    x = "Componente principal", y = "Variable ambiental"
  )

# 3) Mostrar (opcional)
print(p_heat)

# 4) Guardar
ggsave(
  filename = file.path(OUT_DIR, paste0(M_NAME, "_PCA_loadings_heatmap.png")),
  plot = p_heat, width = 7, height = 6, dpi = 300
)

# --- 3. Biplot (relación entre variables y componentes) ---
biplot_path <- file.path(OUT_DIR, paste0(M_NAME, "_PCA_biplot.png"))
png(biplot_path, width = 800, height = 600)
biplot(pca_model, main = paste("Biplot PCA -", toupper(M_NAME)))
dev.off()
cat("??? Gráficos del PCA guardados en:", OUT_DIR, "\n")

################################################################################
### BLOQUE FINAL. EXPORTAR INSUMOS PARA WALLACE + RESUMEN VISUAL ###############
################################################################################

BASE_WALLACE_DIR <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/WALLACE_INPUTS"
dir.create(file.path(BASE_WALLACE_DIR, "presence"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BASE_WALLACE_DIR, "variables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BASE_WALLACE_DIR, "area"),      recursive = TRUE, showWarnings = FALSE)

# ---- Área que deseas exportar (ajusta si trabajas con buffer) ----
M_NAME <- "ecoreg"   # o "buffer"
#M_NAME <- "buffer"

# === 1) PRESENCIAS ===
pres_df <- if (M_NAME == "buffer" && exists("presencias_buffer")) presencias_buffer else presencias_ecoreg
if (!exists("pres_df") || is.null(pres_df)) pres_df <- sp_df  # respaldo mínimo

stopifnot(all(c("scientificName","decimalLongitude","decimalLatitude") %in% names(pres_df)))

pres_export <- pres_df %>%
  st_drop_geometry() %>%
  select(scientificName, decimalLongitude, decimalLatitude) %>%
  distinct() %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude),
         decimalLongitude >= -180, decimalLongitude <= 180,
         decimalLatitude  >= -90,  decimalLatitude  <= 90)

pres_csv <- file.path(BASE_WALLACE_DIR, "presence", paste0("presencias_", M_NAME, ".csv"))
write.csv(pres_export, pres_csv, row.names = FALSE)
message("??? Presencias ??? ", pres_csv)

# === 2) VARIABLES (elige mejor stack disponible, sin bio3/bio7) ===
OUT_DIR <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/SUELO/OUT"
cand_auto <- list.files(OUT_DIR, pattern = paste0("^", M_NAME, "_PCA_components_auto_\\d+PCs\\.tif$"), full.names = TRUE)
cand_full <- list.files(OUT_DIR, pattern = paste0("^", M_NAME, "_PCA_components_\\d+\\.tif$"),       full.names = TRUE)
cand_vif  <- list.files(OUT_DIR, pattern = paste0("^variables_", M_NAME, "_reducidas_VIF\\.tif$"),   full.names = TRUE)

# 1) Elegir fuente priorizando PCA auto > PCA full > VIF
vars_src <- if (length(cand_auto)) cand_auto[which.max(file.mtime(cand_auto))] else
  if (length(cand_full)) cand_full[which.max(file.mtime(cand_full))] else
    if (length(cand_vif))  cand_vif[which.max(file.mtime(cand_vif))]  else NA

if (is.na(vars_src)) stop("No encontré TIFs de variables en OUT_DIR para '", M_NAME, "'. ¿PCA/VIF ya corrió?")

# 2) Cargar raster y descartar bio3/bio7 si existen
preds <- terra::rast(vars_src)
drop_vars <- tolower(c("bio3","bio07","wc2.1_30s_bio_3","wc2.1_30s_bio_7"))
keep_idx <- !tolower(names(preds)) %in% drop_vars
preds <- preds[[keep_idx]]
message("??? Variables mantenidas: ", paste(names(preds), collapse = ", "))

# 3) Guardar multibanda filtrado en WALLACE_INPUTS/variables/
vars_tgt <- file.path(BASE_WALLACE_DIR, "variables", paste0("predictors_", M_NAME, ".tif"))
terra::writeRaster(preds, vars_tgt, overwrite = TRUE)
message("??? Variables (multibanda sin bio3/bio7) ??? ", vars_tgt)

# 4) Exportar también "por capas" (un .tif por variable/PC)
nm <- gsub("[^A-Za-z0-9_]", "_", names(preds))
nm <- make.unique(nm)
names(preds) <- nm

vars_dir_layers <- file.path(BASE_WALLACE_DIR, "variables", paste0("layers_", M_NAME))
dir.create(vars_dir_layers, showWarnings = FALSE, recursive = TRUE)

for (i in 1:terra::nlyr(preds)) {
  terra::writeRaster(
    preds[[i]],
    filename = file.path(vars_dir_layers, paste0(names(preds)[i], ".tif")),
    overwrite = TRUE
  )
}
message("??? Variables por capa (sin bio3/bio7) ??? ", vars_dir_layers)

# === 3) ÁREA (shapefile) ===
area_src <- if (M_NAME == "buffer") {
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_buffer50km.shp"
} else {
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_ecoregiones.shp"
}
stopifnot(file.exists(area_src))
area_sf <- st_read(area_src, quiet = TRUE) |> st_make_valid()
area_wgs <- st_transform(area_sf, 4326)

area_tgt <- file.path(BASE_WALLACE_DIR, "area", paste0("area_", M_NAME, ".shp"))
# limpia shapefiles previos con mismo prefijo
invisible(file.remove(list.files(dirname(area_tgt), pattern = paste0("^area_", M_NAME, "\\."), full.names = TRUE)))
st_write(area_wgs, area_tgt, delete_layer = TRUE, quiet = TRUE)
message("??? Área      ??? ", area_tgt)

# === 4) RESUMEN VISUAL (robusto) ===
# cierra dispositivos gráficos "colgados"
while (dev.cur() > 1) dev.off()

predictors_r <- rast(vars_tgt)
presence_r   <- read.csv(pres_csv)
area_r       <- st_read(area_tgt, quiet = TRUE)

# reproyecta área al CRS del raster si hace falta
if (!identical(st_crs(area_r)$wkt, as.character(crs(predictors_r)))) {
  area_r <- st_transform(area_r, crs(predictors_r))
}

# gráfico rápido
plot(predictors_r[[1]], main = paste("Resumen de insumos Wallace -", toupper(M_NAME)))
plot(st_geometry(area_r), add = TRUE, border = "black", lwd = 1.2)
points(presence_r$decimalLongitude, presence_r$decimalLatitude, col = "red", pch = 20, cex = 0.6)
legend("bottomleft", legend = c("Presencias", "Área M"), col = c("red", "black"),
       pch = c(20, NA), lty = c(NA, 1), bty = "n")

cat("\n??? Insumos listos para Wallace en:", BASE_WALLACE_DIR, "\n")

wallace_vars <- rast("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/WALLACE_INPUTS/variables/predictors_ecoreg.tif")

# Ver información general
wallace_vars
names(wallace_vars)

# Visualizar el primer componente (por ejemplo PC1)
plot(wallace_vars[[1]], main = names(wallace_vars)[1])

##########Solo para pruebas, luego se puede borrar###############################
##Validando cantidad de variables en el raste de variables ambientales###
library(terra)
r <- rast("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/WALLACE_INPUTS/variables/predictors_ecoreg.tif")
nlyr(r)
names(r)

