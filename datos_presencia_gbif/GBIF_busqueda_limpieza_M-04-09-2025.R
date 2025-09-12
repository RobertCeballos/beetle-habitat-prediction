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
plot(col_dem, main = "Modelo de ElevaciÃ³n Digital - Colombia")
col_dem <- rast(col_dem)

# VisualizaciÃ³n de los registros  
ggplot() +
  geom_spatraster(data = col_dem) +
  scale_fill_viridis_c(name = "Elevación") +
  geom_sf(data = col, fill = "lightgrey", color = "black", alpha = 0.1) +
  geom_point(data = sp, aes(x = decimalLongitude, y = decimalLatitude), color = "darkblue")

## limpieza de puntos atÃ­picos basados en literatura (Montoya et al. 2021)
sp <- sp[!(sp$scientificName=="Dichotomius agenor" & sp$decimalLongitude <= -81),] # remove records greater than  81 Long D agenor
sp <- sp[!(sp$scientificName=="Dichotomius agenor" & # Remove D. agenor departments based on Montoya et al. 2021
             sp$stateProvince %in% c("Boyacá",  "Meta","Casanare", "Vichada", "Boyacá")), ] # 1824 records
sp <- sp[!(sp$decimalLongitude == -76.09989 & sp$decimalLatitude == 8.039917), ] # remove UrabÃ¡ recrord, 1823 records

# limpieza de duplicados 745 records
sp$coordinates <- paste(sp$decimalLatitude, sp$decimalLongitude, sep="_")
sp <- sp[!is.na(sp$coordinates) & !duplicated(sp[c("scientificName", "coordinates")]), ]

# 
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

# Leer el archivo de elevaciÃ³n
elev_sf <- readRDS("./elev_ALOS.rds")

# Revisar nombres de columnas
print(names(elev_sf))  # Debe incluir "elev_ALOS"

# Convertir sp a sf
library(sf)
library(dplyr)

sp_sf <- st_as_sf(sp, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Asegurar que elev_sf estÃ¡ en el mismo CRS
elev_sf <- st_transform(elev_sf, crs = 4326)

# Unir elevaciones al dataset original por proximidad espacial
sp_joined <- st_join(sp_sf, elev_sf["elev_ALOS"], join = st_nearest_feature)

# Convertir a AEA
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
sp_aea <- st_transform(sp_joined, crs = AEAstring)

# Filtrar puntos por distancia mÃ­nima de 1 km
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
  stop("elev_ALOS se perdiÃ³ en la conversiÃ³n a data.frame. Revisa antes de filtrar.")
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
# extracciÃ³n del "area de movilidad (M)" basado en regiones biogeogrÃ¡ficas (WWWF ecoregions)
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

################################################################################
#########6. Enmascarar capas ambientales según área de calibración.#############


# --- 6.1. Guardar las áreas de calibración (M) ---
st_write(eco_col, 
         "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_ecoregiones.shp", 
         delete_layer = TRUE)

st_write(M_buffer50k_union, 
         "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_buffer50km.shp", 
         delete_layer = TRUE)

# --- 6.2. Cargar variables bioclimáticas de WorldClim 2.1 ---
variables <- terra::rast(list.files("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/WORLDCLIM",
                                    pattern = "tif$", full.names = TRUE))

# Para ecoregiones
M_ecoreg_vect <- terra::project(M_ecoreg_vect, crs(variables))

# Para buffer
M_buffer_vect <- terra::project(M_buffer_vect, crs(variables))

# --- 6.3. Recorte y enmascarado usando M por ecorregiones ---
M_ecoreg <- st_read("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_ecoregiones.shp")
M_ecoreg_vect <- vect(M_ecoreg) # convertir a SpatVector

variables_ecoreg_crop <- terra::crop(variables, M_ecoreg_vect)
variables_ecoreg_mask <- terra::mask(variables_ecoreg_crop, M_ecoreg_vect)

# --- 6.3. Recorte y enmascarado usando M por buffer ---
M_buffer <- st_read("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/AREAS_CALIBRACION/M_buffer50km.shp")
M_buffer_vect <- vect(M_buffer)

variables_buffer_crop <- terra::crop(variables, M_buffer_vect)
variables_buffer_mask <- terra::mask(variables_buffer_crop, M_buffer_vect)

crs(variables)
crs(M_ecoreg_vect)

# --- 6.4. Guardar resultados en carpetas separadas ---
dir.create("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/variables_M_ecoregiones", showWarnings = FALSE, recursive = TRUE)
dir.create("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/variables_M_buffer", showWarnings = FALSE, recursive = TRUE)

terra::writeRaster(variables_ecoreg_mask,
                   filename = paste0("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/variables_M_ecoregiones/", names(variables_ecoreg_mask), ".tif"),
                   overwrite = TRUE)

terra::writeRaster(variables_buffer_mask,filename = paste0("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/variables_M_buffer/", names(variables_buffer_mask), ".tif"),overwrite = TRUE)

# --- 6.5. (Opcional) Visualizar una variable para comprobar ---
plot(variables_ecoreg_mask[[1]], main = "BIO1 - Recortado por ecorregiones")
plot(M_ecoreg_vect, add = TRUE)

plot(variables_buffer_mask[[1]], main = "BIO1 - Recortado por buffer")
plot(M_buffer_vect, add = TRUE)

################################################################################

# --- 1. Cargar los datos filtrados de presencia ---
# (Asegúrate que sp_df contiene decimalLongitude y decimalLatitude ya filtrados)
# y que 'variables' es tu SpatRaster con BIO1 a BIO19 ya recortado a las áreas de calibración.

# --- 2. Convertir a SpatVector en CRS de variables ---
sp_vect <- vect(sp_df, geom = c("decimalLongitude", "decimalLatitude"), crs = crs(variables))


# Extraer valores para cada punto
puntos_ecoreg <- terra::extract(variables_ecoreg_mask, sp_vect)

# Combinar con el dataframe original de presencias
presencias_ecoreg <- cbind(sp_df, puntos_ecoreg)

# Guardar CSV
write.csv(
  presencias_ecoreg,
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/presencias_variables_ecoregiones.csv",
  row.names = FALSE
)


# Extraer valores para cada punto
puntos_buffer <- terra::extract(variables_buffer_mask, sp_vect)

# Combinar con el dataframe original de presencias
presencias_buffer <- cbind(sp_df, puntos_buffer)

# Guardar CSV
write.csv(
  presencias_buffer,
  "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/beetle-habitat-prediction/datos_presencia_gbif/presencias_variables_buffer50km.csv",
  row.names = FALSE
)

# --- 5. Verificación rápida ---
cat("CSV generados en la carpeta 'datos_presencia_gbif'\n")
cat("Filas ecoregiones:", nrow(presencias_ecoreg), "\n")
cat("Filas buffer 50 km:", nrow(presencias_buffer), "\n")