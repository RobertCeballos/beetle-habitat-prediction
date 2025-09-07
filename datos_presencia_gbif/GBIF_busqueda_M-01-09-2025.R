setwd("F:/Doctorado/asesorias")

# librerias
install.packages("reticulate") # si no lo tienes instalado
library(reticulate)
# 1. Instala miniconda en tu usuario (solo la primera vez)
install_miniconda()

# 2. Crea el entorno r-reticulate
conda_create("r-reticulate")

# 3. Instala la librería de Google Earth Engine en ese entorno
conda_install("r-reticulate", "earthengine-api", pip = TRUE)

# 4. Activa el entorno
use_condaenv("r-reticulate", required = TRUE)
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

################# obtención de registros biológicos ############################

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
#saveRDS(all_synonyms, "all_synonyms.rds") # guardar los sinónimos en RDS
all_synonyms <- readRDS("all_synonyms.rds") # cargar los datos para el siguiente procesamientos
# union de listas de especies (Giraldo et al. 2018) y sinónimos (GBIF)
synonym_list <- unique(sort(all_synonyms$canonicalName))
all_taxa <- sort(unique(c(lista_sp, synonym_list)))
#records <- bind_rows(data_list, .id = "scientificName")

# busqueda de registros geográficos de cada especie (hasta 5000 registros)
# con coordenadas geográficas en GBIF

result <- occ_data(scientificName = all_taxa, limit = 5000, hasCoordinate = TRUE)

# procesamiento del objeto gbif_data para obtener un data frame
data_list <- lapply(result, function(x) if (!is.null(x$data)) x$data else NULL)
data_list <- Filter(Negate(is.null), data_list) 
records <- bind_rows(data_list, .id = "scientificName")
#saveRDS(records, "GBIF_2025-03-23.rds")# activar para guardar el archivo en RDS
records <- readRDS("GBIF_2025-03-23.rds") # cargar el archivo para el siguiente proceso

# Actualización del nombre de las especies según GBIF Backbone Taxonomy
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
#ecoreg <- st_read("F:/Capas/World/wwf/Ecoregions2017/Ecoregions2017.shp") # Capa de ecoregiones WWF 2017 
ecoreg <- st_read("C:/Users/usuario/Documents/Posgrado/Proyecto integrador/dungbeetles/Ecoregions2017/Ecoregions2017.shp") # Capa de ecoregiones WWF 2017 

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
             sp$stateProvince %in% c("Boyacá",  "Meta","Casanare", "Vichada", "BoyacÃ¡")), ] # 1824 records
sp <- sp[!(sp$decimalLongitude == -76.09989 & sp$decimalLatitude == 8.039917), ] # remove Urabá recrord, 1823 records

# limpieza de duplicados 745 records
sp$coordinates <- paste(sp$decimalLatitude, sp$decimalLongitude, sep="_")
sp <- sp[!is.na(sp$coordinates) & !duplicated(sp[c("scientificName", "coordinates")]), ]

# 
library(rgee)

# Inicializar GEE
reticulate::use_condaenv("r-reticulate", required = TRUE) # enlace al entorno especifico de mi pc a conda de python

ee_Initialize()

# Registros a formato GEE
sp_clean <- sp[, c("decimalLongitude", "decimalLatitude", "scientificName")]
sp_sf <- st_as_sf(sp_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
sp_ee <- sf_as_ee(sp_sf)

# Elevación  ALOS PALSAR 30m de GEE
DEM <- ee$Image("JAXA/ALOS/AW3D30/V2_2")$select("AVE_DSM")

# Extraer elevación
elev_ee <- DEM$reduceRegions(
  collection = sp_ee,
  reducer = ee$Reducer$first(),
  scale = 30
)

#  Tranformacion de datos de ee a R
elev_sf <- ee_as_sf(elev_ee)
#saveRDS(elev_sf, "./elev_ALOS.rds")
elev_sf <- read("./elev_ALOS.rds")

sp$elev_ALOS <- elev_sf$first

# limpieza de registros por proximidad de 1km (resolución de las variables predictoras de WorldClim)
# 107 registros

# Crea un sf con registros y marco de coordenadas AEA
sp_sf <- st_as_sf(sp, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
sp_aea <- st_transform(sp_sf, crs = AEAstring)

# Itera sobre cada punto para obtener la distancia entre si  
keep <- rep(TRUE, nrow(sp_aea))
for (i in seq_len(nrow(sp_aea))) {
  if (!keep[i]) next
  
  # Calcula distancias desde el punto i a todos los siguientes
  dists <- st_distance(sp_aea[i, ], sp_aea[(i+1):nrow(sp_aea), ])
  
  # Marca como FALSE a los puntos posteriores que estén a < 1000 m
  close_points <- which(as.numeric(dists) < 1000)
  if (length(close_points) > 0) {
    keep[(i + close_points)] <- FALSE
  }
}
# Filtrar y volver a data.frame
sp_aea_filtered <- sp_aea[keep, ]
sp_aea_filtered <- st_transform(sp_aea_filtered, crs = 4326)
coords <- st_coordinates(sp_aea_filtered)
sp_aea_filtered$decimalLongitude <- coords[, "X"]
sp_aea_filtered$decimalLatitude <- coords[, "Y"]
sp_df <- as.data.frame(sp_aea_filtered)

# filtro por elevación, 105 registros
sp_df <- sp_df %>%
  filter(elev_ALOS >= 0, elev_ALOS <= 1600)

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

################################################################################



