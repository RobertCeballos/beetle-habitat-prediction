setwd("F:/Doctorado/asesorias")

# librerias

library(readxl)
library(rgbif)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
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
records <- bind_rows(data_list, .id = "scientificName")

# busqueda de registros geográficos de cada especie (hasta 5000 registros)
# con coordenadas geográficas en GBIF

result <- occ_data(scientificName = all_taxa, limit = 5000, hasCoordinate = TRUE)

# procesamiento del objeto gbif_data para obtener un data frame
data_list <- lapply(result, function(x) if (!is.null(x$data)) x$data else NULL)
data_list <- Filter(Negate(is.null), data_list) 
records <- bind_rows(data_list, .id = "scientificName")
# saveRDS(records, "GBIF_2025-03-23.rds")# activar para guardar el archivo en RDS
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
# Revisar registros #
#
# subset records de Dichotomius agenor
sp <- data.frame(subset(records, scientificName == "Dichotomius agenor"))
#
# mapas base
ecoreg <- st_read("F:/Capas/World/wwf/Ecoregions2017/Ecoregions2017.shp") # Capa de ecoregiones WWF 2017 
ecoreg <- st_make_valid(ecoreg)
col <- ne_states(country = "colombia", returnclass = "sf")
# gráfica 
ggplot() +
  geom_sf(data = col, fill = "lightgrey", color = "black") +
  geom_point(data = sp, aes(x = decimalLongitude, y = decimalLatitude), color = "darkblue")

## limpieza de puntos atípicos
sp <- sp[!(sp$scientificName=="Dichotomius agenor" & sp$decimalLongitude <= -81),] # remove records greater than  81 Long D agenor
sp <- sp[!(sp$scientificName=="Dichotomius agenor" & sp$stateProvince %in% c("Boyacá",  "Meta","Casanare", "Vichada", "BoyacÃ¡")), ] # Remove D. agenor departments based on Montoya et al. 2021

###
# extracción del "area de movilidad (M)" basado en regiones biogeográficas (WWWF ecoregions)
###
sp_sf <- st_as_sf(sp, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(ecoreg)) 
eco_col <- st_intersects(ecoreg, sp_sf)
eco_col <- ecoreg[lengths(eco_col) > 0, ]

ggplot() +
  geom_sf(data = col, fill = "lightgrey", color = "black") +
  geom_sf(data = eco_col, fill = "darkred", color = "black", alpha=0.7) +
  geom_point(data = sp, aes(x = decimalLongitude, y = decimalLatitude), color = "darkblue", size=3)

###
# Area de movilidad basado en buffer 50km 
###
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
sp_sf_AEA <- st_transform(sp_sf, crs = AEAstring)

M_buffer50k <- st_buffer(sp_sf_AEA, dist = 50000)

ggplot() +
  geom_sf(data = col, fill = "lightgrey", color = "black") +
  geom_sf(data = M_buffer50k, fill = "darkred", color = "black", alpha=0.7) +
  geom_point(data = sp, aes(x = decimalLongitude, y = decimalLatitude), color = "darkblue", size=3)

M_buffer50k_union <- st_union(M_buffer50k)

ggplot() +
  geom_sf(data = col, fill = "lightgrey", color = "black") +
  geom_sf(data = M_buffer50k_union, fill = "darkred", color = "black", alpha=0.7) +
  geom_point(data = sp, aes(x = decimalLongitude, y = decimalLatitude), color = "darkblue", size=3)

################################################################################

# conteo de registros geográficos por especie (Inlcluye coordenadas duplicadas)
canonical_sp_count <- as.data.frame(table(sort(records$species)))
species_count <- as.data.frame(table(sort(records$scientificName)))
names(species_count) <- c("species", "count")

# conteo de registros únicos

records$coordinates <- paste(records$decimalLatitude, records$decimalLongitude, sep="_")
records <- records[!is.na(records$coordinates) & !duplicated(records[c("species", "coordinates")]), ]
unicos <- as.data.frame(table(sort(records$scientificName)))
names(unicos) <- c("species", "count")


ggplot(unicos, aes(y = reorder(species, count), x = count)) +
  geom_col(fill = alpha("darkgreen", 0.6), color = "black") +
  scale_x_continuous(expand = c(0, 0.05)) +
  theme_classic() +
  labs(title = "Número de registros disponibles en GBIF",
       x = "Conteo", y = NULL)

unicos <- unicos[!(unicos$species=="Digitonthophagus gazella"), ]

