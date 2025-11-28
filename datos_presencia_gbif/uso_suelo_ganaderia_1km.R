#getwd()

# Definir la ruta base 
oficina_path <- "D:/Usuario/Desktop/ROBERTO_CEBALLOS/MaestrÃƒÂ­a/SEMESTRE III/PROYECTO INTEGRADOR/"
robert_path <- "C:/Users/rober/Documents/MAESTRIA/PROYECTO INTEGRADOR/REPOSITORIO"
mayra_path <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/"
base_path <- robert_path;

setwd(base_path)

library(terra)
library(sf)
library(rnaturalearth)
library(dplyr)
library(ggplot2)

# --- RUTAS ---
suelo_path <- file.path(base_path, "SUELO", "USO_SUELO")

gpkg_file <- file.path(
  suelo_path,
  "gpw_livestock.animals_gpw.fao.glw3_polygon.samples_20000101_20231231_go_epsg.4326_v1.gpkg"
)

rules_p_2022 <- file.path(
  suelo_path,
  "gpw_livestock.pot.land_grassland.cropland.rules_p_1km_20220101_20221231_go_epsg.4326_v1.tif"
)

# ------------------ 1. Cargar raster ----------------------

uso_pot <- rast(rules_p_2022)

# Importante: reproyectar el raster FAO (Homolosine ??? WGS84)
# Esto evita distorsión y errores de recorte
uso_pot_4326 <- project(uso_pot, "EPSG:4326")

# ------------------ 2. Cargar Colombia --------------------

col_sf <- ne_countries(country = "colombia", returnclass = "sf")
col_v  <- vect(col_sf)   # a terra

# Recortar raster ya en WGS84 a Colombia
uso_pot_col <- crop(uso_pot_4326, col_v)
uso_pot_col <- mask(uso_pot_col, col_v)

# ------------------ 3. Cargar polígonos ganaderos ---------

sf_use_s2(FALSE)

polys <- st_read(gpkg_file, quiet = TRUE) %>%
  st_make_valid()

# Recortar a Colombia (WGS84 ??? WGS84)
polys_col <- st_intersection(polys, col_sf)

# Filtrar el campo cattle_2022
cattle_2022 <- polys_col %>%
  select(gazName, country, cattle_2022, geom) %>%
  filter(!is.na(cattle_2022))

# Reproyectar al mismo CRS del raster final (EPSG:4326)
cattle_2022_4326 <- st_transform(cattle_2022, 4326)

# Convertir a SpatVector de terra
cattle_v <- vect(cattle_2022_4326)

# ------------------ 4. Extraer idoneidad ------------------

vals <- extract(uso_pot_col, cattle_v)

# Agregar valores
cattle_v$idoneidad <- vals[,2]
cattle_v$idoneidad[is.na(cattle_v$idoneidad)] <- 0

# ------------------ 5. Calcular valor final ----------------

cattle_v$valor_final <- cattle_v$cattle_2022 * (cattle_v$idoneidad / 100)

# ------------------ 6. Rasterización final -----------------

raster_final <- rasterize(
  cattle_v,
  uso_pot_col,
  field = "valor_final",
  fun   = "sum",
  background = 0
)

# ------------------ 7. Visualización ----------------------

plot(raster_final, main = "Uso de suelo para ganadería (1 km, WGS84)")
plot(col_v, add = TRUE)

# ------------------ 8. Exportar resultado -----------------

writeRaster(raster_final,
            file.path(suelo_path, "raster_ganaderia_final_1km_WGS84.tif"),
            overwrite = TRUE)

nrow(cattle_2022_4326)
plot(cattle_v)
