# -----------------------------------------------------------
#  SCRIPT COMPLETO: IDONEIDAD CLIMÁTICA × USO POTENCIAL GANADERÍA
#  Incluye mapas profesionales y clasificación
#  Autor: May + Tío Chat
# -----------------------------------------------------------

#getwd()

# Definir la ruta base 
oficina_path <- "D:/Usuario/Desktop/ROBERTO_CEBALLOS/MaestrÃƒÂ­a/SEMESTRE III/PROYECTO INTEGRADOR/"
robert_path <- "C:/Users/rober/Documents/MAESTRIA/PROYECTO INTEGRADOR/REPOSITORIO"
mayra_path <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/"
base_path <- robert_path;

setwd(base_path)

library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(rnaturalearth)

# ------------------ RUTAS -------------------------

idoneidad_file <- file.path(base_path, "beetle-habitat-prediction","datos_presencia_gbif","IDONEIDAD_CLIMATICA", "Dichotomius_agenor_fc.LQ_rm.2_logistic.tif")
uso_pot_file   <- file.path(base_path, "SUELO", "USO_SUELO",
                            "gpw_livestock.pot.land_grassland.cropland.rules_p_1km_20220101_20221231_go_epsg.4326_v1.tif")

# Archivo final
cruce_output <- file.path(base_path, "beetle-habitat-prediction","datos_presencia_gbif","IDONEIDAD_CLIMATICA", "raster_cruce_escarabajo_ganaderia.tif")

# ------------------ 1. Cargar rasters ---------------------
id_clima <- rast(idoneidad_file)
uso_pot  <- rast(uso_pot_file)

# ------------------ 2. Reproyección si es necesario --------
if (crs(id_clima) != crs(uso_pot)) {
  uso_pot <- project(uso_pot, crs(id_clima))
}

# ------------------ 3. Cargar Colombia ---------------------
col_sf <- ne_countries(country = "colombia", returnclass = "sf")
col_v <- vect(col_sf)

# ------------------ 4. Recortar a Colombia -----------------
id_clima_col <- crop(id_clima, col_v) |> mask(col_v)
uso_pot_col  <- crop(uso_pot,  col_v) |> mask(col_v)

# ------------------ 5. Alinear resolución ------------------
uso_pot_resampled <- resample(uso_pot_col, id_clima_col, method="bilinear")

# ------------------ 6. Cruce final -------------------------
cruce_final <- id_clima_col * uso_pot_resampled

writeRaster(cruce_final, cruce_output, overwrite = TRUE)

# ------------------ 7. Para graficar con ggplot ------------
cruce_df <- as.data.frame(cruce_final, xy = TRUE)
names(cruce_df)[3] <- "valor"

# ------------------ 8. Mapa principal ----------------------
ggplot() +
  geom_raster(data = cruce_df, aes(x = x, y = y, fill = valor)) +
  geom_sf(data = col_sf, fill = NA, color = "white", size = 0.4) +
  scale_fill_viridis(
    name = "Idoneidad ×\nUso ganadero",
    option = "magma",
    na.value = "transparent"
  ) +
  labs(
    title = "Cruce: Idoneidad climática del escarabajo × Uso potencial ganadería",
    subtitle = "Resolución: 1 km - Colombia",
    x = "Longitud",
    y = "Latitud"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ------------------ 9. Clasificación por cuantiles ----------
q75 <- quantile(cruce_df$valor, 0.75, na.rm = TRUE)
q90 <- quantile(cruce_df$valor, 0.90, na.rm = TRUE)

cruce_df$zona <- cut(
  cruce_df$valor,
  breaks = c(-Inf, q75, q90, Inf),
  labels = c("Baja", "Media", "Alta"),
  include.lowest = TRUE
)

# ------------------ 10. Mapa categórico --------------------
ggplot() +
  geom_raster(data = cruce_df, aes(x = x, y = y, fill = zona)) +
  geom_sf(data = col_sf, fill = NA, color = "white", size = 0.4) +
  scale_fill_manual(
    name = "Nivel de oportunidad\n(Esc × Ganadería)",
    values = c("Baja" = "#253494",
               "Media" = "#2c7fb8",
               "Alta" = "#fe9929")
  ) +
  labs(
    title = "Zonas de potencial oportunidad ecológico",
    subtitle = "Escarabajo estercolero vs. uso potencial para ganadería",
    x = "Longitud",
    y = "Latitud"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ------------------ 11. Mapa binario (zonas críticas) -------
cruce_df$binario <- ifelse(cruce_df$valor >= q90, 1, 0)

ggplot() +
  geom_raster(data = cruce_df, aes(x = x, y = y, fill = factor(binario))) +
  geom_sf(data = col_sf, fill = NA, color = "white", size = 0.4) +
  scale_fill_manual(
    values = c("0" = "grey90", "1" = "red"),
    labels = c("Baja sinergia", "Alta sinergia ecológica"),
    name = "Clasificación"
  ) +
  labs(
    title = "Zonas de oportunidad donde coincide idoneidad × uso potencial",
    subtitle = "Clasificación binaria (cuantil 90%)",
    x = "Longitud",
    y = "Latitud"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid = element_blank()
  )
