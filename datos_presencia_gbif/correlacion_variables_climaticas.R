#getwd()

# Definir la ruta base 
oficina_path <- "D:/Usuario/Desktop/ROBERTO_CEBALLOS/MaestrÃ­a/SEMESTRE III/PROYECTO INTEGRADOR/"
robert_path <- "C:/Users/rober/Documents/MAESTRIA/PROYECTO INTEGRADOR/REPOSITORIO"
mayra_path <- "C:/Users/usuario/Documents/Posgrado/Proyecto integrador/REPOSITORIO/"
base_path <- robert_path;

setwd(base_path)
# ============================================================
# ANÁLISIS DE CORRELACIÓN + VIF AJUSTADO + HEATMAP BONITO
# ============================================================

library(terra)
library(usdm)
library(ggplot2)
library(reshape2)
library(GGally)


# ------------------------------------------------------------
# 0️⃣ Cargar stack de variables recortadas
# ------------------------------------------------------------

vars_path <- file.path(
  base_path,
  "beetle-habitat-prediction",
  "datos_presencia_gbif",
  "variables_M_ecoregiones"
)

clima <- rast(list.files(vars_path, pattern = "\\.tif$", full.names = TRUE))
cat("Variables cargadas desde:\n", vars_path, "\n")

# ------------------------------------------------------------
# 1️⃣ ELIMINAR bio_3, bio_7, stress y swc_fr
# ------------------------------------------------------------

# Crear patrón general
vars_excluir <- c("bio_3$", "bio_7$", "stress$", "swc_fr$")

# Encontrar índices de las variables a remover
idx_excluir <- grep(paste(vars_excluir, collapse = "|"), names(clima))

# Hacer filtro
clima_cor <- clima[[-idx_excluir]]

cat("\nVariables incluidas en el análisis de correlación (excluyendo bio_3, bio_7, stress, swc_fr):\n")
print(names(clima_cor))

# ------------------------------------------------------------
# 2️⃣ MUESTRA ALEATORIA PARA CORRELACIÓN
# ------------------------------------------------------------

#set.seed(123)
sample_points <- spatSample(clima_cor, size = 10000, method = "random")
df_clima <- as.data.frame(sample_points)

# ------------------------------------------------------------
# 3️⃣ MATRIZ DE CORRELACIÓN
# ------------------------------------------------------------

cor_matrix <- cor(df_clima, use = "complete.obs")
cat("\nMatriz de correlación calculada.\n")

# ------------------------------------------------------------
# 🔥 4️⃣ HEATMAP BONITO
# ------------------------------------------------------------

cor_melt <- melt(cor_matrix)

heatmap_plot <- ggplot(cor_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Correlación"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Heatmap de correlación entre variables climáticas",
    x = "",
    y = ""
  )

print(heatmap_plot)

# ------------------------------------------------------------
# 5️⃣ VIF AJUSTADO
# ------------------------------------------------------------

# === 9. Seleccionar variables con baja colinealidad ===
vif_result <- vifcor(df_clima, th = 0.6)
variables_seleccionadas <- vif_result@results$Variables
print(variables_seleccionadas)

cat("\nVariables seleccionadas con VIF umbral = 0.6:\n")
print(variables_seleccionadas)

clima_sel <- subset(clima_cor, variables_seleccionadas)

# ------------------------------------------------------------
# 6️⃣ GUARDAR RESULTADO
# ------------------------------------------------------------

out_file <- file.path(vars_path, "variables_seleccionadas_vif.tif")
writeRaster(clima_sel, out_file, overwrite = TRUE)

cat("\nArchivo guardado en:\n", out_file, "\n")

# ============================================
# GRÁFICO DE CORRELACIÓN ESTILO GGPAIRS
# ============================================



# Convertimos el raster seleccionado a data frame (solo muestra aleatoria)
set.seed(42)
sample_sel <- spatSample(clima_sel, size = 1000, method = "random")
df_sel <- as.data.frame(sample_sel)

# Gráfico bonito con correlaciones
plot_pairs <- GGally::ggpairs(
  df_sel,
  title = "Relaciones entre variables seleccionadas (no correlacionadas)",
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
  diag   = list(continuous = wrap("densityDiag"))
)

print(plot_pairs)



