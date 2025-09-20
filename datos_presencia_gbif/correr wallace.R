if (!require(remotes)) install.packages("remotes")
remotes::install_local("C:/Users/usuario/Downloads/wallace.zip", force = TRUE)

library(wallace)
wallace::run_wallace()
# Guardar la sesión manualmente
saveRDS(wallace.env, "mi_sesion_wallace18_09_2025.rds")

# Volver a cargarla
wallace.env <- readRDS("mi_sesion_wallace18_09_2025.rds")
