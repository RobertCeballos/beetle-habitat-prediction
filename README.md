# beetle-habitat-prediction
PREDICCIÓN DE LA DISTRIBUCIÓN POTENCIAL DE ESCARABAJOS ESTERCOLEROS MEDIANTE MODELOS DE APRENDIZAJE AUTOMÁTICO Y ANÁLISIS DE VARIABLES CLIMÁTICAS Y GEOGRÁFICAS

Explicación de los 4 archivos principales del proyecto (R)
============================================================================================================================================
GBIF_busqueda_limpieza_M-04-09-2025.R  ======================================================================================================
Búsqueda, depuración y preparación de datos de presencia para el modelado de distribución potencial de Dichotomius agenor

Este script reúne todo el flujo de trabajo necesario para construir el dataset final de presencias y variables ambientales que se usarán en Wallace / MaxEnt. Aquí se integran procesos de descarga, limpieza, validación ecológica, extracción ambiental y reducción de variables. Es básicamente el corazón del preprocesamiento del proyecto.

1. Descarga y consolidación de registros biológicos (GBIF)
Carga la lista de especies base (Giraldo et al. 2018).
Consulta sinónimos taxonómicos usando GBIF Backbone Taxonomy para evitar perder registros válidos.
Realiza una descarga masiva de registros georreferenciados en GBIF (hasta 5000 por especie).
Unifica nombres aceptados y consolida en un solo data frame.
Salida:
GBIF_2025-03-23.rds con todos los registros depurados taxonómicamente.

2. Limpieza espacial y ecológica de los registros (filtros M)
Incluye una depuración avanzada y acorde a la literatura:

🔹 Filtros aplicados:
Eliminación de registros fuera del rango geográfico esperado.
Remoción de departamentos no reportados para D. agenor (Montoya-Molina & Vaz-de-Mello, 2021).
Eliminación de duplicados espaciales (coordenadas idénticas o a <1 km de distancia).
Unión espacial con elevaciones ALOS (Google Earth Engine) y filtro final por rango 0–1600 m.
Salida:
sp_df → Tabla final de presencias limpias y con elevación.

3. Construcción del área de calibración (M)
El script genera dos enfoques distintos:

🔸 Opción A: Áreas biogeográficas (WWF Ecoregions 2017)
Se identifican las ecorregiones donde realmente aparece la especie.

🔸 Opción B: Buffer de 50 km
Se calcula un polígono tipo “M” usando un buffer ecológico alrededor de los registros filtrados.

Ambas áreas se guardan como shapefiles en:
datos_presencia_gbif/AREAS_CALIBRACION/

4. Carga y preparación de variables ambientales
🌦️ Variables bioclimáticas

Carga de WorldClim 2.1 (30 arc-seconds).
Ordenamiento, recorte y enmascaramiento de acuerdo al área M seleccionada.
🟫 Variables de suelo (SoilGrids)
Incluye un proceso robusto hecho por ti:
Recorte de ESRI GRID tiles (swc_fr, alpha) usando tu bounding box.
Reproyección, resample al grid de WC2.1 y alineación exacta.
Cálculo del nuevo predictor: stress = 1 – alpha.
Agregado de capas finales al stack ambiental.
Salida:
swc_fr y stress a 1 km, listos para modelar.

5. Extracción de valores ambientales en puntos de presencia
El script genera dos datasets:
presencias_ecoreg_limpio.csv
presencias_buffer_limpio.csv
Cada uno contiene:
Coordenadas
Variables bioclimáticas WC 30s
Variables de suelo
Registros sin NA en bioclimáticas

6. Reducción de variables (correlación, VIF y PCA)
🔹 Correlación Pearson
Se muestrean 10.000 píxeles y se genera la matriz de correlación.
🔹 VIF (Variance Inflation Factor)
Umbral = 10
Retiene solo variables con baja colinealidad.
🔹 PCA avanzado compatible con terra moderno
Uso de prcomp() + proyección raster por bloques.
Exporta:
mapas de componentes
varianza explicada
loadings
heatmaps, screeplot y biplot
Selección automática de PCs acumulando ≥90% de varianza.
Salida:
Archivos .tif multibanda con los componentes ambientales.

7. Preparación final de insumos para Wallace
Genera automáticamente las carpetas:

WALLACE_INPUTS/
    ├─ presence/
    ├─ variables/
    └─ area/

Y exporta:
✔ Presencias (presencias_M.csv)
Solo columnas necesarias: nombre, lon, lat.
✔ Variables ambientales
Multibanda .tif con PCA o VIF (según tu opción).
Versión por capas individuales.
Remoción automática de bio3 y bio7.
✔ Shapefile del área M
Transformado a WGS84.
Estos son exactamente los insumos que Wallace necesita para correr sin problemas.
🧩 Rol del script dentro del proyecto
Este archivo actúa como todo el pipeline de preprocesamiento del proyecto.
Con él se:
Descargan los datos brutos
Limpian ecológicamente
Construye el área M
Preparan todas las variables
Se reducen correlaciones
Y se deja todo listo para el modelado en Wallace / MaxEnt
Es la base de la fase 3 (Preparación de los datos) y parte de la fase 4 (Modelado) del proyecto integrador.


============================================================================================================================================
correlacion_variables_climaticas.R ==========================================================================================================
Evaluación de colinealidad, selección de variables y visualización climática

Este script se enfoca exclusivamente en evaluar la correlación entre variables bioclimáticas (recortadas previamente al área de calibración por ecorregiones) y seleccionar aquellas con baja colinealidad utilizando VIF. Además, produce visualizaciones limpias y útiles para el informe y para validar qué variables entran al modelado.

1. Carga del stack de variables climáticas (recortadas a M-ecoregiones)
El script inicia cargando todos los .tif generados en la carpeta:
datos_presencia_gbif/variables_M_ecoregiones/
Cada TIF corresponde a una variable bioclimática o derivada.

2. Exclusión de variables problemáticas para el análisis
Se eliminan de manera explícita:
bio_3 (Isotermalidad)
bio_7 (Rango anual de temperatura)
Se agregan:
stress (estrés hídrico)
swc_fr (humedad del suelo)
Esto mejora la estabilidad estadística del análisis y se ajusta al enfoque climático puro.

3. Muestreo aleatorio para análisis de correlación
Se extrae una muestra aleatoria de 10.000 píxeles del stack climático con:
valores completos
buena representatividad espacial
bajo costo computacional
Esto produce el data frame base para todos los análisis.

4. Matriz de correlación + Heatmap elegante
Se calcula una matriz de correlación Pearson y se genera un heatmap limpio con:
gradiente azul → blanco → rojo
nombres de variables rotados
estilo minimalista para informes
Producto visual:
Heatmap de correlación entre variables climáticas.
Sirve para evaluar redundancias entre variables antes del modelado.

5. Selección de variables mediante VIF (umbral = 0.6)
Se aplica vifcor() del paquete usdm:
Umbral estricto = 0.6 para evitar multicolinealidad.
Se devuelven las variables no correlacionadas entre sí.
Se crea un stack reducido con solo las variables aceptadas.
Salida:
variables_seleccionadas_vif.tif
Guardado dentro de:
variables_M_ecoregiones/
Esta es la primera versión limpia de predictores climáticos para modelar.


============================================================================================================================================
uso_suelo_ganaderia_1km.R===================================================================================================================
Construcción del raster de uso potencial de suelo para ganadería a 1 km en Colombia
Este script se encarga de integrar dos tipos de información clave para el proyecto:
Capas FAO-GLW3 relacionadas con uso potencial del suelo (grassland/cropland suitability).
Polígonos de densidad ganadera (cattle 2022) también provenientes de FAO.
El objetivo final es generar un raster de 1 km que represente el uso potencial de suelo para ganadería en Colombia, combinando la idoneidad FAO y la cantidad de ganado reportado por polígono.

1. Carga de insumos FAO GLW3 y estandarización del CRS
Se cargan dos insumos:
Raster de reglas FAO (idoneidad de uso potencial del suelo)
rules_p_2022
Polígonos FAO con densidad de ganado (cattle_2022)
Notas importantes del script:
✔ El raster original viene en una proyección tipo Homolosine, por lo que se reproyecta a WGS84 (EPSG:4326).
✔ Esta corrección evita desplazamientos y errores de recorte posteriores.

2. Recorte del raster a Colombia
Usando rnaturalearth, se carga el polígono de Colombia y se hace:
crop() → recorte espacial del raster
mask() → limpia los bordes y deja solo el país
Esto genera el raster base sobre el que se hará la asignación de valores ganaderos.

3. Procesamiento de los polígonos ganaderos FAO
Se cargan los polígonos provenientes del archivo .gpkg que contienen, entre otros atributos, el dato:
cattle_2022 → cantidad o densidad de ganado para ese polígono
Flujo aplicado:
Validación de geometrías
Recorte a Colombia
Selección de columnas clave
Reproyección a WGS84
Conversión a SpatVector para integrarlo con terra

4. Extracción de idoneidad del suelo FAO por polígono
Para cada polígono, se extrae el valor raster de idoneidad del suelo FAO:
vals <- extract(uso_pot_col, cattle_v)
Luego estos valores se guardan como atributo:
idoneidad
NA → se rellena con 0

5. Cálculo del índice final de uso de suelo ganadero
Se construye un indicador combinado:
valor_final=cattle_2022×(idoneidad/100)

Esto entrega un valor ponderado que refleja:
dónde hay ganado
qué tan idóneo es el suelo según FAO
Perfecto para asociarlo después con tus mapas de idoneidad climática del escarabajo.

6. Rasterización del indicador final (1 km)
Los polígonos FAO, con sus valores ponderados, se convierten en un raster de 1 km:
mismo grid del raster FAO reproyectado
sumando valores por celda
fondo = 0 (zonas sin datos ganaderos)
Resultado → raster_final

7. Visualización rápida
El script incluye una vista:
plot(raster_final)
plot(col_v, add = TRUE)
Sirve para revisar la distribución espacial del uso ganadero.
(Lo vimos juntas: las zonas calientes coinciden con valles interandinos y piedemontes ganaderos).

8. Exportación del raster final
El archivo resultante se guarda como:
SUELO/USO_SUELO/raster_ganaderia_final_1km_WGS84.tif

Este raster es el que se cruza en el siguiente script con el mapa de idoneidad climática del escarabajo para identificar zonas de coincidencia entre ganadería e idoneidad ambiental.


==============================================================================================================================================
cruceespacial_idoneidadclimatica_x_usopotencialganaderia.R====================================================================================
Cruce espacial entre la idoneidad climática del escarabajo y el uso potencial del suelo para ganadería (1 km)
Este script integra los dos grandes productos del proyecto:
El mapa de idoneidad climática modelado para Dichotomius agenor (salida de MaxEnt/Wallace).
El raster FAO de uso potencial del suelo para ganadería (1 km, GLW3).
El objetivo es identificar zonas donde coinciden condiciones ambientales favorables para el escarabajo con áreas aptas para actividad ganadera, generando un insumo clave para el análisis territorial.

1. Carga de insumos principales
Raster de idoneidad climática (formato logistic de MaxEnt).
Raster FAO de uso potencial para ganadería (grassland/cropland suitability).
Ambos se alinean en:
extensión
sistema de referencia (EPSG:4326)
resolución (1 km)
Esto garantiza compatibilidad para el análisis.

2. Recorte espacial a Colombia
Con rnaturalearth, se carga el polígono de Colombia y se realiza:
crop()
mask()
sobre ambos rasters, para trabajar únicamente en el ámbito nacional.

3. Resample del raster ganadero a la resolución del mapa climático
Se ajusta el raster FAO usando:
uso_pot_resampled <- resample(uso_pot_col, id_clima_col, method = "bilinear")
Esto garantiza que ambos grids coincidan celda por celda.

4. Cruce espacial (multiplicación celda a celda)
Se genera el raster final:
Cruce=Idoneidad climatica×Uso potencial ganadero
Este producto identifica áreas donde:
✔ hay buena idoneidad para el escarabajo
✔ hay mayor potencial de uso ganadero
Salida guardada como:
raster_cruce_escarabajo_ganaderia.tif

5. Visualización principal (mapa continuo)
Se transforma el raster en un data frame espacial (cruce_df) y se genera un mapa profesional con:
geom_raster()
color scale tipo magma (viridis)
contorno limpio de Colombia
Este mapa muestra el gradiente de coincidencia Escarabajo × Ganadería.

6. Clasificación por cuantiles (75% y 90%)
El script categoriza el cruce continuo en tres niveles:
Baja
Media
Alta
usando cuantiles para garantizar una clasificación relativa al comportamiento del dataset.
Se genera un mapa categórico con colores intuitivos:
Azul oscuro → Baja
Azul medio → Media
Naranja → Alta
Este mapa permite identificar zonas de oportunidad ecológica.

7. Mapa binario de zonas críticas

Se crea una versión más simple y útil para reportes: