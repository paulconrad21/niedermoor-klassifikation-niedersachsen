# ============================================================
# 03_build_feature_table.R
# Modelltabelle aus Trainingsdaten, GEE-Features und externen Rasterdaten bauen
# ============================================================

library(sf)
library(terra)
library(ggplot2)
library(this.path)

script_dir <- this.path::here()
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(project_dir)

cat("Projektordner:", project_dir, "\n")


# Eingabeordner
processed_dir <- file.path(project_dir, "data", "processed")
gee_dir <- file.path(project_dir, "data", "raw", "gee_exports")
external_dir <- file.path(project_dir, "data", "raw", "external")

# Ausgabeordner
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1 Eingabedateien
# ============================================================

training_polygons_file <- file.path(processed_dir, "training_polygons_clean.gpkg")

s2_file <- file.path(gee_dir, "training_features_s2_2025.csv")
s1_file <- file.path(gee_dir, "training_features_s1_2025.csv")
lst_file <- file.path(gee_dir, "training_features_landsat_lst_2025.csv")
jrc_file <- file.path(gee_dir, "training_features_jrc_water.csv")

dgm_file <- file.path(
  external_dir,
  "dgm200.utm32s.gridascii",
  "dgm200",
  "dgm200_utm32s.asc"
)

ph_0_30_file <- file.path(
  external_dir,
  "thuenen_institut",
  "thuenen_alacarte_pH_0_30.tif"
)

ph_30_100_file <- file.path(
  external_dir,
  "thuenen_institut",
  "thuenen_alacarte_pH_30_100.tif"
)

# ============================================================
# 2 Trainingsdaten laden
# ============================================================

training_polygons = st_read(training_polygons_file, quiet = TRUE)

# Die Attributtabelle wird direkt aus dem GeoPackage erzeugt.
# Dadurch brauchen wir keine separate training_polygons_attributes.csv.
training_attributes <- st_drop_geometry(training_polygons)

cat("\nTrainingspolygone:\n")
print(nrow(training_attributes))

cat("\nKlassenverteilung:\n")
print(table(training_attributes$class_binary, useNA = "ifany"))

# ============================================================
# 3 GEE-Featuretabellen einlesen
# ============================================================

s2_features <- read.csv(s2_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA", "null"))
s1_features <- read.csv(s1_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA", "null"))
lst_features <- read.csv(lst_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA", "null"))
jrc_features <- read.csv(jrc_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA", "null"))

remove_columns <- c("system:index", ".geo", "area_m2", "cls_bin", "cls_det", "is_moor")

s2_features <- s2_features[, !(names(s2_features) %in% remove_columns)]
s1_features <- s1_features[, !(names(s1_features) %in% remove_columns)]
lst_features <- lst_features[, !(names(lst_features) %in% remove_columns)]
jrc_features <- jrc_features[, !(names(jrc_features) %in% remove_columns)]

# Feature-Spalten numerisch machen
for (data_name in c("s2_features", "s1_features", "lst_features", "jrc_features")) {
  data_object <- get(data_name)
  
  feature_columns <- names(data_object)[names(data_object) != "poly_id"]
  
  for (column in feature_columns) {
    data_object[[column]] <- as.numeric(data_object[[column]])
  }
  
  assign(data_name, data_object)
}

cat("\nGEE-Features:\n")
cat("Sentinel-2:", ncol(s2_features) - 1, "Features\n")
cat("Sentinel-1:", ncol(s1_features) - 1, "Features\n")
cat("Landsat LST:", ncol(lst_features) - 1, "Features\n")
cat("JRC Water:", ncol(jrc_features) - 1, "Features\n")

# ============================================================
# 4 GEE-Features mit Trainingsdaten verbinden
# ============================================================

model_table <- merge(training_attributes, s2_features, by ="poly_id", all.x = TRUE)
model_table <- merge(model_table, s1_features, by = "poly_id", all.x = TRUE)
model_table <- merge(model_table, lst_features, by = "poly_id", all.x = TRUE)
model_table <- merge(model_table, jrc_features, by = "poly_id", all.x = TRUE)

cat("\nZwischenstand Modelltabelle:\n")
cat("Zeilen:", nrow(model_table), "\n")
cat("Spalten:", ncol(model_table), "\n")

# ============================================================
# 5 Sentinel-2 Zeitreihenfeatures berechnen
# ============================================================

s2_columns <- grep("^S2_m[0-9]{2}_", names(model_table), value = TRUE)
s2_variables <- sub("^S2_m[0-9]{2}_", "", s2_columns)
s2_variables <- unique(s2_variables)

for (variable in s2_variables) {
  variable_columns <- grep(paste0("^S2_m[0-9]{2}_", variable, "$"), names(model_table), value = TRUE)
  
  values <- as.matrix(model_table[, variable_columns])
  storage.mode(values) <- "numeric"
  
  annual_mean <- rowMeans(values, na.rm = TRUE)
  annual_mean[is.nan(annual_mean)] <- NA
  
  annual_min <- apply(values, 1, min, na.rm = TRUE)
  annual_max <- apply(values, 1, max, na.rm = TRUE)
  annual_sd <- apply(values, 1, sd, na.rm = TRUE)
  
  annual_min[is.infinite(annual_min)] <- NA
  annual_max[is.infinite(annual_max)] <- NA
  annual_sd[is.infinite(annual_sd)] <- NA
  
  model_table[[paste0("S2_", variable, "_annual_mean")]] <- annual_mean
  model_table[[paste0("S2_", variable, "_annual_min")]] <- annual_min
  model_table[[paste0("S2_", variable, "_annual_max")]] <- annual_max
  model_table[[paste0("S2_", variable, "_annual_amplitude")]] <- annual_max - annual_min
  model_table[[paste0("S2_", variable, "_annual_sd")]] <- annual_sd  
}
cat("\nS2-Zeitreihenfeatures ergänzt.\n")

# ============================================================
# 6 Sentinel-1 Zeitreihenfeatures berechnen
# ============================================================

s1_columns <- grep("^S1_m[0-9]{2}_(ASC|DESC)_", names(model_table), value = TRUE)

s1_variables <- sub("^S1_m[0-9]{2}_", "", s1_columns)
s1_variables <- unique(s1_variables)

for (variable in s1_variables) {
  variable_columns <- grep(paste0("^S1_m[0-9]{2}_", variable, "$"), names(model_table), value = TRUE)
  
  values <- as.matrix(model_table[, variable_columns])
  storage.mode(values) <- "numeric"
  
  annual_mean <- rowMeans(values, na.rm = TRUE)
  annual_mean[is.nan(annual_mean)] <- NA
  
  annual_min <- apply(values, 1, min, na.rm = TRUE)
  annual_max <- apply(values, 1, max, na.rm = TRUE)
  annual_sd <- apply(values, 1, sd, na.rm = TRUE)
  
  annual_min[is.infinite(annual_min)] <- NA
  annual_max[is.infinite(annual_max)] <- NA
  annual_sd[is.nan(annual_sd)] <- NA
  
  model_table[[paste0("S1_", variable, "_annual_mean")]] <- annual_mean
  model_table[[paste0("S1_", variable, "_annual_min")]] <- annual_min
  model_table[[paste0("S1_", variable, "_annual_max")]] <- annual_max
  model_table[[paste0("S1_", variable, "_annual_amplitude")]] <- annual_max - annual_min
  model_table[[paste0("S1_", variable, "_annual_sd")]] <- annual_sd
}

cat("\nS1-Zeitreihenfeatures ergänzt.\n")

# ============================================================
# 7 DGM200 Terrainfeatures berechnen
# ============================================================

dgm <- rast(dgm_file)
names(dgm) <- "DGM200_elevation_m"

cat("\nDGM200 Auflösung:\n")
print(res(dgm))

# Trainingspolygone in das CRS des DGM bringen.
# Fuer die Rasterextraktion brauchen wir spaeter ein terra-SpatVector-Objekt.
training_polygons_dgm <- st_transform(training_polygons, crs(dgm))
polygons_dgm <- vect(training_polygons_dgm)

# Nur Umgebung der Trainingsdaten berechnen.
# Dafuer nehmen wir die Bounding Box aus sf und erweitern sie um 4000 m.
# Diese Variante ist robuster als ext(polygons_dgm), das bei manchen Setups zickt.
bbox_dgm <- st_bbox(training_polygons_dgm)

dgm_ext <- terra::ext(c(
  as.numeric(bbox_dgm["xmin"]) - 4000,
  as.numeric(bbox_dgm["xmax"]) + 4000,
  as.numeric(bbox_dgm["ymin"]) - 4000,
  as.numeric(bbox_dgm["ymax"]) + 4000
))

dgm_local <- terra::crop(dgm, dgm_ext)

# Hangneigung
slope <- terrain(dgm_local, v = "slope", unit = "degrees")
names(slope) <- "DGM200_slope_deg"

# Nachbarschaftsfenster:
# DGM200 hat etwa 200 m Rasterweite.
# 5x5 Zellen ≈ 1 km, 15x15 Zellen ≈ 3 km.
window_1km <- matrix(1, nrow = 5, ncol = 5)
window_3km <- matrix(1, nrow = 15, ncol = 15)

mean_1km <- focal(dgm_local, w = window_1km, fun = "mean", na.rm = TRUE)
mean_3km <- focal(dgm_local, w = window_3km, fun = "mean", na.rm = TRUE)

names(mean_1km) <- "DGM200_mean_1km_m"
names(mean_3km) <- "DGM200_mean_3km_m"

tpi_1km <- dgm_local - mean_1km
tpi_3km <- dgm_local - mean_3km

names(tpi_1km) <- "DGM200_TPI_1km_m"
names(tpi_3km) <- "DGM200_TPI_3km_m"

terrain_stack <- c(
  dgm_local,
  slope,
  mean_1km,
  mean_3km,
  tpi_1km,
  tpi_3km
)

terrain_values <- extract(
  terrain_stack,
  polygons_dgm,
  fun = mean,
  na.rm = TRUE,
  exact = TRUE
)

terrain_features <- data.frame(
  poly_id = training_polygons$poly_id[terrain_values$ID],
  terrain_values[, names(terrain_values) != "ID"],
  check.names = FALSE
)

model_table <- merge(model_table, terrain_features, by = "poly_id", all.x = TRUE)

cat("\nDGM200-Features ergänzt.\n")

# ============================================================
# 8 Boden-pH aus Thünen-Rastern extrahieren
# ============================================================

ph_0_30_raster <- rast(ph_0_30_file)
ph_30_100_raster <- rast(ph_30_100_file)

polygons_ph_0_30 <- project(vect(training_polygons), crs(ph_0_30_raster))
polygons_ph_30_100 <- project(vect(training_polygons), crs(ph_30_100_raster))

ph_0_30_values <- extract(
  ph_0_30_raster,
  polygons_ph_0_30,
  fun = mean,
  na.rm = TRUE,
  exact = TRUE
)

ph_30_100_values <- extract(
  ph_30_100_raster,
  polygons_ph_30_100,
  fun = mean,
  na.rm = TRUE,
  exact = TRUE
)

ph_features <- data.frame(
  poly_id = training_polygons$poly_id,
  Soil_pH_0_30cm = ph_0_30_values[[2]],
  Soil_pH_30_100cm = ph_30_100_values[[2]]
)

model_table <- merge(model_table, ph_features, by = "poly_id", all.x = TRUE)

cat("\nBoden-pH-Features ergänzt.\n")

# ============================================================
# 9 Feature-Gruppen erstellen
# ============================================================

id_columns <- c(
  "poly_id",
  "class_binary",
  "is_peatland",
  "class_detail",
  "class_raw",
  "class_spec_raw",
  "area_m2",
  "x_utm",
  "y_utm"
)

predictor_columns <- names(model_table)[!(names(model_table) %in% id_columns)]

feature_groups <- data.frame(
  feature = predictor_columns,
  group = NA_character_,
  stringsAsFactors = FALSE
)

feature_groups$group[grepl("^S2_m[0-9]{2}_", feature_groups$feature)] <- "S2_monthly"
feature_groups$group[grepl("^S2_", feature_groups$feature) & is.na(feature_groups$group)] <- "S2_temporal"

feature_groups$group[grepl("^S1_m[0-9]{2}_", feature_groups$feature)] <- "S1_monthly"
feature_groups$group[grepl("^S1_", feature_groups$feature) & is.na(feature_groups$group)] <- "S1_temporal"

feature_groups$group[grepl("^LST_", feature_groups$feature)] <- "Landsat_LST"
feature_groups$group[grepl("^JRC_", feature_groups$feature)] <- "Hydrology"
feature_groups$group[grepl("^DGM200_", feature_groups$feature)] <- "Terrain"
feature_groups$group[grepl("^Soil_", feature_groups$feature)] <- "Soil"

feature_groups$group[is.na(feature_groups$group)] <- "Other"

cat("\nFeature-Gruppen:\n")
print(table(feature_groups$group))

# ============================================================
# 10 Komplett leere Features entfernen
# ============================================================

empty_features <- character(0)

for (feature in feature_groups$feature) {
  values <- model_table[[feature]]
  
  if (is.numeric(values)) {
    if(all(is.na(values))) {
      empty_features <- c(empty_features,feature)
    }
  }
}
if (length(empty_features) > 0) {
  cat("\nKomplett leere Features werden entfernt:\n")
  print(empty_features)
  
  model_table <- model_table[, !(names(model_table) %in% empty_features)]
  feature_groups <- feature_groups[!(feature_groups$feature %in% empty_features),]
}
cat("\nAnzahl nutzbarer Predictor-Features:\n")
print(nrow(feature_groups))

# ============================================================
# 11 Einfache Kontrollplots
# ============================================================

group_table <- as.data.frame(table(feature_groups$group))
names(group_table) <- c("group", "n_features")

ggplot(group_table, aes(x = reorder(group, n_features), y = n_features, fill = group)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "Anzahl Features je Variablengruppe",
    x = NULL,
    y = "Anzahl Features"
  ) +
  theme_minimal()

# Beispielplot: Höhe nach Klasse
ggplot(model_table, aes(x = class_binary, y = DGM200_elevation_m, fill = class_binary)) +
  geom_boxplot() +
  labs(
    title = "DGM200-Höhe nach Klasse",
    x = NULL,
    y = "Höhe in m"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Beispielplot: Boden-pH nach Klasse
ggplot(model_table, aes(x = class_binary, y = Soil_pH_0_30cm, fill = class_binary)) +
  geom_boxplot() +
  labs(
    title = "Boden-pH 0-30 cm nach Klasse",
    x = NULL,
    y = "pH"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# ============================================================
# 12 Ergebnisse speichern
# ============================================================

model_table_file <- file.path(processed_dir, "model_table.csv")
feature_groups_file <- file.path(processed_dir, "feature_groups.csv")

write.csv(
  model_table,
  model_table_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  feature_groups,
  feature_groups_file,
  row.names = FALSE
)

cat("\nGespeichert:\n")
cat(model_table_file, "\n")
cat(feature_groups_file, "\n")

cat("\nFinale Übersicht:\n")
cat("Polygone:", nrow(model_table), "\n")
cat("Tabellenspalten:", ncol(model_table), "\n")
cat("Predictor-Features:", nrow(feature_groups), "\n")
cat("Niedermoor:", sum(model_table$is_peatland == 1, na.rm = TRUE), "\n")
cat("Kein Niedermoor:", sum(model_table$is_peatland == 0, na.rm = TRUE), "\n")

cat("\nSkript 03 fertig.\n")
