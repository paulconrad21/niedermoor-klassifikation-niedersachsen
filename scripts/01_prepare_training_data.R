# Skript 1: Trainingspolygone vorbereiten

# Libraries laden
library(sf)
library(ggplot2)
library(this.path)

# Projektordner
script_dir <- this.path::here()
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(project_dir)

cat("Projektordner:", project_dir, "\n")

# Eingabeordner: rohe Trainingspolygone
raw_dir <- file.path(project_dir, "data", "raw", "training_polygons")

# Ausgabeordner: bereinigte Daten
processed_dir <- file.path(project_dir, "data", "processed")
gee_dir <- file.path(processed_dir, "training_polygons_for_gee")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gee_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1 Rohdaten einlesen
# ============================================================

training_file <- file.path(raw_dir, "training_polygons.shp")

training_raw <- st_read(training_file, quiet = TRUE)

print("Anzahl Trainingspolygone:")
print(nrow(training_raw))

print("Spalten in den Rohdaten:")
print(names(training_raw))

print("Erste Zeilen der Attributtabelle:")
print(head(st_drop_geometry(training_raw)))

# ============================================================
# 2 In metrisches Koordinatensystem umwandeln
# ============================================================

training_clean <- st_transform(training_raw, crs = 25832)
print(st_crs(training_clean))

# ============================================================
# 3 Klassen und IDs vorbereiten
# ============================================================

# Eindeutige ID pro Polygon
training_clean$poly_id <- sprintf("P%04d", seq_len(nrow(training_clean)))

# Originalklassen als Text bereinigen
training_clean$class_raw <- trimws(as.character(training_clean$Class))
training_clean$class_spec_raw <- trimws(as.character(training_clean$Class_spec))

training_clean$class_raw[training_clean$class_raw == ""] <- NA
training_clean$class_spec_raw[training_clean$class_spec_raw == ""] <- NA

# Binäre Zielklasse erstellen
class_lower <- tolower(training_clean$class_raw)

training_clean$class_binary <- NA_character_
training_clean$class_binary[class_lower == "niedermoor"] <- "Niedermoor"

training_clean$class_binary[
  class_lower %in% c("kein niedermoor", "kein_niedermoor", "keinniedermoor")
] <- "Kein_Niedermoor"

# Numerische Zielvariable für spätere Modelle
# 1 = Niedermoor, 0 = kein Niedermoor
training_clean$is_peatland <- NA_integer_
training_clean$is_peatland[training_clean$class_binary == "Niedermoor"] <- 1L
training_clean$is_peatland[training_clean$class_binary == "Kein_Niedermoor"] <- 0L

# Genauere Landnutzungsklasse
training_clean$class_detail <- training_clean$class_spec_raw
training_clean$class_detail[training_clean$class_detail == "hochmoor"] <- "Hochmoor"

# Für echte Niedermoore ist die Detailklasse einfach Niedermoor
training_clean$class_detail[training_clean$class_binary == "Niedermoor"] <- "Niedermoor"

# Falls bei Nicht-Niedermooren keine Detailklasse angegeben ist
training_clean$class_detail[
  training_clean$class_binary == "Kein_Niedermoor" &
    is.na(training_clean$class_detail)
] <- "Kein_Niedermoor_unspecified"

print("Verteilung der binären Klassen:")
print(table(training_clean$class_binary, useNA = "ifany"))
print("Verteilung der Detailklassen:")
print(table(training_clean$class_detail, useNA = "ifany"))

# ============================================================
# 4 Geometrie, Fläche und Lage berechnen
# ============================================================

# Geometrien reparieren, falls einzelne Polygone kleine Geometrieprobleme haben
training_clean <- st_make_valid(training_clean)

# Fläche in Quadratmetern
training_clean$area_m2 <- as.numeric(st_area(training_clean))

# Punkt innerhalb jedes Polygons
# Wird für räumliche Gruppierung / Cross-Validation genutzt
centroids <- st_point_on_surface(training_clean)
coords <- st_coordinates(centroids)

training_clean$x_utm <- coords[, 1]
training_clean$y_utm <- coords[, 2]

print("Zusammenfassung der Polygonfläche in m2:")
print(summary(training_clean$area_m2))

# ============================================================
# 5 Einfache Auswertung und Plots
# ============================================================

print("Kreuztabelle: Detailklasse nach binärer Klasse:")
print(table(training_clean$class_detail, training_clean$class_binary, useNA = "ifany"))

# Balkendiagramm: binäre Klassen
class_table <- as.data.frame(table(training_clean$class_binary))
names(class_table) <- c("class_binary", "n")

ggplot(class_table, aes(x = class_binary, y = n, fill = class_binary)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Niedermoor" = "#2b8cbe",
                               "Kein_Niedermoor" = "#f03b20")) +
  labs(
    title = "Anzahl Trainingspolygone je Klasse",
    x = "Klasse",
    y = "Anzahl"
  ) +
  theme_minimal()

# Karte der Trainingspolygone
ggplot(training_clean) +
  geom_sf(aes(color = class_binary), fill = NA, linewidth = 0.7) +
  scale_color_manual(values = c("Niedermoor" = "#2b8cbe",
                                "Kein_Niedermoor" = "#f03b20")) +
  labs(
    title = "Räumliche Verteilung der Trainingspolygone",
    color = "Klasse"
  ) +
  theme_minimal()

# ============================================================
# 6 Bereinigte Daten speichern
# ============================================================

# Nur die wichtigen Spalten behalten
training_export <- training_clean[, c(
  "poly_id",
  "class_binary",
  "is_peatland",
  "class_detail",
  "class_raw",
  "class_spec_raw",
  "area_m2",
  "x_utm",
  "y_utm"
)]

# GeoPackage für die weitere lokale Arbeit
gpkg_file <- file.path(processed_dir, "training_polygons_clean.gpkg")
if (file.exists(gpkg_file)) {
  file.remove(gpkg_file)
}

st_write(
  training_export,
  gpkg_file,
  layer = "training_polygons_clean",
  quiet = TRUE
)

# ============================================================
# 7 Version für Google Earth Engine exportieren
# ============================================================

gee_export <- training_export[, c(
  "poly_id",
  "is_peatland",
  "class_binary",
  "class_detail",
  "area_m2"
)]

# Kurze Feldnamen für Shapefile / GEE
gee_export$is_moor <- gee_export$is_peatland
gee_export$cls_bin <- gee_export$class_binary
gee_export$cls_det <- gee_export$class_detail
gee_export$area_m2 <- round(gee_export$area_m2, 1)

gee_export <- gee_export[, c(
  "poly_id",
  "is_moor",
  "cls_bin",
  "cls_det",
  "area_m2"
)]

gee_file <- file.path(gee_dir, "training_polygons_clean.shp")

st_write(
  gee_export,
  gee_file,
  delete_dsn = TRUE,
  quiet = TRUE
)

print("GEE-Shapefile gespeichert:")
print(gee_file)

print("Fertig. Trainingsdaten sind vorbereitet.")
