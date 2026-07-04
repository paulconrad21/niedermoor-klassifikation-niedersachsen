# Niedermoor-Klassifikation Niedersachsen

Dieses Repository enthält den Projektworkflow zur binären Erkennung von **Niedermoor** gegenüber **Kein Niedermoor** in Niedersachsen. Ziel ist es, verschiedene Fernerkundungs- und Umweltvariablen hinsichtlich ihrer Eignung für die Niedermoor-Erkennung zu vergleichen.

## Kurzbeschreibung

Für kuratierte Trainingspolygone wurden satellitenbasierte und externe Umweltvariablen extrahiert und zu einer Modelltabelle zusammengeführt. Anschließend wurden korrelierte und technisch ungeeignete Variablen reduziert. Die finale Modellierung erfolgt mit einem Random-Forest-Modell und nested spatial cross-validation, damit die räumliche Übertragbarkeit realistischer bewertet wird.

## Datenquellen

Verwendete Eingangsdaten:

- Trainingspolygone mit den Klassen `Niedermoor` und `Kein_Niedermoor`
- Sentinel-2 Monatsfeatures für 2025, exportiert über Google Earth Engine
- Sentinel-1 Monatsfeatures für 2025, exportiert über Google Earth Engine
- Landsat-Oberflächentemperaturfeatures, exportiert über Google Earth Engine
- JRC Global Surface Water Features, exportiert über Google Earth Engine
- DGM200 des Bundesamts für Kartographie und Geodäsie
- Boden-pH-Daten des Thünen-Instituts

Die Rohdaten liegen im Ordner `data/raw/`. Abgeleitete Tabellen und vorbereitete Daten liegen in `data/processed/`.

## Externe Rohdaten

Große externe Rohdaten werden **nicht** direkt in diesem Repository versioniert. Sie müssen bei Bedarf von den Originalquellen heruntergeladen und lokal in `data/raw/external/` abgelegt werden.

Für dieses Projekt wurden folgende externe Datensätze verwendet:

### DGM200

Digitales Geländemodell mit 200 m Gitterweite des Bundesamts für Kartographie und Geodäsie.

Verwendete lokale Datei:

```text
data/raw/external/dgm200.utm32s.gridascii/dgm200/dgm200_utm32s.asc
```

Im Projekt wurden daraus Terrainvariablen abgeleitet, unter anderem:

- Höhe
- Hangneigung
- mittlere Umgebungshöhe
- Topographic Position Index, TPI

### Boden-pH des Thünen-Instituts

Boden-pH-Raster aus dem Thünen-ALaCarte-Datenangebot.

Verwendete lokale Dateien:

```text
data/raw/external/thuenen_institut/thuenen_alacarte_pH_0_30.tif
data/raw/external/thuenen_institut/thuenen_alacarte_pH_30_100.tif
```

Im Projekt wurden daraus polygonweise mittlere pH-Werte für zwei Tiefenbereiche abgeleitet:

- `Soil_pH_0_30cm`
- `Soil_pH_30_100cm`

## Projektstruktur

```text
Projekt/
├── data/
│   ├── raw/
│   │   ├── training_polygons/
│   │   ├── gee_exports/
│   │   └── external/
│   └── processed/
├── outputs/
│   ├── figures/
│   ├── tables/
│   └── models/
├── scripts/
└── README.md
```

## Workflow

Die Skripte sind nummeriert und sollten in dieser Reihenfolge ausgeführt werden:

### 01_prepare_training_data.R

Bereitet die Trainingspolygone vor.

- liest die rohen Trainingspolygone ein
- vereinheitlicht Klassen und IDs
- berechnet Fläche und Mittelpunktkoordinaten
- speichert ein bereinigtes GeoPackage
- erzeugt ein Shapefile für Google Earth Engine

Wichtige Outputs:

- `data/processed/training_polygons_clean.gpkg`
- `data/processed/training_polygons_for_gee/`

### 02_export_training_features_gee.js

Google-Earth-Engine-Skript zum Export der Trainingsfeatures.

- extrahiert Sentinel-2-, Sentinel-1-, Landsat- und JRC-Features
- exportiert Tabellen nach Google Drive
- die exportierten CSV-Dateien werden anschließend lokal in `data/raw/gee_exports/` abgelegt

Wichtige Inputs:

- `training_polygons_clean` als GEE Asset

Wichtige Outputs:

- `training_features_s2_2025.csv`
- `training_features_s1_2025.csv`
- `training_features_landsat_lst_2025.csv`
- `training_features_jrc_water.csv`

### 03_build_feature_table.R

Baut die zentrale Modelltabelle.

- liest die bereinigten Trainingspolygone
- liest die GEE-Featuretabellen
- ergänzt externe Rasterdaten wie DGM200 und Boden-pH
- berechnet einfache Zeitreihenkennwerte
- erstellt eine Feature-Gruppentabelle

Wichtige Outputs:

- `data/processed/model_table.csv`
- `data/processed/feature_groups.csv`

### 04_reduce_features.R

Reduziert die Anzahl der Variablen vor der Modellierung.

- entfernt technisch ungeeignete Features
- prüft fehlende Werte und konstante Variablen
- reduziert stark korrelierte Features
- erstellt eine reduzierte Kandidatenliste für die Modellierung

Wichtiger Output:

- `data/processed/feature_candidates_reduced.csv`

### 05_forward_selection_binary_rf.R

Führt die finale Modellierung und Featurebewertung durch.

- nutzt einen binären Random Forest mit `ranger`
- verwendet nested spatial cross-validation
- führt Forward Selection innerhalb der Trainingsfolds durch
- bewertet die äußeren Testfolds unabhängig
- trainiert ein finales Modell mit der Feature-Union
- berechnet Permutation Importance

Wichtige Outputs:

- `outputs/tables/05_nested_selected_features_by_outer_fold.csv`
- `outputs/tables/05_nested_feature_selection_frequency.csv`
- `outputs/tables/05_nested_outer_test_metrics_threshold_05.csv`
- `outputs/tables/05_nested_final_union_permutation_importance.csv`
- `outputs/figures/05_feature_selection_frequency_nested.png`
- `outputs/figures/05_permutation_importance_final_union.png`
- `outputs/models/05_final_binary_ranger_random_forest.rds`

## Methodischer Kern

Das finale Modell nutzt eine **nested spatial cross-validation**:

1. Die Trainingspolygone werden räumlich in fünf Folds eingeteilt.
2. Ein Fold wird jeweils als äußerer Testfold zurückgehalten.
3. Auf den übrigen vier Folds erfolgt die Feature Selection mit innerer räumlicher Cross-Validation.
4. Das ausgewählte Feature-Set wird anschließend auf dem äußeren Testfold bewertet.
5. Nach allen fünf Durchläufen werden die ausgewählten Features zusammengeführt und über Permutation Importance bewertet.

Dadurch wird vermieden, dass der Testfold bereits während der Feature-Auswahl zur Optimierung genutzt wird.

## Software

Die Analyse wurde in R und Google Earth Engine umgesetzt.

Wichtige R-Pakete:

- `sf`
- `terra`
- `ggplot2`
- `ranger`
- `yardstick`
- `this.path`

## Hinweise zur Reproduzierbarkeit

Die Skripte nutzen relative Pfade innerhalb des Projektordners. Sie sollten daher aus dem Ordner `scripts/` heraus oder über RStudio als Skriptdateien ausgeführt werden.

Empfohlene Reihenfolge:

```r
source("scripts/01_prepare_training_data.R")
source("scripts/03_build_feature_table.R")
source("scripts/04_reduce_features.R")
source("scripts/05_forward_selection_binary_rf.R")
```

Das GEE-Skript `02_export_training_features_gee.js` wird separat im Google Earth Engine Code Editor ausgeführt.

## Ergebnisinterpretation

Die bisherige Modellierung zeigt, dass vor allem Sentinel-2-basierte spektrale und saisonale Variablen für die Niedermoor-Erkennung relevant sind. Besonders häufig bzw. wichtig waren Variablen aus den Bereichen NBR/NBR2, Red-Edge, NIR/SWIR und Feuchteindizes. Sentinel-1-Features lieferten ergänzende Informationen. Externe Standortdaten wie Boden-pH oder Temperatur traten vereinzelt auf, waren aber weniger stabil.

Die räumliche Übertragbarkeit ist foldabhängig und damit ein wichtiger Bestandteil der Ergebnisinterpretation.
