# ============================================================
# 04_reduce_features.R
# Feature-Reduktion vor der Modellierung
# ============================================================

library(ggplot2)
library(this.path)

script_dir <- this.path::here()
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(project_dir)

cat("Projektordner:", project_dir, "\n")


# Eingaben
processed_dir <- file.path(project_dir, "data", "processed")

model_table_file <- file.path(processed_dir, "model_table.csv")
feature_groups_file <- file.path(processed_dir, "feature_groups.csv")

# Ausgaben
reduced_features_file <- file.path(processed_dir, "feature_candidates_reduced.csv")
feature_summary_file <- file.path(processed_dir, "feature_reduction_summary.csv")
correlation_removed_file <- file.path(processed_dir, "feature_correlation_removed.csv")

# Einstellungen
max_missing_share <- 0.40
min_unique_values <- 2
correlation_threshold <- 0.90
correlation_method <- "spearman"

cat("\nFeature-Reduktion gestartet.\n")
cat("Maximal erlaubter Anteil fehlender Werte:", max_missing_share, "\n")
cat("Korrelationsschwelle:", correlation_threshold, "\n")
cat("Korrelationsmethode:", correlation_method, "\n")

# ============================================================
# 1 Daten laden
# ============================================================

model_table <- read.csv(
  model_table_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN", "null")
)

feature_groups <- read.csv(
  feature_groups_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nModelltabelle:\n")
cat("Zeilen:", nrow(model_table), "\n")
cat("Spalten:", ncol(model_table), "\n")

cat("\nFeatures laut feature_groups.csv:\n")
print(table(feature_groups$group))

# ============================================================
# 2 Technische Feature-Zusammenfassung
# ============================================================

feature_summary <- data.frame(
  feature = feature_groups$feature,
  group = feature_groups$group,
  is_numeric = NA,
  missing_share = NA,
  unique_values = NA,
  standard_deviation = NA,
  keep_technical = NA,
  remove_reason = NA,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(feature_summary))) {
  feature <- feature_summary$feature[i]
  values <- model_table[[feature]]
  
  feature_summary$is_numeric[i] <- is.numeric(values)
  feature_summary$missing_share[i] <- mean(is.na(values))
  feature_summary$unique_values[i] <- length(unique(values[!is.na(values)]))
  
  if (is.numeric(values)) {
    feature_summary$standard_deviation[i] <- sd(values, na.rm = TRUE)
  }
  
  if (!is.numeric(values)) {
    feature_summary$keep_technical[i] <- FALSE
    feature_summary$remove_reason[i] <- "nicht numerisch"
  } else if (feature_summary$missing_share[i] > max_missing_share) {
    feature_summary$keep_technical[i] <- FALSE
    feature_summary$remove_reason[i] <- "zu viele fehlende Werte"
  } else if (feature_summary$unique_values[i] < min_unique_values) {
    feature_summary$keep_technical[i] <- FALSE
    feature_summary$remove_reason[i] <- "konstant oder fast konstant"
  } else if (is.na(feature_summary$standard_deviation[i]) ||
             feature_summary$standard_deviation[i] == 0) {
    feature_summary$keep_technical[i] <- FALSE
    feature_summary$remove_reason[i] <- "keine Streuung"
  } else {
    feature_summary$keep_technical[i] <- TRUE
    feature_summary$remove_reason[i] <- "behalten"
  }
}

cat("\nTechnische Feature-Prüfung:\n")
print(table(feature_summary$remove_reason))

candidate_features <- feature_summary$feature[feature_summary$keep_technical]

cat("\nFeatures nach technischem Filter:", length(candidate_features), "\n")

# ============================================================
# 3 Einfache Priorität für Features vergeben
# ============================================================

feature_summary$priority_score <- 0

for (i in seq_len(nrow(feature_summary))) {
  feature <- feature_summary$feature[i]
  group <- feature_summary$group[i]
  
  score <- 0
  
  # Zeitliche Zusammenfassungen sind oft stabiler als einzelne Monatswerte
  if (grepl("annual|spring|summer|growing|amplitude", feature, ignore.case = TRUE)) {
    score <- score + 3
  }
  
  # Einzelne Monatswerte sind nützlich, aber oft stärker wetterabhängig
  if (grepl("^S2_m[0-9]{2}_|^S1_m[0-9]{2}_", feature)) {
    score <- score - 1
  }
  
  # Gut interpretierbare Vegetations- und Feuchteindizes
  if (grepl("NDVI|NDWI|NDMI|EVI|MNDWI|NDRE|NBR", feature, ignore.case = TRUE)) {
    score <- score + 2
  }
  
  # Externe Standortdaten sind fachlich gut erklärbar
  if (group %in% c("Terrain", "Hydrology", "Soil", "Landsat_LST")) {
    score <- score + 1
  }
  
  # Weniger fehlende Werte sind besser
  score <- score - feature_summary$missing_share[i]
  
  feature_summary$priority_score[i] <- score
}

# ============================================================
# 4 Korrelationen berechnen
# ============================================================

feature_matrix <- model_table[, candidate_features]

correlation_matrix <- cor(
  feature_matrix,
  use = "pairwise.complete.obs",
  method = correlation_method
)

correlation_matrix[!is.finite(correlation_matrix)] <- NA

cat("\nKorrelationsmatrix berechnet für", ncol(feature_matrix), "Features.\n")

# ============================================================
# 5 Stark korrelierte Features reduzieren
# ============================================================

candidate_summary <- feature_summary[feature_summary$feature %in% candidate_features, ]

# Features sortieren:
# zuerst hohe Priorität, dann wenig fehlende Werte
candidate_summary <- candidate_summary[
  order(
    -candidate_summary$priority_score,
    candidate_summary$missing_share,
    candidate_summary$feature
  ),
]

ordered_features <- candidate_summary$feature

selected_features <- character(0)
removed_features <- data.frame(
  removed_feature = character(0),
  kept_feature = character(0),
  correlation = numeric(0),
  stringsAsFactors = FALSE
)

# ============================================================
# 5 Stark korrelierte Features reduzieren
# ============================================================

candidate_summary <- feature_summary[feature_summary$feature %in% candidate_features, ]

# Features sortieren:
# zuerst hohe Priorität, dann wenig fehlende Werte
candidate_summary <- candidate_summary[
  order(
    -candidate_summary$priority_score,
    candidate_summary$missing_share,
    candidate_summary$feature
  ),
]

ordered_features <- candidate_summary$feature

selected_features <- character(0)
removed_features <- data.frame(
  removed_feature = character(0),
  kept_feature = character(0),
  correlation = numeric(0),
  stringsAsFactors = FALSE
)

for (feature in ordered_features) {
  if (length(selected_features) == 0) {
    selected_features <- c(selected_features, feature)
  } else {
    correlations_to_selected <- correlation_matrix[feature, selected_features]
    max_correlation <- max(abs(correlations_to_selected), na.rm = TRUE)
    
    if (is.na(max_correlation) || max_correlation < correlation_threshold) {
      selected_features <- c(selected_features, feature)
    } else {
      strongest_match <- selected_features[which.max(abs(correlations_to_selected))]
      correlation_value <- correlation_matrix[feature, strongest_match]
      
      removed_features <- rbind(
        removed_features,
        data.frame(
          removed_feature = feature,
          kept_feature = strongest_match,
          correlation = correlation_value,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

cat("\nFeatures vor Korrelationsreduktion:", length(candidate_features), "\n")
cat("Features nach Korrelationsreduktion:", length(selected_features), "\n")
cat("Wegen Korrelation entfernt:", nrow(removed_features), "\n")

# ============================================================
# 6 Reduzierte Feature-Liste bauen
# ============================================================

reduced_features <- feature_summary[feature_summary$feature %in% selected_features, ]

reduced_features <- reduced_features[
  order(reduced_features$group, -reduced_features$priority_score, reduced_features$feature),
]


cat("\nReduzierte Features je Gruppe:\n")
print(table(reduced_features$group))

cat("\nErste reduzierte Features:\n")
print(head(reduced_features[, c("feature", "group", "missing_share", "priority_score")], 20))

# ============================================================
# 7 Einfache Plots
# ============================================================

# Plot 1: Features vor und nach Reduktion je Gruppe
groups <- sort(unique(feature_summary$group))

plot_data <- data.frame(
  group = rep(groups, 2),
  step = c(rep("vor Reduktion", length(groups)), rep("nach Reduktion", length(groups))),
  n_features = NA
)

for (i in seq_along(groups)) {
  group <- groups[i]
  
  plot_data$n_features[plot_data$group == group & plot_data$step == "vor Reduktion"] <-
    sum(feature_summary$group == group & feature_summary$keep_technical)
  
  plot_data$n_features[plot_data$group == group & plot_data$step == "nach Reduktion"] <-
    sum(reduced_features$group == group)
}

ggplot(plot_data, aes(x = group, y = n_features, fill = step)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Feature-Reduktion nach Variablengruppe",
    x = "Feature-Gruppe",
    y = "Anzahl Features",
    fill = "Schritt"
  ) +
  theme_minimal()

# ============================================================
# 8 Ergebnisse speichern
# ============================================================

write.csv(
  reduced_features,
  reduced_features_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  feature_summary,
  feature_summary_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  removed_features,
  correlation_removed_file,
  row.names = FALSE,
  na = ""
)

cat("\nGespeichert:\n")
cat(reduced_features_file, "\n")
cat(feature_summary_file, "\n")
cat(correlation_removed_file, "\n")

cat("\nSkript 04 fertig.\n")
