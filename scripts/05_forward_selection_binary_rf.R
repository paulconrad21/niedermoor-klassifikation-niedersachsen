# ============================================================
# 05_forward_selection_binary_rf.R
# Nested spatial CV + Forward Selection + Ranger Random Forest
# ============================================================
#
# Ziel:
# - Features wissenschaftlich sauberer auswaehlen.
# - Der aeussere Testfold bleibt jeweils komplett unangetastet.
# - Feature Selection passiert nur innerhalb der uebrigen 4 Folds.
# - Am Ende werden alle in den Outer-Folds gewaehlten Features gesammelt.
# - Das finale Random-Forest-Modell nutzt diese Feature-Union.

library(ranger)
library(yardstick)
library(ggplot2)
library(this.path)

# ============================================================
# 0 Pfade und Einstellungen
# ============================================================

script_dir <- this.path::here()
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(project_dir)

processed_dir <- file.path("data", "processed")
table_dir <- file.path("outputs", "tables")
figure_dir <- file.path("outputs", "figures")
model_dir <- file.path("outputs", "models")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

model_table_file <- file.path(processed_dir, "model_table.csv")
reduced_features_file <- file.path(processed_dir, "feature_candidates_reduced.csv")

set.seed(42)

n_folds <- 5
ntree_selection <- 300
ntree_final <- 800
max_selected_features <- 15
min_improvement <- 0.002
selection_metric <- "pr_auc"
target_precision <- 0.80
num_threads <- max(1, parallel::detectCores() - 1)

cat("Projektordner:", project_dir, "\n")
cat("Nested spatial CV mit", n_folds, "Folds\n")
cat("Threads:", num_threads, "\n")

# ============================================================
# 1 Daten laden
# ============================================================

model_data <- read.csv(
  model_table_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN", "null")
)

feature_candidates <- read.csv(
  reduced_features_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_features <- feature_candidates$feature
candidate_features <- candidate_features[candidate_features %in% names(model_data)]

# Verwaltungs-, Ziel- und Lageinformationen duerfen nicht ins Modell.
forbidden_features <- c(
  "poly_id",
  "id_original",
  "class_binary",
  "is_peatland",
  "class_detail",
  "class_raw",
  "class_spec_raw",
  "area_m2",
  "x_utm",
  "y_utm",
  "target",
  "spatial_fold"
)

candidate_features <- candidate_features[!(candidate_features %in% forbidden_features)]

# Nur numerische Features verwenden.
numeric_features <- character(0)

for (feature in candidate_features) {
  if (is.numeric(model_data[[feature]])) {
    numeric_features <- c(numeric_features, feature)
  }
}

candidate_features <- numeric_features

model_data$class_binary <- as.character(model_data$class_binary)
model_data$target <- factor(
  model_data$class_binary,
  levels = c("Kein_Niedermoor", "Niedermoor")
)

cat("\nDaten geladen.\n")
cat("Polygone:", nrow(model_data), "\n")
cat("Kandidaten-Features:", length(candidate_features), "\n")

cat("\nKlassenverteilung:\n")
print(table(model_data$target))

# ============================================================
# 2 Raeumliche Folds erstellen
# ============================================================

coordinates <- data.frame(
  x = model_data$x_utm,
  y = model_data$y_utm
)

coordinates_scaled <- scale(coordinates)

fold_clusters <- kmeans(
  coordinates_scaled,
  centers = n_folds,
  nstart = 50
)

model_data$spatial_fold <- fold_clusters$cluster
fold_ids <- sort(unique(model_data$spatial_fold))

cat("\nSpatial Folds:\n")
print(table(model_data$spatial_fold))

cat("\nKlassen je Fold:\n")
print(table(model_data$spatial_fold, model_data$target))

fold_plot <- ggplot(model_data, aes(x = x_utm, y = y_utm)) +
  geom_point(aes(color = factor(spatial_fold), shape = target), size = 2.5, alpha = 0.85) +
  coord_equal() +
  labs(
    title = "Raeumliche Cross-Validation-Folds",
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Fold",
    shape = "Klasse"
  ) +
  theme_minimal(base_size = 12)

print(fold_plot)

ggsave(
  filename = file.path(figure_dir, "05_spatial_folds.png"),
  plot = fold_plot,
  width = 8,
  height = 6,
  dpi = 220,
  bg = "white"
)

# ============================================================
# 3 Metriken
# ============================================================

calculate_metrics <- function(truth, probability, threshold) {
  prediction <- ifelse(probability >= threshold, "Niedermoor", "Kein_Niedermoor")

  eval_data <- data.frame(
    truth = factor(truth, levels = c("Kein_Niedermoor", "Niedermoor")),
    .pred_Niedermoor = probability,
    .pred_class = factor(prediction, levels = c("Kein_Niedermoor", "Niedermoor"))
  )

  eval_data <- eval_data[is.finite(eval_data$.pred_Niedermoor), ]

  precision_value <- suppressWarnings(yardstick::precision(
    eval_data,
    truth = truth,
    estimate = .pred_class,
    event_level = "second"
  )$.estimate)

  recall_value <- suppressWarnings(yardstick::recall(
    eval_data,
    truth = truth,
    estimate = .pred_class,
    event_level = "second"
  )$.estimate)

  f1_value <- suppressWarnings(yardstick::f_meas(
    eval_data,
    truth = truth,
    estimate = .pred_class,
    event_level = "second"
  )$.estimate)

  accuracy_value <- yardstick::accuracy(
    eval_data,
    truth = truth,
    estimate = .pred_class
  )$.estimate

  roc_auc_value <- suppressWarnings(yardstick::roc_auc(
    eval_data,
    truth = truth,
    .pred_Niedermoor,
    event_level = "second"
  )$.estimate)

  pr_auc_value <- suppressWarnings(yardstick::pr_auc(
    eval_data,
    truth = truth,
    .pred_Niedermoor,
    event_level = "second"
  )$.estimate)

  truth_positive <- eval_data$truth == "Niedermoor"
  prediction_positive <- eval_data$.pred_class == "Niedermoor"

  tp <- sum(truth_positive & prediction_positive)
  fp <- sum(!truth_positive & prediction_positive)
  tn <- sum(!truth_positive & !prediction_positive)
  fn <- sum(truth_positive & !prediction_positive)

  result <- data.frame(
    threshold = threshold,
    precision = precision_value,
    recall = recall_value,
    f1 = f1_value,
    accuracy = accuracy_value,
    roc_auc = roc_auc_value,
    pr_auc = pr_auc_value,
    tp = tp,
    fp = fp,
    tn = tn,
    fn = fn
  )

  return(result)
}

# ============================================================
# 4 Modelltraining und Vorhersage
# ============================================================

fit_predict_ranger <- function(features, train_rows, test_rows, ntree_value, importance_mode) {
  train_x <- model_data[train_rows, features, drop = FALSE]
  test_x <- model_data[test_rows, features, drop = FALSE]
  train_y <- model_data$target[train_rows]

  # Fehlende Werte mit Median aus dem jeweiligen Trainingsset ersetzen.
  for (feature in features) {
    median_value <- median(train_x[[feature]], na.rm = TRUE)

    if (!is.finite(median_value)) {
      median_value <- 0
    }

    train_x[[feature]][is.na(train_x[[feature]])] <- median_value
    test_x[[feature]][is.na(test_x[[feature]])] <- median_value
  }

  train_data <- data.frame(target = train_y, train_x, check.names = FALSE)
  mtry_value <- max(1, floor(sqrt(length(features))))

  rf_model <- ranger(
    dependent.variable.name = "target",
    data = train_data,
    probability = TRUE,
    num.trees = ntree_value,
    mtry = mtry_value,
    importance = importance_mode,
    classification = TRUE,
    seed = 42,
    num.threads = num_threads
  )

  prediction <- predict(rf_model, data = test_x)
  probability <- prediction$predictions[, "Niedermoor"]

  return(list(
    probability = probability,
    model = rf_model
  ))
}

evaluate_split <- function(features, train_rows, test_rows, ntree_value) {
  result <- fit_predict_ranger(
    features = features,
    train_rows = train_rows,
    test_rows = test_rows,
    ntree_value = ntree_value,
    importance_mode = "none"
  )

  metrics <- calculate_metrics(
    truth = model_data$target[test_rows],
    probability = result$probability,
    threshold = 0.5
  )

  return(list(
    metrics = metrics,
    probability = result$probability
  ))
}

evaluate_inner_cv <- function(features, inner_folds) {
  probability <- rep(NA_real_, nrow(model_data))

  for (validation_fold in inner_folds) {
    inner_train_rows <- model_data$spatial_fold %in% inner_folds &
      model_data$spatial_fold != validation_fold

    inner_valid_rows <- model_data$spatial_fold == validation_fold

    result <- fit_predict_ranger(
      features = features,
      train_rows = inner_train_rows,
      test_rows = inner_valid_rows,
      ntree_value = ntree_selection,
      importance_mode = "none"
    )

    probability[inner_valid_rows] <- result$probability
  }

  used_rows <- is.finite(probability)

  metrics <- calculate_metrics(
    truth = model_data$target[used_rows],
    probability = probability[used_rows],
    threshold = 0.5
  )

  return(list(
    metrics = metrics,
    probability = probability
  ))
}

# ============================================================
# 5 Nested spatial Forward Selection
# ============================================================

all_selected_features <- data.frame()
selection_history <- data.frame()
outer_test_metrics <- data.frame()
outer_test_predictions <- data.frame()

for (outer_test_fold in fold_ids) {
  cat("\n############################################################\n")
  cat("Outer Test Fold:", outer_test_fold, "\n")
  cat("Feature Selection nur auf den anderen 4 Folds\n")
  cat("############################################################\n")

  inner_folds <- fold_ids[fold_ids != outer_test_fold]
  outer_train_rows <- model_data$spatial_fold %in% inner_folds
  outer_test_rows <- model_data$spatial_fold == outer_test_fold

  selected_features <- character(0)
  remaining_features <- candidate_features
  best_score <- -Inf

  for (step in seq_len(max_selected_features)) {
    cat("\nOuter Fold", outer_test_fold, "- Forward Schritt", step, "\n")
    cat("Bereits ausgewaehlt:", length(selected_features), "\n")
    cat("Noch uebrig:", length(remaining_features), "\n")

    step_results <- data.frame(
      candidate = character(0),
      precision = numeric(0),
      recall = numeric(0),
      f1 = numeric(0),
      roc_auc = numeric(0),
      pr_auc = numeric(0),
      stringsAsFactors = FALSE
    )

    for (i in seq_along(remaining_features)) {
      candidate <- remaining_features[i]
      test_features <- c(selected_features, candidate)

      if (i == 1 || i %% 25 == 0 || i == length(remaining_features)) {
        cat("  Kandidat", i, "von", length(remaining_features), ":", candidate, "\n")
      }

      result <- evaluate_inner_cv(
        features = test_features,
        inner_folds = inner_folds
      )

      metrics <- result$metrics

      step_results <- rbind(
        step_results,
        data.frame(
          candidate = candidate,
          precision = metrics$precision,
          recall = metrics$recall,
          f1 = metrics$f1,
          roc_auc = metrics$roc_auc,
          pr_auc = metrics$pr_auc,
          stringsAsFactors = FALSE
        )
      )
    }

    metric_values <- step_results[[selection_metric]]
    metric_values[is.na(metric_values)] <- -Inf

    best_id <- which.max(metric_values)
    best_candidate <- step_results$candidate[best_id]
    new_score <- metric_values[best_id]

    if (step == 1) {
      improvement <- NA
    } else {
      improvement <- new_score - best_score
    }

    cat("\nBestes Feature:", best_candidate, "\n")
    cat("Innere", selection_metric, ":", round(new_score, 4), "\n")
    cat("Verbesserung:", ifelse(is.na(improvement), "Startwert", round(improvement, 4)), "\n")

    if (step > 1 && improvement < min_improvement) {
      cat("Abbruch: Verbesserung kleiner als", min_improvement, "\n")
      break
    }

    selected_features <- c(selected_features, best_candidate)
    remaining_features <- remaining_features[remaining_features != best_candidate]
    best_score <- new_score

    selection_history <- rbind(
      selection_history,
      data.frame(
        outer_test_fold = outer_test_fold,
        step = step,
        added_feature = best_candidate,
        n_features = length(selected_features),
        precision = step_results$precision[best_id],
        recall = step_results$recall[best_id],
        f1 = step_results$f1[best_id],
        roc_auc = step_results$roc_auc[best_id],
        pr_auc = step_results$pr_auc[best_id],
        stringsAsFactors = FALSE
      )
    )
  }

  selected_table <- data.frame(
    outer_test_fold = outer_test_fold,
    rank = seq_along(selected_features),
    feature = selected_features,
    stringsAsFactors = FALSE
  )

  all_selected_features <- rbind(all_selected_features, selected_table)

  outer_result <- evaluate_split(
    features = selected_features,
    train_rows = outer_train_rows,
    test_rows = outer_test_rows,
    ntree_value = ntree_final
  )

  outer_metrics <- outer_result$metrics
  outer_metrics$outer_test_fold <- outer_test_fold
  outer_metrics$n_features <- length(selected_features)

  outer_test_metrics <- rbind(outer_test_metrics, outer_metrics)

  fold_predictions <- data.frame(
    poly_id = model_data$poly_id[outer_test_rows],
    class_binary = model_data$class_binary[outer_test_rows],
    outer_test_fold = outer_test_fold,
    probability_niedermoor = outer_result$probability,
    stringsAsFactors = FALSE
  )

  outer_test_predictions <- rbind(outer_test_predictions, fold_predictions)

  cat("\nAusgewaehlte Features fuer Outer Fold", outer_test_fold, ":\n")
  print(selected_table)

  cat("\nTestmetriken fuer Outer Fold", outer_test_fold, ":\n")
  print(outer_metrics)
}

# ============================================================
# 6 Feature-Haeufigkeit und Feature-Union
# ============================================================

feature_counts <- table(all_selected_features$feature)

feature_frequency <- data.frame(
  feature = names(feature_counts),
  selected_in_n_outer_folds = as.integer(feature_counts),
  stringsAsFactors = FALSE
)

mean_rank <- aggregate(
  rank ~ feature,
  data = all_selected_features,
  FUN = mean
)

feature_frequency <- merge(feature_frequency, mean_rank, by = "feature", all.x = TRUE)
names(feature_frequency)[names(feature_frequency) == "rank"] <- "mean_selection_rank"

feature_frequency <- feature_frequency[
  order(-feature_frequency$selected_in_n_outer_folds, feature_frequency$mean_selection_rank, feature_frequency$feature),
]

final_features <- as.character(feature_frequency$feature)
final_features <- final_features[final_features %in% names(model_data)]

cat("\nFeature-Haeufigkeit:\n")
print(feature_frequency)

cat("\nFeature-Union fuer finales Modell:\n")
print(final_features)

# ============================================================
# 7 Aeussere Test-Performance und Threshold
# ============================================================

outer_test_predictions$target <- factor(
  outer_test_predictions$class_binary,
  levels = c("Kein_Niedermoor", "Niedermoor")
)

outer_metrics_05 <- calculate_metrics(
  truth = outer_test_predictions$target,
  probability = outer_test_predictions$probability_niedermoor,
  threshold = 0.5
)

cat("\nGesamtperformance der aeusseren Testfolds bei Threshold 0.5:\n")
print(outer_metrics_05)

thresholds <- seq(0.05, 0.95, by = 0.01)
threshold_metrics <- data.frame()

for (threshold in thresholds) {
  metrics <- calculate_metrics(
    truth = outer_test_predictions$target,
    probability = outer_test_predictions$probability_niedermoor,
    threshold = threshold
  )

  threshold_metrics <- rbind(threshold_metrics, metrics)
}

precision_candidates <- threshold_metrics[
  !is.na(threshold_metrics$precision) &
    threshold_metrics$precision >= target_precision,
]

if (nrow(precision_candidates) > 0) {
  score <- precision_candidates$f1
  score[is.na(score)] <- -Inf
  best_threshold_row <- precision_candidates[which.max(score), ]
  cat("\nThreshold wurde mit Precision-Ziel gewaehlt.\n")
} else {
  score <- threshold_metrics$f1
  score[is.na(score)] <- -Inf
  best_threshold_row <- threshold_metrics[which.max(score), ]
  cat("\nKein Threshold erreicht die Ziel-Precision. Bester F1-Threshold wird genutzt.\n")
}

best_threshold <- best_threshold_row$threshold

cat("\nEmpfohlener Threshold:\n")
print(best_threshold_row)

# ============================================================
# 8 Finales Modell mit Feature-Union auf allen Daten trainieren
# ============================================================

final_x <- model_data[, final_features, drop = FALSE]
final_y <- model_data$target

final_medians <- data.frame(
  feature = final_features,
  median_value = NA,
  stringsAsFactors = FALSE
)

for (feature in final_features) {
  median_value <- median(final_x[[feature]], na.rm = TRUE)

  if (!is.finite(median_value)) {
    median_value <- 0
  }

  final_x[[feature]][is.na(final_x[[feature]])] <- median_value
  final_medians$median_value[final_medians$feature == feature] <- median_value
}

final_train_data <- data.frame(target = final_y, final_x, check.names = FALSE)
final_mtry <- max(1, floor(sqrt(length(final_features))))

final_rf <- ranger(
  dependent.variable.name = "target",
  data = final_train_data,
  probability = TRUE,
  num.trees = ntree_final,
  mtry = final_mtry,
  importance = "permutation",
  classification = TRUE,
  seed = 42,
  num.threads = num_threads
)

cat("\nFinales Ranger-Modell mit Feature-Union:\n")
print(final_rf)

importance_table <- data.frame(
  feature = names(final_rf$variable.importance),
  permutation_importance = as.numeric(final_rf$variable.importance),
  stringsAsFactors = FALSE
)

importance_table <- importance_table[
  order(-importance_table$permutation_importance),
]

# ============================================================
# 9 Plots
# ============================================================

selection_plot <- ggplot(selection_history, aes(x = n_features, y = pr_auc, color = factor(outer_test_fold))) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  labs(
    title = "Nested spatial Forward Selection",
    subtitle = "Feature Selection nur innerhalb der Trainingsfolds",
    x = "Anzahl ausgewaehlter Features",
    y = "Innere PR-AUC",
    color = "Aeusserer Testfold"
  ) +
  theme_minimal(base_size = 12)

print(selection_plot)

ggsave(
  filename = file.path(figure_dir, "05_nested_forward_selection_pr_auc.png"),
  plot = selection_plot,
  width = 8,
  height = 5,
  dpi = 220,
  bg = "white"
)

top_frequency <- head(feature_frequency, 30)

frequency_plot <- ggplot(
  top_frequency,
  aes(x = reorder(feature, selected_in_n_outer_folds), y = selected_in_n_outer_folds)
) +
  geom_col(fill = "#2a9d8f") +
  coord_flip() +
  scale_y_continuous(breaks = 0:n_folds, limits = c(0, n_folds)) +
  labs(
    title = "Feature-Stabilitaet in nested spatial CV",
    subtitle = "Wie oft wurde ein Feature in den aeusseren Trainingslaeufen ausgewaehlt?",
    x = NULL,
    y = "Ausgewaehlt in n Outer Folds"
  ) +
  theme_minimal(base_size = 12)

print(frequency_plot)

ggsave(
  filename = file.path(figure_dir, "05_feature_selection_frequency_nested.png"),
  plot = frequency_plot,
  width = 8,
  height = 6,
  dpi = 220,
  bg = "white"
)

threshold_long <- rbind(
  data.frame(threshold = threshold_metrics$threshold, metric = "Precision", value = threshold_metrics$precision),
  data.frame(threshold = threshold_metrics$threshold, metric = "Recall", value = threshold_metrics$recall),
  data.frame(threshold = threshold_metrics$threshold, metric = "F1", value = threshold_metrics$f1)
)

threshold_plot <- ggplot(threshold_long, aes(x = threshold, y = value, color = metric)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = best_threshold, linetype = "dashed", color = "black") +
  scale_color_manual(values = c(
    "Precision" = "#1b9e77",
    "Recall" = "#d95f02",
    "F1" = "#386cb0"
  )) +
  labs(
    title = "Threshold-Auswahl auf aeusseren Testvorhersagen",
    subtitle = paste("Gestrichelte Linie = empfohlener Threshold", round(best_threshold, 2)),
    x = "Threshold",
    y = "Metrikwert",
    color = "Metrik"
  ) +
  theme_minimal(base_size = 12)

print(threshold_plot)

ggsave(
  filename = file.path(figure_dir, "05_threshold_precision_recall_f1_nested.png"),
  plot = threshold_plot,
  width = 8,
  height = 5,
  dpi = 220,
  bg = "white"
)

importance_plot <- ggplot(
  importance_table,
  aes(x = reorder(feature, permutation_importance), y = permutation_importance)
) +
  geom_col(fill = "#2364aa") +
  coord_flip() +
  labs(
    title = "Permutation Importance im finalen Ranger-Modell",
    subtitle = "Finales Modell mit Feature-Union aus nested spatial CV",
    x = NULL,
    y = "Permutation Importance"
  ) +
  theme_minimal(base_size = 12)

print(importance_plot)

ggsave(
  filename = file.path(figure_dir, "05_permutation_importance_final_union.png"),
  plot = importance_plot,
  width = 8,
  height = 7,
  dpi = 220,
  bg = "white"
)

# ============================================================
# 10 Ergebnisse speichern
# ============================================================

final_features_table <- data.frame(
  rank = seq_along(final_features),
  feature = final_features,
  stringsAsFactors = FALSE
)

write.csv(
  all_selected_features,
  file.path(table_dir, "05_nested_selected_features_by_outer_fold.csv"),
  row.names = FALSE
)

write.csv(
  selection_history,
  file.path(table_dir, "05_nested_forward_selection_history.csv"),
  row.names = FALSE
)

write.csv(
  feature_frequency,
  file.path(table_dir, "05_nested_feature_selection_frequency.csv"),
  row.names = FALSE
)

write.csv(
  final_features_table,
  file.path(processed_dir, "selected_features_union_rf.csv"),
  row.names = FALSE
)

write.csv(
  outer_test_metrics,
  file.path(table_dir, "05_nested_outer_test_metrics_threshold_05.csv"),
  row.names = FALSE
)

write.csv(
  outer_test_predictions,
  file.path(table_dir, "05_nested_outer_test_predictions.csv"),
  row.names = FALSE
)

write.csv(
  threshold_metrics,
  file.path(table_dir, "05_nested_threshold_metrics.csv"),
  row.names = FALSE
)

write.csv(
  best_threshold_row,
  file.path(table_dir, "05_nested_best_threshold.csv"),
  row.names = FALSE
)

write.csv(
  importance_table,
  file.path(table_dir, "05_nested_final_union_permutation_importance.csv"),
  row.names = FALSE
)

write.csv(
  final_medians,
  file.path(processed_dir, "selected_feature_medians.csv"),
  row.names = FALSE
)

saveRDS(
  final_rf,
  file.path(model_dir, "05_final_binary_ranger_random_forest.rds")
)

# ============================================================
# 11 Abschluss
# ============================================================

cat("\n============================================================\n")
cat("Skript 05 abgeschlossen.\n")
cat("Nested Feature-Listen:", file.path(table_dir, "05_nested_selected_features_by_outer_fold.csv"), "\n")
cat("Feature-Haeufigkeit:", file.path(table_dir, "05_nested_feature_selection_frequency.csv"), "\n")
cat("Finale Feature-Union:", file.path(processed_dir, "selected_features_union_rf.csv"), "\n")
cat("Permutation Importance:", file.path(table_dir, "05_nested_final_union_permutation_importance.csv"), "\n")
cat("Empfohlener Threshold:", round(best_threshold, 3), "\n")

cat("\nAeusserer Test bei Threshold 0.5:\n")
print(outer_metrics_05)

cat("\nAeusserer Test beim empfohlenen Threshold:\n")
print(best_threshold_row)

cat("\nFinale Feature-Union:\n")
print(final_features_table)

cat("============================================================\n")
