library("dplyr")
library("ggplot2")
library("tidyr")
library("purrr")
library("dtwclust")

library("xgboost")
library("caTools")
library("caret")

normalize_minmax <- function(x) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

sort_viscode <- function(viscode, months, cdrsb) {
  visit_numeric <- as.numeric(gsub("m", "", gsub("bl", "0", viscode)))
  sorted_indices <- order(visit_numeric)
  list(
    viscode_sorted = viscode[sorted_indices],
    cdrsb_sorted = cdrsb[sorted_indices],
    month_sorted = months[sorted_indices]
  )
}

data <- read.csv("/Users/nathaniel.putera/Projects/UCPH/MRI-Project/ADNIMERGE-Simple.csv")
data <- data %>% filter(!is.na(WholeBrain)) %>%
  mutate(across(where(is.numeric), normalize_minmax))

ptid_counts <- table(data$PTID)

ptid_more_than_once <- names(ptid_counts[ptid_counts > 1])
filtered_data <- data %>%
  filter(data$PTID %in% ptid_more_than_once)

grouped_data <- filtered_data %>%
  drop_na(CDRSB) %>%                
  group_by(PTID) %>%
  summarize(
    code = list(VISCODE),
    months = list(Month),
    cdrsb_scores = list(CDRSB)
  ) %>%
  mutate(
    sorted = pmap(list(code, months, cdrsb_scores), sort_viscode), 
    code = map(sorted, "viscode_sorted"),             
    cdrsb_scores = map(sorted, "cdrsb_sorted"),
    months = map(sorted, "month_sorted")
  ) %>%
  dplyr::select(-sorted) %>%  
  filter(
    !sapply(cdrsb_scores, function(x) all(x == 0))
  )

time_series_data <- grouped_data$cdrsb_scores

dtw_kmeans_result <- tsclust(time_series_data, 
                             type = "partitional",  # Partitional clustering (k-means)
                             k = 4,                 # Number of clusters (adjust k as needed)
                             distance = "dtw",      # Use DTW as distance metric
                             centroid = "pam",      # PAM-based centroid (good for non-Euclidean)
                             seed = 123)            # Seed for reproducibility

cluster_assignments <- dtw_kmeans_result@cluster
dtw_data <- grouped_data
dtw_data$cluster <- as.factor(cluster_assignments)
dtw_data <- dtw_data %>%
  unnest(cols = c(months, cdrsb_scores))

dtw_data <- dtw_data %>%
  mutate(cluster = factor(cluster, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("Mild Progression", 
                                     "Stable", 
                                     "Severe Progression", 
                                     "Moderate Progression")))

ggplot(dtw_data, aes(x = months, y = cdrsb_scores, group = PTID, color = cluster)) +
  geom_line() +
  labs(title = "CDRSB Progression by DTW-Based Clusters",
       x = "Months",
       y = "CDRSB Score") +
  theme_minimal()

filtered_dtw <- dtw_data %>%
  select(PTID, cluster) %>%
  group_by(PTID) %>%
  slice(1)

merged_data <- left_join(filtered_dtw, data, by = "PTID")

# ANOVA cluster and APOE4
anova_apoe4 <- aov(APOE4 ~ cluster, data = merged_data)
summary(anova_apoe4)

# ANOVA cluster and education
anova_edu <- aov(PTEDUCAT ~ cluster, data = merged_data)
summary(anova_edu)

# ANOVA cluster and age
anova_age <- aov(AGE ~ cluster, data = merged_data)
summary(anova_age)

# This to filter to show only 1 data
first_occurrence <- merged_data %>%
  group_by(PTID) %>% 
  dplyr::slice(1) %>%
  ungroup()

first_occurrence <- first_occurrence %>%
  select(c(Month, AGE, PTEDUCAT, APOE4, CDRSB, cluster))

# Try to do XGBoost
# Data Prep
set.seed(123)
sample_split <- sample.split(Y = first_occurrence$cluster, SplitRatio = 0.7)
train_set <- subset(x = first_occurrence, sample_split == TRUE)
test_set <- subset(x = first_occurrence, sample_split == FALSE)

y_train <- as.integer(train_set$cluster) - 1
y_test <- as.integer(test_set$cluster) - 1
X_train <- train_set %>% select(-cluster)
X_test <- test_set %>% select(-cluster)

# Modeling
xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
xgb_params <- list(
  booster = "gbtree",
  eta = 0.01,
  max_depth = 8,
  gamma = 4,
  subsample = 0.75,
  colsample_bytree = 1,
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = length(levels(first_occurrence$cluster))
)

xgb_model <- xgb.train(
  params = xgb_params,
  data = xgb_train,
  nrounds = 5000,
  verbose = 1
)
xgb_model

xgb_preds <- predict(xgb_model, as.matrix(X_test), reshape = TRUE)
xgb_preds <- as.data.frame(xgb_preds)
colnames(xgb_preds) <- levels(first_occurrence$cluster)
xgb_preds

# Accuracy
xgb_preds$PredictedClass <- apply(xgb_preds, 1, function(y) colnames(xgb_preds)[which.max(y)])
xgb_preds$ActualClass <- levels(first_occurrence$cluster)[y_test + 1]
xgb_preds

accuracy <- sum(xgb_preds$PredictedClass == xgb_preds$ActualClass) / nrow(xgb_preds)
accuracy

# Find importance of the variable
importance_matrix <- xgb.importance(model = xgb_model)
print(importance_matrix)
xgb.plot.importance(importance_matrix)