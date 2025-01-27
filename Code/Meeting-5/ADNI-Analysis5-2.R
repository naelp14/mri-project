library("dplyr")
library("ggplot2")
library("tidyr")
library("purrr")
library("dtwclust")

data <- read.csv("/Users/nathaniel.putera/Projects/UCPH/MRI-Project/ADNIMERGE-Simple.csv")

data <- data %>%
  filter(!is.na(WholeBrain))

ptid_counts <- table(data$PTID)

ptid_more_than_once <- names(ptid_counts[ptid_counts > 1])
filtered_data <- data %>%
  filter(data$PTID %in% ptid_more_than_once)

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
grouped_data$cluster <- as.factor(cluster_assignments)
grouped_data <- grouped_data %>%
  mutate(cdrsb_mean = map_dbl(cdrsb_scores, ~ mean(unlist(.))))

grouped_data <- grouped_data %>%
  mutate(cluster = factor(cluster, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("Mild Progression", 
                                     "Stable", 
                                     "Severe Progression", 
                                     "Moderate Progression")))

merged_data <- left_join(grouped_data, data, by = "PTID")

contingency_table <- table(merged_data$cluster, merged_data$APOE4)
# Perform the Chi-square test
chi_square_test <- chisq.test(contingency_table)

# Output the result
print(chi_square_test)