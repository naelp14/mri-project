library("dplyr")
library("ggplot2")
library("stringr")
library("lmtest")
library("stargazer")
library("purrr")
library("tidyr")
library("dtw")
data <- read.csv("/Users/nathaniel.putera/Projects/UCPH/MRI-Project/ADNIMERGE-Simple.csv")

normalize_minmax <- function(x) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

data <- data %>%
  filter(!is.na(WholeBrain)) %>%
  mutate(
    Ventricles_percent = (Ventricles / ICV),
    Hippocampus_percent = (Hippocampus / ICV),
    WholeBrain_percent = (WholeBrain / ICV),
    Entorhinal_percent = (Entorhinal / ICV),
    Fusiform_percent = (Fusiform / ICV),
    MidTemp_percent = (MidTemp / ICV)
  )

ptid_counts <- table(data$PTID)

ptid_more_than_once <- names(ptid_counts[ptid_counts > 1])
filtered_data <- data %>%
  filter(data$PTID %in% ptid_more_than_once)

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
  select(-sorted) %>%  
  filter(
    !sapply(cdrsb_scores, function(x) all(x == 0))
  )

perform_clustering <- function(summary_data, centers = 4, seed = 123) {
  # Clean the data
  cleaned_summary_table <- summary_data %>%
    drop_na()
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Perform k-means clustering
  kmeans_result <- kmeans(cleaned_summary_table %>% select(-PTID), centers = centers)
  
  # Add cluster results to the summary table
  cleaned_summary_table <- cleaned_summary_table %>%
    mutate(cluster = kmeans_result$cluster)
  
  grouped_data_with_clusters <- grouped_data %>%
    left_join(cleaned_summary_table %>% select(PTID, cluster), by = "PTID")
  
  return(list(kmeans_result = kmeans_result, grouped_data_with_clusters = grouped_data_with_clusters))
}

# try to get fluctuating data
summary_table <- grouped_data %>%
  mutate(
    cdrsb_mean = map_dbl(cdrsb_scores, mean, na.rm = TRUE),         
    cdrsb_first = map_dbl(cdrsb_scores, ~ .x[1]),                   
    cdrsb_last = map_dbl(cdrsb_scores, ~ .x[length(.x)]),           
    cdrsb_sd = map_dbl(cdrsb_scores, sd, na.rm = TRUE),
    cdrsb_scores = cdrsb_scores,
    months = months,
    visits = map_int(cdrsb_scores, length)
  ) %>%
  select(PTID, cdrsb_mean, cdrsb_first, cdrsb_last, cdrsb_sd, cdrsb_scores, months, visits)

threshold_diff <- 1
threshold_sd <- 1
fluctuating_data <- summary_table %>%
  filter(abs(cdrsb_first - cdrsb_last) < threshold_diff & cdrsb_sd > threshold_sd)

plot_data <- fluctuating_data %>%
  unnest(c(cdrsb_scores, months))  # Unnest the list columns

ggplot(plot_data, aes(x = months, y = cdrsb_scores, group = PTID, color = PTID)) +
  geom_line() +
  labs(title = "CDRSB Progression Over Time by Patient",
       x = "Months", y = "CDRSB Score") +
  theme_minimal() +
  theme(legend.position = "none")  # Hide legend for individual patients

# Try to normalize the data again for the time part
data <- data %>%
  mutate(across(where(is.numeric), normalize_minmax))

ptid_counts <- table(data$PTID)

ptid_more_than_once <- names(ptid_counts[ptid_counts > 1])
filtered_data <- data %>%
  filter(data$PTID %in% ptid_more_than_once)

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
  select(-sorted) %>%  
  filter(
    !sapply(cdrsb_scores, function(x) all(x == 0))
  )

summary_table <- grouped_data %>%
  mutate(
    cdrsb_mean = map_dbl(cdrsb_scores, mean, na.rm = TRUE),         
    cdrsb_first = map_dbl(cdrsb_scores, ~ .x[1]),                   
    cdrsb_last = map_dbl(cdrsb_scores, ~ .x[length(.x)]),           
    cdrsb_sd = map_dbl(cdrsb_scores, sd, na.rm = TRUE),
    visits_count = map_int(cdrsb_scores, length),
    months_diff = map_dbl(months, ~ .x[length(.x)] - .x[1])
  ) %>%
  select(PTID, cdrsb_mean, cdrsb_first, cdrsb_last, cdrsb_sd)

grouped_data_with_clusters <- perform_clustering(summary_data = summary_table)

long_data <- grouped_data_with_clusters %>%
  unnest(c(months, cdrsb_scores))

long_data <- long_data %>%
  mutate(cluster = factor(cluster, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("Moderate Progression", 
                                     "Severe Progression", 
                                     "Stable",
                                     "Mild Progression")))

ggplot(long_data, aes(x = months, y = cdrsb_scores, group = PTID, color = factor(cluster))) +
  geom_line(alpha = 0.6) +
  labs(title = "CDRSB Progression by Cluster",
       x = "Months of Visit",
       y = "CDRSB Score",
       color = "Cluster") +
  theme_bw()

# Find out where is the fluctuating data located in the cluster
fluctuating_clusters <- fluctuating_data %>%
  inner_join(long_data %>% select(PTID, cluster), by = "PTID")

fluctuating_clusters %>% select(PTID, cdrsb_first, cdrsb_last, cdrsb_mean, cdrsb_sd, cluster) %>% distinct()

# Try to search for best K
data_numeric <- data %>%
  select(where(is.numeric))

# Step 2: Create a function to calculate the WCSS (within-cluster sum of squares) for different k values
wcss <- function(data, k) {
  kmeans_result <- kmeans(data, centers = k, nstart = 25)  # nstart = 25 runs kmeans with different starting points
  return(kmeans_result$tot.withinss)  # Total WCSS for the given k
}

# Step 3: Run the elbow method for k = 2 to 10 and store the WCSS for each k
wcss_values <- sapply(2:10, function(k) wcss(data_numeric, k))

# Step 4: Plot the WCSS vs. the number of clusters (k)
plot(2:10, wcss_values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (k)",
     ylab = "Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Method for Determining Optimal k")