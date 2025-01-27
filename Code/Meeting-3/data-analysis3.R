library("dplyr")
library("ggplot2")
library("stringr")
library("lmtest")
library("stargazer")
library("purrr")
library("tidyr")
library("dtw")
data <- read.csv("/Users/nathaniel.putera/Projects/UCPH/MRI-Project/ADNIMERGE-Simple.csv")
data %>% dplyr::select(PTID)

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

sort_viscode <- function(viscode, months, cdrsb, brainVolumes, hippocampus, ventricles) {
  visit_numeric <- as.numeric(gsub("m", "", gsub("bl", "0", viscode)))
  sorted_indices <- order(visit_numeric)
  list(
    viscode_sorted = viscode[sorted_indices],
    cdrsb_sorted = cdrsb[sorted_indices],
    month_sorted = months[sorted_indices],
    brain_volumes = brainVolumes[sorted_indices],
    hippocampus_volumes = hippocampus[sorted_indices],
    ventricles_volumes = ventricles[sorted_indices]
  )
}

grouped_data <- filtered_data %>%
  drop_na(CDRSB) %>%                
  group_by(PTID) %>%
  summarize(
    code = list(VISCODE),
    months = list(Month),
    cdrsb_scores = list(CDRSB),
    brain_volumes = list(WholeBrain_percent),
    hippocampus_volumes = list(Hippocampus_percent),
    ventricles_volumes = list(Ventricles_percent)
  ) %>%
  mutate(
    sorted = pmap(list(code, months, cdrsb_scores, brain_volumes, hippocampus_volumes, ventricles_volumes), sort_viscode), 
    code = map(sorted, "viscode_sorted"),             
    cdrsb_scores = map(sorted, "cdrsb_sorted"),
    months = map(sorted, "month_sorted"),
    brain_volumes = map(sorted, "brain_volumes"),
    hippocampus_volumes = map(sorted, "hippocampus_volumes"),
    ventricles_volumes = map(sorted, "ventricles_volumes"),
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
  
  return(grouped_data_with_clusters)
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

# Clustering with only first and last cdrsb
summary_table <- grouped_data %>%
  mutate(    
    cdrsb_first = map_dbl(cdrsb_scores, ~ .x[1]),                   
    cdrsb_last = map_dbl(cdrsb_scores, ~ .x[length(.x)])        
  ) %>%
  select(PTID, cdrsb_first, cdrsb_last)

grouped_data_with_clusters <- perform_clustering(summary_data = summary_table)

long_data <- grouped_data_with_clusters %>%
  unnest(c(months, cdrsb_scores, brain_volumes, hippocampus_volumes, ventricles_volumes))

long_data <- long_data %>%
  mutate(cluster = factor(cluster, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("Mild Progression", 
                                     "Moderate Progression", 
                                     "Severe Progression", 
                                     "Stable")))

ggplot(long_data, aes(x = months, y = cdrsb_scores, group = PTID, color = factor(cluster))) +
  geom_line(alpha = 0.6) +
  labs(title = "CDRSB Progression by Cluster",
       x = "Months of Visit",
       y = "CDRSB Score",
       color = "Cluster") +
  theme_bw()

# Clustering with only mean and sd
summary_table <- grouped_data %>%
  mutate(
    cdrsb_mean = map_dbl(cdrsb_scores, mean, na.rm = TRUE),         
    cdrsb_sd = map_dbl(cdrsb_scores, sd, na.rm = TRUE)              
  ) %>%
  select(PTID, cdrsb_mean, cdrsb_sd)

grouped_data_with_clusters <- perform_clustering(summary_data = summary_table)

long_data <- grouped_data_with_clusters %>%
  unnest(c(months, cdrsb_scores, brain_volumes, hippocampus_volumes, ventricles_volumes))

long_data <- long_data %>%
  mutate(cluster = factor(cluster, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("Mild Progression", 
                                     "Moderate Progression", 
                                     "Severe Progression", 
                                     "Stable")))

ggplot(long_data, aes(x = months, y = cdrsb_scores, group = PTID, color = factor(cluster))) +
  geom_line(alpha = 0.6) +
  labs(title = "CDRSB Progression by Cluster",
       x = "Months of Visit",
       y = "CDRSB Score",
       color = "Cluster") +
  theme_bw()

# Cluster with first, last, mean, and sd of cdrsb
summary_table <- grouped_data %>%
  mutate(
    cdrsb_mean = map_dbl(cdrsb_scores, mean, na.rm = TRUE),         
    cdrsb_first = map_dbl(cdrsb_scores, ~ .x[1]),                   
    cdrsb_last = map_dbl(cdrsb_scores, ~ .x[length(.x)]),           
    cdrsb_sd = map_dbl(cdrsb_scores, sd, na.rm = TRUE)              
  ) %>%
  select(PTID, cdrsb_mean, cdrsb_first, cdrsb_last, cdrsb_sd)

grouped_data_with_clusters <- perform_clustering(summary_data = summary_table)

long_data <- grouped_data_with_clusters %>%
  unnest(c(months, cdrsb_scores, brain_volumes, hippocampus_volumes, ventricles_volumes))

long_data <- long_data %>%
  mutate(cluster = factor(cluster, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("Mild Progression", 
                                     "Moderate Progression", 
                                     "Severe Progression", 
                                     "Stable")))

ggplot(long_data, aes(x = months, y = cdrsb_scores, group = PTID, color = factor(cluster))) +
  geom_line(alpha = 0.6) +
  labs(title = "CDRSB Progression by Cluster",
       x = "Months of Visit",
       y = "CDRSB Score",
       color = "Cluster") +
  theme_bw()

# Cluster with first, last, mean, sd, number of visit, month diff
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
  unnest(c(months, cdrsb_scores, brain_volumes, hippocampus_volumes, ventricles_volumes))

long_data <- long_data %>%
  mutate(cluster = factor(cluster, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("Mild Progression", 
                                     "Moderate Progression", 
                                     "Severe Progression", 
                                     "Stable")))

ggplot(long_data, aes(x = months, y = cdrsb_scores, group = PTID, color = factor(cluster))) +
  geom_line(alpha = 0.6) +
  labs(title = "CDRSB Progression by Cluster",
       x = "Months of Visit",
       y = "CDRSB Score",
       color = "Cluster") +
  theme_bw()

# Find out where is the fluctuating data located in the cluster
fluctuating_clusters <- fluctuating_data %>%
  inner_join(long_data %>% select(PTID, cluster), by = "PTID")  # Merge by PTID to fetch clusters

fluctuating_summary <- fluctuating_clusters %>%
  group_by(cluster) %>%
  summarize(fluctuating_count = n())  # Count how many fluctuating patients are in each cluster

print(fluctuating_summary)

# DTW 
cdrsb_scores_list <- grouped_data$cdrsb_scores

# Now calculate pairwise DTW distances between each pair of patients
dtw_distances <- matrix(0, nrow = length(cdrsb_scores_list), ncol = length(cdrsb_scores_list))

for (i in 1:length(cdrsb_scores_list)) {
  for (j in 1:length(cdrsb_scores_list)) {
    alignment <- dtw(cdrsb_scores_list[[i]], cdrsb_scores_list[[j]])
    dtw_distances[i, j] <- alignment$distance
  }
}

hc <- hclust(as.dist(dtw_distances), method = "ward.D2")
plot(hc, main = "Hierarchical Clustering using DTW Distances")

num_clusters <- 4

# Cut the dendrogram to form 'num_clusters' clusters
clusters <- cutree(hc, k = num_clusters)

# Create a data frame with PTID and their corresponding cluster
ptid_cluster <- data.frame(
  PTID = grouped_data$PTID,
  Cluster = clusters
)

merged_data <- grouped_data %>%
  left_join(ptid_cluster, by = "PTID")

ggplot(merged_data, aes(x = months, y = cdrsb_scores, group = PTID, color = as.factor(Cluster))) +
  geom_line() +
  labs(title = "CDRSB Progression by Cluster", x = "Months", y = "CDRSB Score", color = "Cluster") +
  theme_minimal()
