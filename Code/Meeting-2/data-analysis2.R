library("dplyr")
library("ggplot2")
library("stringr")
library("lmtest")
library("stargazer")
library("purrr")
library("tidyr")
data <- read.csv("/Users/nathaniel.putera/Projects/UCPH/MRI-Project/ADNIMERGE-Simple.csv")

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

# Plot for WholeBrain volume
ggplot(filtered_data, aes(x = Month, y = WholeBrain_percent, group = PTID)) +
  geom_line(color = "black", alpha = 0.5) +
  geom_point(color = "black", alpha = 0.7) + 
  geom_smooth(method = "lm", se = FALSE, aes(group = 1, color = "Trend Line")) +
  scale_color_manual(values = rainbow(length(unique(filtered_data$PTID)))) + 
  labs(title = "Whole Brain Volume Across Visits",
       x = "Month of Visit",
       y = "Whole Brain Volume (Percentage)") +
  theme_bw() 

over_100 <- filtered_data %>%
  filter(filtered_data$WholeBrain_percent > 100)

below_50 <- filtered_data %>%
  filter(filtered_data$WholeBrain_percent < 50)

diagnosis_summary <- below_50 %>%
  group_by(DX) %>%
  summarize(count = n(), .groups = 'drop')
diagnosis_summary

# Plot for Hippocampus
ggplot(filtered_data, aes(x = Month, y = Hippocampus_percent, group = PTID)) +
  geom_line() +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, aes(group = 1, color = "Trend Line")) +
  scale_color_manual(values = rainbow(length(unique(filtered_data$PTID)))) + 
  labs(title = "Hippocampus Volume Across Visits",
       x = "Month of Visit",
       y = "Hippocampus Volume (Percentage)") +
  theme_bw()

# Plot for Ventricles
ggplot(filtered_data, aes(x = Month, y = Ventricles_percent, group = PTID)) +
  geom_line() +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, aes(group = 1, color = "Trend Line")) +
  scale_color_manual(values = rainbow(length(unique(filtered_data$PTID)))) + 
  labs(title = "Ventricles Volume Across Visits",
       x = "Month of Visit",
       y = "Ventricles Volume (Percentage)") +
  theme_bw() 

# Plot for Entorhinal
ggplot(filtered_data, aes(x = Month, y = Entorhinal_percent, group = PTID)) +
  geom_line() +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, aes(group = 1, color = "Trend Line")) +
  scale_color_manual(values = rainbow(length(unique(filtered_data$PTID)))) + 
  labs(title = "Entorhinal Volume Across Visits",
       x = "Month of Visit",
       y = "Entorhinal Volume (Percentage)") +
  theme_bw()

# Plot for Fusiform
ggplot(filtered_data, aes(x = Month, y = Fusiform_percent, group = PTID)) +
  geom_line() +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, aes(group = 1, color = "Trend Line")) +
  scale_color_manual(values = rainbow(length(unique(filtered_data$PTID)))) + 
  labs(title = "Fusiform Volume Across Visits",
       x = "Month of Visit",
       y = "Fusiform Volume (Percentage)") +
  theme_bw()

# Plot for Midterm
ggplot(filtered_data, aes(x = Month, y = MidTemp_percent, group = PTID)) +
  geom_line() +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, aes(group = 1, color = "Trend Line")) +
  scale_color_manual(values = rainbow(length(unique(filtered_data$PTID)))) + 
  labs(title = "Midterm Volume Across Visits",
       x = "Month of Visit",
       y = "Midterm Volume (Percentage)") +
  theme_bw()

# Linear regression
model <- lm(WholeBrain ~ AGE + PTGENDER + PTEDUCAT + APOE4, data = data)

stargazer(model, type = "text",
          title = "Regression Results",
          covariate.labels = c("Age", "Gender", "Education", "APOE4"),
          dep.var.labels = "Brain Size")

# Part 3
ggplot(data, aes(x = Month, y = CDRSB, group = PTID, color = PTID)) +
  geom_line(alpha = 0.6) +  # Connect points with lines for each patient
  geom_point(alpha = 0.6) + # Add points for each visit
  labs(title = "Progression of Cognitive Decline Over Time (CDRSB)",
       x = "Month of Visit",
       y = "CDRSB Score") +
  theme_minimal() +
  theme(legend.position = "none") 

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

summary_table <- grouped_data %>%
  mutate(
    cdrsb_mean = map_dbl(cdrsb_scores, mean, na.rm = TRUE),         
    cdrsb_first = map_dbl(cdrsb_scores, ~ .x[1]),                   
    cdrsb_last = map_dbl(cdrsb_scores, ~ .x[length(.x)]),           
    cdrsb_sd = map_dbl(cdrsb_scores, sd, na.rm = TRUE)              
  ) %>%
  select(PTID, cdrsb_mean, cdrsb_first, cdrsb_last, cdrsb_sd)

cleaned_summary_table <- summary_table %>%
  drop_na()

set.seed(123)
kmeans_result <- kmeans(cleaned_summary_table %>% select(-PTID), centers = 4)

cleaned_summary_table <- cleaned_summary_table %>%
  mutate(cluster = kmeans_result$cluster)

grouped_data_with_clusters <- grouped_data %>%
  left_join(cleaned_summary_table %>% select(PTID, cluster), by = "PTID")

long_data <- grouped_data_with_clusters %>%
  unnest(c(months, cdrsb_scores, brain_volumes, hippocampus_volumes, ventricles_volumes))

long_data <- long_data %>%
  mutate(cluster = factor(cluster, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("Mild Progression", 
                                     "Moderate Progression", 
                                     "Severe Progression", 
                                     "Stable")))

# Plot the CDRSB progression for each cluster
ggplot(long_data, aes(x = months, y = cdrsb_scores, group = PTID, color = factor(cluster))) +
  geom_line(alpha = 0.6) +
  labs(title = "CDRSB Progression by Cluster",
       x = "Months of Visit",
       y = "CDRSB Score",
       color = "Cluster") +
  theme_bw()

summary_data <- long_data %>%
  group_by(cluster, months) %>%
  summarize(
    mean_cdrsb = mean(cdrsb_scores, na.rm = TRUE),
    Q1 = quantile(cdrsb_scores, 0.25, na.rm = TRUE),  # 1st quartile (25th percentile)
    Q3 = quantile(cdrsb_scores, 0.75, na.rm = TRUE)   # 3rd quartile (75th percentile)
  )

# Plot the mean CDRSB progression with IQR bands for each cluster
ggplot(summary_data, aes(x = months, y = mean_cdrsb, color = factor(cluster))) +
  geom_line(size = 1.2) +  # Line plot for the mean
  labs(title = "Mean CDRSB Progression by Cluster with IQR",
       x = "Months of Visit",
       y = "Mean CDRSB Score",
       color = "Cluster",
       fill = "Cluster") +
  theme_bw() +
  theme(legend.position = "right")

ggplot(summary_data, aes(x = months, y = mean_cdrsb, color = factor(cluster))) +
  geom_line(size = 1.2) +  # Lines for the mean CDRSB scores of each cluster
  geom_errorbar(aes(ymin = mean_cdrsb - Q1, ymax = mean_cdrsb + Q3), width = 5, size = 1) +  # Error bars
  geom_point(size = 3) +  # Points to emphasize the mean CDRSB score at each time point
  labs(title = "Mean CDRSB Progression by Cluster with Error Bars",
       x = "Months of Visit",
       y = "Mean CDRSB Score",
       color = "Cluster") +
  theme_bw() + 
  theme(legend.position = "right")

# plot mean of each clusters with error bar
summary_data_filtered <- summary_data %>%
  filter(months %% 50 == 0)  # Filter for months divisible by 50

# Plot with error bars using Q1 and Q3
ggplot(summary_data_filtered, aes(x = months, y = mean_cdrsb, color = factor(cluster))) +
  geom_line() +
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 5) +  # Use Q1 and Q3 for error bars
  labs(title = "Mean CDRSB Progression by Cluster with Error Bars Every 50 Months",
       x = "Months of Visit", y = "Mean CDRSB Score")

# Summary of each clusters
summary_by_patient <- long_data %>%
  group_by(PTID, cluster) %>%
  summarize(
    first_cdrsb = first(cdrsb_scores),   # First visit CDRSB
    last_cdrsb = last(cdrsb_scores),     # Last visit CDRSB
    mean_cdrsb = mean(cdrsb_scores, na.rm = TRUE)  # Mean CDRSB across all visits
  )

# Now calculate the mean and standard deviation for each cluster
summary_by_cluster <- summary_by_patient %>%
  group_by(cluster) %>%
  summarize(
    mean_start_cdrsb = mean(first_cdrsb, na.rm = TRUE),  # Mean first visit CDRSB
    mean_last_cdrsb = mean(last_cdrsb, na.rm = TRUE),    # Mean last visit CDRSB
    mean_of_means_cdrsb = mean(mean_cdrsb, na.rm = TRUE), # Mean of all patient means within the cluster
  )

# Whole brain volumes for each clusters
long_data_grouped <- long_data %>%
  group_by(cluster, PTID) %>%
  arrange(months)

# Brain volumes each clusters
ggplot(long_data_grouped, aes(x = months, y = brain_volumes, group = PTID, color = PTID)) +
  geom_line() +
  labs(title = "Brain Volume Change Over Time by Cluster",
       x = "Months", y = "Brain Volume") +
  theme_minimal() +
  theme(legend.position = "none") +  # Hide legend for individual patients
  facet_wrap(~ cluster, scales = "free_y")  # Create one plot per cluster with free y-axis scaling

# Hippocampus each clusters
hippocampus_grouped <- long_data %>%
  filter(!is.na(hippocampus_volumes)) %>%
  group_by(cluster, PTID) %>%
  arrange(months) %>%
  filter(n() > 1)

ggplot(hippocampus_grouped, aes(x = months, y = ventricles_volumes, group = PTID, color = PTID)) +
  geom_line() +
  labs(title = "Hippocampus Volume Change Over Time by Cluster",
       x = "Months", y = "Hippocampus Volume") +
  theme_minimal() +
  theme(legend.position = "none") +  # Hide legend for individual patients
  facet_wrap(~ cluster, scales = "free_y")  # Create one plot per cluster with free y-axis scaling

# Ventricles each clusters
ventricles_grouped <- long_data %>%
  filter(!is.na(ventricles_volumes)) %>%
  group_by(cluster, PTID) %>%
  arrange(months) %>%
  filter(n() > 1)

ggplot(long_data_grouped, aes(x = months, y = ventricles_volumes, group = PTID, color = PTID)) +
  geom_line() +
  labs(title = "Ventricles Volume Change Over Time by Cluster",
       x = "Months", y = "Ventricles Volume") +
  theme_minimal() +
  theme(legend.position = "none") +  # Hide legend for individual patients
  facet_wrap(~ cluster, scales = "free_y")  # Create one plot per cluster with free y-axis scaling

# distribution each cluster
ptid_cluster_dict <- long_data %>%
  distinct(PTID, cluster) 

merged_data <- data %>%
  inner_join(ptid_cluster_dict, by = "PTID") %>%
  distinct(PTID, APOE4, DX, cluster)

# APOE4
apoe4_summary <- merged_data %>%
  group_by(cluster, APOE4) %>%
  summarize(count = n())  # Count how many patients have each APOE4 status in each cluster

apoe4_summary_table <- apoe4_summary %>%
  group_by(cluster) %>%
  summarize(
    apoe4_0 = sum(count[APOE4 == 0], na.rm = TRUE),  # Summing up APOE4=0 count
    apoe4_1 = sum(count[APOE4 == 1], na.rm = TRUE),  # Summing up APOE4=1 count
    apoe4_2 = sum(count[APOE4 == 2], na.rm = TRUE),  # Summing up APOE4=2 count
    apoe4_na = sum(count[is.na(APOE4)], na.rm = TRUE)  # Summing up APOE4 missing values (NA)
  )

print(apoe4_summary_table)

# Diagnosis
dx_summary <- merged_data %>%
  group_by(cluster, DX) %>%
  summarize(count = n())  # Count how many patients have each APOE4 status in each cluster

dx_summary_table <- dx_summary %>%
  group_by(cluster) %>%
  summarize(
    ad = sum(count[DX == "AD"], na.rm = TRUE), 
    cn = sum(count[DX == "CN"], na.rm = TRUE),  
    mci = sum(count[DX == "MCI"], na.rm = TRUE), 
    na = sum(count[is.na(DX)], na.rm = TRUE)  
  )

print(dx_summary_table)
