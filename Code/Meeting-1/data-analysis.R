library("dplyr")
library("ggplot2")
library("stringr")
data <- read.csv("/Users/nathaniel.putera/Projects/UCPH/MRI-Project/ADNIMERGE-Simple.csv")

unique_data <- data[!duplicated(data$PTID), ]

# Gender Distribution
gender_data <- unique_data[c("PTID", "PTGENDER")]

gender_distribution <- table(unique_data$PTGENDER)

ggplot(data = as.data.frame(gender_distribution), aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  labs(title = "Gender Distribution by Unique Patient", x = "Gender", y = "Count") +
  theme_bw()

# Age Distribution
age_distribution_filtered <- unique_data %>%
  select(AGE, DX) %>%
  filter(!is.na(AGE) & !is.na(DX)) %>%
  filter(DX != "")

age_bins <- c(50, 60, 70, 80, 90, 100)
age_labels <- c('50-59', '60-69', '70-79', '80-89', '90-99')

age_distribution_filtered$AgeRange <- cut(age_distribution_filtered$AGE, breaks = age_bins, labels = age_labels, right = FALSE)

age_range_counts <- age_distribution_filtered %>%
  group_by(AgeRange) %>%
  summarise(Count = n())
ggplot(age_range_counts, aes(x = AgeRange, y = Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Patients in Each Age Range",
       x = "Age Range",
       y = "Count") +
  theme_bw()

# Age distribution with diagnosed
ggplot(age_distribution_filtered, aes(x = AgeRange, fill = DX)) +
  geom_bar(position = "dodge") +
  labs(title = "Distribution of Diagnosis by Age",
       x = "Age Range",
       y = "Count",
       fill = "Diagnosis") +
  theme_bw()

# Marital status
total_marital_data <- unique_data %>%
  select(PTID, PTMARRY) %>%
  filter(PTMARRY != "" & PTMARRY != "Unknown")

total_marital_distribution <- table(total_marital_data$PTMARRY)

ad_patients <- unique_data %>%
  filter(DX == "AD")
ad_marital_data <- ad_patients %>%
  select(PTID, PTMARRY)
ad_marital_distribution <- table(ad_marital_data$PTMARRY)
combined_distribution <- rbind(total_marital_distribution, ad_marital_distribution)
rownames(combined_distribution) <- c("Total", "Alzheimer's")

combined_distribution_df <- as.data.frame(t(combined_distribution))

print(combined_distribution_df)

# Race
race_data <- unique_data[c("PTID", "PTRACCAT")]
total_race_distribution <- table(race_data$PTRACCAT)

ad_patients <- unique_data %>%
  filter(DX == "AD" & !is.na(PTRACCAT))
ad_race_distribution <- table(factor(ad_patients$PTRACCAT, levels = names(total_race_distribution)))

combined_race_distribution <- data.frame(
  Race = names(total_race_distribution),
  Total = as.vector(total_race_distribution),
  Alzheimer = as.vector(ad_race_distribution)
)

combined_race_distribution[is.na(combined_race_distribution)] <- 0
combined_race_distribution

# APOE4
apoe4_data <- unique_data[c("PTID", "APOE4")]
apoe4_distribution <- table(apoe4_data$APOE4)

ad_apoe4_data <- ad_patients %>%
  select(PTID, APOE4) %>%
  filter(!is.na(APOE4))
ad_apoe4_distribition <- table(ad_apoe4_data$APOE4)
combined_distribution <- rbind(apoe4_distribution, ad_apoe4_distribition)
rownames(combined_distribution) <- c("Total", "Alzheimer's")
combined_distribution_df <- as.data.frame(t(combined_distribution))

print(combined_distribution_df)
unique(data$VISCODE)

# Brain volume
filtered_volume_data <- data %>%
  filter(!is.na(VISCODE) & !is.na(WholeBrain) & !is.na(PTID))

filtered_volume_data$VISCODE <- factor(filtered_volume_data$VISCODE, levels = unique(filtered_volume_data$VISCODE))

ggplot(filtered_volume_data, aes(x = VISCODE, y = WholeBrain)) +
  geom_boxplot() +
  labs(title = "Box Plot of Whole Brain Volume by Visit",
       x = "Visit Code (VISCODE)",
       y = "Whole Brain Volume") +
  theme_bw()

filtered_volume_data$VISCODE_num <- ifelse(filtered_volume_data$VISCODE == "bl", -1, as.numeric(str_extract(filtered_volume_data$VISCODE, "\\d+")))
unique_viscodes <- filtered_volume_data %>%
  distinct(VISCODE, VISCODE_num) %>%
  arrange(VISCODE_num) %>%
  pull(VISCODE)

# Convert VISCODE to a factor with levels ordered accordingly
filtered_volume_data$VISCODE <- factor(filtered_volume_data$VISCODE, levels = unique_viscodes)

# Group by VISCODE and calculate the median of WHOLEBRAIN for each group
median_volume_by_viscode <- filtered_volume_data %>%
  group_by(VISCODE) %>%
  summarise(median_volume = median(WholeBrain, na.rm = TRUE))

ggplot(median_volume_by_viscode, aes(x = VISCODE, y = median_volume)) +
  geom_line(group = 1) + # group = 1 ensures it is treated as a continuous line
  geom_point() +
  labs(title = "Median Whole Brain Volume by Visit Code",
       x = "Visit Code (VISCODE)",
       y = "Median Whole Brain Volume") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Group by how many times a person visit
filtered_data <- data %>%
  filter(!is.na(VISCODE) & !is.na(PTID))

# Extract numeric part of VISCODE and handle 'bl' as the first level for sorting purposes
filtered_data$VISCODE_num <- ifelse(filtered_data$VISCODE == "bl", 0, as.numeric(str_extract(filtered_data$VISCODE, "\\d+")))

# Create a visit number column that increments for each subsequent visit of the same patient
filtered_data <- filtered_data %>%
  group_by(PTID) %>%
  mutate(VisitNumber = row_number())

# Drop the VISCODE_num column as it is no longer needed
filtered_data <- filtered_data %>%
  select(-VISCODE_num)

median_volume_by_visit_number <- filtered_data %>%
  group_by(VisitNumber) %>%
  summarise(median_volume = median(WholeBrain, na.rm = TRUE))

# Create the plot for median whole brain volume by visit number
ggplot(median_volume_by_visit_number, aes(x = VisitNumber, y = median_volume)) +
  geom_line(group = 1) + # group = 1 ensures it is treated as a continuous line
  geom_point() +
  labs(title = "Median Whole Brain Volume by Visit Number",
       x = "Visit Number",
       y = "Median Whole Brain Volume") +
  theme_minimal()

# CDRSB
median_cdrsb_by_visit_number <- filtered_data %>%
  group_by(VisitNumber) %>%
  summarise(median_cdrsb = median(CDRSB, na.rm = TRUE))

# Print the median CDRSB values by visit number
print(median_cdrsb_by_visit_number)

# Create the plot for median CDRSB by visit number
ggplot(median_cdrsb_by_visit_number, aes(x = VisitNumber, y = median_cdrsb)) +
  geom_line(group = 1) +
  geom_point() +
  labs(title = "Median CDRSB by Visit Number",
       x = "Visit Number",
       y = "Median CDRSB") +
  theme_bw()

median_cdrsb_by_viscode <- filtered_data %>%
  group_by(VISCODE) %>%
  summarise(median_cdrsb = median(CDRSB, na.rm = TRUE))

ggplot(median_cdrsb_by_viscode, aes(x = VISCODE, y = median_cdrsb)) +
  geom_line(group = 1) +
  geom_point() +
  labs(title = "Median CDRSB by VISCODE",
       x = "VISCODE",
       y = "Median CDRSB") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# CDRSB Value
filtered_data <- unique_data %>%
  select(DX, CDRSB) %>%
  filter(!is.na(DX) & !is.na(CDRSB) & DX != "")

# Plotting boxplot or violin plot
ggplot(filtered_data, aes(x = DX, y = CDRSB, fill = DX)) +
  geom_boxplot() + 
  labs(title = "Relationship Between CDRSB and Diagnosis",
       x = "Diagnosis",
       y = "CDRSB") +
  theme_bw()

# Age to CDRSB
age_bins <- c(50, 60, 70, 80, 90, 100)
age_labels <- c('50-59', '60-69', '70-79', '80-89', '90-99')

# Create age range column
data <- data %>%
  mutate(AgeRange = cut(AGE, breaks = age_bins, labels = age_labels, right = FALSE))

summary_stats <- data %>%
  filter(AgeRange != "NA") %>% 
  group_by(AgeRange) %>%
  summarise(
    Mean_CDRSB = mean(CDRSB, na.rm = TRUE),
    Median_CDRSB = median(CDRSB, na.rm = TRUE),
    SD_CDRSB = sd(CDRSB, na.rm = TRUE)
  )

ggplot(data, aes(x = AgeRange, y = CDRSB, fill = AgeRange)) +
  geom_boxplot() +
  labs(title = "CDRSB Values Across Age Ranges",
       x = "Age Range",
       y = "CDRSB") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

correlation <- cor(data$AGE, data$CDRSB, use = "complete.obs")

# Percentage brain volume
# Calculate percentage change in WholeBrain volume for each visit compared to the previous visit
percentage_data <- data %>%
  group_by(PTID) %>%
  mutate(
    PreviousWholeBrain = lag(WholeBrain), # Get the WholeBrain volume of the previous visit
    PercentageChange = ((WholeBrain - PreviousWholeBrain) / PreviousWholeBrain) * 100
  ) %>%
  ungroup()

percentage_data$PercentageChange[is.na(percentage_data$PercentageChange)] <- 0

summary_counts <- percentage_data %>%
  filter(!is.na(PercentageChange) & PercentageChange != 0) %>%
  mutate(ChangeType = ifelse(PercentageChange < 0, "Negative", "Positive")) %>%
  group_by(ChangeType) %>%
  summarise(Count = n())

print(summary_counts)

# apoe4 with cdrsb
apoe4_cdrsb_data <- data %>%
  select(PTID, APOE4, CDRSB) %>%
  filter(!is.na(APOE4) & !is.na(CDRSB) & APOE4 != 'NA')

# Summary statistics for CDRSB by APOE4 alleles
summary_stats <- apoe4_cdrsb_data %>%
  group_by(APOE4) %>%
  summarise(
    Mean_CDRSB = mean(CDRSB, na.rm = TRUE),
    Median_CDRSB = median(CDRSB, na.rm = TRUE),
    SD_CDRSB = sd(CDRSB, na.rm = TRUE),
    Count = n()
  )
# Violin plot for a detailed distribution view
ggplot(apoe4_cdrsb_data, aes(x = as.factor(APOE4), y = CDRSB, fill = as.factor(APOE4))) +
  geom_violin(trim = FALSE) +
  labs(title = "Distribution of CDRSB by APOE4 Alleles",
       x = "Number of APOE4 Alleles",
       y = "CDRSB",
       fill = "APOE4 Alleles") +
  theme_bw()
