# Load required libraries ----
library(ggpubr)
library(ggplot2)
library(dplyr)
library(readxl)
library(stats)
library(rstatix)
library(tidyr)

# Import datasets ----
GMV_cluster_df <- read_excel(file.choose(), sheet = "corrected")
View(GMV_cluster_df)
str(GMV_cluster_df)

GMV_cluster_df <- GMV_cluster_df %>%
  rename(CODE = Subjects)
GMV_cluster_df <- as.data.frame(GMV_cluster_df)
str(GMV_cluster_df)

VBM_df <- read_excel(file.choose())
VBM_df <- as.data.frame(VBM_df)
View(VBM_df)
str(VBM_df)

VBM_df$SEX <- as.factor(VBM_df$SEX)
VBM_df$GROUP <- as.factor(as.character(VBM_df$GROUP))

VBM_df <- VBM_df %>%
  mutate(GROUP = case_when(
    GROUP %in% c("1.1", "1.2") ~ "HO",
    GROUP == "2" ~ "HY",
    GROUP == "3" ~ "MCI",
    TRUE ~ as.character(GROUP)
  ) %>% factor()) 

VBM_df <- VBM_df %>%
  mutate(CODE = gsub(" \\(T\\d+\\)", "", CODE))

VBM_df <- VBM_df %>%
  select(-`TIV_GM (%)`, -`TIV_WM (%)`, -`TIV_CSF (%)`) %>%
  rename(TIV = `TIV (cm3)`)

str(VBM_df)

# Merge both dataframes
merged_df <- inner_join(GMV_cluster_df, VBM_df, by = "CODE")
View(merged_df)
str(merged_df)

# Filter to keep only HO and MCI 
final_df <- merged_df %>%
  filter(GROUP != "HY")

final_df <- final_df %>%
  distinct()

View(final_df)
str(final_df)

path <- dirname(file.choose())
filename <- "final_df.csv"
full_path <- file.path(path, filename)
write.csv(final_df, full_path, row.names = FALSE)

# Perform normalization to HO controls ----

# Method: Calculate percentage change from control group mean
# Formula: ((individual_value - control_mean) / control_mean) * 100

normalize_to_control_mean <- function(volume, group) {
  # Calculate mean of control group
  control_mean <- mean(volume[group == "HO"])
  
  # Calculate percentage change for all subjects
  pct_change <- ((volume - control_mean) / control_mean) * 100
  
  return(pct_change)
}

# Apply normalization to all brain regions
# Normalize all columns from column 2 to the column before TIV
tiv_position <- which(names(final_df) == "TIV")

# Get all brain region columns (between CODE and TIV)
brain_cols <- names(final_df)[2:(tiv_position - 1)]

for(col in brain_cols) {
  final_df[[paste0(col, "_norm")]] <- normalize_to_control_mean(final_df[[col]], final_df$GROUP)
}

View(final_df)

# Verify the normalization ----

cat("\n=== VERIFICATION ===\n")
cat("Control group mean for normalized hippocampus should be ~0:\n")
cat("Mean:", mean(final_df$L_Amygdala_norm[final_df$GROUP == "HO"]), "\n\n")
cat("Mean:", mean(final_df$L_Amygdala_norm[final_df$GROUP == "MCI"]), "\n\n")

# PREPARE DATA FOR PLOTTING ----
### Calculate group means and standard errors ----

prepare_plot_data <- function(data, normalized_col, group_col = "GROUP") {
  data %>%
    group_by(!!sym(group_col)) %>%
    summarise(
      mean = mean(!!sym(normalized_col)),
      sd = sd(!!sym(normalized_col)),
      n = n(),
      se = sd / sqrt(n),
      .groups = 'drop'
    )
}

# Create summary data for each region
# Get all columns that end with "_norm"
norm_cols <- names(final_df)[grep("_norm$", names(final_df))]

# Apply prepare_plot_data to each normalized column and print summaries
for(col in norm_cols) {
  # Extract region name (remove "_norm" suffix)
  region_name <- sub("_norm$", "", col)
  # Print region name
  cat("\n", region_name, "\n", sep = "")
  # Prepare data and print summary
  region_summary <- prepare_plot_data(final_df, col)
  print(region_summary)
}

# Create an empty list to store summaries
summary_clusters <- list()

# Loop through each normalized column
for(col in norm_cols) {
  # Extract region name (remove "_norm" suffix)
  region_name <- sub("_norm$", "", col)
  
  # Print region name
  cat("\n", region_name, "\n", sep = "")
  
  # Prepare data and print summary
  region_summary <- prepare_plot_data(final_df, col)
  print(region_summary)
  
  # Store the summary
  summary_clusters[[region_name]] <- region_summary
}

mean_values <- list()

for(col in norm_cols) {
  region_name <- sub("_norm$", "", col)
  region_data <- summary_clusters[[region_name]]
  mean_values[[region_name]] <- region_data$mean
}

mean_values <- as.data.frame(mean_values)

mean_df <- t(mean_values)
colnames(mean_df) <- c("HO", "MCI")

# Show result
View(mean_df)

### Reshape data to long format for plotting ----
plot_data_long <- final_df %>%
  select(CODE, GROUP, !ends_with("_norm")) %>%
  pivot_longer(
    cols = c(-CODE, -GROUP, -AGE, -SEX),
    names_to = "region",
    values_to = "normalized_volume",
    values_transform = list(normalized_volume = as.character)
  )

View(plot_data_long)

### Calculate summary statistics for the plot ----
summary_data <- plot_data_long %>%
  mutate(normalized_volume = as.numeric(normalized_volume)) %>%
  group_by(region, GROUP) %>%
  summarise(
    mean = mean(normalized_volume),
    sd = sd(normalized_volume),
    n = n(),
    se = sd / sqrt(n),
    .groups = 'drop'
  )

View(summary_data)

### Extract only patient data for cleaner visualization ----
patient_summary <- as.data.frame(mean_df) %>%
  select(-HO)

patient_summary$region <- rownames(patient_summary)
rownames(patient_summary) <- NULL  # Optional: remove row names

View(patient_summary)
patient_summary <- as.data.frame(patient_summary)
str(patient_summary)

max_region <- patient_summary$region[which.min(patient_summary$MCI)]

# Create color palette: steel grey for all, except a highlight color for max
steel_grey <- "#71797E"  # or use any steel grey you prefer
highlight_color <- "#1E90FF"  # red, or choose your preferred highlight color

color_palette <- setNames(
  ifelse(patient_summary$region == max_region, highlight_color, steel_grey),
  patient_summary$region
)

# Create formatted labels with exactly 1 decimal place
patient_summary$label <- sprintf("%.2f", patient_summary$MCI)

### Plot the actual data ----
plt <- ggbarplot(patient_summary, x = "region", y = "MCI", fill = "region",
                 color = "region",
                 palette = color_palette,
                 order = c("L_SFG_Medial_Orbital", "L_NucleusAccumbens", "L_Amygdala", "R_Hippocampus", "R_Mid_Cingulate", "L_Hippocampus"), 
                 orientation = "horiz",
                 width = 0.3,
                 label = patient_summary$label,
                 lab.pos = "in", 
                 lab.col = "white",
                 lab.vjust = 0.5,
                 lab.hjust = 1.4,
                 lab.size = 7,
                 #lab.nb.digits = 1,
                 legend = "none",
                 xlab = "",
                 ylab = "Mean Volume change from controls (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_reverse() +
  scale_fill_manual(values = color_palette) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 22)
    )

print(plt)

# Export plot ----
ggsave("brain_volume_comparison_all_TIV_age_sex.png", plt, width = 10, height = 14, dpi = 900)
