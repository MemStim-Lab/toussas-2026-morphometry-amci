# Load required libraries
library(dplyr)
library(ggpubr)
library(ggplot2)
library(readxl)

# Load the GMV cluster dataframe

df <- read.csv(file.choose(), sep = ",", header = T)
df <- df %>%
  rename(names = CODE)
View(df)

# Load the NPSY tests dataframe

df2 <- read.csv(file.choose(), sep = ",", header = T)
View(df2)

combined_df <- inner_join(df, df2, by = "names")
View(combined_df)

combined_df <- combined_df %>% distinct(names, .keep_all = TRUE)
View(combined_df)

# Create bilateral hippocampus variable
combined_df <- combined_df %>%
  mutate(B_Hippocampus = L_Hippocampus + R_Hippocampus)
View(combined_df)

# Filter data: keep only MCI patients and select relevant columns
mci_data <- combined_df %>%
  select(names, L_Hippocampus, R_Hippocampus, B_Hippocampus,
         LMSB.immediate, LMSB.delayed, X3.Object.3.Places, ROCF.copy, ROCF.recall, TMT.A,
         Age, Gender, Years.of.education)
View(mci_data)
str(mci_data)

mci_data <- mci_data %>%
  rename(Sex = Gender) %>%
  mutate(Sex = case_when(
    Sex == "F" ~ "Female",
    Sex == "M" ~ "Male",
    TRUE ~ Sex
  ))
View(mci_data)

mci_data <- as.data.frame(mci_data)
mci_data <- mci_data %>%
  mutate(
    LMSB.immediate = as.numeric(LMSB.immediate),
    LMSB.delayed = as.numeric(LMSB.delayed),
    X3.Object.3.Places = as.numeric(X3.Object.3.Places),
    ROCF.copy = as.numeric(ROCF.copy),
    ROCF.recall = as.numeric(ROCF.recall),
    TMT.A = as.numeric(TMT.A)
  )


# Display filtered data summary
cat("\nFiltered MCI data summary:\n")
summary(mci_data)

# Function to create correlation plot with statistics
create_correlation_plot <- function(data, x_var, y_var, title, 
                                    manual_r = NULL, manual_p = NULL,
                                    cor_method = "pearson",  # "pearson" or "kendall"
                                    stats_position = "top-left") {  # position options
  
  # Calculate correlation and p-value (only if manual values not provided)
  cor_test <- cor.test(data[[x_var]], data[[y_var]], method = cor_method)
  
  # Determine correlation symbol based on method
  cor_symbol <- if(cor_method == "kendall") "τ" else "r"
  
  # Use manual values if provided, otherwise use calculated values
  if (!is.null(manual_r) && !is.null(manual_p)) {
    r_text <- paste0(cor_symbol, " = ", manual_r)
    p_text <- paste0("p = ", manual_p)
  } else {
    # Format calculated values
    p_text <- if (cor_test$p.value < 0.001) {
      "p < 0.001"
    } else {
      paste0("p = ", sprintf("%.3f", cor_test$p.value))
    }
    r_text <- paste0(cor_symbol, " = ", sprintf("%.3f", cor_test$estimate))
  }
  
  # Combine correlation coefficient and p values
  stats_text <- paste(r_text, p_text, sep = "\n")
  
  # Get data ranges for positioning
  x_range <- range(data[[x_var]], na.rm = TRUE)
  y_range <- range(data[[y_var]], na.rm = TRUE)
  
  # Calculate position based on stats_position parameter
  position_settings <- switch(stats_position,
                              "top-left" = list(
                                x = x_range[1] + 0.05 * diff(x_range),
                                y = y_range[1] + 0.95 * diff(y_range),
                                hjust = 0, vjust = 1
                              ),
                              "top-right" = list(
                                x = x_range[1] + 0.95 * diff(x_range),
                                y = y_range[1] + 0.95 * diff(y_range),
                                hjust = 1, vjust = 1
                              ),
                              "bottom-left" = list(
                                x = x_range[1] + 0.05 * diff(x_range),
                                y = y_range[1] + 0.05 * diff(y_range),
                                hjust = 0, vjust = 0
                              ),
                              "bottom-right" = list(
                                x = x_range[1] + 0.95 * diff(x_range),
                                y = y_range[1] + 0.05 * diff(y_range),
                                hjust = 1, vjust = 0
                              ),
                              # Default to top-left if invalid position
                              list(
                                x = x_range[1] + 0.05 * diff(x_range),
                                y = y_range[1] + 0.95 * diff(y_range),
                                hjust = 0, vjust = 1
                              )
  )
  
  # Create the plot
  p <- ggscatter(data, x = x_var, y = y_var,
                 color = "black",
                 #color = "Sex",
                 #shape = "Sex",
                 size = 5,
                 add = NULL,
                 conf.int = TRUE,
                 cor.coef = FALSE,
                 ylab = "Left Hippocampal Volume (cm3)",
                 xlab = gsub("[._]", " ", x_var),
                 title = "")
  
  p <- p + geom_smooth(aes(x = .data[[x_var]], y = .data[[y_var]]), 
                       method = "lm", se = TRUE, color = "black",
                       inherit.aes = FALSE, data = data)
  
  p <- p +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title.x = element_text(size = 25, face = "bold", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 0)),  # Increased spacing
      axis.title.y = element_text(size = 25, face = "bold",
                                  margin = margin(t = 0, r = 15, b = 0, l = 0)),  # Increased spacing
      axis.text.x = element_text(size = 20,  # Increased size
                                 margin = margin(t = 8, r = 0, b = 0, l = 0)),  # Increased spacing
      axis.text.y = element_text(size = 20,  # Increased size
                                 margin = margin(t = 0, r = 8, b = 0, l = 0)),  # Increased spacing
      axis.ticks.length = unit(0.3, "cm"),  # Increased tick length
      legend.title = element_text(size = 22, face = "bold"),  # Increased size
      legend.text = element_text(size = 20),  # Increased size
      legend.key.size = unit(1.5, "cm")  # Increased size
    ) +
    # Add correlation coefficient and p-value annotation
    annotate("text", 
             x = position_settings$x, 
             y = position_settings$y, 
             hjust = position_settings$hjust, 
             vjust = position_settings$vjust,
             label = stats_text,
             size = 10, 
             color = "black")
  
  return(p)
}
# Add a linear regression for Males
# p <- p + geom_smooth(
#   data = subset(data, Gender == "F"), 
#   aes(x = .data[[x_var]], y = .data[[y_var]], color = "F"),
#   method = "lm", 
#   se = F, 
#   linetype = "solid",
#   linewidth = 2
# )

# Add a loess smoothing for Females
# p <- p + geom_smooth(
#   data = subset(data, Gender == "M"),
#   aes(x = .data[[x_var]], y = .data[[y_var]], color = "M"),
#   method = "lm", 
#   se = F, 
#   linetype = "solid",
#   linewidth = 2
# )

# --- Left Hippocampus ---
plot_L_immediate <- create_correlation_plot(mci_data, "LMSB.immediate", "L_Hippocampus",
                                            "L_Hippocampus vs LMSB Immediate",
                                            manual_r = "0.609", manual_p = "0.025",
                                            cor_method = "pearson",
                                            stats_position = "top-left")

plot_L_delayed <- create_correlation_plot(mci_data, "LMSB.delayed", "L_Hippocampus",
                                          "L_Hippocampus vs LMSB Delayed",
                                          manual_r = "0.537", manual_p = "0.042",
                                          cor_method = "pearson",
                                          stats_position = "top-left")

plot_L_TMTA <- create_correlation_plot(mci_data, "TMT.A", "L_Hippocampus",
                                       "L_Hippocampus vs TMT-A",
                                       manual_r = "-0.368", manual_p = "0.042",
                                       cor_method = "kendall",
                                       stats_position = "top-right")

# --- Right Hippocampus ---
plot_R_3O3P <- create_correlation_plot(mci_data, "X3.Object.3.Places", "R_Hippocampus",
                                       "R_Hippocampus vs 3O3P",
                                       manual_r = "0.427", manual_p = "0.042",
                                       cor_method = "kendall",
                                       stats_position = "top-left")

plot_R_ROCFcopy <- create_correlation_plot(mci_data, "ROCF.copy", "R_Hippocampus",
                                           "R_Hippocampus vs ROCF Copy",
                                           manual_r = "0.499", manual_p = "0.018",
                                           cor_method = "kendall",
                                           stats_position = "top-left")

plot_R_ROCFrecall <- create_correlation_plot(mci_data, "ROCF.recall", "R_Hippocampus",
                                             "R_Hippocampus vs ROCF Recall",
                                             manual_r = "0.351", manual_p = "0.042",
                                             cor_method = "kendall",
                                             stats_position = "top-left")

plot_R_TMTA <- create_correlation_plot(mci_data, "TMT.A", "R_Hippocampus",
                                       "R_Hippocampus vs TMT-A",
                                       manual_r = "-0.353", manual_p = "0.042",
                                       cor_method = "kendall",
                                       stats_position = "top-right")

# --- Bilateral Hippocampus ---
plot_B_3O3P <- create_correlation_plot(mci_data, "X3.Object.3.Places", "B_Hippocampus",
                                       "Bilateral Hippocampus vs 3O3P",
                                       manual_r = "0.419", manual_p = "0.042",
                                       cor_method = "kendall",
                                       stats_position = "top-left")

plot_B_ROCFcopy <- create_correlation_plot(mci_data, "ROCF.copy", "B_Hippocampus",
                                           "Bilateral Hippocampus vs ROCF Copy",
                                           manual_r = "0.494", manual_p = "0.018",
                                           cor_method = "kendall",
                                           stats_position = "top-left")

plot_B_TMTA <- create_correlation_plot(mci_data, "TMT.A", "B_Hippocampus",
                                       "Bilateral Hippocampus vs TMT-A",
                                       manual_r = "-0.345", manual_p = "0.042",
                                       cor_method = "kendall",
                                       stats_position = "top-right")

# --- Save all plots at 1050x788 ---
plots_list <- list(
  plot_L_immediate, plot_L_delayed, plot_L_TMTA,
  plot_R_3O3P, plot_R_ROCFcopy, plot_R_ROCFrecall, plot_R_TMTA,
  plot_B_3O3P, plot_B_ROCFcopy, plot_B_TMTA
)

filenames <- c(
  "L_HPC_LMSB_Imm", "L_HPC_LMSB_Del", "L_HPC_TMTA",
  "R_HPC_3O3P", "R_HPC_ROCF_Copy", "R_HPC_ROCF_Recall", "R_HPC_TMTA",
  "B_HPC_3O3P", "B_HPC_ROCF_Copy", "B_HPC_TMTA"
)

for (i in seq_along(plots_list)) {
  ggsave(paste0(filenames[i], ".png"), plot = plots_list[[i]], width = 1050, height = 788, units = "px", dpi = 120)
  cat("Saved:", filenames[i], "\n")
}

# plot_immediate <- plot_immediate + 
#   theme(plot.margin = margin(t = 10, r = 30, b = 10, l = 10, unit = "pt"))
# 
# 
# plot_TMTA <- plot_TMTA + 
#   theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt"))
# 
# # Then combine
# combined_plot <- ggarrange(plot_immediate, plot_TMTA,
#                            #labels = c("A", "B"),
#                            ncol = 2, nrow = 1,
#                            common.legend = FALSE)
# 
# # Add main title to combined plot
# combined_plot <- annotate_figure(combined_plot,
#                                  top = text_grob("Hippocampal Volume Correlations with LMSB and TMT-A scores in MCI Patients", 
#                                                  color = "black", face = "bold", size = 14))
# 
# # Display combined plot
# print("Displaying combined correlation plot...")
# print(combined_plot)
# 
# # Save plots
# ggsave("L_HPC_LSMB_immediate.png", plot_immediate, width = 15, height = 9, dpi = 900)
# ggsave("L_HPC_combined.png", combined_plot, width = 18, height = 9, dpi = 900)
# 
# # Print detailed correlation statistics
# cat("\n", "="*50, "\n")
# cat("DETAILED CORRELATION STATISTICS\n")
# cat("="*50, "\n")
# 
# # L_HPC vs LSMB_immediate
# cor_immediate <- cor.test(mci_data$L_HPC, mci_data$LSMB_immediate, method = "pearson")
# cat("\nL_HPC vs LSMB_immediate:\n")
# cat("Pearson's r =", round(cor_immediate$estimate, 3), "\n")
# cat("95% CI: [", round(cor_immediate$conf.int[1], 3), ", ", round(cor_immediate$conf.int[2], 3), "]\n")
# cat("p-value =", format.pval(cor_immediate$p.value, digits = 3), "\n")
# cat("Degrees of freedom =", cor_immediate$parameter, "\n")
# 
# # L_HPC vs LSMB_delayed  
# cor_delayed <- cor.test(mci_data$L_HPC, mci_data$LSMB_delayed, method = "pearson")
# cat("\nL_HPC vs LSMB_delayed:\n")
# cat("Pearson's r =", round(cor_delayed$estimate, 3), "\n")
# cat("95% CI: [", round(cor_delayed$conf.int[1], 3), ", ", round(cor_delayed$conf.int[2], 3), "]\n")
# cat("p-value =", format.pval(cor_delayed$p.value, digits = 3), "\n")
# cat("Degrees of freedom =", cor_delayed$parameter, "\n")
