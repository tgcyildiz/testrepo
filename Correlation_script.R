# Load Libraries & Options
rm(list=ls())
library(OpenMx)
library(psych); library(polycor)
library(mets)
library(dplyr)
library(gridExtra)
source("miFunctions.R")

setwd('C:\\Users\\tuyild\\OneDrive - Karolinska Institutet\\Documents\\Twins_OpenMX\\Cerebellum\\New_correlation_plots')
# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA
# Load Data
data_file <- "Cerebellum_volumes.csv"
# Import data
dat <- read.csv(file=data_file, header=TRUE, sep = ";")
variable_input1 <- "age"

# Remove NAN values
complete_cases <- complete.cases(dat[[variable_input1]])
rows_with_nan <- which(!complete_cases)
pair_ids_with_nan <- dat$pair_id[rows_with_nan]
dat <- dat[!(dat$pair_id %in% pair_ids_with_nan), ]

#Split the MZ and DZ groups
MZ_var1 <- subset(dat, zyg=="1")
DZ_var1 <- subset(dat, zyg=="2")

twin_variable_one <- variable_input1
twin_variable_one1 <- sprintf("%s%s", twin_variable_one, "1")
twin_variable_one2 <- sprintf("%s%s", twin_variable_one, "2")


# Split on MZ and DZ for first variable input.
twin_variables <- fast.reshape(dat, id="pair_id",varying=c(twin_variable_one))
mz_var1 <- subset(twin_variables, zyg=="1")[,c(twin_variable_one1,twin_variable_one2)]
dz_var1 <- subset(twin_variables, zyg=="2")[,c(twin_variable_one1,twin_variable_one2)]

#All MZ
mean_var1_MZ <- mean(MZ_var1[[variable_input1]], na.rm = TRUE)
variance_var1_MZ <- var(MZ_var1[[variable_input1]], na.rm = TRUE)
# MZ T1 and T2
mean_var1_MZ_T1 <- mean(mz_var1[[twin_variable_one1]], na.rm = TRUE)
variance_var1_MZ_T1 <- var(mz_var1[[twin_variable_one1]], na.rm = TRUE)
mean_var1_MZ_T2 <- mean(mz_var1[[twin_variable_one2]], na.rm = TRUE)
variance_var1_MZ_T2 <- var(mz_var1[[twin_variable_one2]], na.rm = TRUE)

#All DZ
mean_var1_DZ <- mean(DZ_var1[[variable_input1]], na.rm = TRUE)
variance_var1_DZ <- var(DZ_var1[[variable_input1]], na.rm = TRUE)
# MZ T1 and T2
mean_var1_DZ_T1 <- mean(dz_var1[[twin_variable_one1]], na.rm = TRUE)
variance_var1_DZ_T1 <- var(dz_var1[[twin_variable_one1]], na.rm = TRUE)
mean_var1_DZ_T2 <- mean(dz_var1[[twin_variable_one2]], na.rm = TRUE)
variance_var1_DZ_T2 <- var(dz_var1[[twin_variable_one2]], na.rm = TRUE)

#Covariances
covariance_MZ <- cov(mz_var1[[twin_variable_one1]], mz_var1[[twin_variable_one2]], use = "complete.obs")
covariance_DZ <- cov(dz_var1[[twin_variable_one1]], dz_var1[[twin_variable_one2]], use = "complete.obs")

# Create a table including means, variances, and covariances
results_table <- data.frame(
  Group = c("All MZ", "All MZ T1", "All MZ T2", "All DZ", "All DZ T1", "All DZ T2"),
  Mean = c(mean_var1_MZ, mean_var1_MZ_T1, mean_var1_MZ_T2, mean_var1_DZ, mean_var1_DZ_T1, mean_var1_DZ_T2),
  Variance = c(variance_var1_MZ, variance_var1_MZ_T1, variance_var1_MZ_T2, variance_var1_DZ, variance_var1_DZ_T1, variance_var1_DZ_T2),
  Covariance = c(covariance_MZ, NA, NA, covariance_DZ, NA, NA)
)

# Print the table
knitr::kable(results_table, format = "markdown")

# Calculate correlations for MZ and DZ twin pairs
correlation_MZ <- cor(mz_var1[[twin_variable_one1]], mz_var1[[twin_variable_one2]], use = "complete.obs")
correlation_DZ <- cor(dz_var1[[twin_variable_one1]], dz_var1[[twin_variable_one2]], use = "complete.obs")

# Plot correlations
par(mfrow = c(1, 2)) # Arrange plots in 1 row, 2 columns

# Plot for MZ twin pairs
plot(mz_var1[[twin_variable_one1]], mz_var1[[twin_variable_one2]],
     main = paste("Correlation for MZ Twin Pairs\nCorrelation:", round(correlation_MZ, 2)),
     xlab = twin_variable_one1,
     ylab = twin_variable_one2)
abline(lm(mz_var1[[twin_variable_one2]] ~ mz_var1[[twin_variable_one1]]), col = "red") # Add regression line

# Plot for DZ twin pairs
plot(dz_var1[[twin_variable_one1]], dz_var1[[twin_variable_one2]],
     main = paste("Correlation for DZ Twin Pairs\nCorrelation:", round(correlation_DZ, 2)),
     xlab = twin_variable_one1,
     ylab = twin_variable_one2)
abline(lm(dz_var1[[twin_variable_one2]] ~ dz_var1[[twin_variable_one1]]), col = "red") # Add regression line

pair_ids_with_nan


# Function to calculate proportions of males and females
calculate_proportions <- function(data) {
  num_male <- sum(data == 0, na.rm = TRUE)
  num_female <- sum(data == 1, na.rm = TRUE)
  proportion_male <- num_male / length(data)
  proportion_female <- num_female / length(data)
  return(c(Male = proportion_male, Female = proportion_female))
}

# Calculate proportions for each group
proportions_mz_twin1 <- calculate_proportions(mz_var1[[twin_variable_one1]])
proportions_mz_twin2 <- calculate_proportions(mz_var1[[twin_variable_one2]])
proportions_dz_twin1 <- calculate_proportions(dz_var1[[twin_variable_one1]])
proportions_dz_twin2 <- calculate_proportions(dz_var1[[twin_variable_one2]])

# Combine proportions into a data frame
proportions_df <- data.frame(
  Group = c("MZ Twin 1", "MZ Twin 2", "DZ Twin 1", "DZ Twin 2"),
  Male = c(proportions_mz_twin1["Male"], proportions_mz_twin2["Male"], proportions_dz_twin1["Male"], proportions_dz_twin2["Male"]),
  Female = c(proportions_mz_twin1["Female"], proportions_mz_twin2["Female"], proportions_dz_twin1["Female"], proportions_dz_twin2["Female"])
)

proportions_plot <- ggplot(proportions_df, aes(x = Group)) +
  geom_bar(aes(y = Male), stat = "identity", fill = "blue", position = "dodge") +
  geom_bar(aes(y = Female), stat = "identity", fill = "purple", position = "dodge") +
  labs(title = "Proportion of Males and Females in Twin Pairs", y = "Proportion", fill = "Sex") +
  scale_y_continuous(labels = scales::percent_format()) + # Convert y-axis to percentage
  theme_minimal() + # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  guides(fill = guide_legend(title = "Sex", override.aes = list(fill = c("blue", "purple"), labels = c("Male", "Female")))) # Add legend with custom labels

# Display plot
print(proportions_plot)


# Combine age1 and age2 data into a single data frame for each group
mz_age <- rbind(
  cbind(data.frame(age = mz_var1[[twin_variable_one1]], group = "MZ Twin 1")),
  cbind(data.frame(age = mz_var1[[twin_variable_one2]], group = "MZ Twin 2"))
)

dz_age <- rbind(
  cbind(data.frame(age = dz_var1[[twin_variable_one1]], group = "DZ Twin 1")),
  cbind(data.frame(age = dz_var1[[twin_variable_one2]], group = "DZ Twin 2"))
)

# Create histograms for MZ twins
hist_mz <- ggplot(mapping = aes(x = c(mz_var1[[twin_variable_one1]], mz_var1[[twin_variable_one2]]))) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "Histogram of Age for MZ Twins", x = "Age", y = "Count") +
  theme_minimal()

# Create histograms for DZ twins
hist_dz <- ggplot(mapping = aes(x = c(dz_var1[[twin_variable_one1]], dz_var1[[twin_variable_one2]]))) +
  geom_histogram(binwidth = 1, fill = "green", color = "black") +
  labs(title = "Histogram of Age for DZ Twins", x = "Age", y = "Count") +
  theme_minimal()

# Display plots side by side
grid.arrange(hist_mz, hist_dz, nrow = 1)

# Combine age data for MZ twins
age_mz <- c(mz_var1[[twin_variable_one1]], mz_var1[[twin_variable_one2]])

# Combine age data for DZ twins
age_dz <- c(dz_var1[[twin_variable_one1]], dz_var1[[twin_variable_one2]])

# Summary statistics for MZ twins
summary_mz <- summary(age_mz)
cat("Summary Statistics for Ages of MZ Twins:\n")
print(summary_mz)

# Summary statistics for DZ twins
summary_dz <- summary(age_dz)
cat("\nSummary Statistics for Ages of DZ Twins:\n")
print(summary_dz)

# Perform t-test
t_test_result <- t.test(age_mz, age_dz)
cat("\nT-Test Result:\n")
print(t_test_result)