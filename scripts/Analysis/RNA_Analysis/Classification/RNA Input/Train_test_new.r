# This script processes RNA data files by:

# 	1.	Setting up the directory paths.
# 	2.	Listing the RNA data files to be processed.
# 	3.	Loading sample IDs from each file and storing them.
# 	4.	Identifying common sample IDs across all datasets.
# 	5.	Randomly selecting 75% of the common sample IDs for training.
# 	6.	Preparing directories for training, testing, and full sets.
# 	7.	Splitting each data file into training, testing, and full sets based on the common IDs.
# 	8.	Saving the processed datasets to their respective directories.

library(dplyr)

# Set the path to the directory containing RNA data files
path <- ""
setwd(path)

full_path <- paste0(path, "All/")

rna_list <- c(
  "log_scld_tpm_soi.csv",
  "scld_tpm_soi.csv",
  "log_tpm_soi.csv",
  "tpm_soi.csv",
  "exp_soi.csv",
  "log_exp_soi.csv",
  "log_scld_exp_soi.csv",
  "scld_exp_soi.csv"
)

# Create a list to store sample IDs from each file
sampleIDs_list <- list()

# Load sample IDs from all files and store in the list
for (file in rna_list) {
  dataset <- read.csv(paste0(full_path,"full_", file))
  sampleIDs_list[[file]] <- dataset[, 1]
  cat("\n Loaded Sample IDs from: ", file, "\n Number of IDs: ", length(dataset[, 1]), "\n")
}

# # Find the common sample IDs across all datasets
common_sampleIDs <- Reduce(intersect, sampleIDs_list)
cat("\nNumber of common sample IDs: ", length(common_sampleIDs), "\n")

# # Get random 75% of the common sample IDs
set.seed(99)
num_samples_to_select <- round(0.75 * length(common_sampleIDs))
random_sample_ids <- sample(common_sampleIDs, num_samples_to_select, replace = FALSE)

# # Prepare directories for train and test sets and the full common sets
train_folder <- paste0(path, "Train/")
test_folder <- paste0(path, "Test/")
Full_folder <- paste0(path, "Full/")
if (!dir.exists(train_folder)) dir.create(train_folder)
if (!dir.exists(test_folder)) dir.create(test_folder)
if (!dir.exists(Full_folder)) dir.create(Full_folder)

# Process each file and split into training, testing and full sets based on common IDs
for (file in rna_list) {
  dataset <- read.csv(paste0(full_path, "full_", file))
  assign(paste0("full_", file), dataset)
  
  # Filter data for training and testing, then arrange by row names
  train_data <- dataset %>%
    filter(.[, 1] %in% random_sample_ids) %>%
    arrange(.[, 1])
  
  # test_data <- dataset %>%
    filter(!.[, 1] %in% random_sample_ids) %>%
    arrange(.[, 1])
  
  full_data <- dataset %>%
    filter(.[, 1] %in% common_sampleIDs) %>%
    arrange(.[, 1])
  
  # Identify any IDs present in the filtered data that shouldn't be there
  ids_to_remove <- setdiff(full_data$X, common_sampleIDs)
  
  # If there are any such IDs, remove them
  if (length(ids_to_remove) > 0) {
    full_data <- full_data %>%
      filter(!.$X %in% ids_to_remove)
    cat("Removed ", length(ids_to_remove), " incorrect IDs from full_data.\n")
  }
  assign(paste0("full_", file), full_data)
  
  # Print dimensions for verification
  cat("\n\n File: ", file, "\n")
  print(dim(train_data))
  print(dim(test_data))
  print(dim(full_data))
  
  # # Write datasets to respective files
  write.csv(train_data, paste0(train_folder, "train_", file), row.names = FALSE)
  write.csv(test_data, paste0(test_folder, "test_", file), row.names = FALSE)
  write.csv(full_data, paste0(Full_folder, "full_", file), row.names = FALSE)
}