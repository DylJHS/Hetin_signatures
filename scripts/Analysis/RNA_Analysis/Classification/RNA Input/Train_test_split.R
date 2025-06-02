library(dplyr)

path <- "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/RNA_Data/Model_Input/"
setwd(path)



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

# Get the sample IDs from the first data frame
set <- read.csv(paste0(path, "Full/", rna_list[1]))
all_sample_ids <- set[, 1]

# Calculate 25% of the total samples
num_samples_to_select <- round(0.25 * length(all_sample_ids))

# Get a random sample of sample IDs
set.seed(99)
random_sample_ids <- sample(
  all_sample_ids, num_samples_to_select, replace = FALSE
)

for (file in rna_list) {
  set <- read.csv(paste0(path, "Full/", file))
  cat("\n\n File: ", file, "\n")
  print(dim(set))

  train_set <- set[!set[, 1] %in% random_sample_ids, ]
  test_set <- set[set[, 1] %in% random_sample_ids, ]

  test_name <- paste0(path, "test", "_", file)
  train_name <- paste0(path, "train", "_", file)
  cat("Train: ", train_name, "\n")
  print(dim(train_set))

  # write.csv(
  #   test_set,
  #   test_name,
  #   row.names = FALSE
  # )

  # write.csv(
  #   train_set,
  #   train_name,
  #   row.names = FALSE
  # )
}
