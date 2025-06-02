# Purpose: Split the RNA data into a training and testing set for the cancer specific data
# ensuring they all contain the same sample IDs
# and ordering the data identically

library(dplyr)

path2cancersets <- "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/Cancer_specific_data/Model_input/RNA/"
setwd(path2cancersets)



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

for (cancer in list.files(".")) {
  cat("\n\n\n\t\t\t\t\t Cancer: ", cancer, "\n")

  path <- paste0(path2cancersets, cancer, "/")
  full_path <- paste0(path, "Full/")
  if (!dir.exists(full_path)) {
    next
  }

  # Create empty lists for thee sampleIDs in each rna set
  tpm_sampleIDs <- c()
  scld_tpm_sampleIDs <- c()
  log_scld_tpm_sampleIDs <- c()
  log_tpm_sampleIDs <- c()
  exp_sampleIDs <- c()
  scld_exp_sampleIDs <- c()
  log_exp_sampleIDs <- c()
  log_scld_exp_sampleIDs <- c()

  for (file in rna_list) {
    name <- substr(file, 1, nchar(file) - 8)
    cat("\n\t\t\t Name: ", name, "\n")
    ids <- read.csv(
      paste0(full_path, file),
      row.names = 1
    ) %>%
      rownames()

    cat("\n\t\t File: ", file, "\n")
    print(length(ids))

    assign(paste0(name, "_sampleIDs"), ids)
  }

  # Find the common sample IDs in each rna set
  common_sampleIDs <- Reduce(intersect, list(
    tpm_sampleIDs, scld_tpm_sampleIDs,
    log_scld_tpm_sampleIDs, log_tpm_sampleIDs, exp_sampleIDs,
    scld_exp_sampleIDs, log_exp_sampleIDs, log_scld_exp_sampleIDs
  ))

  cat("Number of common sample IDs: ", length(common_sampleIDs), "\n\n")

  # Get random 75% of the common sample IDs
  num_samples_to_select <- round(0.75 * length(common_sampleIDs))
  random_sample_ids <- sample(
    common_sampleIDs, num_samples_to_select,
    replace = FALSE
  )

  train_folder <- paste0(path, "Train")
  if (!dir.exists(train_folder)) {
    dir.create(train_folder)
  }
  test_folder <- paste0(path, "Test")
  if (!dir.exists(test_folder)) {
    dir.create(test_folder)
  }

  # Create new rna sets with only the common sample IDs from the random sample
  for (file in rna_list) {
    name <- substr(file, 1, nchar(file) - 8)
    data <- read.csv(
      paste0(full_path, file)
    )

    train_data <- data[data[, 1] %in% random_sample_ids, ] %>%
      arrange(row.names(.))
    test_set <- data[!data[, 1] %in% random_sample_ids, ] %>%
      arrange(row.names(.))
    
    cat("\n Training & Testing data: ", name, "\n")
    print(dim(train_data))
    print(dim(test_set))

    write.csv(
      train_data,
      paste0(train_folder, "/train_", file),
      row.names = FALSE
    )

    write.csv(
      test_set,
      paste0(test_folder, "/test_", file),
      row.names = FALSE
    )
  }
}
