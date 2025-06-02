# This script is to correct the training sets of the rna seq data
# by ensuring they all contain the same sample IDs
# and ordering the data identically

library(dplyr)


rna_data_path <- "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/RNA_Data/Model_Input/Train/"

rna_list <- list("tpm", "scld_tpm",
  "log_scld_tpm", "log_tpm",
  "exp", "scld_exp",
  "log_exp", "log_scld_exp"
)

# Create empty lists for thee sampleIDs in each rna set
tpm_sampleIDs <- c()
scld_tpm_sampleIDs <- c()
log_scld_tpm_sampleIDs <- c()
log_tpm_sampleIDs <- c()
exp_sampleIDs <- c()
scld_exp_sampleIDs <- c()
log_exp_sampleIDs <- c()
log_scld_exp_sampleIDs <- c()

# Get the sample IDs for each rna set
for (file in rna_list) {
  ids <- read.csv(paste0(rna_data_path, "mismatched/train_", file, "_soi.csv"),
    row.names = 1) %>%
    rownames()

  assign(paste0(file, "_sampleIDs"), ids)
}

# Find the common sample IDs in each rna set
common_sampleIDs <- Reduce(intersect, list(tpm_sampleIDs, scld_tpm_sampleIDs,
  log_scld_tpm_sampleIDs, log_tpm_sampleIDs, exp_sampleIDs,
  scld_exp_sampleIDs, log_exp_sampleIDs, log_scld_exp_sampleIDs)) 

cat("Number of common sample IDs: ", length(common_sampleIDs), "\n")

# create new rna sets with only the common sample IDs
for (file in rna_list) {
  data <- read.csv(paste0(rna_data_path, "mismatched/train_", file, "_soi.csv"),
    row.names = 1)

  new_data <- data[common_sampleIDs, ] %>%
  arrange(row.names(.))

  write.csv(new_data, paste0(rna_data_path, file, "_soi.csv"))

  # print(dim(data))
}