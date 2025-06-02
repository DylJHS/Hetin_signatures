library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

rna_data_path <- "Data/RNA_Data/Model_Input/Train/train_"
index <- 4

# Default parameters


# RNA EC SOI sets
log_scld_exp <- read.csv(
  paste0(
    rna_data_path,
    "log_scld_exp_soi.csv"
  ),
  row.names = 1
)

scld_exp_set <- read.csv(
  paste0(
    rna_data_path,
    "scld_exp_soi.csv"
  ),
  row.names = 1
)

cat("\n\n Log-Scaled Expected Counts df: \n")
print(head(log_scld_exp[, 1:10]))

# RNA TPM sets
tpm_set <- read.csv(
  paste0(
    rna_data_path,
    "tpm_soi.csv"
  ),
  row.names = 1
)
cat("\n\n TPM df: \n")
print(head(tpm_set[, 1:10]))

scld_tpm_set <- read.csv(
  paste0(
    rna_data_path,
    "scld_tpm_soi.csv"
  ),
  row.names = 1
)

log_scld_tpm <- read.csv(
  paste0(
    rna_data_path,
    "log_scld_tpm_soi.csv"
  ),
  row.names = 1
)

log_tpm <- read.csv(
  paste0(
    rna_data_path,
    "log_tpm_soi.csv"
  ),
  row.names = 1
)

# HRD scores
# import HRD data
ori_hrd <- read_tsv("Data/CIN_Features/TCGA.HRD_withSampleID.txt")
t_hrd <- as.data.frame(t(ori_hrd)) # Transpose the data
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ] # Set the first row as the column names
rm(ori_hrd)

# Remove the first row, convert the data to numeric, rename the column, reformat the rownames, and remove the HRD column
hrd <- as.data.frame(first_hrd[-1, ]) %>%
  mutate_all(as.numeric) %>%
  rename(loh_hrd = "hrd-loh") %>%
  mutate(new = str_replace_all(rownames(.), "-", "\\.")) %>%
  select(-"HRD")

rownames(hrd) <- hrd$new # Set the rownames to the new column
hrd <- hrd %>% # Remove the new column
  select(-new)

cat("\n\n hrd data: \n")
print(head(hrd))
rm(t_hrd)
rm(first_hrd)


# Pericentromeric CNVs
# Import pericentromeric CNV data, replace NA values with 0
# and reformat the sampleID column
peri_cnv <- read.csv(
  "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/CIN_Features/CNV_Data/lim_alpha_incl_TCGA_pericentro.csv"
) %>%
  mutate_all(~ replace(., is.na(.), 0)) %>%
  mutate(sampleID = gsub("-", ".", sampleID))
cat("\n\n pericentromeric data: \n")
print(head(peri_cnv[, 1:5]))
cat("\n\n All dfs loaded \n")

# Merge the HRD and pericentromeric CNV data to get full CIN data
full_cin <- merge(
  hrd,
  peri_cnv,
  by.x = "row.names",
  by.y = "sampleID"
) %>%
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")
# Get the target feature list
aneu_reg_feature_list <- colnames(full_cin)

cat("\n\n All feature names: ", colnames(full_cin), "\n")
cat("\n\n full cin data: \n")
print(head(full_cin[, 1:5]))
rm(hrd)
rm(peri_cnv)

# import the hyperparameter df
reg_hyperparam_df <- read.csv(
  "Data/XGB_hyperparameters.csv"
)

print(head(reg_hyperparam_df))

# Select the specific parameters for this array
selected_combination <- reg_hyperparam_df[index, ]
selected_feature <- selected_combination$Feature
selected_rna_set <- selected_combination$RNA_set
selected_trees <- selected_combination$Trees
selected_depth <- selected_combination$Max_depth
selected_eta <- selected_combination$Eta
selected_child_weight <- 1

rna_list <- list(
  transcripts_per_million = tpm_set,
  scaled_transcripts_per_million = scld_tpm_set,
  log_scaled_transcripts_per_million = log_scld_tpm,
  log_transcripts_per_million = log_tpm,
  scaled_expected_counts = scld_exp_set,
  log_scaled_expected_counts = log_scld_exp
)
rna_names <- names(rna_list)

#     MODELLING
cat(paste0(
  "\n\n Running model for feature: ",
  selected_feature,
  " RNA set: ",
  selected_rna_set,
  " Trees: ",
  selected_trees, "\n"
))

# Now select the data based on these choices
rna_data <- rna_list[[selected_rna_set]]
cat("\n\n RNA data: \n")
print(head(rna_data[, 1:5]))
cat("\n\n")

full_df <- merge(rna_data, full_cin, by = "row.names")
cat("\n\n full_df: \n")
print(head(full_df[, 1:5]))

aneu_reg_metrics_df <- data.frame(
  RNA_Set = character(),
  Trees = numeric(),
  Feature = character(),
  Depth = numeric(),
  Child_weight = numeric(),
  Learning_Rate = numeric(),
  Gamma = numeric(),
  Trained_RMSE = numeric(),
  Test_RMSE = numeric()
)

y <- as.numeric(full_df[[selected_feature]])
X <- full_df %>% select(-c("Row.names", colnames(full_cin)))
cat("\n\n Predictors: \n")
print(head(X[, 1:5]))
cat("\n\n ")

xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

grid <- expand.grid(
  gam = seq(0, 0.15, 0.05)
)

set.seed(99)
for (j in 1:nrow(grid)) { # nolint
  cat(paste0(
    "\t\t eta: ", lr,
    "\t\t gamma: ", grid$gam[j],
    "\t\t depth: ", depth,
    "\t\t trees: ", selected_trees,
    "\t\t child_weight: ", min_child,
    "\n"
  ))

  m_xgb_untuned <- xgb.cv(
    data = xgb_data,
    nrounds = selected_trees,
    objective = "reg:squarederror",
    eval_metric = "rmse",
    early_stopping_rounds = 10,
    nfold = 2,
    max_depth = selected_depth,
    min_child_weight = selected_child_weight,
    eta = selected_eta,
    gamma = grid$gam[j],
    verbose = 1
  )

  best_iteration <- 0

  # First, check if best_iteration is valid
  if (is.null(
    m_xgb_untuned$best_iteration
  ) ||
    m_xgb_untuned$best_iteration < 1) {
    cat(paste0(
      "Warning: No valid best_iteration found.",
      " Using last iteration values instead.\n"
    ))
    # Use the last iteration if best_iteration is not valid
    best_iteration <- nrow(m_xgb_untuned$evaluation_log)
  } else {
    # Ensure that the best_iteration does not exceed the number of rows logged
    if (m_xgb_untuned$best_iteration > nrow(m_xgb_untuned$evaluation_log)) {
      cat(paste0(
        "Warning: best_iteration exceeds the number of rows in evaluation_log.",
        " Adjusting to maximum available.\n"
      ))
      best_iteration <- nrow(m_xgb_untuned$evaluation_log)
    } else {
      best_iteration <- m_xgb_untuned$best_iteration
    }
  }

  # Accessing the RMSE values safely
  best_rmse_trained <- if (best_iteration > 0) {
    m_xgb_untuned$evaluation_log$train_rmse_mean[best_iteration]
  } else {
    NA # Or appropriate default/error value
  }

  best_rmse_test <- if (best_iteration > 0) {
    m_xgb_untuned$evaluation_log$test_rmse_mean[best_iteration]
  } else {
    NA # Or appropriate default/error value
  }

  cat(paste0(
    "The best iteration occurs with tree #: ",
    best_iteration, "\n\n"
  ))


  aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df, data.frame(
    RNA_Set = selected_rna_set,
    Trees = selected_trees,
    Feature = selected_feature,
    Depth = selected_depth,
    Child_weight = selected_child_weight,
    Learning_Rate = selected_eta,
    Gamma = grid$gam[j],
    Trained_RMSE = best_rmse_trained,
    Test_RMSE = best_rmse_test
  ))
}


datetime <- Sys.time() %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "_")

name <- paste0(
  "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/Model_output/categorical",
  selected_feature, "_",
  selected_rna_set, "_",
  datetime, ".csv"
) %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "_")


cat("\n Completed processing for index: ", index, "\n")
