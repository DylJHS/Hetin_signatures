library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # This is the SLURM_ARRAY_TASK_ID

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")
rna_data_path <- "Data/RNA_Data/Model_Input/Train/train_"

cat("The index for this run is: ", index, "\n")

selected_depth <- 1
selected_min_child <- 1
selected_lr <- 0.3
selected_gamma <- 0
selected_trees <- 10000

# RNA SOI SETS
# Expected Counts
exp_set <- read.csv(
  paste0(
    rna_data_path,
    "exp_soi.csv"
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

log_exp <- read.csv(
  paste0(
    rna_data_path,
    "log_exp_soi.csv"
  ),
  row.names = 1
)

log_scld_exp <- read.csv(
  paste0(
    rna_data_path,
    "log_scld_exp_soi.csv"
  ),
  row.names = 1
)

# Transcripts Per Million
tpm_set <- read.csv(
  paste0(
    rna_data_path,
    "tpm_soi.csv"
  ),
  row.names = 1
)

cat("\n\n TPM df: \n")
print(head(tpm_set[, 1:5]))

scld_tpm_set <- read.csv(
  paste0(
    rna_data_path,
    "scld_tpm_soi.csv"
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

log_scld_tpm <- read.csv(
  paste0(
    rna_data_path,
    "log_scld_tpm_soi.csv"
  ),
  row.names = 1
)

# Arm Level Aneuploidies
# Load the data
chr_cnv <- read_tsv(
  "Data/CIN_Features/CNV_Data/PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
) %>%
  replace(is.na(.), 0) %>%
  select(-c("Type", "Aneuploidy Score")) %>%
  mutate(Sample = str_replace_all(Sample, "-", "\\.")) %>%
  column_to_rownames("Sample") %>%
  mutate_all(~ replace(., . == 1, 2)) %>%
  mutate_all(~ replace(., . == 0, 1)) %>%
  mutate_all(~ replace(., . == -1, 0))

# Calc the class weights
freq_df <- chr_cnv %>%
  as.data.frame() %>%
  gather(key = "arm", value = "freq") %>%
  group_by(arm) %>%
  count(freq) %>%
  as.data.frame() %>%
  spread(key = freq, value = n) %>%
  replace(is.na(.), 0)

arm_weights <- freq_df %>%
  mutate(total = rowSums(select(., -arm))) %>%
  mutate_at(vars(-arm, -total), list(~ 1 - round(. / total, 2))) %>%
  mutate(total = rowSums(select(., -arm, -total))) %>%
  mutate_at(vars(-arm, -total), list(~ round(. / total, 2))) %>%
  select(-total) %>%
  t() %>%
  as.data.frame() %>%
  setNames(make.unique(unlist(.[1, ]))) %>% # Convert first row to column names
  .[-1, ]

print(arm_weights[, 1:5])
rm(freq_df)

# Define the features to be used
aneu_cat_feature_list <- colnames(chr_cnv)
selected_feature <- aneu_cat_feature_list[[index]]
rm(aneu_cat_feature_list)

# Determine the class weights for the target feature
target_weights <- arm_weights[, index]
selected_weights <- target_weights
print(selected_weights)

target <- as.character(colnames(arm_weights)[index])
cat("\n", target, "weights: ")
print(target_weights)
cat("\n\n")
rm(arm_weights)

# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) target_weights[as.numeric(x)])
}

# MODELLING

aneu_cat_metrics_df <- data.frame(
  RNA_Set = character(),
  Trees = numeric(),
  Feature = character(),
  Depth = numeric(),
  Child_weight = numeric(),
  Learning_Rate = numeric(),
  Gamma = numeric(),
  Weight_loss = numeric(),
  Weight_normal = numeric(),
  Weight_gain = numeric(),
  Trained_mlogloss = numeric(),
  Test_mlogloss = numeric()
)

rna_list <- list(
  transcripts_per_million = tpm_set,
  scaled_transcripts_per_million = scld_tpm_set,
  log_scaled_transcripts_per_million = log_scld_tpm,
  log_transcripts_per_million = log_tpm,
  expected_counts = exp_set,
  scaled_expected_counts = scld_exp_set,
  log_expected_counts = log_exp,
  log_scaled_expected_counts = log_scld_exp
)
rna_names <- names(rna_list)

for (i in 1:length(rna_list)) {
  rna_data <- rna_list[[i]]
  selected_rna_set <- names(rna_list)[i]
  cat(paste0("\t", selected_rna_set, "\n"))

  full_df <- merge(rna_data,
    chr_cnv,
    by = "row.names"
  )

  y <- as.integer(full_df[[selected_feature]])
  X <- full_df %>% select(-c("Row.names", colnames(chr_cnv)))

  train_y_factor <- factor(y, levels = c(0, 1, 2))
  weights <- as.numeric(feature_digit_function(train_y_factor))

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)

  cat(paste0(
    "\t\t Depth: ", selected_depth,
    "\n"
  ))

  m_xgb_untuned <- xgb.cv(
    data = xgb_data,
    nrounds = selected_trees,
    objective = "multi:softmax",
    eval_metric = "mlogloss",
    early_stopping_rounds = 100,
    nfold = 5,
    max_depth = selected_depth,
    eta = selected_lr,
    gamma = selected_gamma,
    num_class = 3,
    print_every_n = 15
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

  best_mlogloss_train <- if (best_iteration > 0) {
    m_xgb_untuned$evaluation_log$train_mlogloss_mean[best_iteration]
  } else {
    NA # Or appropriate default/error value
  }

  best_mlogloss_test <- if (best_iteration > 0) {
    m_xgb_untuned$evaluation_log$test_mlogloss_mean[best_iteration]
  } else {
    NA # Or appropriate default/error value
  }

  cat(paste0(
    "The best iteration occurs with tree #: ",
    best_iteration, "\n\n"
  ))

  aneu_cat_metrics_df <- rbind(aneu_cat_metrics_df, data.frame(
    RNA_Set = selected_rna_set,
    Trees = best_iteration,
    Feature = selected_feature,
    Depth = selected_depth,
    Child_weight = selected_min_child,
    Learning_Rate = selected_lr,
    Gamma = selected_gamma,
    Weight_loss = selected_weights[1],
    Weight_normal = selected_weights[2],
    Weight_gain = selected_weights[3],
    Trained_mlogloss = best_mlogloss_train,
    Test_mlogloss = best_mlogloss_test
  ))
}

datetime <- as.character(Sys.time()) %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "_") %>%
  str_replace_all("\\.", "_")

file_name <- paste0(
  "Data/Model_output/categorical/mlogloss_measures/",
  selected_feature, "_",
  index, "_",
  datetime, ".csv"
) %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "_")

write.csv(
  aneu_cat_metrics_df,
  file = file_name,
  row.names = FALSE
)

cat(paste0("\n Completed processing for index: ", index, "\n"))
