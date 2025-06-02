library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # This is the SLURM_ARRAY_TASK_ID

print(index)

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

rna_data_path <- "Data/RNA_Data/Model_Input/Train/train_"
selected_depth <- 5
selected_min_child <- 1
selected_trees <- 1000
selected_lr <- 0.3

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

cat("\n\n arm level data: \n")
print(head(chr_cnv[, 1:5]))

aneu_cat_feature_list <- colnames(chr_cnv)

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
  Trained_Logloss = numeric(),
  Test_Logloss = numeric()
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

combinations <- expand.grid(
  feature = aneu_cat_feature_list,
  RNA_Set = rna_names,
  stringsAsFactors = FALSE
)

cat("\n\n All combinations: ")
print(combinations)

total_combinations <- nrow(combinations)
cat("\n\n Number fo total combinations: ", total_combinations)

# Select the specific feature and RNA set based on the SLURM task ID
selected_combination <- combinations[index, ]
selected_feature <- selected_combination$feature
selected_rna_set <- selected_combination$RNA_Set

cat(paste0(
  "\n\n Running model for feature: ",
  selected_feature,
  " and RNA set: ",
  selected_rna_set, "\n"
))

# Now select the data based on these choices
rna_data <- rna_list[[selected_rna_set]]

cat("\n\n RNA data: \n")
print(head(rna_data[, 1:5]))
cat("\n\n R data: \n")
print(head(rna_data[, 1:5]))
cat("\n\n")

full_df <- merge(rna_data,
  chr_cnv,
  by = "row.names"
)
cat("\n\n full_df: \n")
print(head(full_df[, 1:5]))
cat("\n\n")

y <- as.integer(full_df[[selected_feature]])
X <- full_df %>% select(-c("Row.names", colnames(chr_cnv)))
cat("\n\n Predicotrs: \n")
# print(head(X[, 1:5]))

grid <- expand.grid(
  w1 = seq(0.50, 0.70, 0.05),
  w2 = seq(0.05, 0.06, 0.05),
  w3 = seq(0.50, 0.70, 0.05)
)
set.seed(101)
for (j in 1:nrow(grid)) { # nolint
  selected_weights <- c(grid$w1[j], grid$w2[j], grid$w3[j])

  cat(paste0(
    # "\t\t eta: ", selected_lr,
    # "\t\t gamma: ", grid$gam[j],
    # "\t\t depth: ", selected_depth,
    # "\t\t trees: ", selected_trees,
    # "\t\t child_weight: ", selected_min_child,
    "\t\t Loss Weight: ", selected_weights[1],
    "\t\t Norm Weight: ", selected_weights[2],
    "\t\t Gain Weight: ", selected_weights[3],
    "\n"
  ))

  # Function to map factor levels to weights
  feature_digit_function <- function(factors) {
    sapply(factors, function(x) selected_weights[as.numeric(x)])
  }

  train_y_factor <- factor(y, levels = c(0, 1, 2))
  weights <- as.numeric(feature_digit_function(train_y_factor))

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)

  m_xgb_untuned <- xgb.cv(
    data = xgb_data,
    nrounds = selected_trees,
    objective = "multi:softprob",
    eval_metric = "auc",
    early_stopping_rounds = 100,
    nfold = 5,
    max_depth = selected_depth,
    min_child_weight = selected_min_child,
    eta = selected_lr,
    gamma = 0,
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

  best_auc_train <- if (best_iteration > 0) {
    m_xgb_untuned$evaluation_log$train_auc_mean[best_iteration]
  } else {
    NA # Or appropriate default/error value
  }

  best_auc_test <- if (best_iteration > 0) {
    m_xgb_untuned$evaluation_log$test_auc_mean[best_iteration]
  } else {
    NA # Or appropriate default/error value
  }


  cat(paste0(
    "The best iteration occurs with tree #: ",
    best_iteration, "\n\n"
  ))


  aneu_cat_metrics_df <- rbind(aneu_cat_metrics_df, data.frame(
    RNA_Set = selected_rna_set,
    Trees = selected_trees,
    Feature = selected_feature,
    Depth = selected_depth,
    Child_weight = selected_min_child,
    Learning_Rate = selected_lr,
    Gamma = 0,
    Weight_loss = selected_weights[1],
    Weight_normal = selected_weights[2],
    Weight_gain = selected_weights[3],
    Trained_AUC = best_auc_train,
    Test_AUC = best_auc_test
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

write.csv(
  aneu_cat_metrics_df,
  file = name,
  row.names = FALSE
)

cat(paste0("\n Completed processing for index: ", index, "\n"))
