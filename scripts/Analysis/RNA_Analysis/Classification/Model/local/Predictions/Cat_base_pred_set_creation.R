# This script is inteneded to test out the stacking of the models and the generated predictions for the categorical models only


# Load the packages
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

# Load the functions
# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) selected_weights[as.numeric(x)])
}

# Set the constant variables
# The RNA list
rna_list <- list(
  transcripts_per_million = "tpm",
  scaled_transcripts_per_million = "scld_tpm",
  log_scaled_transcripts_per_million = "log_scld_tpm",
  log_transcripts_per_million = "log_tpm",
  expected_counts = "exp",
  scaled_expected_counts = "scld",
  log_expected_counts = "log",
  log_scaled_expected_counts = "log_scld"
)

# The path to the RNA sets
rna_data_path <- "Data/RNA_Data/Model_Input/Train/"

# The CIN response features
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

# Create the folds to be used

full_data <- merge(read.csv(
    paste0(
      rna_data_path,
      "tpm",
      "_soi.csv"
    ),
    row.names = 1
  ),
  chr_cnv,
  by = "row.names") 

folds <- createFolds(full_data[["1p"]], k = 5, list = TRUE, returnTrain = TRUE)

# Get the parameters from stored Parameters file
parameters <- read.csv("Data/CIN_Features/XGB_cat_parameters.csv")

for (index in 1:length(rownames(parameters))) {
  # select the parameters and weights corrsponding to the index
  selected_parameters <- parameters[index, ]

  selected_feature <- selected_parameters$Feature
  selected_rna_set <- selected_parameters$RNA_set
  selected_trees <- as.numeric(selected_parameters$Trees)
  selected_eta <- selected_parameters$Eta
  selected_gamma <- selected_parameters$Gamma
  selected_max_depth <- selected_parameters$Max_depth
  selected_weights <- as.numeric(selected_parameters[c("Weight_loss", "Weight_normal", "Weight_gain")])
  selected_min_child <- 1
  selected_seed <- 99

  cat(
    "\n Selected Parameters: \n",
    "Feature: ", selected_feature, "\n",
    "RNA Set: ", selected_rna_set, "\n",
    "Trees: ", selected_trees, "\n",
    "Eta: ", selected_eta, "\n",
    "Gamma: ", selected_gamma, "\n",
    "Max Depth: ", selected_max_depth, "\n",
    "Weights: ", selected_weights, "\n",
    "Min Child: ", selected_min_child, "\n",
    "Seed: ", selected_seed, "\n"
  )

  rna_selection_name <- rna_list[[selected_rna_set]]
  rna_set <- read.csv(
    paste0(
      rna_data_path,
      rna_selection_name,
      "_soi.csv"
    ),
    row.names = 1
  )

  print("")
  print(head(rna_set[, 1:5]))

  full_df <- merge(rna_set,
    chr_cnv,
    by = "row.names"
  )

  y <- as.integer(full_df[[selected_feature]])
  X <- full_df %>% select(-c("Row.names", colnames(chr_cnv)))

  train_y_factor <- factor(y, levels = c(0, 1, 2))
  weights <- as.numeric(feature_digit_function(train_y_factor))

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)

  xgb_model <- xgb.cv(
    data = xgb_data,
    nrounds = selected_trees,
    objective = "multi:softmax",
    eval_metric = "mlogloss",
    early_stopping_rounds = 50,
    nfold = 2,
    folds = folds,
    max_depth = selected_max_depth,
    min_child_weight = selected_min_child,
    eta = selected_eta,
    gamma = selected_gamma,
    num_class = 3,
    stratified = TRUE,
    print_every_n = 10,
    prediction = TRUE
  )
}
