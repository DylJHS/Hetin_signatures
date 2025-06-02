# General template script for the hyperparameter tuning for the categorical learners of the base layer
# Uses a predefined parameters tabel with the estimated parameters that will be used by the feauture learner
# But also helps refine the parameters table which is the output in order to best select the hyperparamters in terms of the resulting evaluation metrics 
# in case the predefined parameter table does not exist, the intial parameters can be set manually
# The index should be an integer from 1 - 39
# But can be run from the R console by setting the index manually

library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

# Set the index for the learner to tune for from the command line (should be an integer from 1 - 39)
args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

# Set the working directory
setwd("")

# Get the parameters from stored Parameters file
parameters <- read.csv("XGB_cat_parameters.csv")
cat("\n Parameters: \n")
print(dim(parameters))

# select the parameters and weights corrsponding to the index' feature 
selected_parameters <- parameters[index, ]
rm(parameters)
cat("\n Parameters: \n")
print(selected_parameters)
cat("\n\n")

# The parameters are initialised based on the hyperparameter file but can be modified manually here
selected_feature <- selected_parameters$Feature
selected_rna_set <- selected_parameters$RNA_set
selected_trees <- as.numeric(selected_parameters$Trees)
selected_eta <- selected_parameters$Eta
selected_gamma <- selected_parameters$Gamma
selected_depth <- selected_parameters$Max_depth
selected_weights <- as.numeric(selected_parameters[c("Weight_loss", "Weight_normal", "Weight_gain")])

rm(selected_parameters)

# Define the extra parameters
selected_min_child <- 1
selected_seed <- 99

# Get the corresonding rna set
rna_list <- list(
  transcripts_per_million = "tpm",
  scaled_transcripts_per_million = "scld_tpm",
  log_scaled_transcripts_per_million = "log_scld_tpm",
  log_transcripts_per_million = "log_tpm",
  expected_counts = "exp",
  scaled_expected_counts = "scld_exp",
  log_expected_counts = "log_exp",
  log_scaled_expected_counts = "log_scld_exp"
)

# Get the corresonding rna set for predictor data
rna_selection_name <- rna_list[[selected_rna_set]]

rna_data_path <- "train_"

rna_set <- read.csv(
  paste0(
    rna_data_path,
    rna_selection_name,
    "_soi.csv"
  ),
  row.names = 1
)

# Arm Level Aneuploidies
# Load the data
chr_cnv <- read_tsv(
  "PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
) %>%
  replace(is.na(.), 0) %>%
  select(-c("Type", "Aneuploidy Score")) %>%
  mutate(Sample = str_replace_all(Sample, "-", "\\.")) %>%
  column_to_rownames("Sample") %>%
  mutate_all(~ replace(., . == 1, 2)) %>%
  mutate_all(~ replace(., . == 0, 1)) %>%
  mutate_all(~ replace(., . == -1, 0))

cat("\n\n CNV size: \n", dim(chr_cnv))

# MODELLING

aneu_cat_metrics_df <- data.frame(
  RNA_Set = character(),
  Trees = numeric(),
  Feature = character(),
  Depth = numeric(),
  Child_weight = numeric(),
  Eta = numeric(),
  Gamma = numeric(),
  Weight_loss = numeric(),
  Weight_normal = numeric(),
  Weight_gain = numeric(),
  Trained_mlogloss = numeric(),
  Test_mlogloss = numeric(),
  Seed = numeric()
)

# Determine the class weights for the target feature
cat("\n selected weights: ")
print(selected_weights)
cat("\n\n")

# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) selected_weights[as.numeric(x)])
}

full_df <- merge(rna_set,
  chr_cnv,
  by = "row.names"
)
rm(rna_set)
cat("\n\n full_df: \n")
print(head(full_df[, 1:5]))
cat("\n\n")

# Define the target labels
y <- as.integer(full_df[[selected_feature]])
# Define the predictors
X <- full_df %>% select(-c("Row.names", colnames(chr_cnv)))
cat("\n\n Predicotrs: \n")


# Set the weights for the classes
train_y_factor <- factor(y, levels = c(0, 1, 2))
weights <- as.numeric(feature_digit_function(train_y_factor))
rm(train_y_factor)

# Create the DMatrix
xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)
rm(weights)
rm(X)
rm(y)

# Define the grid for the hyperparameter tuning
# The extent of the grid can be modified to include more values
grid <- expand.grid(
  min_child = seq(1, 10, 1),
  eta = seq(0.01, 0.3, 0.05),
  gamma = seq(0, 0.5, 0.1),
  depth = seq(3, 10, 1)
)

# Iterate over the grid
for (j in 1:nrow(grid)) {
  for (param in names(grid)) {
    assign(paste0("selected_", param), grid[j, param])
  }

  set.seed(selected_seed)

  cat(paste0(
    "\t\t eta: ", selected_eta,
    "\t\t gamma: ", selected_gamma,
    "\t\t depth: ", selected_depth,
    "\t\t trees: ", selected_trees,
    "\t\t child_weight: ", selected_min_child,
    "\n"
  ))

  # Train the model with the selected hyperparameters and evaluate
  m_xgb_untuned <- xgb.cv(
    data = xgb_data,
    nrounds = selected_trees,
    objective = "multi:softmax",
    eval_metric = "mlogloss",
    early_stopping_rounds = 150,
    nfold = 10,
    max_depth = selected_depth,
    min_child_weight = selected_min_child,
    eta = selected_eta,
    gamma = selected_gamma,
    num_class = 3,
    stratified = TRUE,
    print_every_n = 20
  )

  # initialise the best iteration
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

    # Add the hyperparameters and metrics from the best iteration for the current grid iteration to the dataframe
    aneu_cat_metrics_df <- rbind(aneu_cat_metrics_df, data.frame(
      RNA_Set = selected_rna_set,
      Trees = best_iteration,
      Feature = selected_feature,
      Depth = selected_depth,
      Child_weight = selected_min_child,
      Eta = selected_eta,
      Gamma = selected_gamma,
      Weight_loss = selected_weights[1],
      Weight_normal = selected_weights[2],
      Weight_gain = selected_weights[3],
      Trained_mlogloss = best_mlogloss_train,
      Test_mlogloss = best_mlogloss_test,
      Seed = selected_seed
    ))
  }
}

name <- paste0(
  selected_feature, "_",
  index, "_",
   ".csv"
) %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "_")

# Save the metrics to a file
write.csv(
  aneu_cat_metrics_df,
  file = name,
  row.names = FALSE
)

cat(paste0("\n Completed processing for index: ", index, "\n"))
