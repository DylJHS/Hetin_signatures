# This script is intended to train the base and meta models for the CIN prediction and save the predictions in their original scale

start_time <- Sys.time()

# Load the packages
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)

print_every <- 150
early_stop <- 150
fold_n <- 5

cat("\n Loading the data:  \n")
input_path <- "/hpc/shared/prekovic/dhaynessimmons/data/"

# Load the functions
# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) selected_weights[as.numeric(x)])
}

# Function to extract the fold indices from the folds
combine_all_folds <- function(folds) {
  all_indices <- c()
  for (fold in folds) {
    all_indices <- c(all_indices, fold)
  }
  return(all_indices)
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

# The CIN response features
# Categorical features
cat_cin <- read_tsv(
  paste0(
    input_path,
    "CIN/PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
  ),
  show_col_types = FALSE
) %>%
  replace(is.na(.), 0) %>%
  select(-c("Type", "Aneuploidy Score")) %>%
  mutate(Sample = str_replace_all(Sample, "-", "\\.")) %>%
  column_to_rownames("Sample") %>%
  mutate_all(~ replace(., . == 1, 2)) %>%
  mutate_all(~ replace(., . == 0, 1)) %>%
  mutate_all(~ replace(., . == -1, 0)) %>%
  rename(
    "13q" = "13 (13q)",
    "14q" = "14 (14q)",
    "15q" = "15 (15q)",
    "21q" = "21 (21q)",
    "22q" = "22 (22q)"
  )

cat("\n Categorical CIN:  \n")
print(head(cat_cin[20:30]))

# Numerical features
# HRD scores
ori_hrd <- read_tsv(
  paste0(input_path, "CIN/TCGA.HRD_withSampleID.txt"),
  show_col_types = FALSE
)

t_hrd <- as.data.frame(t(ori_hrd))
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ]
hrd <- as.data.frame(first_hrd[-1, ]) %>%
  mutate_all(as.numeric) %>%
  rename(loh_hrd = "hrd-loh") %>%
  select(-HRD) %>%
  mutate(new = str_replace_all(rownames(.), "-", "\\."))

rm(first_hrd)
rm(ori_hrd)
rm(t_hrd)

rownames(hrd) <- hrd$new
hrd <- hrd %>%
  select(-new)

# print(head(hrd))

# Pericentromeric CNVs
peri_cnv <- read.csv(
  paste0(
    input_path,
    "CIN/lim_alpha_incl_TCGA_pericentro_cnv.csv"
  )
) %>%
  mutate_all(~ replace(., is.na(.), 0)) %>%
  mutate(sampleID = gsub("-", ".", sampleID)) %>%
  column_to_rownames("sampleID")

# print(dim(peri_cnv))

base_cat_parameters <- read.csv(
  paste0(
    input_path,
    "hyperparameters/base_cat_hyperparams.csv"
  )
)
base_reg_parameters <- read.csv(
  paste0(
    input_path,
    "hyperparameters/base_reg_hyperparams.csv"
  )
)

reg_cin <- merge(
  hrd,
  peri_cnv,
  by = "row.names"
) %>%
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")
cat("\n Regression CIN:  \n")
print(head(reg_cin[1:5]))

rm(peri_cnv)
rm(hrd)

full_cin <- merge(
  cat_cin,
  reg_cin,
  by = "row.names"
) %>%
  column_to_rownames("Row.names")
cat("\n Full CIN:  \n")
print(head(full_cin[1:5]))

# All response feature names
response_features <- colnames(full_cin)
reg_features <- colnames(reg_cin)
cat_features <- colnames(cat_cin)
cat(
  "\n Response Features # :", length(response_features), "\n",
  "Regression Features # :", length(reg_features), "\n",
  "Categorical Features # :", length(cat_features), "\n"
)

rm(reg_cin)
rm(cat_cin)

# Create the folds to be used
# Create the full data by merging the soi & cancerous RNA data with the CIN data
full_data <- merge(
  read.csv(
    paste0(
      input_path,
      "mRNA/gbm_input/Train/train_tpm_soi.csv"
    ),
    row.names = 1
  ),
  full_cin,
  by = "row.names"
)

cat("\n Full Data:  \n")
print(head(full_data[1:5]))

# Creating the folds and returning the indices for the out-of-fold predictions only
folds <- createFolds(full_data[["1p"]], k = fold_n, list = TRUE, returnTrain = FALSE)
cat("\n Folds # :", length(folds), "\n")

# Create the empty dataframe that will store the out-of-fold predictions for each model
base_oof_predictions <- data.frame(
  matrix(
    ncol = length(response_features),
    nrow = length(combine_all_folds(folds))
  )
)
colnames(base_oof_predictions) <- paste0("pred_", response_features)

cat("\n Predictions Dataframe:  \n")
print(head(base_oof_predictions[1:5]))

# Get the actual (true) fold labels
labels <- full_data %>%
  select(all_of(intersect(response_features, colnames(.)))) %>%
  mutate(
    index = as.factor(row_number())
  ) %>%
  rename_with(~ paste0("act_", .))

# merge the fold labels with the empty predictions dataframe
# making sure to keep the order of the original indices
base_oof_predictions <- cbind(base_oof_predictions, labels) %>%
  arrange(act_index)

for (feature in response_features) {
  cat(
    "Feature: ", feature, "\n\n"
  )
  if (feature %in% cat_features) {
    # Get the parameters from stored Parameters file
    parameters <- base_cat_parameters

    # select the parameters and weights corresponding to the index
    selected_parameters <- parameters[parameters$Feature == feature, ]

    selected_feature <- selected_parameters$Feature
    selected_rna_set <- selected_parameters$RNA_set
    selected_trees <- as.numeric(selected_parameters$Trees)
    selected_eta <- selected_parameters$Eta
    selected_gamma <- selected_parameters$Gamma
    selected_max_depth <- selected_parameters$Max_depth
    selected_min_child <- selected_parameters$Child_weight
    selected_weights <- as.numeric(
      selected_parameters[c(
        "Weight_loss",
        "Weight_normal",
        "Weight_gain"
      )]
    )

    cat(
      "\n Selected Parameters: \n",
      "Feature: ", selected_feature, "\t",
      "RNA Set: ", selected_rna_set, "\t",
      "Trees: ", selected_trees, "\n",
      "Eta: ", selected_eta, "\t",
      "Gamma: ", selected_gamma, "\t",
      "Max Depth: ", selected_max_depth, "\t",
      "Weights: ", selected_weights, "\t",
      "Min Child: ", selected_min_child, "\n"
    )

    rna_selection_name <- rna_list[[selected_rna_set]]
    rna_set <- read.csv(
      paste0(
        input_path,
        "mRNA/gbm_input/Train/train_",
        rna_selection_name,
        "_soi.csv"
      ),
      row.names = 1
    )

    full_df <- merge(rna_set,
      full_cin,
      by = "row.names"
    )

    y <- as.integer(full_df[[selected_feature]])
    X <- full_df %>% select(-c("Row.names", all_of(response_features)))

    train_y_factor <- factor(y, levels = c(0, 1, 2))
    weights <- as.numeric(feature_digit_function(train_y_factor))

    xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)

    xgb_model <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      objective = "multi:softmax",
      eval_metric = "mlogloss",
      folds = folds,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      early_stopping_rounds = early_stop,
      eta = selected_eta,
      gamma = selected_gamma,
      num_class = 3,
      print_every_n = print_every,
      prediction = TRUE
    )

    # Store the predictions in the corresponding column
    base_oof_predictions[[paste0("pred_", selected_feature)]] <- xgb_model$pred[, 1]
  } else if (feature %in% reg_features) {
    # Get the parameters from stored Parameters file
    parameters <- base_reg_parameters

    # select the parameters and weights corresponding to the index
    selected_parameters <- parameters[parameters$Feature == feature, ]

    selected_feature <- selected_parameters$Feature
    selected_rna_set <- selected_parameters$RNA_set
    selected_trees <- as.numeric(selected_parameters$Trees)
    selected_eta <- selected_parameters$Eta
    selected_gamma <- selected_parameters$Gamma
    selected_max_depth <- selected_parameters$Max_depth
    selected_min_child <- selected_parameters$Child_weight

    cat(
      "\n Selected Parameters: \n",
      "Feature: ", selected_feature, "\t",
      "RNA Set: ", selected_rna_set, "\t",
      "Trees: ", selected_trees, "\n",
      "Eta: ", selected_eta, "\t",
      "Gamma: ", selected_gamma, "\t",
      "Max Depth: ", selected_max_depth, "\t",
      "Min Child: ", selected_min_child, "\n"
    )

    rna_selection_name <- rna_list[[selected_rna_set]]
    rna_set <- read.csv(
      paste0(
        input_path,
        "mRNA/gbm_input/Train/train_",
        rna_selection_name,
        "_soi.csv"
      ),
      row.names = 1
    )

    full_df <- merge(rna_set,
      full_cin,
      by = "row.names"
    )

    y <- as.numeric(full_df[[selected_feature]])
    X <- full_df %>% select(-c("Row.names", all_of(response_features)))

    xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

    xgb_model <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      objective = "reg:squarederror",
      eval_metric = "rmse",
      early_stopping_rounds = early_stop,
      folds = folds,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      eta = selected_eta,
      gamma = selected_gamma,
      print_every_n = print_every,
      prediction = TRUE
    )

    # Store the predictions in the corresponding column
    base_oof_predictions[[paste0("pred_", selected_feature)]] <- xgb_model$pred
  } else {
    cat("Feature not found")
  }
}
rm(base_cat_parameters)
rm(base_reg_parameters)
rm(rna_selection_name)

base_time <- Sys.time() - start_time
cat("Training the base models completed in :" base_time)

write.csv(
  base_oof_predictions,
  "/hpc/shared/prekovic/dhaynessimmons/data/base_predictions/full_base_predictions.csv"
)

# Reconvert the data into the original scale
final_base_oof_predictions <- base_oof_predictions %>%
  mutate_all(as.numeric) %>%
  mutate_all(~ replace(., . == 0, -1)) %>%
  mutate_all(~ replace(., . == 1, 0)) %>%
  mutate_all(~ replace(., . == 2, 1))

# Save the predictions
write.csv(
  final_base_oof_predictions,
  "hpc/shared/prekovic/dhaynessimmons/data/model_output/Final_base_predictions.csv"
)
rm(final_base_oof_predictions)

# Load the meta parameters and remove the base ones
parameters <- read.csv(
  paste0(
    input_path,
    "hyperparameters/meta_hyperparams.csv"
  )
)

# Create the training folds
meta_folds <- createFolds(base_oof_predictions[["act_1p"]], k = fold_n, list = TRUE, returnTrain = FALSE)

# Create the empty dataframe that will store the out-of-fold predictions for each model of the meta learners
meta_oof_predictions <- data.frame(
  matrix(
    ncol = length(response_features),
    nrow = length(combine_all_folds(meta_folds))
  )
)
colnames(meta_oof_predictions) <- paste0("pred_", response_features)

cat("\n Predictions Dataframe:  \n")
print(head(meta_oof_predictions[1:5]))

meta_oof_predictions <- cbind(
  meta_oof_predictions,
  labels
) %>%
  arrange(act_index) %>%
  column_to_rownames("act_index")
  
rm(labels)

# Create the predictor data for the meta learners
meta_input_data <- base_oof_predictions %>%
  select(starts_with("pred_"))

for (feature in response_features) {
  cat(
    "Feature: ", feature, "\n\n"
  )

  # select the parameters and weights corresponding to the index
  selected_parameters <- parameters[parameters$Feature == feature, ]

  selected_feature <- selected_parameters$Feature
  selected_trees <- as.numeric(selected_parameters$Trees)
  selected_eta <- selected_parameters$Eta
  selected_gamma <- selected_parameters$Gamma
  selected_max_depth <- selected_parameters$Max_depth
  selected_min_child <- selected_parameters$Child_weight

  cat(
    "\n Selected Parameters: \n",
    "Feature: ", selected_feature, "\t",
    "Trees: ", selected_trees, "\t",
    "Eta: ", selected_eta, "\n",
    "Gamma: ", selected_gamma, "\t",
    "Max Depth: ", selected_max_depth, "\t",
    "Min Child: ", selected_min_child, "\n"
  )

  y <- as.integer(base_oof_predictions[[paste0("act_", feature)]])

  xgb_data <- xgb.DMatrix(data = as.matrix(meta_input_data), label = y)

  if (feature %in% cat_features) {
    xgb_meta_model <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      objective = "multi:softmax",
      eval_metric = "mlogloss",
      folds = meta_folds,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      early_stopping_rounds = early_stop,
      eta = selected_eta,
      gamma = selected_gamma,
      num_class = 3,
      print_every_n = print_every,
      prediction = TRUE
    )

    # Store the predictions in the corresponding column
    meta_oof_predictions[[paste0("pred_", selected_feature)]] <- xgb_meta_model$pred[, 1]
  } else if (feature %in% reg_features) {
    xgb_meta_model <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      objective = "reg:squarederror",
      eval_metric = "rmse",
      early_stopping_rounds = early_stop,
      folds = meta_folds,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      eta = selected_eta,
      gamma = selected_gamma,
      print_every_n = print_every,
      prediction = TRUE
    )

    # Store the predictions in the corresponding column
    meta_oof_predictions[[paste0("pred_", selected_feature)]] <- xgb_meta_model$pred
  } else {
    cat("Feature not found")
  }
}

meta_time <- Sys.time() - start_time
cat("Training the meta models completed in :" meta_time)

final_meta_oof_predictions <- meta_oof_predictions %>%
  mutate_all(as.numeric) %>%
  mutate_all(~ replace(., . == 0, -1)) %>%
  mutate_all(~ replace(., . == 1, 0)) %>%
  mutate_all(~ replace(., . == 2, 1))

write.csv(
  final_meta_oof_predictions,
  "/hpc/shared/prekovic/dhaynessimmons/data/model_output/Final_meta_predictions.csv"
)
