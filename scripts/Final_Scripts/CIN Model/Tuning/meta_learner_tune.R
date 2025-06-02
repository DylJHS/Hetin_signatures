# This script is the third part of the meta model prediction pipeline.
# Its goal is to tune the hyperparameters for the meta learner models.
# Uses a predefined parameters table with the estimated parameters that will be used by the feauture learner
# But also helps refine the parameters table which is the output in order to best select the hyperparamters
# in case the predefined parameter table does not exist, the intial parameters can be set manually
# This file should be run from the command line using parallel processing specifying the index of the features to tune for
# The index should be an integer from 1 - 56
# But can be run from the R console by setting the index manually

# Load the packages
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)

# Set the index for the learner to tune for from the command line
args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

# Set the working directory
setwd("")

# Get the parameters from stored Parameters file
hyperparameters <- read.csv("meta_hyperparams.csv")

# Define the constant hyperparameters
selected_parameters <- hyperparameters[index, ]
cat("\n\n Selected Hyperparameters: ")
print(selected_parameters)

# The parameters are initialised based on the hyperparameter file but can be modified manually here
selected_feature <- selected_parameters$Feature
selected_trees <- selected_parameters$Trees + 500
selected_eta <- selected_parameters$Eta
selected_gamma <- selected_parameters$Gamma
selected_max_depth <- selected_parameters$Max_depth
selected_min_child <- selected_parameters$Child_weight
selected_seed <- 99

print(selected_max_depth)

# Load the full base prediction data for predictor data for the meta learners
predictions <- read.csv("Full_base_predictions.csv",
  row.names = 1
)

# Define the feature labels
feature_labels <- grep("act", names(predictions), value = TRUE)

# Define the feature names
response_features <- gsub("act_", "", feature_labels)

# Define the categorical features
cat_features <- response_features[grep("q$|p$", response_features)]

# Define the continuous regression features
reg_features <- setdiff(response_features, cat_features)

# Create the empty df that will contain the hyperparameters and the associated results
meta_model_metrics <- data.frame(
  "Type" = character(),
  "Feature" = character(),
  "Trees" = integer(),
  "Eta" = numeric(),
  "Gamma" = numeric(),
  "Max_depth" = integer(),
  "Train_result" = numeric(),
  "Test_result" = numeric(),
  "Seed" = integer()
)

set.seed(selected_seed)

# Define the hyperparameter grid
# The extent of the grid can be modified to include/exclude more values
hyper_grid <- expand.grid(
  depth = seq(1, 6, 1),
  min_child = seq(selected_min_child-1, selected_min_child+2, 1)
  eta = seq(0.02, 0.2, 0.03),
  gamma = seq(0, 0.5, 0.05)
)

# Start the model tuning based on the selected feature
# Set the predictor data
X <- predictions %>%
  select(-c("SampleID", feature_labels))

# Check if the selected feature learner is optimised for categorical or continuous tasks
if (selected_feature %in% cat_features) { # if the feature is categorical
  # Define the type
  selected_type <- "Categorical"

  # Loop over the hyperparameter grid
  for (i in 1:nrow(hyper_grid)) {
    # Get the hyperparameters
    for (param in names(hyper_grid)) {
      assign(paste0("selected_", param), hyper_grid[i, param])
    }

    # Print the selected hyperparameters
    cat(
      "\n\n\t Selected Parameters: \n",
      "Type: ", selected_type, "\n",
      "Feature: ", selected_feature, "\n",
      "Trees: ", selected_trees, "\n",
      "Eta: ", selected_eta, "\n",
      "Gamma: ", selected_gamma, "\n",
      "Max Depth: ", selected_max_depth, "\n",
      "Min Child Weight: ", selected_min_child, "\n",
      "Seed: ", selected_seed, "\n\n"
    )

    # Set the label target data
    y <- as.numeric(
      predictions[[paste0("act_", selected_feature)]]
    )

    # Create the xgb data
    xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

    # Set common parameters
    xgb_params <- list(
      data = xgb_data,
      nrounds = selected_trees,
      nfold = 10,
      early_stopping_rounds = 150,
      objective = if (selected_feature %in% cat_features) "multi:softmax" else "reg:squarederror", # Set the objective based on the learner type
      eval_metric = if (selected_feature %in% cat_features) "mlogloss" else "rmse", # Set the evaluation metric based on the learner type
      eta = selected_eta,
      gamma = selected_gamma,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      stratified = TRUE,
      print_every_n = 50,
      prediction = TRUE
    )
    
    # Add num_class only if it's a classification task
    if (selected_feature %in% cat_features) {
      xgb_params$num_class <- 3
    }
    
    # Train model using the parameters list
    xgb_model <- do.call(xgb.cv, xgb_params)

    # initialise the best iteration
    best_iteration <- 0

    # First, check if best_iteration is valid
    if (is.null(
      xgb_model$best_iteration
    ) ||
      xgb_model$best_iteration < 1) {
      cat(paste0(
        "Warning: No valid best_iteration found.",
        " Using last iteration values instead.\n"
      ))
      # Use the last iteration if best_iteration is not valid
      best_iteration <- nrow(xgb_model$evaluation_log)
    } else {
      # Ensure that the best_iteration does not exceed the number of rows logged
      if (xgb_model$best_iteration > nrow(xgb_model$evaluation_log)) {
        cat(paste0(
          "Warning: best_iteration exceeds the number of rows in evaluation_log.",
          " Adjusting to maximum available.\n"
        ))
        best_iteration <- nrow(xgb_model$evaluation_log)
      } else {
        best_iteration <- xgb_model$best_iteration
      }
    }

    best_mlogloss_train <- if (best_iteration > 0) {
      xgb_model$evaluation_log$train_mlogloss_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    best_mlogloss_test <- if (best_iteration > 0) {
      xgb_model$evaluation_log$test_mlogloss_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    cat(paste0(
      "The best iteration occurs with tree #: ",
      best_iteration, "\n\n"
    ))

    # Append the results to the parameters df
    meta_model_metrics <- rbind(
      meta_model_metrics,
      data.frame(
        "Type" = selected_type,
        "Feature" = selected_feature,
        "Trees" = best_iteration,
        "Eta" = selected_eta,
        "Gamma" = selected_gamma,
        "Max_depth" = selected_max_depth,
        "Child_weight" = selected_min_child,
        "Train_result" = best_mlogloss_train,
        "Test_result" = best_mlogloss_test,
        "Seed" = selected_seed
      )
    )
  }

  # Else if the learner task is continuous
} else if (selected_feature %in% reg_features) {
  # Define the type
  selected_type <- "Regression"

  # Loop over the hyperparameter grid
  for (i in 1:nrow(hyper_grid)) {
    # Get the hyperparameters
    for (param in names(hyper_grid)) {
      assign(paste0("selected_", param), hyper_grid[i, param])
    }

    # Print the selected hyperparameters
    cat(
      "\n\n\t Selected Parameters: \n",
      "Type: ", selected_type, "\n",
      "Feature: ", selected_feature, "\n",
      "Trees: ", selected_trees, "\n",
      "Eta: ", selected_eta, "\n",
      "Gamma: ", selected_gamma, "\n",
      "Max Depth: ", selected_max_depth, "\n",
      "Min Child Weight: ", selected_min_child, "\n",
      "Seed: ", selected_seed, "\n\n"
    )

    # Set the label
    y <- as.numeric(
      predictions[[paste0("act_", selected_feature)]]
    )

    # Create the xgb data
    xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

    # Train the model
    xgb_model <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      nfold = 5,
      early_stopping_rounds = 100,
      objective = "reg:squarederror",
      eval_metric = "rmse",
      eta = selected_eta,
      gamma = selected_gamma,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      print_every_n = 25,
      prediction = TRUE
    )

    # initialise the best iteration
    best_iteration <- 0

    # First, check if best_iteration is valid
    if (is.null(
      xgb_model$best_iteration
    ) ||
      xgb_model$best_iteration < 1) {
      cat(paste0(
        "Warning: No valid best_iteration found.",
        " Using last iteration values instead.\n"
      ))
      # Use the last iteration if best_iteration is not valid
      best_iteration <- nrow(xgb_model$evaluation_log)
    } else {
      # Ensure that the best_iteration does not exceed the number of rows logged
      if (xgb_model$best_iteration > nrow(xgb_model$evaluation_log)) {
        cat(paste0(
          "Warning: best_iteration exceeds the number of rows in evaluation_log.",
          " Adjusting to maximum available.\n"
        ))
        best_iteration <- nrow(xgb_model$evaluation_log)
      } else {
        best_iteration <- xgb_model$best_iteration
      }
    }

    best_rmse_train <- if (best_iteration > 0) {
      xgb_model$evaluation_log$train_rmse_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    best_rmse_test <- if (best_iteration > 0) {
      xgb_model$evaluation_log$test_rmse_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    cat(paste0(
      "The best iteration occurs with tree #: ",
      best_iteration, "\n\n"
    ))

    # Append the results to the parameters df
    meta_model_metrics <- rbind(
      meta_model_metrics,
      data.frame(
        "Type" = selected_type,
        "Feature" = selected_feature,
        "Trees" = best_iteration,
        "Eta" = selected_eta,
        "Gamma" = selected_gamma,
        "Max_depth" = selected_max_depth,
        "Child_weight" = selected_min_child,
        "Train_result" = best_rmse_train,
        "Test_result" = best_rmse_test,
        "Seed" = selected_seed
      )
    )
  }
}

# Save the results
write.csv(
  meta_model_metrics,
  paste0(
    "Hyperparams_",
    selected_feature, "_", ".csv"
  ),
  row.names = FALSE
)


cat(
  "\n\n The hyperparameter tuning for the meta learner is complete.\n\n"
)
