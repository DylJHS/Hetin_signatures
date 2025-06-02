# This script is the third part of the cancer-specific model prediction pipeline.
# Its goal is to tune the hyperparameters for the meta learner models.

# Load the packages
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # for the feature index


cat("\n\n\t\t Index: ", index, "\n")

all_features <- c(
  "1p", "1q", "2p", "2q", "3p",
  "3q", "4p", "4q", "5p", "5q",
  "6p", "6q", "7p", "7q", "8p",
  "8q", "9p", "9q", "10p", "10q",
  "11p", "11q", "12p", "12q",
  "13q", "14q", "15q", "16p",
  "16q", "17p", "17q", "18p",
  "18q", "19p", "19q", "20p",
  "20q", "21q", "22q", "ai1",
  "lst1", "loh_hrd", "peri_22_3",
  "peri_1", "peri_2", "peri_3",
  "peri_4", "peri_5", "peri_7",
  "peri_8", "peri_9", "peri_10",
  "peri_16_1", "peri_16_2", "peri_17",
  "peri_20"
)

cancer_types <- c(
  "BLCA", "BRCA", "CESC",
  "HNSC", "LGG", "LIHC",
  "LUAD", "LUSC", "OV",
  "PRAD", "STAD", "THCA"
)

feature <- all_features[index]

# Load the functions
# Function to map factor levels to weights
feature_digit_function <- function(factors, reference) {
  sapply(factors, function(x) reference[as.numeric(x)])
}

# Get the parameters from stored Parameters file
hyper_param_folder <- file.path("/hpc/shared/prekovic/dhaynessimmons/data/Cancer_specific_data/Model_input/Parameters/Hyperparameters/")

for (cancer in cancer_types) {
  selected_cancer <- cancer
  cat("\n\n\t\t\t\t\t\t\t\t\t\t\t\t Cancer: ", selected_cancer, "\n",
      "Feature: ", feature, "\n")

  # Define the hyperparameters to tune
  hyperparameter_file <- paste0(
    hyper_param_folder,
    selected_cancer,
    "/",
    "meta_hyperparameters.csv"
  )

  # Load the hyperparameters if they exist
  if (!file.exists(hyperparameter_file)) {
    selected_feature <- feature
    selected_trees <- 10000
    selected_eta <- 0.3
    selected_gamma <- 0
    selected_max_depth <- 5
    selected_min_child <- 2
    selected_seed <- 99
  } else {
    hyperparameters <- read.csv(hyperparameter_file, header = TRUE, stringsAsFactors = FALSE)

    # Check that the feature is part of the hyperparameters
    if (!(feature %in% hyperparameters$Feature)) {
      cat("\n\n The feature is not in the hyperparameters. \n Using default intial params \n")
      selected_feature <- feature
      selected_trees <- 10000
      selected_eta <- 0.3
      selected_gamma <- 0
      selected_max_depth <- 5
      selected_min_child <- 2
      selected_seed <- 99
    } else {
      # Define the constant hyperparameters
      selected_parameters <- hyperparameters[hyperparameters$Feature == feature, ]

      selected_feature <- selected_parameters$Feature
      selected_trees <- selected_parameters$Trees + 500
      selected_eta <- ifelse(selected_trees > 88000, selected_parameters$Treesselected_parameters$Trees * 10,
                             ifelse(selected_trees > 44000,  selected_parameters$Treesselected_parameters$Trees * 5,
                                    ifelse(selected_trees > 32000,  selected_parameters$Treesselected_parameters$Trees * 3, 0.3)))
      selected_gamma <- selected_parameters$Gamma
      selected_max_depth <- selected_parameters$Max_depth
      selected_min_child <- selected_parameters$Child_weight
      selected_seed <- 99
    }
  }

  cat("Loaded the hyperparameters.\n")

  # Load the arm class weights
  arm_weights <- read.csv(
    paste0(
      "/hpc/shared/prekovic/dhaynessimmons/data/Cancer_specific_data/Model_input/Parameters/Arm_class_weights/",
      selected_cancer,
      "_arm_weights.csv"
    ),
    row.names = 1
  )
  cat("\n\n The arm weights have been loaded.\n")

  # Load the full base prediction data
  input_predictions_file <- file.path(
    "/hpc/shared/prekovic/dhaynessimmons/data/Cancer_specific_data/Model_output",
    selected_cancer,
    "Results/Predictions/Full_base_predictions.csv"
  )

  predictions <- read.csv(
    input_predictions_file,
    row.names = 1
  ) %>%
    arrange(act_index) %>%
    column_to_rownames("act_index")
  cat("\n\n The predictions data has been loaded.\n\n")
  print(head(predictions[1:3], 3))

  # Define the feature labels
  labels <- predictions %>%
    select(starts_with("act_")) %>%
    colnames(.)

  cat("\n\n The labels have been defined.\n")
  print(head(labels, 3))

  # Define the feature names
  response_features <- labels %>%
    gsub("act_", "", .)
  cat("\n\n The response features have been defined.\n")
  print(head(response_features, 3))

  # skip if the feature is not in the response features
  if (!(feature %in% response_features)) {
    cat("\n\n The feature is not in the response features. \n Skipping to next Cancer type \n")
    next
  }

  # Define the categorical features
  cat_features <- response_features[grep("q$|p$", response_features)]
  cat("\n\n The categorical features have been defined.\n")
  print(head(cat_features, 3))

  # Define the continuous features
  reg_features <- setdiff(response_features, cat_features)
  cat("\n\n The continuous features have been defined.\n")
  print(head(reg_features, 3))

  cat("\n\n Response Features: ")
  print(response_features)


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
  hyper_grid <- expand.grid(
    max_depth = seq(selected_max_depth - 1, selected_max_depth + 1, 1),
    min_child = seq(selected_min_child - 3, selected_min_child + 3, 1),
    eta = seq(selected_eta - 0.05, selected_eta, 0.05)
  )

  # Start the model tuning based on the selected feature
  # Set the data
  X <- predictions %>%
    select(-all_of(labels))

  # Loop over the hyperparameter grid
  for (i in 1:nrow(hyper_grid)) {
    # Get the hyperparameters
    for (param in names(hyper_grid)) {
      param_value <- hyper_grid[i, param]
      assign(paste0("selected_", param), param_value)

      # Ensure that the hyperparameters are not negative
      # If they are, set them to random value between 1 and 12
      if (param_value <= 0) {
        assign(paste0("selected_", param), sample(1:12, 1))
      }
    }

    # Check that the tuning has not been done before for the exact same parameters
    if (nrow(meta_model_metrics) > 0) {
      if (any(
        meta_model_metrics$Max_depth == selected_max_depth &
          meta_model_metrics$Child_weight == selected_min_child &
          meta_model_metrics$Eta == selected_eta &
          meta_model_metrics$Gamma == selected_gamma
      )) {
        cat("\n\t\t\t Hyperparameters already tuned. Skipping.\n\n")
        next
      }
    }

    # Print the selected hyperparameters
    cat(
      "\n\n\t Selected Parameters: \n",
      "Feature: ", selected_feature, "\n",
      "Trees: ", selected_trees, "\n",
      "Eta: ", selected_eta, "\n",
      "Gamma: ", selected_gamma, "\n",
      "Max Depth: ", selected_max_depth, "\n",
      "Min Child Weight: ", selected_min_child, "\n",
      "Seed: ", selected_seed, "\n\n"
    )

    # Set the label
    y <- ifelse(feature %in% cat_features, as.integer, as.numeric)(predictions[[paste0("act_", selected_feature)]])
    cat("\n\n Target: ", selected_feature, "\n")
    print(head(y))

    # Adjust weights for categorical features
    if (selected_feature %in% cat_features) {
      selected_weights <- arm_weights[selected_feature, ]
      cat("\n\n Selected Weights: \n")
      print(selected_weights)

      selected_ref <- as.numeric(
        selected_weights[c(
          "Weight_loss",
          "Weight_normal",
          "Weight_gain"
        )]
      )
      weights <- as.numeric(feature_digit_function(factor(y, levels = c(0, 1, 2)), selected_ref))
      # Create the xgb data
      xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)
      cat("\n Loaded xgb data with weights.\n")
    } else {
      # Create the xgb data
      cat("\n Loaded xgb data without weights.\n")
      xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)
    }

    # Set common parameters
    xgb_params <- list(
      data = xgb_data,
      nrounds = selected_trees,
      nfold = 10,
      early_stopping_rounds = 150,
      objective = if (selected_feature %in% cat_features) "multi:softmax" else "reg:squarederror",
      eval_metric = if (selected_feature %in% cat_features) "mlogloss" else "rmse",
      eta = selected_eta,
      gamma = selected_gamma,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      stratified = TRUE,
      print_every_n = 50,
      prediction = TRUE
    )

    # Add num_class only if it's a classification problem
    if (selected_feature %in% cat_features) {
      xgb_params$num_class <- 3
    }

    # Train model using the parameters list
    xgb_model <- do.call(xgb.cv, xgb_params)

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

    best_train <- if (best_iteration > 0) {
      if (selected_feature %in% cat_features) {
        xgb_model$evaluation_log$train_mlogloss_mean[best_iteration]
      } else {
        xgb_model$evaluation_log$train_rmse_mean[best_iteration]
      }
    } else {
      NA # Or appropriate default/error value
    }

    best_test <- if (best_iteration > 0) {
      if (selected_feature %in% cat_features) {
        xgb_model$evaluation_log$test_mlogloss_mean[best_iteration]
      } else {
        xgb_model$evaluation_log$test_rmse_mean[best_iteration]
      }
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
        "Feature" = selected_feature,
        "Trees" = best_iteration,
        "Eta" = selected_eta,
        "Gamma" = selected_gamma,
        "Max_depth" = selected_max_depth,
        "Child_weight" = selected_min_child,
        "Train_result" = best_train,
        "Test_result" = best_test,
        "Seed" = selected_seed
      )
    ) %>%
      arrange(Test_result)
  }

  # Save the hyperparameters
  output_file <- file.path(
    hyper_param_folder,
    selected_cancer,
    "individual_features"
  )

  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
  }

  unique_id <- Sys.time() %>%
    as.numeric() %>%
    round()

  write.csv(
    meta_model_metrics,
    paste0(
      output_file,
      "/",
      selected_feature,
      "_",
      unique_id,
      "_meta_hyperparams.csv"
    ),
    row.names = FALSE
  )
}

cat("Task Completed")
# This script is the third part of the cancer-specific model prediction pipeline.
# Its goal is to tune the hyperparameters for the meta learner models.

# Load the packages
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # for the feature index


cat("\n\n\t\t Index: ", index, "\n")

all_features <- c(
  "1p", "1q", "2p", "2q", "3p",
  "3q", "4p", "4q", "5p", "5q",
  "6p", "6q", "7p", "7q", "8p",
  "8q", "9p", "9q", "10p", "10q",
  "11p", "11q", "12p", "12q",
  "13q", "14q", "15q", "16p",
  "16q", "17p", "17q", "18p",
  "18q", "19p", "19q", "20p",
  "20q", "21q", "22q", "ai1",
  "lst1", "loh_hrd", "peri_22_3",
  "peri_1", "peri_2", "peri_3",
  "peri_4", "peri_5", "peri_7",
  "peri_8", "peri_9", "peri_10",
  "peri_16_1", "peri_16_2", "peri_17",
  "peri_20"
)

cancer_types <- c(
  "BLCA", "BRCA", "CESC",
  "HNSC", "LGG", "LIHC",
  "LUAD", "LUSC", "OV",
  "PRAD", "STAD", "THCA"
)

feature <- all_features[index]

# Load the functions
# Function to map factor levels to weights
feature_digit_function <- function(factors, reference) {
  sapply(factors, function(x) reference[as.numeric(x)])
}

# Get the parameters from stored Parameters file
hyper_param_folder <- file.path("/hpc/shared/prekovic/dhaynessimmons/data/Cancer_specific_data/Model_input/Parameters/Hyperparameters/")

for (cancer in cancer_types) {
  selected_cancer <- cancer
  cat("\n\n\t\t\t\t\t\t Cancer: ", selected_cancer, "\t",
      "Feature: ", feature, "\n")

  # Define the hyperparameters to tune
  hyperparameter_file <- paste0(
    hyper_param_folder,
    selected_cancer,
    "/",
    "meta_hyperparameters.csv"
  )

  # Load the hyperparameters if they exist
  if (!file.exists(hyperparameter_file)) {
    selected_feature <- feature
    selected_trees <- 10000
    selected_eta <- 0.3
    selected_gamma <- 0
    selected_max_depth <- 5
    selected_min_child <- 2
    selected_seed <- 99
  } else {
    hyperparameters <- read.csv(hyperparameter_file, header = TRUE, stringsAsFactors = FALSE)

    # Check that the feature is part of the hyperparameters
    if (!(feature %in% hyperparameters$Feature)) {
      cat("\n\n The feature is not in the hyperparameters. \n Using default intial params \n")
      selected_feature <- feature
      selected_trees <- 10000
      selected_eta <- 0.3
      selected_gamma <- 0
      selected_max_depth <- 5
      selected_min_child <- 2
      selected_seed <- 99
    } else {
      # Define the constant hyperparameters
      selected_parameters <- hyperparameters[hyperparameters$Feature == feature, ]

      selected_feature <- selected_parameters$Feature
      selected_trees <- selected_parameters$Trees + 500
      selected_eta <- ifelse(selected_trees > 88000, 3,
                             ifelse(selected_trees > 44000, 2,
                                    ifelse(selected_trees > 32000, 1, 0.3)))
      selected_gamma <- selected_parameters$Gamma
      selected_max_depth <- selected_parameters$Max_depth
      selected_min_child <- selected_parameters$Child_weight
      selected_seed <- 99
    }
  }

  cat("Loaded the hyperparameters.\n")

  # Load the arm class weights
  arm_weights <- read.csv(
    paste0(
      "/hpc/shared/prekovic/dhaynessimmons/data/Cancer_specific_data/Model_input/Parameters/Arm_class_weights/",
      selected_cancer,
      "_arm_weights.csv"
    ),
    row.names = 1
  )
  cat("\n\n The arm weights have been loaded.\n")

  # Load the full base prediction data
  input_predictions_file <- file.path(
    "/hpc/shared/prekovic/dhaynessimmons/data/Cancer_specific_data/Model_output",
    selected_cancer,
    "Results/Predictions/Full_base_predictions.csv"
  )

  predictions <- read.csv(
    input_predictions_file,
    row.names = 1
  ) %>%
    arrange(act_index) %>%
    column_to_rownames("act_index")
  cat("\n\n The predictions data has been loaded.\n\n")
  print(head(predictions[1:3], 3))

  # Define the feature labels
  labels <- predictions %>%
    select(starts_with("act_")) %>%
    colnames(.)

  cat("\n\n The labels have been defined.\n")
  print(head(labels, 3))

  # Define the feature names
  response_features <- labels %>%
    gsub("act_", "", .)
  cat("\n\n The response features have been defined.\n")
  print(head(response_features, 3))

  # skip if the feature is not in the response features
  if (!(feature %in% response_features)) {
    cat("\n\n The feature is not in the response features. \n Skipping to next Cancer type \n")
    next
  }

  # Define the categorical features
  cat_features <- response_features[grep("q$|p$", response_features)]
  cat("\n\n The categorical features have been defined.\n")
  print(head(cat_features, 3))

  # Define the continuous features
  reg_features <- setdiff(response_features, cat_features)
  cat("\n\n The continuous features have been defined.\n")
  print(head(reg_features, 3))

  cat("\n\n Response Features: ")
  print(response_features)


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
  hyper_grid <- expand.grid(
    max_depth = seq(selected_max_depth - 2, selected_max_depth + 2, 1),
    min_child = seq(selected_min_child - 4, selected_min_child + 4, 1),
    eta = seq(selected_eta - 0.1, selected_eta + 0.1, 0.1)
  )

  # Start the model tuning based on the selected feature
  # Set the data
  X <- predictions %>%
    select(-all_of(labels))

  # Loop over the hyperparameter grid
  for (i in 1:nrow(hyper_grid)) {
    # Get the hyperparameters
    for (param in names(hyper_grid)) {
      param_value <- hyper_grid[i, param]
      assign(paste0("selected_", param), param_value)

      # Ensure that the hyperparameters are not negative
      # If they are, set them to random value between 1 and 12
      if (param_value <= 0) {
        assign(paste0("selected_", param), sample(1:12, 1))
      }
    }

    # Check that the tuning has not been done before for the exact same parameters
    if (nrow(meta_model_metrics) > 0) {
      if (any(
        meta_model_metrics$Max_depth == selected_max_depth &
          meta_model_metrics$Child_weight == selected_min_child &
          meta_model_metrics$Eta == selected_eta &
          meta_model_metrics$Gamma == selected_gamma
      )) {
        cat("\n\t\t\t Hyperparameters already tuned. Skipping.\n\n")
        next
      }
    }

    # Print the selected hyperparameters
    cat(
      "\n\n\t Selected Parameters: \n",
      "Feature: ", selected_feature, "\n",
      "Trees: ", selected_trees, "\n",
      "Eta: ", selected_eta, "\n",
      "Gamma: ", selected_gamma, "\n",
      "Max Depth: ", selected_max_depth, "\n",
      "Min Child Weight: ", selected_min_child, "\n",
      "Seed: ", selected_seed, "\n\n"
    )

    # Set the label
    y <- ifelse(feature %in% cat_features, as.integer, as.numeric)(predictions[[paste0("act_", selected_feature)]])
    cat("\n\n Target: ", selected_feature, "\n")
    print(head(y))

    # Adjust weights for categorical features
    if (selected_feature %in% cat_features) {
      selected_weights <- arm_weights[selected_feature, ]
      cat("\n\n Selected Weights: \n")
      print(selected_weights)

      selected_ref <- as.numeric(
        selected_weights[c(
          "Weight_loss",
          "Weight_normal",
          "Weight_gain"
        )]
      )
      weights <- as.numeric(feature_digit_function(factor(y, levels = c(0, 1, 2)), selected_ref))
      # Create the xgb data
      xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)
      cat("\n Loaded xgb data with weights.\n")
    } else {
      # Create the xgb data
      cat("\n Loaded xgb data without weights.\n")
      xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)
    }

    # Set common parameters
    xgb_params <- list(
      data = xgb_data,
      nrounds = selected_trees,
      nfold = 10,
      early_stopping_rounds = 150,
      objective = if (selected_feature %in% cat_features) "multi:softmax" else "reg:squarederror",
      eval_metric = if (selected_feature %in% cat_features) "mlogloss" else "rmse",
      eta = selected_eta,
      gamma = selected_gamma,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      stratified = TRUE,
      print_every_n = 50,
      prediction = TRUE
    )

    # Add num_class only if it's a classification problem
    if (selected_feature %in% cat_features) {
      xgb_params$num_class <- 3
    }

    # Train model using the parameters list
    xgb_model <- do.call(xgb.cv, xgb_params)

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

    best_train <- if (best_iteration > 0) {
      if (selected_feature %in% cat_features) {
        xgb_model$evaluation_log$train_mlogloss_mean[best_iteration]
      } else {
        xgb_model$evaluation_log$train_rmse_mean[best_iteration]
      }
    } else {
      NA # Or appropriate default/error value
    }

    best_test <- if (best_iteration > 0) {
      if (selected_feature %in% cat_features) {
        xgb_model$evaluation_log$test_mlogloss_mean[best_iteration]
      } else {
        xgb_model$evaluation_log$test_rmse_mean[best_iteration]
      }
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
        "Feature" = selected_feature,
        "Trees" = best_iteration,
        "Eta" = selected_eta,
        "Gamma" = selected_gamma,
        "Max_depth" = selected_max_depth,
        "Child_weight" = selected_min_child,
        "Train_result" = best_train,
        "Test_result" = best_test,
        "Seed" = selected_seed
      )
    ) %>%
      arrange(Test_result)
  }

  # Save the hyperparameters
  output_file <- file.path(
    hyper_param_folder,
    selected_cancer,
    "individual_features"
  )

  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
  }

  unique_id <- Sys.time() %>%
    as.numeric() %>%
    round()

  write.csv(
    meta_model_metrics,
    paste0(
      output_file,
      "/",
      selected_feature,
      "_",
      unique_id,
      "_meta_hyperparams.csv"
    ),
    row.names = FALSE
  )
}

stop("Task Completed")
