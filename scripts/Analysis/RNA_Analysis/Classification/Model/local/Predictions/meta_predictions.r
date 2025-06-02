# Create the predictor data for the meta learners
setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)

# Function to map factor levels to weights
feature_digit_function <- function(factors, weight_ref) {
  sapply(factors, function(x) weight_ref[as.numeric(x)])
}

feat_imp <- function(imp_df, top_gene_num = 10, basis = "Gain", Type) {
  # Create the feature importance matrix
  created_imp <- imp_df %>%
    group_by(Feature) %>%
    summarise_all(mean) %>%
    # create the combined normalised importance
    mutate(
      Combined = rowSums(across(.cols = -Feature))
    ) %>%
    mutate(
      Normalised = Combined / sum(Combined)
    ) %>%
    select(-Combined) %>%
    arrange(desc(basis))
  
  type_file_folder <- ifelse(Type == "Meta", "Meta_feat_imp/", "Base_feat_imp/")
  csv_filefold <- file.path("Data/Gen_model_output/Results/Feature_importance/",
                            type_file_folder)
  if (!dir.exists(csv_filefold)) {
    csv_filefold  }
  # Save the feature importance df
  write.csv(
    created_imp,
    paste0(
      csv_filefold,
      feature,
      "_feature_importance.csv"
    )
  )
  
  
  # Create the plot for the top 10 features
  top_genes <- created_imp %>%
    select(Feature, basis) %>%
    top_n(top_gene_num, basis) %>%
    arrange(desc(basis))
  
  # Create the plot
  top_genes_plot <- ggplot(top_genes[1:top_gene_num,],
                           aes(x = reorder(Feature, .data[[basis]]), 
                               y = .data[[basis]])
  ) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(
      title = paste0("Top ", top_gene_num, " Features for ", feature),
      x = "Feature",
      y = basis
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save the plot
  plt_filefold <- file.path("Plots/Model_Plots/General_model/Results/Feature_importance/",
                            type_file_folder)
  if (!dir.exists(plt_filefold)) {
    dir.create(plt_filefold, recursive = TRUE)
  }
  
  ggsave(
    filename = paste0(
      plt_filefold,
      feature,
      "_Top_",
      top_gene_num,
      "_Features.pdf"
    ),
    plot = top_genes_plot,
    width = 10,
    height = 10
  )
  
  return(list(feature_imp_df = created_imp, top_genes = top_genes))
}

input_path <- "Data/Gen_model_input/"

# Upload the base predictions
base_oof_predictions <- read.csv(
  paste0(
    input_path,
    "Base_predictions/Full_base_predictions.csv"
  ),
  row.names = 1
)

# Load the meta parameters
meta_parameters <- read.csv(
  paste0(
    input_path,
    "Hyperparameters/meta_hyperparams.csv"
  )
)

# Create the training folds
meta_folds <- createFolds(base_oof_predictions[["act_1p"]], k = 10, list = TRUE, returnTrain = FALSE)

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

# Create the empty dataframe that will store the importance matrix for each model of the meta learners
full_meta_importance <- data.frame()

# Get list of regression and categorical features
cat_features <- response_features[str_detect(response_features, regex("^\\d"))]
reg_features <- response_features[!response_features %in% cat_features]
cat("\n\n Regression Features: \n")
print(reg_features)
cat("\n\n Categorical Features: \n")
print(cat_features)

# Create the predictor data for the meta learners
meta_input_data <- base_oof_predictions %>%
  select(starts_with("pred_")) %>%
  rename_all(~ gsub("pred_", "", .)) 

cat("\n\n Meta Input Data: \n")
print(head(meta_input_data[, 1:5], 3))

for (feature in response_features) {
  cat(
    "\n\n\t\t\t\t\t\t Feature: ", feature, "\n\n"
  )

  # Create df to store the importance matrix
  feature_imp_df <- data.frame()

  # Create the out-of-fold predictions column
  oof_predictions <- numeric(nrow(base_oof_predictions))

  # Determine parameter set based on feature type
  parameters <- if (feature %in% meta_parameters$Feature) {
    meta_parameters
  } else {
    cat("Feature not found")
    next
  }

  # Select the parameters and weights corresponding to the feature
  selected_parameters <- parameters[parameters$Feature == feature, ]

  if (nrow(selected_parameters) == 0) {
    cat("No parameters found for feature: ", feature, "\n")
    next
  }

  # Print selected parameters
  cat("\n Selected Parameters:\n")
  print(selected_parameters)

  y <- ifelse(feature %in% cat_features, as.integer, as.numeric)(meta_input_data[[feature]])
  cat("\n Y: ")
  print(head(y))

  # Fold-wise training
  for (fold_index in seq_along(folds)) {
    cat(
      "\n\t\t\t\t Fold: ", fold_index, "\n"
    )
    train_indices <- setdiff(seq_len(nrow(meta_input_data)), folds[[fold_index]])
    valid_indices <- folds[[fold_index]]

    train_data <- xgb.DMatrix(data = as.matrix(meta_input_data[train_indices, ]), label = y[train_indices])
    cat("\n Train Data: \n")
    print(head(y[train_indices], 10))

    valid_data <- xgb.DMatrix(data = as.matrix(meta_input_data[valid_indices, ]), label = y[valid_indices])

    if (feature %in% cat_features) {
      selected_ref <- as.numeric(
        selected_parameters[c(
          "Weight_loss",
          "Weight_normal",
          "Weight_gain"
        )]
      )
      weights <- as.numeric(feature_digit_function(factor(y[train_indices], levels = c(0, 1, 2)), selected_ref))
      setinfo(train_data, "weight", weights)

      cat(
        "\n\n Weight Loss: ", selected_ref[1], "\n",
        "Weight Normal: ", selected_ref[2], "\n",
        "Weight Gain: ", selected_ref[3], "\n\n")

      cat("\n Weights: \n")
      print(head(weights, 10))
      cat("\n\n")
  
    }

    # Adjust weights for categorical features

    # Train model
    params_list <- list(
      data = train_data,
      nrounds = selected_parameters$Trees,
      objective = ifelse(feature %in% cat_features, "multi:softmax", "reg:squarederror"),
      eval_metric = ifelse(feature %in% cat_features, "mlogloss", "rmse"),
      max_depth = selected_parameters$Max_depth,
      min_child_weight = selected_parameters$Child_weight,
      early_stopping_rounds = early_stop,
      eta = selected_parameters$Eta,
      gamma = selected_parameters$Gamma,
      watchlist = list(eval = valid_data),
      print_every_n = print_every
    )

    # Only add num_class for categorical features
    if (feature %in% cat_features) {
      params_list$num_class <- 3
    }

    xgb_model <- do.call(xgb.train, params_list)

    # Store OOF predictions
    oof_predictions[valid_indices] <- predict(xgb_model, valid_data)
    cat(
      "\n OOF Predictions: \n"
    )
    print(head(oof_predictions, 10))

    # Create the feature importance matrix
    importance_matrix <- xgb.importance(
      feature_names = colnames(meta_input_data),
      model = xgb_model
    )

    # Store the importance matrix in the dataframe
    feature_imp_df <- rbind(feature_imp_df, as.data.frame(importance_matrix))
  }

  # Store the predictions in the corresponding column
  meta_oof_predictions[[paste0("pred_", selected_feature)]] <- oof_predictions

  # get the feature imp df
  imp <- feat_imp(imp_df = feature_imp_df, top_gene_num = 3)

  top_genes <- imp$top_genes

  full_meta_importance <- rbind(
    full_meta_importance,
    imp$feature_imp_df %>%
      mutate(Target = feature)
  ) %>%
    arrange(desc(Gain))
}

cat("Training of the meta models complete")

# Save the out-of-fold predictions
write.csv(
  meta_oof_predictions,
  "Data/Gen_model_output/Predictions/Meta_predictions/Full_meta_predictions.csv"
)
write.csv(
  meta_oof_predictions,
  "Data/Gen_model_input/Meta_predictions/Full_meta_predictions.csv"
)
