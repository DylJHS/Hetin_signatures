# This script is the first part of the cancer-specific model prediction pipeline.
# Its goal is to generate predictions for the all base models in parallel

# Load the packages
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # Cancer folder index

early_stop <- 250
print_every <- 50
cv <- 10

cancer_types <- c(
  "BLCA", "BRCA", "CESC",
  "HNSC", "LGG", "LIHC",
  "LUAD", "LUSC", "OV",
  "PRAD", "STAD", "THCA"
)

input_path <- "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/Cancer_specific_data/Model_input/"
output_path <- "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/Cancer_specific_data/Model_output"

# Load the functions
# Function to map factor levels to weights
feature_digit_function <- function(factors, reference) {
  sapply(factors, function(x) reference[as.numeric(x)])
}

# Feature importance function to create the df and sort based on the Gain and plot the top 10 features
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
  csv_filefold <- file.path(
    output_path,
    selected_cancer, "/",
    "Results/Feature_importance/",
    type_file_folder
  )

  if (!dir.exists(csv_filefold)) {
    dir.create(csv_filefold, recursive = TRUE)
  }
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
  top_genes_plot <- ggplot(
    top_genes[1:top_gene_num, ],
    aes(
      x = reorder(Feature, .data[[basis]]),
      y = .data[[basis]]
    )
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
  plt_filefold <- file.path(
    output_path,
    selected_cancer, "/",
    "Plots/",
    "Feature_importance/",
    type_file_folder
  )

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
  scaled_expected_counts = "scld_exp",
  log_expected_counts = "log_exp",
  log_scaled_expected_counts = "log_scld_exp"
)

# print(rna_list)

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
print(head(cat_cin[1:5], 3))

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
    "CIN/lim_alpha_incl_TCGA_pericentro.csv"
  )
) %>%
  mutate_all(~ replace(., is.na(.), 0)) %>%
  mutate(sampleID = gsub("-", ".", sampleID)) %>%
  column_to_rownames("sampleID")

# print(dim(peri_cnv))

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
# print(head(full_cin[1:5]))

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

# Loop through the cancer types
selected_cancer <- cancer_types[index]

# Create the folds to be used
full_data <- merge(
  read.csv(
    paste0(
      input_path,
      "RNA/",
      selected_cancer,
      "/Train/train_tpm_soi.csv"
    ),
    row.names = 1
  ),
  full_cin,
  by = "row.names"
)
cat("\n Full Data:  \n")
print(head(full_data[1:5]))

# Creating the folds and returning the indices for the out-of-fold predictions only
folds <- createFolds(full_data[["1p"]], k = cv, list = TRUE, returnTrain = FALSE)
cat("\n Folds # :", length(folds), "\n")


# Get the actual fold labels
fold_labels <- full_data %>%
  select(all_of(intersect(response_features, colnames(.)))) %>%
  mutate(
    index = as.factor(row_number())
  ) %>%
  rename_with(~ paste0("act_", .))

# merge the fold labels with the predictions dataframe
base_oof_predictions <- data.frame(
  matrix(
    ncol = length(response_features),
    nrow = length(combine_all_folds(folds))
  )
)

colnames(base_oof_predictions) <- paste0("pred_", response_features)
base_oof_predictions <- cbind(base_oof_predictions, fold_labels) %>%
arrange(act_index)

cat("\n Predictions Dataframe:  \n")
print(head(base_oof_predictions[1:5]))

base_cat_parameters <- read.csv(
  paste0(
    input_path,
    "Parameters/Hyperparameters/",
    selected_cancer,
    "/base_cat_hyperparams.csv"
  )
)
base_reg_parameters <- read.csv(
  paste0(
    input_path,
    "Parameters/Hyperparameters/",
    selected_cancer,
    "/base_reg_hyperparams.csv"
  )
)


# Create the full feature importance dataframe
full_feature_imp <- data.frame()

for (feature in response_features) {
  cat( 
    "Feature: ", feature, "\n\n"
  )

  # Create df to store the importance matrix
  feature_imp_df <- data.frame()

  # Create the out-of-fold predictions column
  oof_predictions <- numeric(nrow(base_oof_predictions))

  # Determine parameter set based on feature type
  parameters <- if (feature %in% cat_features) {
    base_cat_parameters
  } else if (feature %in% reg_features) {
    base_reg_parameters
  } else {
    cat("Feature not found")
    # remove the feature from the base_oof_predictions
    base_oof_predictions <- base_oof_predictions %>%
      select(-c(
        paste0("pred_", feature),
        paste0("act_", feature)
      ))
    next
  }

  # Select the parameters and weights corresponding to the feature
  selected_parameters <- parameters[parameters$Feature == feature, ]

  if (nrow(selected_parameters) == 0 | is.na(selected_parameters$Trees)) {
    cat("No parameters found for feature: ", feature, "\n")
    base_oof_predictions <- base_oof_predictions %>%
      select(-c(
        paste0("pred_", feature),
        paste0("act_", feature)
      ))
    next
  }

  # Print selected parameters
  cat("\n Selected Parameters:\n")
  print(selected_parameters)

  # Load and prepare data
  rna_selection_name <- rna_list[[selected_parameters$RNA_set]]
  rna_set <- read.csv(
    paste0(
      input_path,
      "RNA/",
      selected_cancer,
      "/Train/train_",
      rna_selection_name,
      "_soi.csv"
    ),
    row.names = 1
  )

  full_df <- merge(rna_set, full_cin, by = "row.names")

  y <- ifelse(feature %in% cat_features, as.integer, as.numeric)(full_df[[feature]])
  # cat("\n y: ")
  # print(head(y))

  X <- full_df %>% select(-c("Row.names", all_of(response_features)))

  # Initialize containers for predictions and feature importance
  oof_predictions <- numeric(nrow(X))
  feature_imp_df <- data.frame()

  # Fold-wise training
  for (fold_index in seq_along(folds)) {
    cat(
      "\n Fold: ", fold_index, "\n"
    )
    train_indices <- setdiff(seq_len(nrow(full_df)), folds[[fold_index]])
    valid_indices <- folds[[fold_index]]

    train_data <- xgb.DMatrix(data = as.matrix(X[train_indices, ]), label = y[train_indices])
    valid_data <- xgb.DMatrix(data = as.matrix(X[valid_indices, ]), label = y[valid_indices])

    # Adjust weights for categorical features
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
    }

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
    # cat(
    #   "\n OOF Predictions: \n",
    #   dim(oof_predictions), "\n"
    # )
    # print(head(oof_predictions, 10))

    # Create the feature importance matrix
    importance_matrix <- xgb.importance(feature_names = colnames(X), model = xgb_model)

    # Store the importance matrix in the dataframe
    feature_imp_df <- rbind(feature_imp_df, as.data.frame(importance_matrix))
  }

  # Store the predictions in the corresponding column
  base_oof_predictions[[paste0("pred_", feature)]] <- oof_predictions
  # print(base_oof_predictions[[paste0("pred_", feature)]][1:10])

  # get the feature imp df
  imp <- feat_imp(imp_df = feature_imp_df, top_gene_num = 10, basis = "Gain", Type = "Base")

  top_genes <- imp$top_genes

  full_feature_imp <- rbind(
    full_feature_imp,
    imp$feature_imp_df %>%
      mutate(Target = feature)
  ) %>% 
    arrange(desc(Gain))
}

# Save the feature importance df
Feature_imp_folder <- file.path(
  output_path,
  selected_cancer,
  "Results/Feature_importance/"
)
if (!dir.exists(Feature_imp_folder)) {
  dir.create(Feature_imp_folder, recursive = TRUE)
}
write.csv(
  full_feature_imp,
  paste0(
    Feature_imp_folder,
    "Full_feature_importance.csv"
  )
)

# Save the predictions
Predictions_folder <- file.path(
  output_path,
  selected_cancer,
  "Results/Predictions/"
)
if (!dir.exists(Predictions_folder)) {
  dir.create(Predictions_folder, recursive = TRUE)
}
write.csv(
  base_oof_predictions,
  paste0(
    Predictions_folder,
    "Full_base_predictions.csv"
  )
)
cat("\n\n\t\t\t\t Training the base models for ", selected_cancer, " complete")
