# This script is intended to generate the predictions for the all base models 
# for the training/tuning of the meta model

# Needs to be updated to reflect the changes that are made to the pred_p1 and pred_p2 scripts


# Load the packages
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")
input_path <- "Data/Model_input/"

# Set the order for the response features
response_features_order <- c(
  "1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q","6p", "6q", "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p", "12q", "13q", "14q", "15q", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p", "20q", "21q", "22q", "ai1", "lst1", "loh_hrd", "peri_1",
  "peri_2", "peri_3", "peri_4", "peri_5", "peri_7", "peri_8", "peri_9", "peri_10", "peri_16_1", "peri_16_2", "peri_17", "peri_20", "peri_22_3"
)

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

# The CIN response features
# Categorical features
cat_cin <- read_tsv(
  paste0(input_path, "CIN_Features/PANCAN_ArmCallsAndAneuploidyScore_092817.txt")
) %>%
  replace(is.na(.), 0) %>%
  select(-c("Type", "Aneuploidy Score")) %>%
  mutate(Sample = str_replace_all(Sample, "-", "\\.")) %>%
  column_to_rownames("Sample") %>%
  mutate_all(~ replace(., . == 1, 2)) %>%
  mutate_all(~ replace(., . == 0, 1)) %>%
  mutate_all(~ replace(., . == -1, 0))

# Numerical features
# HRD scores
ori_hrd <- read_tsv(paste0(input_path, "CIN_Features/TCGA.HRD_withSampleID.txt"))

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

# Pericentromeric CNVs
peri_cnv <- read.csv(paste0(input_path, "/CIN_Features/lim_alpha_incl_TCGA_pericentro.csv")) %>% 
  mutate_all(~ replace(., is.na(.), 0)) %>%
  mutate(sampleID = gsub("-", ".", sampleID)) %>% 
  column_to_rownames("sampleID")

reg_cin <- merge(
  hrd,
  peri_cnv,
  by = "row.names"
) %>%
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")

rm(peri_cnv)
rm(hrd)

full_cin <- merge(
  cat_cin,
  reg_cin,
  by = "row.names"
) %>%
  column_to_rownames("Row.names")

# All response feature names
response_features <- colnames(full_cin)
reg_features <- colnames(reg_cin)
cat_features <- colnames(cat_cin)

rm(reg_cin)
rm(cat_cin)

# Create the folds to be used
full_data <- merge(read.csv(
    paste0(
      input_path,
      "RNA/Train/train_tpm_soi.csv"
    ),
    row.names = 1
  ),
  full_cin,
  by = "row.names")

# Creating the folds and returning the indices for the out-of-fold predictions only
folds <- createFolds(full_data[["1p"]], k = 2, list = TRUE, returnTrain = FALSE)

# Get the indices of the folds
fold_indices <- combine_all_folds(folds)

# Create the empty dataframe that will store the predictions from each model
predictions_df <- data.frame(
  matrix(
    ncol = length(response_features),
    nrow = length(folds$Fold1)*length(folds)
  )
)
colnames(predictions_df) <- paste0("pred_", response_features)
rownames(predictions_df) <- as.factor(fold_indices)

# Get the actual fold labels
fold_labels <- full_data %>% 
  select(all_of(intersect(response_features, colnames(.)))) %>% 
  mutate(index = as.factor(row_number()),
         sort_order = match(rownames(.), fold_indices)) %>%
  rename_with(~ paste0("act_", .)) %>%
  column_to_rownames("act_index") %>% 
  arrange(act_sort_order)

# merge the fold labels with the predictions dataframe
predictions_df <- merge(predictions_df, fold_labels, by = "row.names") %>% 
  arrange(act_sort_order) %>% 
  column_to_rownames("Row.names") %>% 
  select(-act_sort_order)
s
cat_parameters <- read.csv(paste0(input_path, "/CIN_Features/XGB_cat_parameters.csv"))
reg_parameters <- read.csv(paste0(input_path, "/CIN_Features/XGB_reg_parameters.csv"))

for (index in 1:length(response_features)) {
  feature <- response_features[index]
  cat("\n index: ", index, "\n",
      "Feature: ", feature, "\n\n")
  
  if (feature %in% cat_features) {
    # Get the parameters from stored Parameters file
    parameters <- cat_parameters
    
    # select the parameters and weights corresponding to the index
    selected_parameters <- parameters[parameters$Feature == feature, ]
    
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
      "Seed: ", selected_seed, "\n\n"
    )
    
    rna_selection_name <- rna_list[[selected_rna_set]]
    rna_set <- read.csv(
      paste0(input_path,
        "RNA/Train/train_",
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
      X <- full_df %>% select(-c("Row.names", response_features))
      
      train_y_factor <- factor(y, levels = c(0, 1, 2))
      weights <- as.numeric(feature_digit_function(train_y_factor))
      
      xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)
      
      xgb_model <- xgb.cv(
        data = xgb_data,
        nrounds = selected_trees,
        objective = "multi:softmax",
        eval_metric = "mlogloss",
        early_stopping_rounds = 50,
        folds = folds,
        max_depth = selected_max_depth,
        min_child_weight = selected_min_child,
        eta = selected_eta,
        gamma = selected_gamma,
        num_class = 3,
        print_every_n = 10,
        prediction = TRUE
      )
      
      # Store the predictions in the corresponding column
      predictions_df[[paste0("pred_", selected_feature)]] <- xgb_model$pred[, 1]
      
    } else if (feature %in% reg_features) {
      # Get the parameters from stored Parameters file
      parameters <- reg_parameters
      
      # select the parameters and weights corresponding to the index
      selected_parameters <- parameters[parameters$Feature == feature, ]
      
      selected_feature <- selected_parameters$Feature
      selected_rna_set <- selected_parameters$RNA_set
      selected_trees <- as.numeric(selected_parameters$Trees)
      selected_eta <- selected_parameters$Eta
      selected_gamma <- selected_parameters$Gamma
      selected_max_depth <- selected_parameters$Max_depth
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
        "Min Child: ", selected_min_child, "\n",
        "Seed: ", selected_seed, "\n\n"
      )
      
      rna_selection_name <- rna_list[[selected_rna_set]]
      rna_set <- read.csv(
        paste0(input_path,
          "RNA/Train/train_",
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
        X <- full_df %>% select(-c("Row.names", response_features))
        
        xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)
        
        xgb_model <- xgb.cv(
          data = xgb_data,
          nrounds = selected_trees,
          objective = "reg:squarederror",
          eval_metric = "rmse",
          early_stopping_rounds = 50,
          folds = folds,
          max_depth = selected_max_depth,
          min_child_weight = selected_min_child,
          eta = selected_eta,
          gamma = selected_gamma,
          print_every = 15,
          prediction = TRUE
          )
        
        # Store the predictions in the corresponding column
        predictions_df[[paste0("pred_", selected_feature)]] <- xgb_model$pred
        
    } else {
      cat("Feature not found")
    }
}


# Reorder the columns of the predictions dataframe
suffixes <- sapply(strsplit(names(predictions_df), "_"), tail, 1)
order <- order(suffixes)
df_ordered <- predictions_df[, order]
  
write.csv(predictions_df, paste0(input_path,"/Base_predictions/full_XGB_base_predictions.csv"))
write.csv(predictions_df, "Data/Model_output/XGB_base_predictions.csv")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
