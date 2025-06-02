library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # This is the SLURM_ARRAY_TASK_ID
index <- 1

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")
rna_data_path <- "Data/Cancer_specific_data/Model_input/RNA/"

cat("The index for this run is: ", index, "\n")

selected_max_depth <- 5
selected_min_child <- 1
selected_eta <- 0.3
selected_gamma <- 0
selected_trees <- 10000

cancer_types <- c("BLCA", "BRCA", "CESC", "HNSC", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PRAD", "STAD", "THCA")

# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) selected_weights[as.numeric(x)])
}

# Load the constant data
# Arm Level Aneuploidies
chr_cnv <- read_tsv(
  "Data/Cancer_specific_data/Model_input/CIN/PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
) %>%
  replace(is.na(.), 0) %>%
  select(-c("Type", "Aneuploidy Score")) %>%
  mutate(
    Sample = str_replace_all(Sample, "-", "\\."),
    across(-Sample, ~ case_when(
      . == 1 ~ 2,
      . == 0 ~ 1,
      . == -1 ~ 0,
      is.na(.) ~ 1,
      TRUE ~ .
    ))
  ) %>%
  rename( "13q" = "13 (13q)",
          "14q" = "14 (14q)",
          "15q" = "15 (15q)",
          "21q" = "21 (21q)",
          "22q" = "22 (22q)"
  ) %>%
  column_to_rownames("Sample")

cat("\n\n CNV df: \n")
print(head(chr_cnv[, 1:5], 3))

# Define the features to be used
aneu_cat_feature_list <- colnames(chr_cnv)
selected_feature <- aneu_cat_feature_list[[index]]
cat(paste0("\n\n Selected feature: ", selected_feature, "\n"))
rm(aneu_cat_feature_list)

# Loop over the cancer types
for (cancer in cancer_types){
  cat(paste0("\n\n Processing for cancer type: ", cancer, "\n"))
  rna_folder <- paste0(rna_data_path, cancer, "/")

  # RNA SOI SETS
  # Expected Counts
  exp_set <- read.csv(
    paste0(
      rna_folder,
      "exp_soi.csv"
    ),
    row.names = 1
  )

  scld_exp_set <- read.csv(
    paste0(
      rna_folder,
      "scld_exp_soi.csv"
    ),
    row.names = 1
  )

  log_exp <- read.csv(
    paste0(
      rna_folder,
      "log_exp_soi.csv"
    ),
    row.names = 1
  )

  log_scld_exp <- read.csv(
    paste0(
      rna_folder,
      "log_scld_exp_soi.csv"
    ),
    row.names = 1
  )

  # Transcripts Per Million
  tpm_set <- read.csv(
    paste0(
      rna_folder,
      "tpm_soi.csv"
    ),
    row.names = 1
  )

  scld_tpm_set <- read.csv(
    paste0(
      rna_folder,
      "scld_tpm_soi.csv"
    ),
    row.names = 1
  )

  log_tpm <- read.csv(
    paste0(
      rna_folder,
      "log_tpm_soi.csv"
    ),
    row.names = 1
  )

  log_scld_tpm <- read.csv(
    paste0(
      rna_folder,
      "log_scld_tpm_soi.csv"
    ),
    row.names = 1
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

  # Load the cancer-specific class weights for the arms
  arm_weights <- read.csv(
    paste0(
      "Data/Cancer_specific_data/Model_input/Parameters/Arm_class_weights/",
      cancer,
      "_arm_weights.csv"
    ),
    row.names = 1
  )
  cat("\n\n Arm Weights: \n")
  print(head(arm_weights, 3))

  # Determine the class weights for the target feature
  selected_weights <- arm_weights[selected_feature, ]
 
  cat("\n\n Target Weights: \n")
  print(selected_weights)
  cat("\n\n")
  rm(arm_weights)


  # MODELLING

  aneu_cat_metrics_df <- data.frame(
    RNA_set = character(),
    Trees = numeric(),
    Feature = character(),
    Max_depth = numeric(),
    Child_weight = numeric(),
    Eta = numeric(),
    Gamma = numeric(),
    Weight_loss = numeric(),
    Weight_normal = numeric(),
    Weight_gain = numeric(),
    Trained_mlogloss = numeric(),
    Test_mlogloss = numeric()
  )

  for (i in 1:length(rna_list)) {
    rna_data <- rna_list[[i]]
    selected_rna_set <- names(rna_list)[i]
    cat(paste0("\t", selected_rna_set, "\n"))
    print(head(rna_data[, 1:5], 3))

    full_df <- merge(rna_data,
      chr_cnv,
      by = "row.names"
    )
    cat("\n\n Full df: \n")
    print(head(full_df[, 1:5], 3))

    y <- as.integer(full_df[[selected_feature]])
    cat("\n\n Y: \n")
    print(head(y))

    X <- full_df %>% select(-c("Row.names", colnames(chr_cnv)))

    train_y_factor <- factor(y, levels = c(0, 1, 2))
    weights <- as.numeric(feature_digit_function(train_y_factor))
    cat("\n\n Weights: \n")
    print(head(weights))

    xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)

    cat(paste0(
      "\t\t Max_depth: ", selected_max_depth,
      "\n"
    ))

    m_xgb_untuned <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      objective = "multi:softmax",
      eval_metric = "mlogloss",
      early_stopping_rounds = 250,
      nfold = 10,
      max_depth = selected_max_depth,
      eta = selected_eta,
      gamma = selected_gamma,
      num_class = 3,
      print_every_n = 10
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
      RNA_set = selected_rna_set,
      Trees = best_iteration,
      Feature = selected_feature,
      Max_depth = selected_max_depth,
      Child_weight = selected_min_child,
      Eta = selected_eta,
      Gamma = selected_gamma,
      Weight_loss = selected_weights[1],
      Weight_normal = selected_weights[2],
      Weight_gain = selected_weights[3],
      Trained_mlogloss = best_mlogloss_train,
      Test_mlogloss = best_mlogloss_test
    ))
  }

  new_dir <- paste0(
    "Data/Cancer_specific_data/Model_output/Hyperparameters/",
    cancer,
    "/"
  )

  if (!dir.exists(new_dir)) {
    dir.create(new_dir)
  }

  datetime <- as.character(Sys.time()) %>%
    str_replace_all(" ", "_") %>%
    str_replace_all(":", "_") %>%
    str_replace_all("\\.", "_")

  file_name <- paste0(
    new_dir,
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
  cat(paste0("\n Completed processing for cancer type: ", cancer, "\n"))
}

cat(paste0("\n Completed processing for feature: ", selected_feature, "\n"))
