library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # This is the SLURM_ARRAY_TASK_ID

cat("\n\n The index for this run is: ", index, "\n")

selected_max_depth <- 5
selected_min_child <- 1
selected_eta <- 0.3
selected_gamma <- 0
selected_trees <- 10000

cancer_types <- c("BLCA", "BRCA", "CESC", "HNSC", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PRAD", "STAD", "THCA")

# Load the constant data
## HRD scores
ori_hrd <- read_tsv("/hpc/shared/prekovic/dhaynessimmons/data/CIN/TCGA.HRD_withSampleID.txt")

# Pericentromeric CNVs
peri_cnv <- read.csv("/hpc/shared/prekovic/dhaynessimmons/data/CIN/lim_alpha_incl_TCGA_pericentro_cnv.csv") %>%
  mutate_all(~ replace(., is.na(.), 0)) %>%
  mutate(sampleID = gsub("-", ".", sampleID))

cat("\n\n pericentromeric data: \n")
print(head(peri_cnv[, 1:7], 3))

cat("\n\n All dfs loaded \n")

t_hrd <- as.data.frame(t(ori_hrd))
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ]
hrd <- as.data.frame(first_hrd[-1, ]) %>%
  mutate_all(as.numeric) %>%
  rename(loh_hrd = "hrd-loh") %>%
  mutate(new = str_replace_all(rownames(.), "-", "\\."))

rownames(hrd) <- hrd$new
hrd <- hrd %>%
  select(-new)

cat("\n\n hrd data: \n")
print(head(hrd), 3)

rm(t_hrd)
rm(first_hrd)

full_cin <- merge(
  hrd,
  peri_cnv,
  by.x = "row.names",
  by.y = "sampleID"
) %>%
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")

cat("\n\n full cin data: \n")
print(head(full_cin[, 1:7], 3))
cat("\n\n Full cin features: \n")
print(colnames(full_cin))

rm(hrd)
rm(peri_cnv)

# Define the features to be used
aneu_reg_feature_list <- colnames(full_cin)
selected_feature <- aneu_reg_feature_list[[index]]
cat(paste0("\n\n Selected feature: ", selected_feature, "\n"))
rm(aneu_reg_feature_list)

# Loop over the cancer types
for (cancer in cancer_types){
  cat(paste0("\n\n\t\t\t\t\t\t Processing for cancer type: ", cancer, "\n"))
  rna_folder <- paste0("/hpc/shared/prekovic/dhaynessimmons/data/mRNA/cancer_specific/", cancer, "/")

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

  # MODELLING

    aneu_reg_metrics_df <- data.frame(
    RNA_set = character(),
    Trees = numeric(),
    Feature = character(),
    Max_depth = numeric(),
    Child_weight = numeric(),
    Eta = numeric(),
    Gamma = numeric(),
    Trained_RMSE = numeric(),
    Test_RMSE = numeric()
  )

  for (i in 1:length(rna_list)) {
    rna_data <- rna_list[[i]]
    selected_rna_set <- names(rna_list)[i]
    cat(paste0("\n\t", selected_rna_set, "\n"))
    print(head(rna_data[, 1:5], 3))

    full_df <- merge(rna_data,
      full_cin,
      by = "row.names"
    )
    cat("\n\n Full df: \n")
    print(head(full_df[, 1:5], 3))

    y <- as.integer(full_df[[selected_feature]])

    # Check that not all the values are the same
    if (length(unique(y)) == 1) {
      cat(paste0(
        "Warning: All values for ", selected_feature,
        " are the same. Skipping this feature.\n"
      ))

      # fill the aneu_reg_metrics_df with NA values
      aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df, data.frame(
        RNA_set = selected_rna_set,
        Trees = NA,
        Feature = selected_feature,
        Max_depth = NA,
        Child_weight = NA,
        Eta = NA,
        Gamma = NA,
        Trained_RMSE = NA,
        Test_RMSE = NA
      ))
      next
    }
    
    cat("\n\n Y: \n")
    print(head(y))

    X <- full_df %>% select(-c("Row.names", colnames(full_cin)))

    xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

    cat(paste0(
      "\t\t Max_depth: ", selected_max_depth,
      "\n"
    ))

    m_xgb_untuned <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      objective = "reg:squarederror",
      eval_metric = "rmse",
      early_stopping_rounds = 250,
      nfold = 10,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      eta = selected_eta,
      gamma = selected_gamma,
      print_every_n = 100
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

    best_rmse_trained <- if (best_iteration > 0) {
      m_xgb_untuned$evaluation_log$train_rmse_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    best_rmse_test <- if (best_iteration > 0) {
      m_xgb_untuned$evaluation_log$test_rmse_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    cat(paste0(
      "The best iteration occurs with tree #: ",
      best_iteration, "\n\n"
    ))

    aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df, data.frame(
      RNA_set = selected_rna_set,
      Trees = best_iteration,
      Feature = selected_feature,
      Max_depth = selected_max_depth,
      Child_weight = selected_min_child,
      Eta = selected_eta,
      Gamma = selected_gamma,
      Trained_RMSE = best_rmse_trained,
      Test_RMSE = best_rmse_test
    ))
  }

  new_dir <- paste0(
    "/hpc/shared/prekovic/dhaynessimmons/data/model_output/cancer_specific/",
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
    aneu_reg_metrics_df,
    file = file_name,
    row.names = FALSE
  )
  cat(paste0("\n Completed processing for cancer type: ", cancer, "\n"))
}

cat(paste0("\n Completed processing for feature: ", selected_feature, "\n"))
