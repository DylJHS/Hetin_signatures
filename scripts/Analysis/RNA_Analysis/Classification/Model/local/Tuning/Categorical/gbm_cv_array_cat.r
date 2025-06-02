library(dplyr)
library(readr)
library(mltools)
library(knitr)
library(factoextra)
library(data.table)
library(psych)
library(limma)
library(edgeR)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project") # nolint

extract_element <- function(strings, index) {
  # Split each string by "." and extract the third element
  element_list <- sapply(strsplit(strings, "\\."), function(x) x[index])
  return(element_list)
}

untransform_exp <- function(x) { # Function to convert the log transformed counts back into original counts # nolint
  return(ceiling((2^x) - 1))
}

transform_exp <- function(x) { # Function to convert the log transformed counts back into original counts # nolint
  return(log2(x + 1))
}

untransform_tpm <- function(x) { # Function to convert the log transformed counts back into original counts # nolint: line_length_linter.
  return(ceiling((2^x) - 0.001))
}

transform_tpm <- function(x) { # Function to convert the log transformed counts back into original counts # nolint
  return(log2(x + 0.001))
}


args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # This is the SLURM_ARRAY_TASK_ID

# SOI genes
soi <- read.csv("Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_SOI.csv")
soi_genes <-soi[,2]
rm(soi)

# Metadata

# Sample meta
tss_meta <- read.csv("Data/Other/TCGA_meta/tissueSourceSite.tsv", sep = "\t")
abbrv_meta <- read.csv("Data/Other/TCGA_meta/bcrBatchCode.tsv", sep = "\t")

# Gene Metadata
gene_ids <- read.delim("Data/Other/TCGA_meta/TCGA_PanCan_TPM_Gene_Annotations.txt")

meta <- left_join(
  tss_meta %>%
    dplyr::select(c("TSS.Code", "Study.Name")) %>%
    distinct() %>%
    sapply(trimws) %>%
    as.data.frame(),
  abbrv_meta %>%
    dplyr::select(c("Study.Abbreviation", "Study.Name")) %>%
    distinct() %>%
    sapply(trimws) %>%
    as.data.frame(),
  by = "Study.Name"
)
rm(tss_meta)
rm(abbrv_meta)

# Predictor variables

# TPM counts
ori_tpm <- read.csv("Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_Full.csv")
order_tpm1 <- ori_tpm[,order(colnames(ori_tpm))] %>% 
  dplyr::select(-"id")
rm(ori_tpm)
cat("\n loaded tpm \n")

# Expected Counts
ori_exp <- read.csv("Data/RNA_Data/TCGA_Norm/tcga_gene_expected_count.csv")
order_exp0 <- ori_exp[,order(colnames(ori_exp))]
rm(ori_exp)

# Convert the Gene Ids into names
order_exp1 <- right_join(
  gene_ids %>%
    dplyr::select(c("id", "gene")) %>%
    sapply(trimws) %>%
    as.data.frame(),
  order_exp0,
  by = c("id" = "sample")
) %>%
  dplyr::select(-"id") %>%
  rename(Gene = "gene")
rm(gene_ids)
rm(order_exp0)
cat("\n loaded exp \n")

# Number of random columns to select

n_samples <- 8000

order_tpm <- order_tpm1[, 0:n_samples]
order_exp <- order_exp1[, 0:n_samples]
rm(order_tpm1)
rm(order_exp1)
rm(n_samples)

# Response data

# HRD scores
ori_hrd <- read_tsv("Data/CIN_Features/TCGA.HRD_withSampleID.txt")

# Arm level aneuploidies
ori_arm_cnv <- read_tsv("Data/CIN_Features/CNV_Data/PANCAN_ArmCallsAndAneuploidyScore_092817.txt") # nolint

trans_exp <- order_exp # Log transformed Expected Counts
rownames(trans_exp) <- NULL
trans_tpm <- order_tpm # Log transformed TPM Counts
rownames(trans_tpm) <- NULL

rm(order_exp)
rm(order_tpm)

count_exp <- trans_exp %>%
  mutate_at(vars(-1), untransform_exp)
count_tpm <- trans_tpm %>%
  mutate_at(vars(-1), untransform_tpm)

rm(trans_exp)
rm(trans_tpm)

# Remove the non-cancerous sample types from the set
codes_to_use <- c("01", "02", "03", "04", "05", "08", "09")

exp_samples_to_use <- count_exp %>%
  dplyr::select(c("Gene", ends_with(codes_to_use)))

tpm_samples_to_use <- count_tpm %>%
  dplyr::select(c("Gene", ends_with(codes_to_use)))
rm(count_exp)
rm(count_tpm)

counts_exp <- exp_samples_to_use %>%
  mutate(Gene = trimws(Gene))

# Convert data to a data.table for faster processing during the grouping of the duplicate genes # nolint
setDT(counts_exp)
setDT(tpm_samples_to_use)

# Combine duplicate genes together using the median of the expression
grouped_exp <- counts_exp[, lapply(.SD, function(x) if (length(x) > 1) median(x, na.rm = TRUE) else x), by = Gene, .SDcols = -"Gene"] # nolint
grouped_tpm <- tpm_samples_to_use[, lapply(.SD, function(x) if (length(x) > 1) median(x, na.rm = TRUE) else x), by = Gene, .SDcols = -"Gene"] # nolint
rm(counts_exp)
rm(tpm_samples_to_use)
cat("\n grouped by median \n")

# Reset the row names to the Gene names
groups_exp <- as.data.frame(grouped_exp)
groups_tpm <- as.data.frame(grouped_tpm)
rm(grouped_exp)
rm(grouped_tpm)

exp_data <- distinct(groups_exp) %>%
  column_to_rownames(var = "Gene")
tpm_data <- distinct(groups_tpm) %>%
  column_to_rownames(var = "Gene")
rm(groups_exp)
rm(groups_tpm)

# Remove samples (cols) with 0 values throughout
cols_zeros_tpm <- which(apply(tpm_data, 2, function(x) all(x == 0)))
cols_zeros_exp <- which(apply(exp_data, 2, function(x) all(x == 0)))

if (length(cols_zeros_tpm) > 0) {
  data_complete_tpm <- tpm_data[, -cols_zeros_tpm]
} else {
  data_complete_tpm <- tpm_data
}

if (length(cols_zeros_exp) > 0) {
  data_complete_exp <- exp_data[, -cols_zeros_exp]
} else {
  data_complete_exp <- exp_data
}

# Identify genes with 0 expression across all samples
zero_genes_exp <- rowSums(data_complete_exp == 0) == ncol(data_complete_exp)

# Remove genes with 0 expression across all samples
filt_dt_exp <- data_complete_exp[!zero_genes_exp, ]
rm(zero_genes_exp)

# Turn into a DGE object
matrix_exp <- filt_dt_exp %>% as.matrix()
rm(filt_dt_exp)

d_exp <- DGEList(matrix_exp)
rm(matrix_exp)

# Calculate the normalisation factor
Normfact_exp <- calcNormFactors(d_exp, method = "TMM")
rm(d_exp)
Normfactors <- as.data.frame(Normfact_exp$samples) %>% select("norm.factors")
rm(Normfact_exp)
cat("\n calced norms \n")

# match the column names from the normalisation factor df with that of the other dfs in order to divide the count values by the library size factors for the correct samples  # nolint
matching_cols_tpm <- colnames(tpm_data)[match(rownames(Normfactors), colnames(tpm_data))]

matching_cols_tpm <- matching_cols_tpm[!is.na(matching_cols_tpm)]

scld_cnts_tpm <- round(sweep(tpm_data[, matching_cols_tpm], 2, Normfactors$norm.factors, "/"), 2)
rm(tpm_data)

matching_cols_exp <- colnames(exp_data)[match(rownames(Normfactors), colnames(exp_data))]
matching_cols_exp <- matching_cols_exp[!is.na(matching_cols_exp)]

scld_cnts_exp <- round(sweep(exp_data[, matching_cols_exp], 2, Normfactors$norm.factors, "/"), 2)
rm(exp_data)

# Expected Counts
exp_set <- data_complete_exp %>% # Raw expected counts
  filter(rownames(.) %in% soi_genes) %>%
  t() %>%
  as.data.frame()
rm(data_complete_exp)

scld_exp_set <- scld_cnts_exp %>% # Library size normalised Expected counts
  filter(rownames(.) %in% soi_genes) %>%
  t() %>%
  as.data.frame()
rm(scld_cnts_exp)

log_exp <- exp_set %>% # Log transformed expected counts
  mutate_at(vars(everything()), transform_exp)

log_scld_exp <- scld_exp_set %>% # Log transformed Library size normalised expected counts
  mutate_at(vars(everything()), transform_exp)

# Transcript per million (TPM)
tpm_set <- data_complete_tpm %>% # Raw TPM counts
  filter(rownames(.) %in% soi_genes) %>%
  t() %>%
  as.data.frame()
rm(data_complete_tpm)

scld_tpm_set <- scld_cnts_tpm %>% # Library size normalised TPM counts
  filter(rownames(.) %in% soi_genes) %>%
  t() %>%
  as.data.frame()
rm(scld_cnts_tpm)

log_tpm <- tpm_set %>% # Log transformed TPM counts
  mutate_at(vars(everything()), transform_tpm)

log_scld_tpm <- scld_tpm_set %>% # Log transformed Library size normalised TPM counts
  mutate_at(vars(everything()), transform_tpm)

t_hrd <- as.data.frame(t(ori_hrd))
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ]
hrd <- as.data.frame(first_hrd[-1, ]) %>%
  mutate_all(as.numeric) %>%
  rename(loh_hrd = "hrd-loh")
rm(t_hrd)
rm(first_hrd)
cat("\n created all sets \n")

# ARM-LEVEL ANEUPLOIDIES
# Replace the NAs with 0
cvn_arm <- ori_arm_cnv %>% replace(is.na(.), 0)
rm(ori_arm_cnv)

# Re-organise the columns
cnvs_arm <- cvn_arm %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-"Type") %>%
  mutate_all(as.numeric)
rm(cvn_arm)

full_cin <- merge(
  hrd,
  cnvs_arm,
  by = "row.names"
) %>%
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")
rm(hrd)
rm(cnvs_arm)

aneu_feature_list <- colnames(full_cin[1, 6:length(full_cin)])

feature <- aneu_feature_list[[index]]
cat(paste0("\n feature:", feature, "\n"))

aneu_cat_metrics_df <- data.frame(
  RNA_Set = character(),
  Feature = character(),
  Depth = numeric(),
  Learning_Rate = numeric(),
  Gamma = numeric(),
  Logloss = numeric()
)

rna_list <- list(
  transcripts_per_million = tpm_set, # Seems to be the best performing
  scaled_transcripts_per_million = scld_tpm_set, # not too useful (scaled)
  log_scaled_transcripts_per_million = log_scld_tpm,
  log_transcripts_per_million = log_tpm,
  expected_counts = exp_set,
  scaled_expected_counts = scld_exp_set,
  log_expected_counts = log_exp,
  log_scaled_expected_counts = log_scld_exp
)


for (i in 1:length(rna_list)) {
  rna <- rna_list[[i]]
  name <- names(rna_list)[i]
  cat(paste0("\t", name, "\n"))

  full_df <- merge(rna, full_cin, by = "row.names")
  y <- as.integer(full_df[[feature]])
  X <- full_df %>% select(-c("Row.names", colnames(full_cin)))

  y[y == -1] <- 0
  y[y == 1] <- 2
  y[y == 0] <- 1

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

  grid <- expand.grid(
    max_depth = seq(1, 12, 3),
    gamma = seq(0, 2, 0.5),
    eta = seq(0.01, 0.1, 0.3)
  )

  for (j in 1:nrow(grid)) { # nolint
    cat(paste0(
      "\t\t Depth: ", grid$max_depth[j], 
      "\n")
      )

    m_xgb_untuned <- xgb.cv(
      data = xgb_data,
      nrounds = 5000,
      objective = "multi:softmax",
      early_stopping_rounds = 50,
      nfold = 3,
      max_depth = grid$max_depth[j],
      eta = grid$eta[j],
      gamma = grid$gamma[j],
      num_class = 3,
      verbose = 0
    )

    best_loss <- m_xgb_untuned$evaluation_log$test_mlogloss_mean[
      m_xgb_untuned$best_iteration
    ]

    aneu_cat_metrics_df <- rbind(aneu_cat_metrics_df, data.frame(
      Feature = feature,
      RNA_Set = name,
      Depth = grid$max_depth[j],
      Learning_Rate = grid$eta[j],
      Gamma = grid$gamma[j],
      Logloss = best_loss
    ))
  }
}
write.csv(
  aneu_cat_metrics_df,
  paste0(
    "Data/aneu_cat_xgb_metrics_params_", feature, "_", Sys.Date(), ".csv"
  )
)
cat("\n Completed processing for index: ", index, "\n")
