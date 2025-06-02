# This script is used to extract the aneuploidies of the arms from the CNV data,
# calculate the cancer-specific class weights for the arms and save them to a csv file\

# Load the necessary libraries
library(dplyr)
library(tidyverse) 

# Folder 
data_folder <- "Data/Cancer_specific_data/Model_input/Hyperparams/Arm_class_weights/"

# TCGA metadata
meta_data_folder <- "Data/Other/TCGA_meta"
tss_code_file <- read_tsv(file.path(meta_data_folder, "tissueSourceSite.tsv"),show_col_types = FALSE, na = character())
abbreviations_file <- read_tsv(file.path(meta_data_folder, "diseaseStudy.tsv"), show_col_types = FALSE)

# Arm Level Aneuploidies
chr_cnv <- read_tsv(
  "Data/CIN_Features/CNV_Data/PANCAN_ArmCallsAndAneuploidyScore_092817.txt",
  na = c("", "NA")
) %>% 
  select(-"Aneuploidy Score", -"Type") %>%
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
  mutate(
    Code = str_extract(Sample, "(?<=\\.)[^\\.]+"),
    .before = 1
  ) %>% 
  left_join(
    tss_code_file %>% 
      select("TSS Code", "Study Name"), 
    by = c("Code" = "TSS Code")
  ) %>% 
  left_join(abbreviations_file, by = "Study Name") %>%
  rename("Cancer" = "Study Abbreviation") %>%
  select(Cancer, everything(), -c(Code, "Study Name")) %>% 
  column_to_rownames("Sample")

# Loop over the cancer types
for (cancer in unique(chr_cnv$Cancer)) {
  cat("\n\n Cancer: ", cancer, "\n")
  
  # Filter the data for the cancer type
  cancer_cnv <- chr_cnv %>% 
    filter(Cancer == cancer) %>%
    select(-Cancer)
  
  # Calculate the frequency of the aneuploidies
  freq <- cancer_cnv %>%
    as.data.frame() %>%
    gather(key = "arm", value = "freq") %>%
    group_by(arm) %>%
    count(freq)%>%
    as.data.frame() %>%
    spread(key = freq, value = n) %>%
    replace(is.na(.), 0)
  
  # Calculate the weights of each class for the arms
  arm_weights <- freq %>%
    mutate(total = rowSums(select(., -arm))) %>%
    mutate_at(vars(-arm, -total), list(~1- round(. / total,2))) %>%
    mutate(total = rowSums(select(., -arm,-total))) %>%
    mutate_at(vars(-arm, -total), list(~ round(./total, 2))) %>%
    select(-total) %>%
    t() %>%
    as.data.frame() %>%
    setNames(make.unique(unlist(.[1, ]))) %>%   # Convert first row to column names
    .[-1, ]
  
  arm_order <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q",
                 "5p", "5q", "6p", "6q", "7p", "7q", "8p", "8q",
                 "9p", "9q", "10p", "10q", "11p", "11q", "12p",
                 "12q", "13q", "14q", "15q",
                 "16p", "16q", "17p", "17q", "18p", "18q",
                 "19p", "19q", "20p", "20q", "21q", "22q")
  
  # Reorder the arms and save the weights as data frame
  n_arm_weights <- arm_weights[, arm_order] %>%
    as.data.frame() %>%
    rownames_to_column("arm") %>%
    mutate(arm = case_when(
      arm == "0" ~ "Weight_loss",
      arm == "1" ~ "Weight_normal",
      arm == "2" ~ "Weight_gain",
      TRUE ~ arm  
    ),
    # Replace the 0.00 with 0.01 so as not to generate errors during the training
    across(-arm, ~ if_else(as.numeric(as.character(.)) == 0.00, 0.01, as.numeric(as.character(.))))
    )%>%
    column_to_rownames("arm") %>%
    t() %>%
    as.data.frame()
  
  # Save the weights to a csv file
  write.csv(n_arm_weights,
            paste0(data_folder,cancer, "_arm_weights.csv"),
            row.names = TRUE)
}

