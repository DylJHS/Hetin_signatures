# This script is used to extract the aneuploidies of the arms from the CNV data,
# calculate the class weights for the arms and save them to a csv file

#Set working directory
setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

# Load the required libraries
library(tidyverse)

# Arm Level Aneuploidies
# Load the data
chr_cnv <- read_tsv(
  "Data/CIN_Features/CNV_Data/PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
)%>% 
  replace(is.na(.), 0) %>%
  select(-"Aneuploidy Score", -"Type") %>%
  mutate(Sample = str_replace_all(Sample, "-", "\\.")) %>% 
  column_to_rownames("Sample") %>%
  mutate_all(~ replace(., . == 1, 2)) %>% 
  mutate_all(~ replace(., . == 0, 1)) %>% 
  mutate_all(~ replace(., . == -1, 0))

# Calculate the frequency of the aneuploidies
freq <- chr_cnv %>% 
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

# replace any values under 0.1 with 0.1
arm_weights[arm_weights < 0.1] <- 0.1

arm_order <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", 
               "5p", "5q", "6p", "6q", "7p", "7q", "8p", "8q", 
               "9p", "9q", "10p", "10q", "11p", "11q", "12p", 
               "12q", "13 (13q)", "14 (14q)", "15 (15q)", 
               "16p", "16q", "17p", "17q", "18p", "18q", 
               "19p", "19q", "20p", "20q", "21 (21q)", "22 (22q)")

# Reorder the arms and save the weights as data frame
n_arm_weights <- arm_weights[, arm_order] %>%
  rownames_to_column("arm") %>%
  mutate(arm = as.integer(arm)-1) %>%
  column_to_rownames("arm") %>% 
  t() %>% 
  as.data.frame() 

#save the weights to a csv file
# write.csv(arm_weights, "Data/CIN_Features/CNV_Data/arm_weights.csv", row.names = TRUE)



# PER CANCER TYPE
