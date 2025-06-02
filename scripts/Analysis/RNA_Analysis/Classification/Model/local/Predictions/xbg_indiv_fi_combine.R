# This script is the second part of the base model prediction pipeline.
# Its goal is to gather the base layer individual feature importances and merge them into
# a single data frame

# Load the packages
library(dplyr)
library(readr)
library(tidyverse)


setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")


fi_folder <- "Data/Gen_model_output/Results/Feature_importance"
indiv_folder <- file.path(fi_folder, "Base_feat_imp/")

# Create the full empty feat_imp data frame
feat_imp <- data.frame()

# Loop over the individual feat_imp csv files
for (file in list.files(indiv_folder)[1:length(list.files(indiv_folder))]) {
  cat("\n\n file \n", file, "\n\n")
  # Read the file
  feat_file <- read.csv(paste0(indiv_folder, file),
                        row.names = 1) %>% 
    mutate(Target = strsplit(file, "_feature_importance.csv")[[1]][1])
  
  # Add the feat_imp to the full data frame
  feat_imp <- rbind(
    feat_imp,
    feat_file
  )
}


print(dim(feat_imp))
cat("feat_imp \n\n")
print(head(feat_imp))

# Save the feat_imp with the row names since they represent the actual index
write.csv(feat_imp, file.path(fi_folder, "Full_base_feature_importance.csv"),
          row.names = TRUE
)