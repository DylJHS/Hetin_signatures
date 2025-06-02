# This script is the second part of the base model prediction pipeline.
# Its goal is to gather the base layer predictions and merge them into
# a single data frame

# Load the packages
library(dplyr)
library(readr)
library(tidyverse)


setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")


predictions_folder <- "Data/Gen_model_output/Predictions/Base_predictions/indiv/"

# Create the full empty predictions data frame
predictions <- read.csv(
  paste0(predictions_folder, list.files(predictions_folder)[1])
) %>%
  rename("act_index" = X)

# Loop over the individual predictions csv files
for (file in list.files(predictions_folder)[2:length(list.files(predictions_folder))]) {
  cat("\n\n file \n", file, "\n\n")
  # Read the file
  preds <- read.csv(paste0(predictions_folder, file), ) %>%
    rename("act_index" = X)

  # Add the predictions to the full data frame
  predictions <- merge(
    predictions,
    preds,
    by = "act_index",
  )
}

predictions <- predictions %>% 
  column_to_rownames("act_index")

print(dim(predictions))
cat("predictions \n\n")
print(head(predictions[, 1:4]))

# Save the predictions with the row names since they represent the actual index
write.csv(predictions, "Data/Gen_model_output/Predictions/Base_predictions/Full_base_predictions.csv",
  row.names = TRUE
)

write.csv(predictions, "Data/Gen_model_input/Base_predictions/Full_base_predictions.csv",
  row.names = TRUE
)
