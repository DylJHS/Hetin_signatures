#!/bin/bash

# Specify the directory path here
base_dir="/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/Cancer_specific_data/Model_input/Parameters/Hyperparameters"

# Arrays for old and new filenames
old_names=("Best_cat_loss.csv" "Best_reg_loss.csv" "cat_full_df.csv" "reg_full_df.csv")
new_names=("base_best_cat_loss.csv" "base_best_reg_loss.csv" "base_cat_full_df.csv" "base_reg_full_df.csv")

# Loop through each subdirectory in the base directory
for dir in "$base_dir"/*/; do
  if [ -d "$dir" ]; then
    echo "Processing directory: $dir"
    
    # Loop through each file pair
    for i in "${!old_names[@]}"; do
      old_name="${old_names[$i]}"
      new_name="${new_names[$i]}"
      
      # Check if the old file exists and rename it
      if [ -f "$dir/$old_name" ]; then
        mv "$dir/$old_name" "$dir/$new_name"
        echo "  Renamed $old_name to $new_name"
      else
        echo "  $old_name not found in $dir"
      fi
    done
  fi
done

echo "Processing complete."