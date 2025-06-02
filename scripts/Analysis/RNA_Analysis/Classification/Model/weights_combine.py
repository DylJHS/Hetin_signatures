import os 
import sys
import pandas as pd

directory = "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/Model_output/local_cat"

full_df = pd.DataFrame()
for file in os.listdir(directory):
    if file.endswith(".csv"):
        file_path = os.path.join(directory, file)
        df = pd.read_csv(file_path, index_col=0)
        print(df)
        full_df = pd.concat([full_df, df], axis=0)
    else:
        print("No files found")

print(full_df)
        

        
        