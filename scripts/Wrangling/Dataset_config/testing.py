import random
import pandas as pd
import time

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)




with open('/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/RNA_Data/TCGA_Norm/tcga_RSEM_Hugo_norm_count', 'r') as file:
	raw_data = file.readlines()

# print(raw_data)




# print(raw_data[0])
# print(raw_data[1])
print(raw_data[2])




# Process raw data to separate sample headers 
# and reformat the data into a usable structure.

columns = raw_data[0].strip('\n').split('\t')
print(columns)
# values = raw_data[1:]
# reformated_values = [item.strip('\n').split('\t') for item in values]

# GTEx_tpm_df.iloc[:, 2:] = np.log2(raw_GTEx_tpm_df.iloc[:, 2:].replace(0, 0.001).apply(pd.to_numeric, errors='coerce'))