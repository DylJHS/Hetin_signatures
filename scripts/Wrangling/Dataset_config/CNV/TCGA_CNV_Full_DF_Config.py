"""
This script performs a Dataset Dataframe creation. 
Taking the raw CNV dataset and configuring it to contain 
all the genes as rows and the samples as columns.

"""

import pandas as pd
from io import StringIO


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


# Read raw CNV data from a file and store each line in raw_data.

with open('../../Data/CNV_Data/TCGA_PanCan_CNV.txt', 'r') as file:
	raw_data = file.readlines()


# Process raw data to separate sample headers 
# and reformat the data into a usable structure.

samples = raw_data[0].strip('\n').split('\t')
data = raw_data[1:]
reformated_data = [item.strip('\n').split('\t') for item in data]

# print(len(reformated_data))
print(samples[:3])


# # Construct a dataframe from genes of interest, 
# # with genes as rows and samples as columns.

CNV_df = (pd.DataFrame(reformated_data, columns=samples)
	.rename(columns = {'Sample':'Gene'})
	.set_index("Gene")
	.apply(pd.to_numeric, errors='coerce')
	)

print(CNV_df.iloc[:5,:5])


# Save the formatted gene set and control set dataframes as CSV files.

# CNV_df.to_csv('../../Data/CNV_Data/CNV_Full.csv', index = True)



