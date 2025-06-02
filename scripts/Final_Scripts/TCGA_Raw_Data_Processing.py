"""
This script performs a Dataframe creation. 
Taking the raw TCGA dataset (TPM or Expected Counts) and configuring it to contain 
all the genes as rows and the samples as columns.

"""

import pandas as pd
from io import StringIO


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


# Read raw data from a file and store each line in raw_data.

with open('', 'r') as file:
	raw_data = file.readlines()


# Process raw data to separate sample headers 
# and reformat the data into a usable structure.

columns = raw_data[0].strip('\n').split('\t')
print(columns)
values = raw_data[1:]
reformated_values = [item.strip('\n').split('\t') for item in values]

print(len(reformated_values),'\n')


# # # Construct a dataframe from genes of interest, 
# # # with genes as rows and samples as columns.
Df = pd.DataFrame(reformated_values, columns=columns)

print(Df.iloc[:10,:10])


# # Save the formatted gene set and control set dataframes as CSV files.
Df.to_csv('', index = False)


