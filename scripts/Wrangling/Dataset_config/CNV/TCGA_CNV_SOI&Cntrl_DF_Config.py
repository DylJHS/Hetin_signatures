"""
XENABROWSER TCGA CNV DATASET IMPORT AND TRANSFORMATION

This script transforms raw TCGA CNV data into a readable CSV format,
where Genes are represented as columns and Samples as rows. It relies 
on the mRNA datasets to retrieve the same Genes of Interest (GoIs) and 
Random Control Genes (RCGs) in order to initialise the dataframe without
having to query the STRING DB.


The script follows the steps below:
1. Obtain the pre-defined Set of Interest (SoI) list from the 
	corresponding mRNA dataset.
2. Read and reformat raw CNV TCGA data.
3. Create the SoI CNV dataset.
4. Obtain the control set lists from the RNA control 
	datasets that were created.
4. Create each CNV control dataset with the same RNA control set genes, 
	and some additional random control genes in case the CNV SoI 
	and Control set differ in gene number due to differences 
	in the CNV and RNA datasets.
5. Convert data to numeric format and save the resulting datasets 
as CSV files.

Note: 
- The script assumes the access to the raw CNV TCGA file as well as
 the two additional Set-of-interest, and control set RNA datasets
- Output data is saved in CSV format with specific naming conventions.
"""

import random
import time 
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from io import StringIO


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# Load the mRNA Set-of-Interest dataset. 
# list containing the same set of genes 
# as the mRNA dataframe.

SoI_list =  pd.read_csv('../../Data/RNA_Data/Genes_of_Interest_RNA_Data_df.csv', index_col=0, nrows=0).columns.tolist()


# Read the raw CNV dataset.

with open('../../Data/CNV_Data/TCGA_PanCan_CNV.txt', 'r') as file:
	raw_data = file.readlines()


# Process raw data to separate sample headers 
# and reformat the data into a usable structure.

samples = raw_data[0].strip('\n').split('\t')
gene_data = raw_data[1:]
reformated_data = [gene_datum.strip('\n').split('\t') for gene_datum in gene_data]


# List for sample CNV data from the set-of-interest.

SoI = []

for gene_datum in reformated_data:
    if gene_datum[0] in SoI_list:
        SoI.append(gene_datum)


# Construct the dataframe based on the genes of interest, 
# with genes as rows and samples as columns.
# Values are the CNVs

SoI_CNV_df = (
    pd.DataFrame(SoI, columns=samples)
	.set_index('Sample')
	.T
	.reset_index()
	.rename_axis(None, axis=1)
	.rename(columns = {'index':'Sample'})
	.drop_duplicates('Sample')
)


# Convert all values in dataframe to float, handling non-numeric values as NaN.

for col in SoI_CNV_df.columns[1:]:
	SoI_CNV_df[col] = pd.to_numeric(SoI_CNV_df[col], errors='coerce')


# Save the formatted gene set and control set dataframes as CSV files.

# SoI_CNV_df.to_csv('../../Data/CNV_Data/Genes_of_Interest_CNV_Data_df.csv', index = False)

print(SoI_CNV_df.shape)





# Loop over the mRNA control datasets.


# Set the first RNA control df to obtain

n = 411

while True:

	try:
		start_time = time.time()

		# Load the mRNA control dataset. 
		# list containing the same genes as the RNA control df.

		cntrl_list =  pd.read_csv(f'../../Data/RNA_Data/Control_Genes_RNA_Data_df{n}.csv', index_col=0, nrows=0).columns.tolist()


		# List for sample CNV data from the random control genes

		ctrl_set = []
		leftover_ctrl_set = []


		for gene_datum in reformated_data:
		    if gene_datum[0] in cntrl_list:
		        ctrl_set.append(gene_datum)
		    elif gene_datum[0] not in SoI_list:
		    	leftover_ctrl_set.append(gene_datum)


		# If the number of genes in the set of interest exceeds that of
		# the control-set then add the difference from the leftover control
		# list using random selection.

		if len(SoI) > len(ctrl_set):
			gene_dif = len(SoI) - len(ctrl_set)
			random_cntrls = random.sample(leftover_ctrl_set, k=gene_dif)
			ctrl_set = ctrl_set + random_cntrls


		# Construct the dataframe based on the random control genes, with genes 
		# as rows and samples as columns. Values are the CNVs.

		ctrl_set_CNV_df = (
		    pd.DataFrame(ctrl_set, columns=samples)
		    .drop_duplicates()
		    .set_index('Sample')
		    .T
		    .reset_index()
		    .rename_axis(None, axis=1)
		    .rename(columns = {'index':'Sample'})
		    .drop_duplicates('Sample')
		)


		# Convert all values in both dataframes to float, handling non-numeric values as NaN.

		for col in ctrl_set_CNV_df.columns[1:]:
			ctrl_set_CNV_df[col] = pd.to_numeric(ctrl_set_CNV_df[col], errors='coerce')


		# Save the formatted gene set and control set dataframes as CSV files.

		ctrl_set_CNV_df.to_csv(f'../../Data/CNV_Data/Control_Genes_CNV_Data_df{n}.csv', index = False)

		end_time = time.time()
		elapsed_time = end_time - start_time
		print(f"\nGenerated control dataset {n} in {elapsed_time:.2f} seconds\n")
		n += 1

	except FileNotFoundError:
		print("\n\n\n Last dataset created \n\n\n")
		break






