"""
XENABROWSER TCGA CNV DATASET IMPORT AND TRANSFORMATION

This script transforms raw TCGA CNV data into a readable CSV format,
where Genes are represented as columns and Samples as rows. It relies 
on the mRNA datasets to retrieve the same Genes of Interest (GoIs) and 
Random Control Genes (RCGs) in order to initialise the dataframe without
having to query the STRING DB.


The script follows these main steps:
1. Obtain pre-defined Set of Interest (SoI) and Control set
	from the RNA datasets that were created
2. Read and reformat raw CNV TCGA data.
3. Create the dataset focusing only on genes of interest.
4. Create the control dataset with the control set, and some additional
	random control genes in case the SoI and Control set differ in 
	gene number due to differences in the CNV and RNA datasets.
5. Convert data to numeric format and save the resulting datasets 
as CSV files.

Note: 
- The script assumes the access to the raw CNV TCGA file as well as
 the two additional Set-of-interest, and control set RNA datasets
- Output data is saved in CSV format with specific naming conventions.
"""

import random
import multiprocessing
import time 
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from io import StringIO


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


def process_control_data(n):
	''' Loop over the mRNA control datasets.
	'''

	start_time = time.time()


	# Load the mRNA control dataset for the current iteration
	cntrl_RNA_df = pd.read_csv(f'../../Data/RNA_Data/Control_Genes_RNA_Data_df{n}.csv')

	# list containing the same genes as the RNA control df
	cntrl_list = list(cntrl_RNA_df.columns[1:])

	# List for sample CNV data from the random control genes
	ctrl_set = []
	leftover_ctrl_set = []

	for gene_datum in reformatted_data:
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
	    .rename(columns={'index': 'Sample'})
	    .drop_duplicates('Sample')
	)

	# Convert all values in both dataframes to float, handling non-numeric values as NaN.

	for col in ctrl_set_CNV_df.columns[1:]:
	    ctrl_set_CNV_df[col] = pd.to_numeric(ctrl_set_CNV_df[col], errors='coerce')

	# Save the formatted gene set and control set dataframes as CSV files.
	ctrl_set_CNV_df.to_csv(f'../../Data/CNV_Data/Control_Genes_CNV_Data_df{n}_parall.csv', index=False)
	print(f"Processed control dataset {n}")

	end_time = time.time()
	elapsed_time = end_time - start_time
	print(f"Generated control dataset {n} in {elapsed_time:.2f} seconds")




# Load the mRNA Set-of-Interest dataset. 
SoI_RNA_df = pd.read_csv('../../Data/RNA_Data/Genes_of_Interest_RNA_Data_df.csv')


# list containing the same set of genes 
# as the mRNA dataframe.

SoI_list = list(SoI_RNA_df.columns[1:])


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






if __name__ == '__main__':

	# Define the number of processes you want to run concurrently
	num_processes = 4  # You can adjust this based on your CPU's capabilities
	pool = multiprocessing.Pool(num_processes)

	# Iterate over control datasets and process them in parallel
	for n in range(190, 193 + 1):
	    pool.apply_async(process_control_data, args=(n,))

	# Close the pool and wait for all processes to finish
	pool.close()
	pool.join()

