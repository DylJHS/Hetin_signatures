import random
import pandas as pd



pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

"""
TCGA RNA TPM Configuration Script

The purpose of the script is to clean the TCGA RNA TPM dataset so that
it contains the genes that are used in the original script which used the batch 
corrected dataset. 
This includes:
- A dataset based solely on the Genes of Interest
- Multiple controls whose number is the same as the Cancerous ones
and whose composition is as close to the TCGA controls as possible. 

There are 60500 identifiers (genes) in the TPM file but only
20500 in the batch effect one. 
So will take longer too 
"""

# Reading the full RNA TPM dataset
TPM_Samp_df = pd.read_csv("../../../Data/RNA_Data/TCGA_mRNA_TPM.csv")

# Reading the  RNA TPM metadata dataset that contains the gene names
TPM_annotations_df = pd.read_csv("../../../Data/Other/TCGA_PanCan_TPM_Annotations.csv")

# Load the Genes that are being used in the original TCGA analysis 
Original_Genes =  pd.read_csv("../../../Data/RNA_Data/TCGA_mRNA_Norm.csv", usecols = [0])


Full_df = (TPM_annotations_df[['id','gene']]
	.merge(TPM_Samp_df, left_on = "id", right_on = "sample")
	.drop(columns = ['id', 'sample'])
	.rename(columns = {'gene':'Gene'})
	)


# Save the full TCGA TPM dataset with the gene names
Full_df.to_csv('../../../Data/RNA_Data/TCGA_mRNA_TPM_Gened_Full.csv', index = False)


Original_Genes_list = sorted(list(Original_Genes['sample']))

Filtered_df = Full_df[Full_df['Gene'].isin(Original_Genes_list)].reset_index(drop = True)

# print(Filtered_df.iloc[:5,:5])


# Save the Filtered TCGA TPM dataset with the gene names
# Filtered_df.to_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_Gened_Filtered.csv', index = False)



# Load the Genes that are being used in the original TCGA analysis 
Original_Genes =  pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_Gened_Filtered.csv')
# print(Original_Genes.shape)


Original_Genes = Original_Genes.groupby('Gene').mean().reset_index().set_index('Gene')
# print(Original_Genes.shape)
# print(Original_Genes.iloc[:5,:5])


# Get the Genes of Interest that are used in the TCGA PanCan 

SoI_list =  pd.read_csv('../../../Data/RNA_Data/TCGA_Batch/Genes_of_Interest_RNA_Data_df.csv', index_col=0, nrows=0).columns.tolist()
SoI_list = [word.upper() for word in SoI_list]

print(len([item for item in list(Original_Genes.index) if item in SoI_list]))
## 307 Genes are common to both the TCGA SoI and the TCGA TPM gene-set ###


# Create the dataframe for the GTEx using the genes of interest.

TPM_SOI_df = (Original_Genes[Original_Genes.index.isin(SoI_list)]
	.T
	.reset_index(drop = False)
	.rename(columns = {'index':'SampleID'})
	.rename_axis(None, axis = 1)
	.set_index('SampleID')
	)

print(TPM_SOI_df.shape)

for col in TPM_SOI_df.columns[1:]:
	TPM_SOI_df[col] = pd.to_numeric(TPM_SOI_df[col], errors='coerce')

# TPM_SOI_df.to_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_RNA_TPM_SOI.csv')



# Create the dataframe for the GTEx using the genes that are not of interest.

TPM_ctrl_df = (Original_Genes[~Original_Genes.index.isin(SoI_list)]
	.T
	.reset_index(drop = False)
	.rename(columns = {'index':'SampleID'})
	.rename_axis(None, axis = 1)
	.set_index('SampleID')
	)


for col in TPM_ctrl_df.columns[1:]:
	TPM_ctrl_df[col] = pd.to_numeric(TPM_ctrl_df[col], errors='coerce')


# TPM_ctrl_df.to_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_RNA_TPM_Ctrl_Full.csv')

print(TPM_ctrl_df.iloc[:5,:5])
print(TPM_ctrl_df.shape)


# Create the control sets that match up as closely as possible to the GTEx controls

If the Full control set has already been created
TPM_ctrl_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_RNA_TPM_Ctrl_Full.csv')

n = 1

while n < 413:
	print('\n',n)

	try:

		# Load the mRNA control dataset. 
		# list containing the same genes as the RNA control df.
		GTEx_ctrl =  pd.read_csv(f'../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM_Healthy/Controls/GTEx_Control_RNA_Data_df{n}.csv', 
			index_col=0, nrows=0).columns.tolist()
		GTEx_ctrl_list =  ['SampleID'] + GTEx_ctrl

		TPM_ctrl_set = TPM_ctrl_df[[col for col in GTEx_ctrl_list if col in TPM_ctrl_df]]
		TPM_ctrl_set_list = list(TPM_ctrl_set.columns)

		leftover_ctrl_set = TPM_ctrl_df[[col for col in TPM_ctrl_df if col not in GTEx_ctrl_list]]
		leftover_ctrl_set_list = list(leftover_ctrl_set.columns)

		# If the number of genes in the GTEx control exceeds that of
		# the TPM then add the difference from the leftover control
		# list using random selection.

		gene_dif = len(GTEx_ctrl_list) - len(TPM_ctrl_set_list)
		if gene_dif > 0:	

			random_genes = random.sample(leftover_ctrl_set_list, k=gene_dif)
			all_ctrls = list(TPM_ctrl_set_list + random_genes)

			TPM_ctrl_set_list_new = TPM_ctrl_df[all_ctrls]

			print(TPM_ctrl_set_list_new.iloc[:5,:5], '\n', TPM_ctrl_set_list_new.shape, '\n')


			# Save the formatted gene set and control set dataframes as CSV files.

			# TPM_ctrl_set_list_new.to_csv(f'../../../Data/RNA_Data/TCGA_TPM/TPM_RNA_Control_Data_df{n}.csv', index = False)

		n += 1

	except FileNotFoundError:
		print("\n\n\n Last dataset created \n\n\n")
		break



