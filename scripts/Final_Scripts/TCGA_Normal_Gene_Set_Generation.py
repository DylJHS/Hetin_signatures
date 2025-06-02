import random
import pandas as pd
import numpy as np
import subprocess

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

"""
This script generates the normal (non-cancerous) gene expression datasets from the TCGA Pan-Cancer (PANCAN) project for the exploratory analysis.
"""

nrows = None

####### TCGA Healthy SUBSET ALL GENES


## If the full TCGA Healthy Subset All Genes tpm df is already created and saved, read it in.
TCGA_Normal_full = pd.read_csv('TCGA_Normal_mRNA_TPM_Full.csv', nrows = nrows)
print("TCGA Healthy TPM DF",'\n',TCGA_Normal_full.iloc[:5,:5],'\n', TCGA_Normal_full.shape,'\n\n\n')


###Else 
### get the samples from the full data
column_names = pd.read_csv('TCGA_mRNA_TPM_Full.csv', nrows = 1)
print(column_names.iloc[:5,:5])


sampledf = pd.DataFrame(column_names.columns[2:], columns = ["Samples"])
sampledf["Participant"] = sampledf.Samples.str.split('-').str[2]
sampledf["Type"] = sampledf.Samples.str.split('-').str[3]

### list of the sample type codes that refer to a Healthy sample form sampletype
healthy_codes = ['10','12', '14' ,'11']
metastatic = ['06','07']

### select the participants whose ids types are associated with the Healthy types but not metastatic types
## First select those participants with codes in the Healthy types
HealthyParticipants = sampledf[sampledf.Type.isin(healthy_codes)].Participant.to_list()
print(len(HealthyParticipants))

## Then select those participants with codes in the metastatic types
MetastaticParticipants = sampledf[sampledf.Type.isin(metastatic)].Participant.to_list()
# print(len(MetastaticParticipants))

## Then select those Samples with participants in the first but not the second list above
UsableSamples = sampledf[(sampledf.Participant.isin(HealthyParticipants)) & (~sampledf.Participant.isin(MetastaticParticipants))]
# print(UsableSamples.sort_values('Participant').head())


participants2keep = UsableSamples[UsableSamples.Type.isin(healthy_codes)].Samples.to_list()
print(len(participants2keep))

cols_to_use = column_names.columns[:2].to_list() + participants2keep

### tpm dataset with all genes but subset of samples
TCGA_Normal_full = pd.read_csv('TCGA_mRNA_TPM_Full.csv', usecols = cols_to_use)
print(TCGA_Normal_full.shape)


### save the Healthy Subset with all Genes
## TCGA_Normal_full.to_csv('TCGA_mRNA_TPM_Normal_Full.csv', index = False)



###### TCGA Healthy SOI SUBSET

# If the TCGA Healthy Set of Interest tpm df is already created and saved, read it in.
TCGA_Normal_SOI_tpm_df = pd.read_csv('TCGA_Normal_mRNA_TPM_SOI.csv', nrows = nrows)
print("TCGA Healthy SOI SUBSET",'\n',TCGA_Normal_SOI_tpm_df.iloc[:5,:5],'\n', TCGA_Normal_SOI_tpm_df.shape,'\n\n\n')


###Else Create it

TCGA_Normal_SOI_tpm_df = TCGA_Normal_full[
	(TCGA_Normal_full['id'].isin(TCGA_SOI_tpm_df.id.to_list())) | 
	(TCGA_Normal_full['Gene'].isin(TCGA_SOI_tpm_df.Gene.to_list()))
	].reset_index(drop = True)

print("TCGA Healthy SOI SUBSET",'\n',TCGA_Normal_SOI_tpm_df.iloc[:5,:5],'\n', TCGA_Normal_SOI_tpm_df.shape,'\n\n\n')

# # Save the TCGA Healthy SOI Tpm df.
#TCGA_Normal_SOI_tpm_df.to_csv('TCGA_Normal_mRNA_TPM_SOI.csv', index = False)



###### TCGA Healthy ALL CONTROL GENES SUBSET


# If the TCGA Healthy All Controls tpm df is already created and saved, read it in.
TCGA_Normal_All_Ctrl_df = pd.read_csv('TCGA_Normal_mRNA_TPM_CTRL.csv', nrows = nrows)
print("TCGA Healthy ALL CONTROL GENES",'\n',TCGA_Normal_SOI_tpm_df.iloc[:5,:5],'\n', TCGA_Normal_SOI_tpm_df.shape,'\n\n\n')


###Else Create it

# Create the df with the genes which are not in the SOI set
TCGA_Normal_All_Ctrl_df = TCGA_Normal_full[
	(~TCGA_Normal_full['id'].isin(TCGA_SOI_tpm_df.id.to_list())) & 
	(~TCGA_Normal_full['Gene'].isin(TCGA_SOI_tpm_df.Gene.to_list()))
	].reset_index(drop = True)

print("TCGA Healthy ALL CONTROL GENES",'\n',TCGA_Normal_All_Ctrl_df.iloc[:5,:5],'\n', TCGA_Normal_All_Ctrl_df.shape,'\n\n\n')

## Save the TCGA Healthy All Controls Tpm df.
# TCGA_Normal_All_Ctrl_df.to_csv('TCGA_Normal_mRNA_TPM_CTRL.csv', index = False)



###### INDIVIDUAL SMALLER CONTROL SETS

# Need to read through each of the control sets that were created for the TCGA Cancerous controls using the SOI_and_Control_Set_Configuration and for each of them
# create a Healthy set using the same exact genes based on the full Healthy Control samples
# repeat until each of the small control sets has been viewed and their corresponding Healthy types have been created


df_number = 1


while True:

	# Find and read the genes that are part of the corresponding Tumour control set.
	Cancer_ctrl_geneids = pd.read_csv(f'TCGA_TPM_RNA_Control_df{df_number}.csv', usecols = ["id"]).id.to_list()

	Normal_ctrl_set = TCGA_Normal_All_Ctrl_df[TCGA_Normal_All_Ctrl_df.id.isin(Cancer_ctrl_geneids)].reset_index(drop = True)
	print(f"TCGA NORMAL CONTROL SET {df_number}",'\n',Normal_ctrl_set.iloc[:5,:5],'\n', Normal_ctrl_set.shape,'\n\n\n')

	## Save the TCGA Healthy Control Sets.
	Normal_ctrl_set.to_csv(f'TCGA_Normal_mRNA_TPM_CTRL_Set{df_number}.csv', index = False)


	df_number += 1
























