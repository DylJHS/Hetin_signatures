

'''

This Script is a python script which serves to clean and preprocess the
TCGA TPM RNA data that contains the cancerous samples that are to be used in 
the Part one of the Analysis. 

'''



import pandas as pd
import dask.dataframe as dd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.stats.multitest as smm


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


# Loading of the full RNA TPM dataset.
Data = dd.read_csv(
	'../../../Data/RNA_Data/TCGA_mRNA_TPM.csv',
    storage_options={'anon': True, 'use_ssl': False},
    sample = 100000000
    )

# Loading of the full RNA TCGA Annotation dataset.
Annotations = pd.read_csv('../../../Data/Other/TCGA_PanCan_TPM_Annotations.csv')


print(Annotations.head())
print(Data.columns)
# The Curent number of samples in the TCGA Dataset stands at 10536


# Get the Gene names from the Annotation set 
# that correspond to the sample ids in the Data

# Healthy_Subjects = GTEx_meta_df[GTEx_meta_df.DTHHRDY.isin([ 1, 2, 3, 0])].SUBJID.unique()

# # Get the column names from the full RNA GTEX dataset.
# cols = list(Data.columns[3:])

# # Get the columsn that should be dropped from the dataset.
# cols_to_drop = [subj for subj in cols if subj[:10] not in Healthy_Subjects]

# # Drop the columns.
# GTEx_Samp_df = GTEx_Samp_df.drop(columns = cols_to_drop)
# GTEx_Samp_df = GTEx_Samp_df.drop(columns = ['Unnamed: 0'])

# # Save the new dataset
# # GTEx_Samp_df.to_csv('../../../Data/RNA_Data/GTEx_RNA_Healthy/GTEx_RNA_Healthy*.csv',
# #     index = False)



