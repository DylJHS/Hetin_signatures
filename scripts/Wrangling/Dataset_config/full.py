import random
import pandas as pd
import numpy as np
import dask.dataframe as dd
import subprocess
import time


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

"""
Creating a full TPM dataset with all genes and samples
"""

Full_tpm_df = pd.read_csv("../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_Full.csv")


# raw_tpm_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/tcga_RSEM_gene_tpm.txt', sep = "\t")
# TCGA_TPM_annotations_df = pd.read_csv("../../../Data/Other/TCGA_meta/TCGA_PanCan_TPM_Annotations.csv")



# Full_tpm_df = (TCGA_TPM_annotations_df[['id','gene']]
# 	.merge(raw_tpm_df, left_on = "id", right_on = "sample")
# 	.drop(columns = ['sample'])
# 	.drop_duplicates('sample')
# 	.rename(columns = {"gene": "Gene"})
# 	.reset_index(drop = True)
# 	)



# print(Full_tpm_df.iloc[:5,:5])


print(Full_tpm_df.iloc[:5,:5])
