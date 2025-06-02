import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
import numpy as np

sns.set_theme(style="whitegrid")

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


"""INITIAL DATASET CONFIGURATION
"""

# Loading the datasets with Pandas. 
cntrl_RNA_df = pd.read_csv('../../Data/RNA/Control_Genes_RNA_df.csv')
GoI_RNA_df = pd.read_csv('../../Data/RNA/Genes_of_Interest_RNA_df.csv')
codes = pd.read_csv('../../Data/Codes.csv')
hist_df = pd.read_csv('../../Data/Cancer_Groups.csv')



codes['TSS Code'] = codes['TSS Code'].astype(str)
codes['TSS Code'] = codes['TSS Code'].apply(lambda x: '0' + x if len(x) == 1 else x)
codes['TSS Code'] = codes['TSS Code'].replace('nan', 'NA')


TSS_mapping_dict = pd.Series(codes['Study Name'].values, index=codes['TSS Code']).to_dict()
Hist_mapping_dict = pd.Series(hist_df['Type'].values, index=hist_df['Name']).to_dict()



"""PART 1: INITIAL DISTRIBUTION ANALYSIS

PURPOSE: Add broad metric variables (Mean, Median, Standard deviation) 
to the data and analyse the general distribution of the data based on 
group type.

"""


# Create the columns for analysis metrics in the
# set and control dataframes.

Stats_Set_RNA_df = GoI_RNA_df.iloc[:,:]
Stats_Cntrl_RNA_df = cntrl_RNA_df.iloc[:,:]
print(Stats_Set_RNA_df.iloc[:5,-5:])

Stats_Set_RNA_df['Mean'] = Stats_Set_RNA_df.iloc[:, 1:].mean(axis=1)
Stats_Set_RNA_df['Median'] = Stats_Set_RNA_df.iloc[:, 1:-1].median(axis=1)
Stats_Set_RNA_df['Std'] = Stats_Set_RNA_df.iloc[:, 1:-2].std(axis=1)

Stats_Cntrl_RNA_df['Mean'] = Stats_Cntrl_RNA_df.iloc[:, 1:].mean(axis=1)
Stats_Cntrl_RNA_df['Median'] = Stats_Cntrl_RNA_df.iloc[:, 1:-1].median(axis=1)
Stats_Cntrl_RNA_df['Std'] = Stats_Cntrl_RNA_df.iloc[:, 1:-2].std(axis=1)
 
# Create column defining the groups.

Stats_Set_RNA_df['Group'] = 'Set'
Stats_Cntrl_RNA_df['Group'] = 'Control'


# Reduce the dataframes to the necessary columns.

reduced_set = Stats_Set_RNA_df[['Sample', 'Group', 'Mean','Median', 'Std']]
cntrl_RNA_df_red = Stats_Cntrl_RNA_df[['Sample', 'Group', 'Mean','Median', 'Std']]


## FULL SAMPLE CREATION


# Merge the reduced dataframes for further analysis in R 

Reduced_full = pd.merge(reduced_set, cntrl_RNA_df_red, how = 'inner', on= 'Sample', suffixes=('_set', '_cntrl'))

print(Reduced_full.head())


# Add the tissue source type to the dataframe.

Reduced_full['TSS'] = Reduced_full.Sample.apply(lambda x: x.split('-')[1])
Reduced_full['TSS'] = Reduced_full['TSS'].map(TSS_mapping_dict)


# Add the Histology Group to the dataframe.

Reduced_full['Hist_group'] = Reduced_full['TSS'].map(Hist_mapping_dict)
TSS_abbr_mapping_dict = pd.Series(hist_df['Abbrvs'].values, index=hist_df['Name']).to_dict()
Reduced_full['TSS_Abbrvs'] = Reduced_full['TSS'].map(TSS_abbr_mapping_dict)

# print(Reduced_full.head())
# Reduced_full.to_csv('../../Data/RNA/Full_RNA_Sample_metrics.csv', index=False)


### SET OF INTEREST MODIFICATIONS


# Add the tissue source type to the dataframe.

# GoI_RNA_df['TSS'] = GoI_RNA_df.Sample.apply(lambda x: x.split('-')[1])
# GoI_RNA_df['TSS'] = GoI_RNA_df['TSS'].map(TSS_mapping_dict)


# # Add the Histology Group to the dataframe.

# GoI_RNA_df['Hist_group'] = GoI_RNA_df['TSS'].map(Hist_mapping_dict)
# TSS_abbr_mapping_dict = pd.Series(hist_df['Abbrvs'].values, index=hist_df['Name']).to_dict()
# GoI_RNA_df['TSS_Abbrvs'] = GoI_RNA_df['TSS'].map(TSS_abbr_mapping_dict)


# print(GoI_RNA_df.iloc[:5,-9:])

