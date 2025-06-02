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


# Loadidng the static datasets.

GoI_RNA_df = pd.read_csv('../../Data/RNA_Data/Genes_of_Interest_RNA_df.csv')
GoI_CNV_df = pd.read_csv('../../Data/CNV_Data/Genes_of_Interest_CNV_df.csv')

codes = pd.read_csv('../../Data/Other/Codes.csv')
hist_df = pd.read_csv('../../Data/Other/Cancer_Groups.csv')


# Define the TSS and Histology groups and dataframes 

codes['TSS Code'] = codes['TSS Code'].astype(str)
codes['TSS Code'] = codes['TSS Code'].apply(lambda x: '0' + x if len(x) == 1 else x)
codes['TSS Code'] = codes['TSS Code'].replace('nan', 'NA')

TSS_mapping_dict = pd.Series(codes['Study Name'].values, index=codes['TSS Code']).to_dict()
Hist_mapping_dict = pd.Series(hist_df['Type'].values, index=hist_df['Name']).to_dict()
TSS_abbr_mapping_dict = pd.Series(hist_df['Abbrvs'].values, index=hist_df['Name']).to_dict()


for df in [GoI_CNV_df,GoI_RNA_df]:
	df['Mean'] = df.iloc[:, 1:].mean(axis=1)
	df['Median'] = df.iloc[:, 1:-1].median(axis=1)
	df['Std'] = df.iloc[:, 1:-2].std(axis=1)

	df['Group'] = 'Set'

	reduced_df = df[['Sample', 'Group', 'Mean','Median', 'Std']]


	# Commence the iteration over the different control sets.

	# Set the intial df number.

	n = 1


	while True:

		# Loading the control datasets. 

		if df is GoI_RNA_df:
			ctrl_df = pd.read_csv(f'../../Data/RNA_Data/Control_Genes_RNA_Data_df{n}.csv')
			sort = "RNA"
		else:
			ctrl_df = pd.read_csv(f'../../Data/CNV_Data/Control_Genes_CNV_Data_df{n}.csv')
			sort = "CNV"

		# Create the columns for analysis metrics in the
		# set and control dataframes.

		ctrl_df['Mean'] = ctrl_df.iloc[:, 1:].mean(axis=1)
		ctrl_df['Median'] = ctrl_df.iloc[:, 1:-1].median(axis=1)
		ctrl_df['Std'] = ctrl_df.iloc[:, 1:-2].std(axis=1)


		# Create column defining the groups.

		ctrl_df['Group'] = 'Control'


		# Reduce the dataframes to the necessary columns.

		reduced_cntrl = ctrl_df[['Sample', 'Group', 'Mean','Median', 'Std']]


		## FULL SAMPLE CREATION


		# Merge the reduced dataframes for further analysis in R 

		reduced_full = pd.merge(reduced_df, reduced_cntrl, how = 'inner', on= 'Sample', suffixes=('_set', '_cntrl'))


		# Add the tissue source type to the dataframe.

		reduced_full['TSS'] = reduced_full.Sample.apply(lambda x: x.split('-')[1])
		reduced_full['TSS'] = reduced_full['TSS'].map(TSS_mapping_dict)


		# Add the Histology Group to the dataframe.

		reduced_full['Hist_group'] = reduced_full['TSS'].map(Hist_mapping_dict)
		reduced_full['TSS_Abbrvs'] = reduced_full['TSS'].map(TSS_abbr_mapping_dict)

		print('\n', n)


		# Save the datasets.

		reduced_full.to_csv(f'../../Data/{sort}_Data/Full_{sort}_metrics{n}.csv', index=False)

		n +=1 


