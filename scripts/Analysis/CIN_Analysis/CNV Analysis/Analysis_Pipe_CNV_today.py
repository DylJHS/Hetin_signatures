import pandas as pd
import seaborn as sns
import matplotlib
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

def wilcoxon_test(column):
	stat, p_value = stats.wilcoxon(column )  # Assuming comparison against 0
	return stat, p_value
 
def sample_code(x):
	x.split('-')

def plot_sign_level(p_value):
	'''Determine the significance level for annotating the plots'''
	if p_value < 0.001:
	    sig = '***'  # Highly significant
	elif p_value < 0.01:
	    sig = '**'   # Very significant
	elif p_value < 0.05:
	    sig = '*'    # Significant
	else:
	    sig = 'n.s.'  # Not significant
	return sig

"""INITIAL DATASET CONFIGURATION
"""

# Loading the datasets with Pandas. 
cntrl_CNV_df = pd.read_csv('../../../Data/CNV/Control_Genes_CNV_df.csv')
GoI_CNV_df = pd.read_csv('../../../Data/CNV/Genes_of_Interest_CNV_df.csv')
codes = pd.read_csv('../../../Data/Codes.csv')
hist_df = pd.read_csv('../../../Data/Cancer_Groups.csv')


# Rectify the NA error in the codes dataframe.
codes['TSS Code'] = codes['TSS Code'].astype(str)
codes['TSS Code'] = codes['TSS Code'].apply(lambda x: '0' + x if len(x) == 1 else x)
codes['TSS Code'] = codes['TSS Code'].replace('nan', 'NA')


# create the dictionaries for the mappings/

TSS_mapping_dict = pd.Series(codes['Study Name'].values, index=codes['TSS Code']).to_dict()
Hist_mapping_dict = pd.Series(hist_df['Type'].values, index=hist_df['Name']).to_dict()



SoI_sample = GoI_CNV_df.iloc[:,:]
cntrl_sample = cntrl_CNV_df.iloc[:,:]




""" PART I: 1. INTER-GENOMIC EXPRESSION ANALYSIS 

Questions: 

1. Does the expression of the genes in the set of interest differ
 from the expression of other genes not in the list? 

2. If so, does this differential expression remain true for 
all Cancer types (TSS)?

"""




"""PART I.1: INITIAL DISTRIBUTION ANALYSIS

PURPOSE: Add broad metric variables (Mean, Median, Standard deviation) 
to the data and analyse the general distribution of the data based on 
group type.

"""



# Create the columns for analysis metrics in the
# set and control dataframes.

SoI_sample['Mean'] = SoI_sample.iloc[:, 1:].mean(axis=1)
SoI_sample['Median'] = SoI_sample.iloc[:, 1:-1].median(axis=1)
SoI_sample['Std'] = SoI_sample.iloc[:, 1:-2].std(axis=1)

cntrl_sample['Mean'] = cntrl_sample.iloc[:, 1:].mean(axis=1)
cntrl_sample['Median'] = cntrl_sample.iloc[:, 1:-1].median(axis=1)
cntrl_sample['Std'] = cntrl_sample.iloc[:, 1:-2].std(axis=1)


# Create column defining the groups.

SoI_sample['Group'] = 'Set'
cntrl_sample['Group'] = 'Control'


# Reduce the dataframes to the necessary columns.

reduced_SoI = SoI_sample[['Sample', 'Group', 'Mean','Median', 'Std']]
reduced_cntrl = cntrl_sample[['Sample', 'Group', 'Mean','Median', 'Std']]


print(reduced_SoI.iloc[:5,-5:])
print(reduced_cntrl.iloc[:5,-5:])

# # Plot the distribution of the the metric columns for each Dataframe.
# # This is done in order to show the underlying distribution of these 
# # metrics and hence what type of statistical test can be used when 
# # comparing the expression levels between the selected gene set and 
# # the control group(s)


	# fig = plt.figure(layout = 'constrained', figsize=(12, 5))
	# subfigs = fig.subfigures(1, 3, wspace=0.1, hspace = 0.1)



	# # Subfigure for Mean
	# axs_mean = subfigs[0].subplots(2, 1)
	# subfigs[0].set_facecolor('0.75')
	# sns.histplot(SoI_sample['Mean'], kde=True, ax=axs_mean[0])
	# axs_mean[0].set_title('Gene Set')
	# axs_mean[0].set_xlabel('')
	# sns.histplot(cntrl_sample['Mean'], kde=True, ax=axs_mean[1])
	# axs_mean[1].set_title('Control')
	# axs_mean[1].set_xlabel('log2(norm_value+1)')
	# subfigs[0].suptitle('Mean RNA Distribution', fontsize=16)

	# # Subfigure for Median
	# axs_median = subfigs[1].subplots(2, 1)
	# subfigs[1].set_facecolor('0.6')
	# sns.histplot(SoI_sample['Median'], kde=True, ax=axs_median[0])
	# axs_median[0].set_title('Gene Set')
	# axs_median[0].set_xlabel('')
	# axs_median[0].set_ylim(0,500)
	# sns.histplot(cntrl_sample['Median'], kde=True, ax=axs_median[1])
	# axs_median[1].set_title('Control')
	# axs_median[1].set_xlabel('log2(norm_value+1)')
	# axs_median[1].set_ylim(0,500)
	# subfigs[1].suptitle('Median RNA Distribution', fontsize=16)

	# # Subfigure for Std
	# axs_std = subfigs[2].subplots(2, 1)
	# subfigs[2].set_facecolor('0.75')
	# sns.histplot(SoI_sample['Std'], kde=True, ax=axs_std[0])
	# axs_std[0].set_title('Gene Set')
	# axs_std[0].set_xlabel('')
	# sns.histplot(cntrl_sample['Std'], kde=True, ax=axs_std[1])
	# axs_std[1].set_title('Control')
	# axs_std[1].set_xlabel('log2(norm_value+1)')
	# subfigs[2].suptitle('Std RNA Distribution', fontsize=16)

	# plt.show()




"""PART 1.2: INITIAL COMPARATIVE ANALYSIS

PURPOSE: Compare the groups using the broad metrics that were defined
in Part 1.

"""


# Combine the dataframes into a single df.

full_sample = pd.concat([reduced_SoI, reduced_cntrl]).reset_index(drop = True)
# print(full_sample.head())


# Reshape the data into a suitable format

melted_full = pd.melt(
	full_sample, 
	id_vars=['Sample', 'Group'], 
	value_vars=['Mean'], 
	var_name='Variable', 
	value_name='CNV'
	)


# Perform Paired t-test between the control and gene set groups.

set_means = reduced_SoI['Mean']
cntrl_means = reduced_cntrl['Mean']
ttest_stat, ttest_pvalue = stats.ttest_rel(set_means, cntrl_means)


wilcoxon_results = stats.wilcoxon(set_means, cntrl_means)
wilc_t_statistic = wilcoxon_results.statistic 
wilcoxon_pvalue = wilcoxon_results.pvalue 


print(wilc_t_statistic, wilcoxon_pvalue)


# Plot Boxplot of the Mean and Median distributions for the 
# control and gene set data.


# Create a box plot
plt.figure(figsize=(8, 8))
sns.boxplot(hue='Group', y='CNV', data=melted_full, palette='Set2')

plt.title('''Comparative Distribution of Mean CNV \n in Target vs. Random Control Gene Sets''', fontsize=14)
plt.xlabel('Metric')
plt.ylabel('Copy Number')
plt.xticks(fontweight = "bold")

# Coordinates for the line and text
mean_x1, mean_x2 = -0.2, 0.2  # x-coordinates for the line endpoints
mean_y, mean_h = full_sample['Mean'].max() + full_sample['Mean'].max()*0.2, 0.1  # y-coordinate and height for the line and text
plt.plot([mean_x1, mean_x1, mean_x2, mean_x2], [mean_y, mean_y+mean_h, mean_y+mean_h, mean_y], lw=1.5, c='black')
plt.text((mean_x1+mean_x2)*0.5, mean_y+mean_h, plot_sign_level(wilcoxon_pvalue), ha='center', va='bottom', color='black')


plt.show()



""" PART 2: INTRA GENE-SET & INTER-TUMOUR TYPE EXPRESSION ANALYSIS

Questions: 

1. Does the expression of the genes in my set list differ across/between 
tumour types? i.e. For each tumour type does the general expression of the
 genes take on a different pattern from the overall trend?

2. Are these expression differences characteristics of the tumour types 
themselves? i.e. Can the the expression levels across the genes be used to 
define individual cancer types? And if not what type so groups can be formed 
based on these gene expression levels
"""



# print(GoI_RNA_df.iloc[:5,-5:])








