import pandas as pd
import dask.dataframe as dd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.dataset as ds
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.stats.multitest as smm
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


## OPEN DFs WITH PANDAS
df = pd.read_csv('../Data/TCGA_PanCan_CNV_3.csv')
codes = pd.read_csv('../Data/Codes.csv')
# print(codes.head())


df = df.T
headers = df.iloc[0]
df = df.iloc[1:]
df.columns = headers
df.columns.name = None
df = df.convert_dtypes()

for col in df.columns:
	df[col] = pd.to_numeric(df[col], errors='coerce')

# Convert the series to a float data type
df = df.astype(float)
print(df.iloc[:5,:5])



###													CNV TYPE (GAIN/LOSS) ANALYSIS      									####
##GET THE VALUE COUNTS OF THE CNV TYPE AND PLOT W BARPLOT
new_df = df.reset_index()
G_L_df = pd.melt(new_df, id_vars='index', var_name='Gene', value_name='Value') 
G_L_df['direction'] = G_L_df.apply(lambda x: "Losses" if x.Value < 0 else "Gains", axis = 1)
G_L_counts = G_L_df.direction.value_counts().reset_index()
print(G_L_df)
print(G_L_counts)

f, ax = plt.subplots(figsize=(6, 8))
palette = sns.color_palette("coolwarm", len(G_L_counts))
sns.barplot(x="direction", y="count", data=G_L_counts, palette=palette)
ax.set(ylabel="Counts",
       xlabel="CNV type")
plt.setp(ax.get_xticklabels(), fontweight='bold')
ax.set_title('Number of CNVs by Type', fontsize=17)
ax.grid(True, linestyle='--', alpha=0.5)
# plt.savefig("../Plots/Barplot_CNV_counts_by_type.png")

plt.show()


###              									TUMOUR TYPE ANALYSIS       								###

## dict to map the codes to the tumour type
mapping_dict = pd.Series(codes['Study Name'].values, index=codes['TSS Code']).to_dict()
# print(mapping_dict)


type_df = df.reset_index(names = 'Sample')
# print(type_df.iloc[:5,:5])

type_df['Sample'] = type_df.Sample.apply(lambda x: x.split('-')[1])
type_df['Sample'] = type_df['Sample'].map(mapping_dict)



##Median Fitlered CNV to Tumour type clustermap

	type_df_grouped = type_df.groupby('Sample').median()
	type_df_grouped.columns.name = None
	# print(type_df_grouped.iloc[:,:4])
	# print(type_df_grouped.shape)


	gene_descrp = type_df_grouped.describe().T
	T_gene_descrp = gene_descrp[['min','max']]
	# print(T_gene_descrp)

	non_sign_genes = T_gene_descrp[(T_gene_descrp['min'] > -0.49) & (T_gene_descrp['max'] < 0.99)]
	nsg_index_list = non_sign_genes.index.tolist()

	# print(index_list)


	canc_descrp = type_df_grouped.T.describe().T
	T_canc_descrp = canc_descrp[['min','max']]
	# print(T_canc_descrp)

	non_sign_canc = T_canc_descrp[(T_canc_descrp['min'] > -0.49) & (T_canc_descrp['max'] < 0.99)]
	nsc_index_list = non_sign_canc.index.tolist()
	# print(nsc_index_list)

	new_type_grp = type_df_grouped.drop(columns=nsg_index_list, axis=1)
	new_type_grp = new_type_grp.drop(nsc_index_list, axis=0)
	print(new_type_grp.T)
	print(new_type_grp.shape)

	grid = sns.clustermap(new_type_grp, cmap='coolwarm', figsize=(15, 10), dendrogram_ratio=(.06, .09))

	# # # # # Additional formatting
	grid.fig.suptitle('Clustermap of median CNV by Gene and Tumour Type (conservative approach)', size= 20, y=0.98)
	grid.fig.subplots_adjust(top=0.95, left = 0.05) 

	plt.setp(grid.ax_heatmap.get_yticklabels(), fontweight='bold', size = 10)
	plt.setp(grid.ax_heatmap.get_xticklabels(), fontweight='bold', size = 8)
	plt.setp(grid.ax_heatmap.set_ylabel('Tissue Source'), size = 15)
	plt.setp(grid.ax_heatmap.set_xlabel('Gene'), size = 15)
	plt.ylabel('CNV')
	cbar = grid.ax_cbar
	cbar.set_position([.01, .03, .01, .2])  # [left, bottom, width, height]
	plt.setp(grid.ax_heatmap.get_yticklabels(), fontweight='bold', size = 10)

	# grid.ax_heatmap.set_position([0.05, 0.05, 0.7, 0.7])
	plt.savefig("../Plots/Clustermap of median CNV by Gene and Tumour Type (conservative approach)jpg")
	# plt.show()

#Testicular canc gene identification
	testicular_gains = new_type_grp.T['Testicular Germ Cell Tumors'].reset_index()
	testicular_gains = testicular_gains.rename(columns = {'index':"Gene", "Testicular Germ Cell Tumors":"Median CNV"})
	testicular_gains = testicular_gains[testicular_gains['Median CNV'] > 1]
	# print(testicular_gains)

	testicular_genes = testicular_gains.Gene.tolist()
	# print(testicular_genes)

	testicular_genes_corr = df[testicular_genes].corr()
	print(testicular_genes_corr)

	mask = np.triu(np.ones_like(testicular_genes_corr, dtype=bool))
	f, ax = plt.subplots(figsize=(11, 9))

	heatmap = sns.heatmap(testicular_genes_corr, mask=mask, cmap='coolwarm', vmax=1, vmin=0.6, center=0.5,
	            square=True, linewidths=.2, annot = True, cbar_kws={'label':'R (Pearson\'s)',"shrink": .5})
	plt.title('Gene Subset Correlations', size = 17)
	cbar = heatmap.collections[0].colorbar
	cbar.ax.set_position([0.87, 0.15, 0.03, 0.7])
	cbar.ax.yaxis.set_label_position('left')
	plt.savefig("../Plots/Testic_Canc_Gene_Corr.png")
	plt.show()

##Median CNV to Tumour type clustermap

	# type_df_grouped = type_df.groupby('Sample').median()
	# type_df_grouped.columns.name = None

	# print(type_df_grouped.iloc[:,:5])



	# grid = sns.clustermap(type_df_grouped, cmap='coolwarm', figsize=(25, 12))

	# # # # # Set gene names as x-axis labels
	# # # # ax.set_xticklabels(df.columns, rotation=45)

	# # # # # If you want to add sample names as y-axis labels (assuming they are the index of 'sample')
	# # # # ax.set_yticklabels(df.index, rotation=0)

	# # # # # # Additional formatting
	# grid.fig.suptitle('Clustermap of median CNV by Gene and Tumour Type', size= 20)
	# plt.setp(grid.ax_heatmap.get_yticklabels(), fontweight='bold')
	# plt.ylabel('CNV')
	# plt.show()

##Average CNV to Tumour type clustermap

	# type_df_grouped = type_df.groupby('Sample').mean()
	# type_df_grouped.columns.name = None

	# print(type_df_grouped.iloc[:,:5])



	# grid = sns.clustermap(type_df_grouped, cmap='coolwarm', figsize=(25, 12))

	# # # # Set gene names as x-axis labels
	# # # ax.set_xticklabels(df.columns, rotation=45)

	# # # # If you want to add sample names as y-axis labels (assuming they are the index of 'sample')
	# # # ax.set_yticklabels(df.index, rotation=0)

	# # # # # Additional formatting
	# grid.fig.suptitle('Clustermap of average CNV by Gene and Tumour Type', size= 20)
	# plt.setp(grid.ax_heatmap.get_yticklabels(), fontweight='bold')
	# plt.ylabel('CNV')



	# plt.show()


### 										ANALYSIS OF CNV GENES BY AVERAGE CNV 									###
''' 
p_values = [wilcoxon_test(df[col])[1] for col in df.columns]

##Benjamini-Hochberg Correction
corrected_p_values = smm.multipletests(p_values, alpha=0.05, method='fdr_bh')[1]

# Create a new DataFrame with corrected p-values
results_df = pd.DataFrame({
    'P_value': p_values,
    'Corrected_P_Value': corrected_p_values
}, 
index = df.columns)


mean = df.mean()

full_CNV_mean = pd.concat([mean, results_df], axis = 1)
full_CNV_mean = full_CNV_mean.rename(columns = {0:'mean_CNV'})
full_CNV_mean.index.name = 'Genes'

# print('\n', full_CNV_mean.head(), full_CNV_mean.shape)

signif_cnv_mean = full_CNV_mean[(full_CNV_mean.Corrected_P_Value <= 0.05) & 
	((full_CNV_mean.mean_CNV > 0.07) | (full_CNV_mean.mean_CNV < -0.05))].sort_values('mean_CNV', ascending = False).reset_index(drop = False)

# print('\n', signif_cnv_mean.head(),'\n')

cnv_mean_top10 = list(signif_cnv_mean.Genes.head(10))
cnv_mean_bottom10 = list(signif_cnv_mean.Genes.tail(10))
cnv_topbottom_10 =  cnv_mean_top10 + cnv_mean_bottom10 

# print('\n', cnv_topbottom_10)

cnv_mean_topbottom_10 = df[cnv_topbottom_10].head()
# print(cnv_mean_topbottom_10)

metled_cnv_mean_topbottom_10 = cnv_mean_topbottom_10.reset_index().melt(id_vars='index', var_name='GENES', value_name='CNVs')
print(metled_cnv_mean_topbottom_10)

f, ax = plt.subplots(figsize=(12, 8))

sns.boxplot(metled_cnv_mean_topbottom_10, x="CNVs", y="GENES", width=.6, palette="vlag", linewidth=.75)
sns.stripplot(metled_cnv_mean_topbottom_10, x="CNVs", y="GENES", size=1, color=".3")

plt.title('CNV distribution of top and bottom 10 genes by average CNV', size = 18)
ax.set_xlim(-1, 1) 

for label in ax.get_yticklabels():
    label.set_fontweight('bold')

plt.show()
'''








