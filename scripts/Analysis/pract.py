import pandas as pd



pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

RNA = pd.read_csv('/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/RNA_Data/TCGA_Norm/Ex_Counts/tcga_gene_expected_count.csv')

print(RNA.shape)

print("Done!")