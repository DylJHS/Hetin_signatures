import pandas as pd
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
import matplotlib.pyplot as plt
import numpy as np 
import random 
import sys
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)



data = pd.read_csv('/Users/Dyll/Downloads/GSE120795_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '	')



print(data.iloc[:5,:5])