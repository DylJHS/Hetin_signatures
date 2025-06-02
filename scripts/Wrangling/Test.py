import random
import time 
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from io import StringIO



n= 5

file_path = f'../../Data/RNA_Data/Control_Genes_RNA_Data_df{2}.csv'


# Load the mRNA control dataset. 
cntrl_RNA_df = pd.read_csv(file_path)


# list containing the same genes as the RNA control df.
cntrl_list = list(cntrl_RNA_df.columns[1:])

print(cntrl_list[33])



column_headers = pd.read_csv(file_path, index_col=0, nrows=0).columns.tolist()

print(column_headers[33])


