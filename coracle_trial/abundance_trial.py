import numpy
import pandas as pd
import numpy as np
import csv
dir2 = '/home/tally/Documents/Thesis_Part_2/Newt_TTX_Project/coracle_files/ttx_2/'
md = pd.read_csv(dir2 + 'ttx_metadata.tsv', sep = '\t', header = 0)
ttx_2_df = pd.read_csv(dir2 + 'feature-table.tsv', sep = '\t', header = 0)
print(ttx_2_df)
#ttx_2_df.rename(columns={'#OTU ID':'SampleID'}, inplace = True)
print(ttx_2_df)
x_file = md['Whole TTX'].to_frame()
print(x_file)


merged_df = pd.merge(ttx_2_df, md, on='SampleID', how = 'outer')
output_file = (dir2 + 'combined_file.csv')
merged_df.to_csv(output_file, sep = '\t', index=False)
new_df = pd.read_csv(dir2 + 'combined_file.csv', sep= '\t', index_col=0)
#print(new_df.head())
#actual coracle start
y_file = new_df['Whole TTX'].to_frame()
print(y_file)
output_file2 = (dir2 + 'y_file.csv')
y_file.to_csv(output_file2, sep = '\t', index = True)