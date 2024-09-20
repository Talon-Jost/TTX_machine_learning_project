import pandas as pd
import numpy as np
import csv
#_______________________Initial data setup____________________________________________________-
#metadata
dir = '/home/tally/Documents/Thesis_Part_2/Newt_TTX_Project/coracle_files/ttx_coracle/'
md_df = pd.read_csv(dir + 'ttx_metadata.tsv', sep = '\t', header=0)

#matrix import
df_asv = pd.read_csv(dir + 'asv_tabulated.csv', index_col=0)
print(df_asv.head(10))
#transposed_asv = np.transpose(ndf_asv)
#merging matrix and metadata based on column

merged_df = pd.concat([md_df, df_asv])
output_file = (dir + 'combined_metadata_frequency.csv')
merged_df.to_csv(output_file, sep = '\t', index=False)
new_df = pd.read_csv(dir + 'combined_metadata_frequency.csv', sep= '\t', index_col=0)
#print(new_df.head())
#actual coracle start
true_y = new_df['Whole TTX'].to_frame()
print(true_y)
output_file2 = (dir + 'y_ture.csv')
true_y.to_csv(output_file2, sep = '\t', index = True)
