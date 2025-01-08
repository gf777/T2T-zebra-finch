import pandas as pd
import sys

# Usage: eg. python combine_annotations.py assembly.colors.csv utg.ls utg_to_chr.nosign.csv hapmers.csv telo_utgs.csv

# Load the tables
data = []
for arg in sys.argv[2:]:
	data.append(pd.read_csv(arg, index_col='utg'))

hic_colors = pd.read_csv(sys.argv[1], sep='\t')
hic_colors.drop(columns=['mat', 'pat', 'mat:pat'], axis=1, inplace=True)
hic_colors = hic_colors.rename(columns={'node': 'utg', 'color': 'hic_color'})
hic_colors.set_index('utg', inplace=True)
print(hic_colors)
# Combine the tables
combined_df = pd.merge(data[0], data[1], how='left', on='utg')
combined_df = pd.merge(combined_df, data[2], how='left', on='utg')
combined_df = pd.merge(combined_df, hic_colors, how='left', on='utg')
combined_df.update(data[3])

combined_df.to_csv('combined_df.csv')
