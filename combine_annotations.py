import pandas as pd
import sys

# Usage: eg. python combine_annotations.py utg.ls utg_to_chr.nosign.csv hapmers.csv telo_utgs.csv

# Load the tables
data = []
for arg in sys.argv[1:]:
	data.append(pd.read_csv(arg, index_col='utg'))

# Combine the tables
combined_df = pd.merge(data[0], data[1], how='left', on='utg')
combined_df = pd.merge(combined_df, data[2], how='left', on='utg')
combined_df.update(data[3])

combined_df.to_csv('combined_df.csv')
