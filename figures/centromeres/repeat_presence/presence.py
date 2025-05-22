import matplotlib.pyplot as plt
import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Scatter plot of "Tgut716" and "Tgut191A" repeats from two CSV files.')
parser.add_argument('mat_csv', type=str, help='Path to the CSV file for "Tgut716" repeats')
parser.add_argument('pat_csv', type=str, help='Path to the CSV file for "Tgut191A" repeats')

# Parse the arguments
args = parser.parse_args()

# Read the "Tgut716" and "Tgut191A" CSV files
mat_repeats = pd.read_csv(args.mat_csv, header=None, names=['Chromosome', 'Value'])
pat_repeats = pd.read_csv(args.pat_csv, header=None, names=['Chromosome', 'Value'])

# Sort both "Tgut716" and "Tgut191A" by chromosome name to match corresponding repeats
mat_repeats = mat_repeats.sort_values(by='Chromosome')
pat_repeats = pat_repeats.sort_values(by='Chromosome')

# Create the plot with a narrow figure size
fig, ax = plt.subplots(figsize=(2, 8))  # Adjusted size (width, height)

# Scatter plot for 'Tgut716' repeats
ax.scatter([1] * len(mat_repeats), mat_repeats['Value'], label='Tgut716 repeats', color='blue')

# Scatter plot for 'Tgut191A' repeats
ax.scatter([2] * len(pat_repeats), pat_repeats['Value'], label='Tgut191A repeats', color='red')

# Add a dashed line at 10 kbp
ax.axhline(y=10000, color='black', linestyle='--', label='10 kbp')

# Customize plot
ax.set_xticks([1, 2])
ax.set_xticklabels(['Tgut716 repeats', 'Tgut191A repeats'])
ax.set_yscale('log')
ax.set_ylabel('Repeat Values (bp)')

# Add a grid and a legend, and adjust layout to make sure everything fits
ax.grid(True)
ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)  # Position the legend outside the plot

# Adjust the layout to ensure labels and legend fit
plt.tight_layout()

# Save the plot as a PNG file
plt.savefig('scatter_repeats.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()
