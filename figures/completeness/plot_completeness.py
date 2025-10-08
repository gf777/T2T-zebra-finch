#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from io import StringIO

# Ensure SVG text remains as text (not paths)
plt.rcParams.update({
    "svg.fonttype": "none",   # keep text as <text>
    "text.usetex": False,     # avoid TeX -> paths
    "font.family": "sans-serif"
})

# Input data
data = """Release date\tCompleteness (%)
02/08/2013\t86.63
02/21/2020\t90.79
05/04/2021\t92.12
05/09/2022\t96.31
03/21/2025\t100.00"""

# Load into DataFrame
df = pd.read_csv(StringIO(data), sep='\t')
df['Release date'] = pd.to_datetime(df['Release date'], format='%m/%d/%Y')

# Plot
plt.figure(figsize=(3, 5))
plt.plot(df['Release date'], df['Completeness (%)'],
         marker='o', linestyle='--')

# Format axes and labels
plt.xticks(rotation=45, ha='right')
plt.xlabel("Release date")
plt.ylabel("Completeness (%)")
plt.title("Assembly completeness")
plt.grid(True)

# Date formatting for X axis
date_fmt = DateFormatter('%m/%d/%Y')
plt.gca().xaxis.set_major_formatter(date_fmt)

plt.tight_layout()

# Save as SVG with text preserved
plt.savefig("completeness_plot.svg", format="svg")
plt.close()

