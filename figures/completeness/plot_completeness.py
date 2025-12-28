#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from io import StringIO

# Plot parameters
plt.rcParams.update({
    "svg.fonttype": "none",   # keep text as <text>
    "text.usetex": False,     # avoid TeX -> paths
    "font.family": "Arial",
    "font.size": 9,
    "axes.linewidth": 0.5,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "lines.linewidth": 1.5,
    "lines.markersize": 6,
})

# Input data
data = """Release date\tCompleteness (%)
02/08/2013\t86.63
02/21/2020\t90.79
05/04/2021\t92.12
05/09/2022\t96.31
03/21/2025\t100.00"""

df = pd.read_csv(StringIO(data), sep='\t')
df['Release date'] = pd.to_datetime(df['Release date'], format='%m/%d/%Y')

# Plot
plt.figure(figsize=(3, 3.5))
plt.plot(df['Release date'], df['Completeness (%)'],
        marker='o', linestyle='--')

# Format axes and labels
plt.xticks(rotation=45, ha='right')
plt.xlabel("Release date")
plt.ylabel("Sequence completeness (%)")
plt.title("Assembly completeness")

# X-axis ticks
plt.gca().set_xticks(pd.to_datetime(['2013-02-08', '2020-02-21', '2021-05-04', '2022-05-09', '2025-03-21']))
plt.gca().set_xticklabels(['2013', '2020', '2021', '2022', '2025'])

# Y-axis ticks and grid lines
ax = plt.gca()
ax.set_yticks([86, 88, 90, 92, 94, 96, 98, 100])
ax.set_yticklabels(['', '88', '', '92', '', '96', '', '100'])
ax.set_axisbelow(True)
ax.yaxis.grid(True, linewidth=0.5, color='#D3D3D3')

# Set vertical grid lines every 2 years (without tick marks)
ax.set_xticks(pd.to_datetime(['2014-01-01', '2016-01-01', '2018-01-01', '2020-01-01', '2022-01-01', '2024-01-01']), minor=True)
ax.tick_params(axis='x', which='minor', length=0)
ax.xaxis.grid(True, which='minor', linewidth=0.5, color='#D3D3D3')

# Date formatting for X axis
date_fmt = DateFormatter('%Y')
plt.gca().xaxis.set_major_formatter(date_fmt)

plt.tight_layout()

# Save as SVG with text preserved
plt.savefig("completeness_plot.svg", format="svg")
plt.savefig("completeness_plot.png", format="png", dpi=600)
plt.savefig("completeness_plot.pdf", format="pdf")
plt.close()