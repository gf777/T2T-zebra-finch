import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# Input data
data = """Release date\tCompleteness (%)
02/08/2013\t86.63
02/21/2020\t90.79
05/04/2021\t92.12
05/09/2022\t96.31
03/21/2025\t100.00"""

# Load into DataFrame
from io import StringIO
df = pd.read_csv(StringIO(data), sep='\t')
df['Release date'] = pd.to_datetime(df['Release date'])

# Plot
plt.figure(figsize=(3, 5))
plt.plot(df['Release date'], df['Completeness (%)'],
         marker='o', linestyle='--')
plt.xticks(df['Release date'], rotation=45, ha='right')
plt.xlabel("Release date")
plt.ylabel("Completeness (%)")
plt.title("Assembly completeness")
plt.grid(True)
plt.tight_layout()

# Save to file
plt.savefig("completeness_plot.png", dpi=300)
plt.close()