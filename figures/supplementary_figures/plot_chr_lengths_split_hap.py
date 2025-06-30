import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def load_fai_with_haplotypes(fai_path):
    df = pd.read_csv(fai_path, sep='\t', header=None, usecols=[0, 1], names=['chr', 'length'])
    df['hap'] = df['chr'].str.extract(r'_(mat|pat)$')[0]
    df['chr_base'] = df['chr'].str.replace(r'_(mat|pat)$', '', regex=True)
    return df[['chr_base', 'hap', 'length']]

def main(fai_path, output):
    df = load_fai_with_haplotypes(fai_path)

    # Pivot to get one row per chromosome with mat/pat columns
    pivot = df.pivot_table(index='chr_base', columns='hap', values='length', aggfunc='sum').fillna(0)

    # Ensure both columns exist
    for hap in ['mat', 'pat']:
        if hap not in pivot.columns:
            pivot[hap] = 0

    # For sorting: use diploid-equivalent length
    def compute_diploid_length(row):
        if row['mat'] > 0 and row['pat'] > 0:
            return row['mat'] + row['pat']
        else:
            return 2 * max(row['mat'], row['pat'])

    pivot['sort_total'] = pivot.apply(compute_diploid_length, axis=1)
    pivot = pivot.sort_values('sort_total', ascending=False).drop(columns='sort_total')

    # Extract sorted labels and values
    chroms = pivot.index.tolist()
    x = np.arange(len(chroms))
    bar_width = 0.4
    mat_vals = pivot['mat'].values
    pat_vals = pivot['pat'].values

    # Plot
    fig, ax = plt.subplots(figsize=(16, 5))
    ax.bar(x - bar_width/2, mat_vals, width=bar_width, label='Maternal', color='skyblue')
    ax.bar(x + bar_width/2, pat_vals, width=bar_width, label='Paternal', color='steelblue')

    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=90)
    ax.set_xlim(x[0] - 0.6, x[-1] + 0.6)
    ax.set_ylabel('Chromosome Length (Gbp)')
    ax.set_title('Chromosome Lengths in Gbp (Maternal and Paternal)')
    ax.legend()

    # Format y-axis to Gbp
    yticks = ax.get_yticks()
    ax.set_yticklabels([f'{y/1e9:.2f}' for y in yticks])

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    print(f"[INFO] Saved plot to {output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot haplotype-resolved chromosome lengths from single .fai")
    parser.add_argument("fai", help=".fai file with _mat and _pat chromosomes")
    parser.add_argument("-o", "--output", default="chr_lengths_split.png", help="Output PNG file")
    args = parser.parse_args()
    main(args.fai, args.output)
