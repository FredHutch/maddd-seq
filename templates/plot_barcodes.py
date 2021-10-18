#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read in the raw barcodes
raw_df = pd.read_csv("barcode_counts.csv.gz")
print(f"Read in {raw_df.shape[0]:,} lines from barcode_counts.csv.gz")

# Read in the corrected barcodes
corr_df = pd.read_csv("barcode_corrections.csv.gz")
print(f"Read in {corr_df.shape[0]:,} lines from barcode_corrections.csv.gz")

# Subset the original table to those barcodes which were filtered out
filtered_barcodes = set(raw_df['barcode'].tolist()) - set(corr_df['barcode'].tolist())
print(f"Number of barcodes filtered out: {len(filtered_barcodes):,}")
filt_df = raw_df.loc[raw_df['barcode'].isin(filtered_barcodes)]

plot_df = pd.concat(
    [
        filt_df['count'].value_counts().reset_index().assign(label='Filtered Out'),
        corr_df.groupby('corrected')['count'].sum().value_counts().reset_index().assign(label='After Correction'),
        raw_df['count'].value_counts().reset_index().assign(label='Raw Barcodes'),
        
    ]
).rename(
    columns=dict(index="nreads", count="nbarcodes")
).reset_index(drop=True)

g = sns.lineplot(
    data=plot_df,
    x="nreads",
    y='nbarcodes',
    hue='label'
)
g.legend_.set_title(None)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of reads per barcode")
plt.ylabel("Number of barcode per bin")
plt.title("Barcode Frequency Histogram")
plt.savefig("${specimen}.barcodes.pdf")
print("DONE")