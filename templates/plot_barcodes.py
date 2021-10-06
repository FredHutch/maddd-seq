#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read in the raw barcodes
raw_df = pd.read_csv("barcode_counts.csv.gz")

def plot_hist(v, l, log_scale=(True, True), bins=100):
    """Plot a histogram following a similar pattern"""
    sns.histplot(
        x=v,
        log_scale=log_scale,
        bins=bins
    )
    plt.xlabel(f'Number of reads per {l} barcode')
    plt.ylabel('Number of barcodes')
    plt.title(f'Barcodes - {l}')
    plt.savefig(f'barcodes.{l}.pdf', bbox_inches='tight')
    plt.close()

# Make a plot showing the frequency distribution of RAW barcodes
plot_hist(
    raw_df['count'],
    l='raw'
)

# Read in the corrected barcodes
corr_df = pd.read_csv("barcode_corrections.csv.gz")

# Subset the original table to those barcodes which were filtered out
filtered_barcodes = set(raw_df['barcode'].tolist()) - set(corr_df['barcode'].tolist())
filt_df = raw_df.loc[raw_df['barcode'].isin(filtered_barcodes)]

# Make a plot showing the frequency distribution of barcodes which were filtered out
plot_hist(filt_df['count'], 'removed')

# Make a plot showing the frequency distribution of CORRECTED barcodes
plot_hist(corr_df.groupby('corrected')['count'].sum(), 'corrected')
