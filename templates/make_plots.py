#!/usr/bin/env python3

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns

# This folder contains a set of many CSV files with the following file name syntax
# {specimen}.SSC.csv.gz
# {specimen}.by_read_position.csv.gz
# {specimen}.by_chr.csv.gz
# {specimen}.snps_by_base.csv.gz
# {specimen}.adducts_by_base.csv.gz

# Generic function to read and merge the CSV files with the same extension
def read_files(suffix, melt=False):
    
    # Populate a list
    output = []

    # Iterate over every file
    for fp in os.listdir('.'):

        # Check for the file extension
        if not fp.endswith(suffix):
            continue

        # Parse the specimen name
        specimen = fp.replace(suffix, '')

        # Read the file
        df = pd.read_csv(fp)

        # If the melt flag is set
        if melt:

            # Then melt into long format
            df = df.melt(id_vars=df.columns.values[0])

        # Add the specimen name
        df = df.assign(specimen=specimen)

        # Add to the list
        output.append(df)

    # If there were no files with the extension
    if len(output) == 0:
        return None

    # Otherwise, join the DataFrames and return
    return pd.concat(output).reset_index(drop=True)


def plot_distribution(
    df,
    col_name,
    y=None,
    pdf=None,
    kind='kde',
    hue='specimen',
    title=None,
    xlabel=None,
    ylabel=None,
):

    # Try to make the plot
    try:
        sns.displot(data=df, x=col_name, y=y, hue=hue, kind=kind)
        annotate_and_save(xlabel=xlabel, ylabel=ylabel, pdf=pdf, title=title)

    # If there is an index error, then the underlying data may not have
    # had enough variation to plot the distribution
    except IndexError:
        print(f"Could not make plot: x='{col_name}', y='{y}', hue='{hue}', title='{title}'")


def plot_lines(
    df,
    x,
    y,
    hue='specimen',
    pdf=None,
    title=None,
    xlabel=None,
    ylabel=None,
    alpha=0.5,
):
    # Make the plot
    sns.lineplot(data=df, x=x, y=y, hue=hue, alpha=alpha)
    plt.legend(bbox_to_anchor=[1.1, 0.9])
    annotate_and_save(xlabel=xlabel, ylabel=ylabel, pdf=pdf, title=title)


def plot_bars(
    v,
    pdf=None,
    title=None,
    xlabel=None,
    ylabel=None,
):
    # Make the plot - sorted by values
    v.sort_values(ascending=False).plot(kind='barh')
    annotate_and_save(xlabel=xlabel, ylabel=ylabel, pdf=pdf, title=title)


def annotate_and_save(xlabel=None, ylabel=None, pdf=None, title=None, logx=False, logy=False):

    # Log-transform x-axis
    if logx:
        plt.xscale('log')

    # Log-transform y-axis
    if logy:
        plt.yscale('log')

    # Set the title
    if title is not None:
        plt.title(title)

    # Set the xlabel
    if xlabel is not None:
        plt.xlabel(xlabel)

    # Set the ylabel
    if ylabel is not None:
        plt.ylabel(ylabel)

    # Save to PDF
    if pdf is not None:
        pdf.savefig(bbox_inches="tight")

    # Close the plot
    plt.close()


# Specimen-level metrics
def plot_specimen_summary(df, pdf):

    # If there is no data
    if df is None:

        # Do not make any plot
        return

    # Plot the SNP rate per specimen
    plot_bars(
        df.snp_rate,
        pdf=pdf,
        title="SNP Rate",
        xlabel="# of SNPs / # of sequenced bases"
    )

    # Plot the adduct rate per specimen
    plot_bars(
        df.adduct_rate,
        pdf=pdf,
        title="Adduct Rate",
        xlabel="# of adducts / # of sequenced bases"
    )

# SSC summary metrics
# {specimen}.SSC.csv.gz
def plot_ssc_summary(pdf):

    # Read the table
    df = read_files('.SSC.csv.gz')

    # If there is no data
    if df is None:

        # Do not make any plot
        return

    # Calculate the number of reads per DSC as the sum of
    # the number of reads on both strands
    df = df.assign(
        nreads_dsc = df.nreads_pos + df.nreads_neg
    )

    # Plot the distribution of read lengths
    plot_distribution(
        df,
        "rlen_fwd",
        pdf=pdf,
        title="Single Read Length Distribution",
        xlabel="Length of Sequence Reads"
    )

    # Plot the number of reads per DSC
    plot_distribution(
        df,
        "nreads_dsc",
        pdf=pdf,
        title="DSC Sequencing Depth",
        xlabel="Number of reads per DSC"
    )
    
# Barcode saturation - number of reads per barcode
# {specimen}.corrected_barcodes.csv
def plot_barcode_counts(pdf):

    # Read the table
    df = read_files('.corrected_barcodes.csv')

    # If there is no data
    if df is None:

        # Do not make any plot
        return

    # Plot the number of reads per barcode, per specimen
    sns.lineplot(
        data=df,
        x='nreads',
        y='nbarcodes',
        hue='specimen'
    )
    annotate_and_save(
        xlabel="Number of reads per barcode",
        ylabel="Number of barcodes",
        title="Barcode Saturation",
        pdf=pdf,
        logx=True,
        logy=True
    )

# Summary of mutations by read position
# {specimen}.by_read_position.csv.gz
def plot_read_position(pdf):
    
    # Read the table
    df = read_files('.by_read_position.csv.gz')

    # If there is no data
    if df is None:

        # Do not make any plot
        return

    # Calculate the proportion of SNPs and adducts
    df = df.assign(
        snp_prop=df.snps / df.nreads,
        adduct_prop=df.adducts / df.nreads,
    )

    plot_lines(
        df,
        'pos',
        'nreads',
        pdf=pdf,
        xlabel="Position in Read",
        ylabel="Number of Reads",
        title="Sequencing Depth"
    )

    plot_lines(
        df,
        'pos',
        'snp_prop',
        pdf=pdf,
        xlabel="Position in Read",
        ylabel="SNP Proportion",
        title="SNP Rate by Read Position"
    )

    plot_lines(
        df,
        'pos',
        'adduct_prop',
        pdf=pdf,
        xlabel="Position in Read",
        ylabel="Adduct Proportion",
        title="Adduct Rate by Read Position"
    )


def plot_heatmap(
    suffix=None,
    norm=None,
    csv_fp=None,
    pdf=None,
    title=None
):

    print("Plotting files with the suffix: " + suffix)
    
    # Read the tables
    df = read_files(suffix, melt=True)

    # If there is no data
    if df is None:

        # Do not make any plot
        return

    # Divide the value by the normalization factor by specimen
    df = df.assign(
        prop=df.value / df.specimen.apply(norm.get)
    )

    # Format the base change as a string
    df = df.assign(
        base_change=df.apply(
            lambda r: f"{r['base']} -> {r['variable']}",
            axis=1
        )
    ).pivot(
        columns="specimen",
        index='base_change',
        values="prop"
    )

    # Save the table to CSV
    print(f"Saving to {csv_fp}")
    df.to_csv(csv_fp)

    # Remove any rows which lack observations
    df = df.loc[df.sum(axis=1) > 0]

    # If there is no data to plot
    if df.shape[0] == 0:

        # Don't make the plot
        print(f"No data found for plotting")
        return

    # Otherwise, if there is data to plot

    # Make a plot with the rate of counts per change
    sns.heatmap(df, cmap="Blues")
    plt.yticks(rotation=0)
    plt.ylabel("Base Change")
    plt.xlabel("")
    if title is not None:
        plt.title(title)
    pdf.savefig(bbox_inches="tight")
    plt.close()


def read_specimen_summary(suffix='.SSC.csv.gz'):
    """Parse the SSC summary tables to get summary metrics."""

    # Get the list of files with the suffix
    fp_list = [
        fp
        for fp in os.listdir('.')
        if fp.endswith(suffix)
    ]

    # If there are no files with the suffix
    if len(fp_list) == 0:

        # Return None
        return None

    # Otherwise

    # Join together data from all of those files
    df = pd.DataFrame(
        [
            parse_specimen_summary(fp, fp[:-(len(suffix))])
            for fp in fp_list
        ]
    ).set_index('specimen')

    # Calculate the rate of SNPs and adducts
    df = df.assign(
        adduct_rate=df.adducts / df.bases,
        snp_rate=df.snps / df.bases,
    )

    return df

def parse_specimen_summary(fp, specimen):
    """Get specimen summary metrics from a single SSC.csv.gz file."""
    df = pd.read_csv(fp)

    return dict(
        specimen=specimen,
        dsc=df.shape[0],
        bases=df.merged_len.sum(),
        adducts=df.n_adducts.sum(),
        snps=df.n_mutations.sum()
    )

# Get the total number of sequenced bases per specimen
specimen_summary = read_specimen_summary()

# If there is data
if specimen_summary is not None:

    # Save the specimen summary
    specimen_summary.to_csv("summary.csv")

# Open up a PDF for plotting
with PdfPages("report.pdf") as pdf:

    # Plot the barcode saturation
    plot_barcode_counts(pdf)

    # Plot a summary of SSC metrics
    plot_ssc_summary(pdf)

    # Plot a summary of each specimen
    plot_specimen_summary(specimen_summary, pdf)

    # Only plot heatmaps if the specimen summary is available
    if specimen_summary is not None:
    
        # Summary of mutations by base -> base
        # {specimen}.snps_by_base.csv.gz
        plot_heatmap(
            suffix=".snps_by_base.csv.gz",
            csv_fp="snps_by_base.csv",
            norm=specimen_summary.bases,
            pdf=pdf,
            title="SNPs by Base"
        )
        # Summary of adducts by base -> base
        # {specimen}.adducts_by_base.csv.gz
        plot_heatmap(
            suffix=".adducts_by_base.csv.gz",
            csv_fp="adducts_by_base.csv",
            norm=specimen_summary.bases,
            pdf=pdf,
            title="Adducts by Base"
        )

    plot_read_position(pdf)
