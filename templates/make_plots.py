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

    # Join and return
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
    # Make the plot
    sns.displot(data=df, x=col_name, y=y, hue=hue, kind=kind)
    annotate_and_save(xlabel=xlabel, ylabel=ylabel, pdf=pdf, title=title)


def plot_lines(
    df,
    x,
    y,
    hue='specimen',
    pdf=None,
    title=None,
    xlabel=None,
    ylabel=None,
):
    # Make the plot
    sns.lineplot(data=df, x=x, y=y, hue=hue)
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

    # Make the plot - sorted by specimen names
    v.sort_index().plot(kind='barh')
    annotate_and_save(xlabel=xlabel, ylabel=ylabel, pdf=pdf, title=title)


def annotate_and_save(xlabel=None, ylabel=None, pdf=None, title=None):

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

# SSC summary metrics
# {specimen}.SSC.csv.gz
def plot_ssc_summary():

    # Read the table
    df = read_files('.SSC.csv.gz')

    # Open the output
    with PdfPages("DSC.summary.pdf") as pdf:

        # Plot the distribution of read lengths
        plot_distribution(
            df,
            "rlen_fwd",
            pdf=pdf,
            title="Single Read Length Distribution",
            xlabel="Length of Sequence Reads"
        )

        # Plot the distribution of merged sequences
        plot_distribution(
            df,
            "merged_len",
            pdf=pdf,
            title="Merged Read Length Distribution",
            xlabel="Length of Merged Reads"
        )

        # Plot the distribution of sequencing depth per DSC
        plot_distribution(
            df,
            "nreads_pos",
            pdf=pdf,
            title="Sequencing Depth per DSC",
            xlabel="Number of reads per strand"
        )

        # Plot the number of DSCs per specimen
        plot_bars(
            df.specimen.value_counts(),
            pdf=pdf,
            title="Unique Sequencing Depth",
            xlabel="Number of DSCs per Specimen"
        )

        # Compare the sequencing depth to the number of adducts
        plot_distribution(
            df,
            "nreads_neg",
            y="n_adducts",
            pdf=pdf,
            title="Adducts by Sequencing Depth",
            xlabel="# of Reads per DSC",
            ylabel="# of Adducts per DSC"
        )

        # Compare the sequencing depth to the number of SNPs
        plot_distribution(
            df,
            "nreads_neg",
            y="n_mutations",
            pdf=pdf,
            title="SNPs by Sequencing Depth",
            xlabel="# of Reads per DSC",
            ylabel="# of SNPs per DSC"
        )


# Summary of mutations by read position
# {specimen}.by_read_position.csv.gz
def plot_read_position():
    # Read the table
    df = read_files('.by_read_position.csv.gz')

    # Calculate the proportion of SNPs and adducts
    df = df.assign(
        snp_prop=df.snps / df.nreads,
        adduct_prop=df.adducts / df.nreads,
    )
    
    # Open the output
    with PdfPages("read_position.pdf") as pdf:

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


def plot_heatmap(suffix=None, pdf_fp=None):
    
    # Read the tables
    df = read_files(suffix, melt=True)

    # Format the base change as a string
    df = df.assign(
        base_change=df.apply(
            lambda r: f"{r['base']} -> {r['variable']}",
            axis=1
        )
    ).pivot(
        columns="specimen",
        index='base_change',
        values="value"
    )
    df = df.loc[df.sum(axis=1) > 0]

    # Open the output
    with PdfPages(pdf_fp) as pdf:

        sns.heatmap(df, cmap="Blues")
        plt.yticks(rotation=0)
        plt.ylabel("Base Change")
        plt.xlabel("")
        pdf.savefig(bbox_inches="tight")

# Summary of mutations by base -> base
# {specimen}.snps_by_base.csv.gz
plot_heatmap(suffix=".snps_by_base.csv.gz", pdf_fp="snps_by_base.pdf")
# Summary of adducts by base -> base
# {specimen}.adducts_by_base.csv.gz
plot_heatmap(suffix=".adducts_by_base.csv.gz", pdf_fp="adducts_by_base.pdf")

plot_ssc_summary()
plot_read_position()
