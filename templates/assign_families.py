#!/usr/bin/env python3

import pandas as pd

# Read in the input data
df = pd.read_csv("read_positions.csv.gz")

def assign_families(bc_seq, bc_df):
    """Assign the family label for the reads with the same barcode."""

    # Count up the number of alignments for each read on each strand (R1/R2)
    vc = bc_df.groupby(['id', 'strand']).apply(len)

    # Only keep reads which align uniquely in a single pair
    reads_to_keep = set([
        read_id
        for read_id in bc_df['id'].unique()
        if vc[read_id].get('R1') == 1 and vc[read_id].get('R2') == 1
    ])
    if len(reads_to_keep) == 0:
        return

    bc_df = bc_df.loc[bc_df['id'].isin(reads_to_keep)]

    # Count up the number of alignments for each read in each direction (fwd/rev)
    vc = bc_df.groupby(['id', 'direction']).apply(len)

    # Only keep reads which align uniquely in a single pair
    reads_to_keep = set([
        read_id
        for read_id in bc_df['id'].unique()
        if vc[read_id].get('fwd') == 1 and vc[read_id].get('rev') == 1
    ])
    if len(reads_to_keep) == 0:
        return

    bc_df = bc_df.loc[bc_df['id'].isin(reads_to_keep)]

    # For each read, get the start positions of the fwd and reverse reads
    read_df = bc_df.pivot(
        index=["id", "chr"],
        columns="direction",
        values="pos"
    ).dropna(
    ).applymap(
        int
    ).sort_values(
        by="fwd"
    ).reset_index()

    # Assign a unique family ID to each set of reads which align
    # to the same fwd and rev positions
    read_family = {
        r['id']: f"{bc_seq}-{r['chr']}-{r['fwd']}-{r['rev']}"
        for _, r in read_df.iterrows()
    }

    read_df = read_df.assign(
        family=read_df['id'].apply(read_family.get)
    ).reindex(
        columns=['id', 'family']
    )

    return read_df

# Make families for barcode independently
df = pd.concat(
    [
        d
        for d in [
            assign_families(bc_seq, bc_df)
            for bc_seq, bc_df in df.groupby("barcode")
        ]
        if d is not None
    ]
)

# Write out
df.to_csv("families.csv.gz", index=None)
