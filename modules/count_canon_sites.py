import argparse
import pandas as pd
import pybedtools


def modify_coordinates(cleavage_df):
    """
    Modify the cleavage site coordinates based on the strand:
    [-25, +25] for '+' strand and [-25, +25] for '-' strand.

    Args:
        cleavage_df (pandas.DataFrame): The cleavage sites DataFrame.

    Returns:
        pandas.DataFrame: A modified copy of the cleavage_df with adjusted coordinates.
    """
    modified = cleavage_df.copy()

    # Adjust coordinates based on strand
    plus_strand = modified[5] == "+"
    minus_strand = modified[5] == "-"

    # For '+' strand
    modified.loc[plus_strand, 1] = modified.loc[plus_strand, 1].astype(int) - 25
    modified.loc[plus_strand, 2] = modified.loc[plus_strand, 2].astype(int) + 25

    # For '-' strand
    modified.loc[minus_strand, 1] = modified.loc[minus_strand, 1].astype(int) - 25
    modified.loc[minus_strand, 2] = modified.loc[minus_strand, 2].astype(int) + 25

    return modified


def intersect_and_count(cleavage_bed, canonical_bed):
    """
    Intersect the modified cleavage BED with the canonical BED to count overlaps.

    Args:
        cleavage_bed (pybedtools.BedTool): Cleavage sites BED file as a BedTool object.
        canonical_bed (pybedtools.BedTool): Canonical sites BED file as a BedTool object.

    Returns:
        pandas.DataFrame: A DataFrame with the overlap counts for each cleavage site.
    """
    # Perform intersection with strandedness
    intersect_result = cleavage_bed.intersect(canonical_bed, s=True, c=True)

    # Convert intersection result back to DataFrame
    intersect_df = pd.read_csv(
        intersect_result.fn, 
        sep="\t", 
        header=None, 
        names=list(range(cleavage_bed.field_count())) + ["overlap_count"]  # Fixed field_count()
    )

    return intersect_df



def main(input_bed, canonical_bed, output_bed):
    # Read BED files into pandas DataFrames
    cleavage_df = pd.read_csv(input_bed, sep="\t", header=None)
    canonical_df = pd.read_csv(canonical_bed, sep="\t", header=None)

    # Modify coordinates based on strand
    modified_cleavage_df = modify_coordinates(cleavage_df)

    # Convert DataFrames to BedTool objects
    cleavage_bed = pybedtools.BedTool.from_dataframe(modified_cleavage_df)
    canonical_bed = pybedtools.BedTool.from_dataframe(canonical_df)

    # Intersect and count overlaps
    intersect_df = intersect_and_count(cleavage_bed, canonical_bed)

    # Merge the overlap counts with the original cleavage DataFrame on the ID column (assumed col 3)
    merged_df = cleavage_df.copy()
    merged_df = pd.merge(
        merged_df,
        intersect_df[[3, "overlap_count"]],
        left_on=3,
        right_on=3,
        how="left",
    )

    # Fill any missing overlap counts with 0
    merged_df["overlap_count"] = merged_df["overlap_count"].fillna(0).astype(int)

    # Write the final result to a BED file
    merged_df.to_csv(output_bed, sep="\t", header=False, index=False)
    print(f"Output written to {output_bed}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Intersect cleavage sites with canonical sites.")
    parser.add_argument("-i", "--input-bed", required=True, help="Input PAS BED file")
    parser.add_argument("-c", "--canonical-bed", required=True, help="BED file with canonical sites")
    parser.add_argument("-o", "--output-bed", required=True, help="Output BED file")

    args = parser.parse_args()
    main(args.input_bed, args.canonical_bed, args.output_bed)
