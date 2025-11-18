import argparse
import pandas as pd
import pysam

# Constants: your RNA‐motifs
MOTIFS = [
    "AAUAAA",
    "AUUAAA",
    "UAUAAA",
    "AGUAAA",
    "AAUGAA",
    "CAUAAA",
    "AAUACA",
    "AAUAUA",
    "GAUAAA",
    "ACUAAA",
    "AAUAGA",
    "AAGAAA",
]


def write_to_bed(df, out_file):
    """Write the original PAS‐atlas columns plus new motif‐presence columns."""
    base_cols = [
        "chr",
        "start",
        "end",
        "cs_id",
        "score",
        "strand",
        "support_percent",
        "support_count",
        "score2",
        "cluster_annotation",
        "orig_motifs",
    ]
    df.to_csv(
        out_file,
        sep="\t",
        columns=base_cols + MOTIFS,
        index=False,
        header=False,
    )


def process_motifs(pas_df, fasta_path):
    """
    Add columns for motif presence to the DataFrame.
    """
    # Name the input columns so we can refer by name:
    pas_df.columns = [
        "chr",
        "start",
        "end",
        "cs_id",
        "score",
        "strand",
        "support_percent",
        "support_count",
        "score2",
        "cluster_annotation",
        "orig_motifs",
    ]

    fasta_f = pysam.FastaFile(fasta_path)
    motif_results = {motif: [] for motif in MOTIFS}

    for _, row in pas_df.iterrows():
        chrom = row["chr"]
        strand = row["strand"]

        # ——— NEW: extract the exact PAS position from the cs_id
        # cs_id looks like "1:16450:-"
        _, pos_str, _ = row["cs_id"].split(":")
        cleavage_site = int(pos_str)

        # Define the scanning window exactly as before
        if strand == "-":
            sequence_start = cleavage_site + 10
            sequence_end = cleavage_site + 35
        else:
            sequence_start = cleavage_site - 35
            sequence_end = cleavage_site - 10

        dna_sequence = fasta_f.fetch(chrom, sequence_start, sequence_end + 1)
        rna_sequence = change_dna_to_rna(dna_sequence, strand)

        for motif in MOTIFS:
            motif_results[motif].append(1 if motif in rna_sequence else 0)

    # Attach one column per motif
    for motif, values in motif_results.items():
        pas_df[motif] = values

    return pas_df


def change_dna_to_rna(sequence, strand):
    """
    Convert a DNA sequence to RNA based on the strand direction.
    """
    if strand == "-":
        sequence = sequence[::-1]
        conversion = {"A": "U", "T": "A", "G": "C", "C": "G"}
    else:
        conversion = {"A": "A", "T": "U", "G": "G", "C": "C"}
    return "".join(conversion.get(base, base) for base in sequence)


def get_args():
    parser = argparse.ArgumentParser(
        description="Analyze motifs in a PAS‐atlas BED file."
    )
    parser.add_argument(
        "-i", "--input-bed", required=True, help="Input PAS atlas BED file"
    )
    parser.add_argument(
        "-f", "--fasta-genome", required=True, help="Path to FASTA file"
    )
    parser.add_argument(
        "-o",
        "--output-bed",
        required=True,
        help="Output BED file with motif‐presence flags",
    )
    return parser.parse_args()


def run_process():
    args = get_args()
    pas_df = pd.read_csv(args.input_bed, delimiter="\t", header=None)

    print("Processing motifs...")
    final_df = process_motifs(pas_df, args.fasta_genome)

    print("Writing output...")
    write_to_bed(final_df, args.output_bed)
    print("Done!")


if __name__ == "__main__":
    run_process()
