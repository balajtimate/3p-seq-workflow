import argparse
import pandas as pd
import pysam


def write_to_bed(df, out_file):
    """Write the DataFrame to a BED file."""
    df.to_csv(out_file, sep='\t', header=False, index=False)


def contains_motif(sequence, motif="AAUAAAA"):
    """Check if the given RNA sequence contains the specified motif."""
    return motif in sequence


def process_cleavage_sites(pas_df, fasta_path):
    """
    Add a column indicating the presence of the "AAUAAAA" motif.
    """
    fasta_f = pysam.FastaFile(fasta_path)
    motif_results = []

    for _, row in pas_df.iterrows():
        chrom = row[0]               # Chromosome
        start = int(row[1])          # Cleavage site position
        strand = row[5]              # Strand

        # Define sequence boundaries
        if strand == '-':
            sequence_start = start + 10
            sequence_end = start + 35
        else:
            sequence_start = start - 35
            sequence_end = start - 10

        # Fetch the DNA sequence
        dna_sequence = fasta_f.fetch(reference=chrom, start=sequence_start, end=sequence_end + 1)
        rna_sequence = change_dna_to_rna(dna_sequence, strand)

        # Check for the "AAUAAAA" motif
        motif_results.append(1 if contains_motif(rna_sequence) else 0)

    # Add the motif presence column
    pas_df["AAUAAAA"] = motif_results

    return pas_df


def change_dna_to_rna(sequence, strand):
    """
    Convert a DNA sequence to RNA based on the strand direction.
    """
    if strand == '-':
        sequence = sequence[::-1]  # Reverse for the negative strand
        conversion = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
    else:
        conversion = {'A': 'A', 'T': 'U', 'G': 'G', 'C': 'C'}

    return ''.join(conversion.get(base, base) for base in sequence)


def get_args():
    parser = argparse.ArgumentParser(
        description='Add an "AAUAAAA" motif column to a cleavage site BED file'
    )
    parser.add_argument(
        '-i', '--input-bed', required=True, help='Input PAS BED file'
    )
    parser.add_argument(
        '-f', '--fasta-genome', required=True, help='Path to FASTA file'
    )
    parser.add_argument(
        '-o', '--output-bed', required=True, help='Output BED file'
    )

    args = parser.parse_args()
    return args.input_bed, args.fasta_genome, args.output_bed


def run_process():
    input_bed, fasta_genome, output_bed = get_args()
    pas_df = pd.read_csv(input_bed, delimiter='\t', header=None)

    print('Processing cleavage sites...')
    final_df = process_cleavage_sites(pas_df, fasta_genome)

    print('Writing output...')
    write_to_bed(final_df, output_bed)
    print('Done!')


if __name__ == '__main__':
    run_process()
