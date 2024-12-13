import argparse
import pandas as pd
import pysam


def write_to_bed(df, out_file):
    """Write the DataFrame to a BED file."""
    df.to_csv(out_file, sep='\t', header=False, index=False)


def count_stretches(sequence, stretch="AAAA", window=(-5, 5)):
    """
    Count the number of non-overlapping occurrences of a stretch in the sequence.
    """
    count = 0
    for i in range(len(sequence) - len(stretch) + 1):
        if sequence[i:i + len(stretch)] == stretch:
            count += 1
    return count


def process_AAAA_stretches(pas_df, fasta_path):
    """
    Add columns for AAAA stretches to the DataFrame.
    """
    fasta_f = pysam.FastaFile(fasta_path)
    stretch_results = {f"AAAA_{window[0]}:{window[1]}": [] for window in STRETCH_WINDOWS}

    for _, row in pas_df.iterrows():
        chrom = row[0]
        cleavage_site = int(row[1])
        strand = row[5]

        # Process AAAA stretches
        for window_name, window in zip(stretch_results.keys(), STRETCH_WINDOWS):
            sequence_start = cleavage_site + window[0]
            sequence_end = cleavage_site + window[1]
            dna_sequence = fasta_f.fetch(reference=chrom, start=sequence_start, end=sequence_end + 1)
            rna_sequence = change_dna_to_rna(dna_sequence, strand)
            stretch_results[window_name].append(count_stretches(rna_sequence, stretch="AAAA"))

    # Add results to the DataFrame
    for column, values in stretch_results.items():
        pas_df[column] = values

    return pas_df


def change_dna_to_rna(sequence, strand):
    """
    Convert a DNA sequence to RNA based on the strand direction.
    """
    if strand == '-':
        sequence = sequence[::-1]
        conversion = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
    else:
        conversion = {'A': 'A', 'T': 'U', 'G': 'G', 'C': 'C'}

    return ''.join(conversion.get(base, base) for base in sequence)


def get_args():
    parser = argparse.ArgumentParser(description='Analyze AAAA stretches in a cleavage site BED file.')
    parser.add_argument('-i', '--input-bed', required=True, help='Input PAS BED file')
    parser.add_argument('-f', '--fasta-genome', required=True, help='Path to FASTA file')
    parser.add_argument('-o', '--output-bed', required=True, help='Output BED file')

    return parser.parse_args()


def run_process():
    args = get_args()
    pas_df = pd.read_csv(args.input_bed, delimiter='\t', header=None)

    print('Processing AAAA stretches...')
    final_df = process_AAAA_stretches(pas_df, args.fasta_genome)

    print('Writing output...')
    write_to_bed(final_df, args.output_bed)
    print('Done!')


# Constants
STRETCH_WINDOWS = [(-5, 5), (-10, 10), (-20, 20)]

if __name__ == '__main__':
    run_process()
