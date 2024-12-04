import argparse
import pandas as pd
import pysam


def write_to_bed(df, out_file):
    """Write the DataFrame to a BED file."""
    df.to_csv(out_file, sep='\t', header=False, index=False)


def count_stretches(sequence, stretch="AAAA", window=(-5, 5)):
    """
    Count the number of non-overlapping occurrences of a stretch in the sequence.
    Overlaps are counted separately.
    """
    count = 0
    for i in range(len(sequence) - len(stretch) + 1):
        if sequence[i:i + len(stretch)] == stretch:
            count += 1
    return count


def contains_motif(sequence, motifs):
    """Check for the presence of multiple motifs in the given sequence."""
    return {motif: (1 if motif in sequence else 0) for motif in motifs}


def process_cleavage_sites(pas_df, fasta_path):
    """
    Add columns for AAAA stretches and motif presence.
    """
    fasta_f = pysam.FastaFile(fasta_path)
    motif_results = {motif: [] for motif in MOTIFS}
    stretch_results = {f"AAAA_{window[0]}:{window[1]}": [] for window in STRETCH_WINDOWS}

    for _, row in pas_df.iterrows():
        chrom = row[0]               # Chromosome
        cleavage_site = int(row[1])  # Cleavage site position
        strand = row[5]              # Strand

        # Calculate stretch windows
        for window_name, window in zip(stretch_results.keys(), STRETCH_WINDOWS):
            sequence_start = cleavage_site + window[0]
            sequence_end = cleavage_site + window[1]
            dna_sequence = fasta_f.fetch(reference=chrom, start=sequence_start, end=sequence_end + 1)
            rna_sequence = change_dna_to_rna(dna_sequence, strand)
            stretch_results[window_name].append(count_stretches(rna_sequence, stretch="AAAA"))

        # Define motif detection window
        if strand == '-':
            sequence_start = cleavage_site + 10
            sequence_end = cleavage_site + 35
        else:
            sequence_start = cleavage_site - 35
            sequence_end = cleavage_site - 10

        # Fetch DNA and convert to RNA
        dna_sequence = fasta_f.fetch(reference=chrom, start=sequence_start, end=sequence_end + 1)
        rna_sequence = change_dna_to_rna(dna_sequence, strand)

        # Check for motifs
        motif_presence = contains_motif(rna_sequence, MOTIFS)
        for motif, presence in motif_presence.items():
            motif_results[motif].append(presence)

    # Add results to the DataFrame
    for column, values in stretch_results.items():
        pas_df[column] = values

    for motif, values in motif_results.items():
        pas_df[motif] = values

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
        description='Add motif and AAAA stretch columns to a cleavage site BED file'
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


# Motif list and stretch windows
MOTIFS = [
    "AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATGAA",
    "CATAAA", "AATACA", "AATATA", "GATAAA", "ACTAAA",
    "AATAGA", "AAGAAA"
]
STRETCH_WINDOWS = [(-5, 5), (-10, 10), (-20, 20)]

if __name__ == '__main__':
    run_process()
