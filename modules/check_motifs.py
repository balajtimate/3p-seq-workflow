import argparse
import pandas as pd
import pysam


def write_to_bed(df, out_file):
    """Write the DataFrame to a BED file."""
    header=['chr', 'start', 'end', 'cs_id', 'count', 'strand', 'gene_id', 'iqr', 'AAAA_5', 'AAAA_10', 'AAAA_20',
            "AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATGAA", "CATAAA", "AATACA", "AATATA", "GATAAA", "ACTAAA",
            "AATAGA", "AAGAAA"]
    df.to_csv(out_file, sep='\t', header=header, index=False)


def process_motifs(pas_df, fasta_path):
    """
    Add columns for motif presence to the DataFrame.
    """
    fasta_f = pysam.FastaFile(fasta_path)
    motif_results = {motif: [] for motif in MOTIFS}

    for _, row in pas_df.iterrows():
        chrom = row[0]
        cleavage_site = int(row[1])
        strand = row[5]

        # Process motif presence
        if strand == '-':
            sequence_start = cleavage_site + 10
            sequence_end = cleavage_site + 35
        else:
            sequence_start = cleavage_site - 35
            sequence_end = cleavage_site - 10

        dna_sequence = fasta_f.fetch(reference=chrom, start=sequence_start, end=sequence_end + 1)
        rna_sequence = change_dna_to_rna(dna_sequence, strand)

        for motif in MOTIFS:
            motif_results[motif].append(1 if motif in rna_sequence else 0)

    # Add results to the DataFrame
    for motif, values in motif_results.items():
        pas_df[motif] = values

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
    parser = argparse.ArgumentParser(description='Analyze motifs in a cleavage site BED file.')
    parser.add_argument('-i', '--input-bed', required=True, help='Input PAS BED file')
    parser.add_argument('-f', '--fasta-genome', required=True, help='Path to FASTA file')
    parser.add_argument('-o', '--output-bed', required=True, help='Output BED file')

    return parser.parse_args()


def run_process():
    args = get_args()
    pas_df = pd.read_csv(args.input_bed, delimiter='\t', header=None)

    print('Processing motifs...')
    final_df = process_motifs(pas_df, args.fasta_genome)

    print('Writing output...')
    write_to_bed(final_df, args.output_bed)
    print('Done!')


# Constants
MOTIFS = [
    "AAUAAA", "AUUAAA", "UAUAAA", "AGUAAA", "AAUGAA",
    "CAUAAA", "AAUACA", "AAUAUA", "GAUAAA", "ACUAAA",
    "AAUAGA", "AAGAAA"
]

if __name__ == '__main__':
    run_process()
