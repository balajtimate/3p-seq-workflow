import pysam
import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(
        description='Filter BAM file based on read IDs from BED file and modify BED file.'
    )
    parser.add_argument(
        '-b', '--bam', dest='bam_file', required=True, help='Input BAM file'
    )
    parser.add_argument(
        '-i', '--input-bed', dest='input_bed', required=True, help='Input BED file'
    )
    parser.add_argument(
        '-o', '--output-bam', dest='output_bam', required=True, help='Output BAM file'
    )
    parser.add_argument(
        '-ob', '--output-bed', dest='output_bed', required=True, help='Output BED file'
    )

    args = parser.parse_args()
    return args.bam_file, args.input_bed, args.output_bam, args.output_bed


def filter_bam(bam_file, input_bed, output_bam, output_bed):
    # Read the BED file into a DataFrame
    bed_df = pd.read_csv(input_bed, sep='\t', header=None)

    # Extract read IDs from the 6th column and flatten into a set
    read_ids_column = 6
    read_ids = set()
    for ids in bed_df[read_ids_column]:
        read_ids.update(ids.split(','))

    # Filter BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    filtered_bam = pysam.AlignmentFile(output_bam, "wb", template=bam)

    for read in bam:
        if read.query_name in read_ids:
            filtered_bam.write(read)

    bam.close()
    filtered_bam.close()

    # Remove the ID column from the BED DataFrame
    bed_df.drop(columns=[read_ids_column], inplace=True)

    # Write the updated BED file
    bed_df.to_csv(output_bed, sep='\t', header=False, index=False)


if __name__ == '__main__':
    bam_file, input_bed, output_bam, output_bed = get_args()
    filter_bam(bam_file, input_bed, output_bam, output_bed)
