import pandas as pd
import os
import argparse

# Function to extract unique gene IDs from motifs.bed files
def extract_gene_ids(motifs_files_dir):
    unique_gene_ids = set()
    for motifs_file in os.listdir(motifs_files_dir):
        if motifs_file.endswith(".motifs.bed"):
            file_path = os.path.join(motifs_files_dir, motifs_file)
            df = pd.read_csv(file_path, sep="\t", header=0)
            unique_gene_ids.update(df['gene_id'].unique())
    return unique_gene_ids

# Function to filter the non-overlapping genes BED file
def filter_non_overlapping_bed(non_overlapping_file, gene_ids):
    # Load the non-overlapping BED file
    col_names = ["chr", "start", "end", "gene_id", "score", "strand", "source", "feature", "frame", "attributes"]
    non_overlapping_df = pd.read_csv(non_overlapping_file, sep="\t", names=col_names, header=None)

    # Filter for genes present in gene_ids
    filtered_df = non_overlapping_df[non_overlapping_df['gene_id'].isin(gene_ids)]

    # Extract relevant information from attributes column
    def parse_attributes(attr):
        gene_name = ""
        gene_type = ""
        for item in attr.split(';'):
            if "gene_name" in item:
                gene_name = item.split('"')[1]
            if "gene_type" in item:
                gene_type = item.split('"')[1]
        return gene_name, gene_type

    filtered_df[['gene_name', 'gene_type']] = filtered_df['attributes'].apply(lambda x: pd.Series(parse_attributes(x)))

    # Create output BED format
    output_df = filtered_df[["chr", "start", "end", "gene_id"]]
    output_df['gene'] = filtered_df['gene_id'] + "," + filtered_df['gene_name'] + "," + filtered_df['gene_type']
    output_df = output_df[["chr", "start", "end", "gene"]]

    return output_df

# Main function to process all files and produce one output BED file
def process_all_files(motifs_files_dir, non_overlapping_file, output_file):
    # Extract unique gene IDs from all motifs files
    gene_ids = extract_gene_ids(motifs_files_dir)

    # Filter the non-overlapping genes BED file
    filtered_bed = filter_non_overlapping_bed(non_overlapping_file, gene_ids)

    # Save the combined filtered BED file
    filtered_bed.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Filtered BED file saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter BED files based on motifs.")
    parser.add_argument("--motifs_dir", "-m", required=True, help="Directory containing motifs.bed files.")
    parser.add_argument("--non_overlapping_file", "-n", required=True, help="Non-overlapping genes BED file.")
    parser.add_argument("--output_file", "-o", required=True, help="Path to save the combined filtered BED file.")

    args = parser.parse_args()

    process_all_files(args.motifs_dir, args.non_overlapping_file, args.output_file)
