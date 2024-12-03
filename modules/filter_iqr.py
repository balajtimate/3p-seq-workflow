import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def calculate_iqr(read_lengths):
    """Calculates the Interquartile Range (IQR) of a list of read lengths."""
    read_lengths = [int(length) for length in read_lengths]
    return np.percentile(read_lengths, 75) - np.percentile(read_lengths, 25)


def process_bed_file(filename):
    """Reads a BED file, calculates IQR, and adds it as a column."""
    df = pd.read_csv(filename, sep='\t', header=None)
    df.columns = ['chr', 'start', 'end', 'id', 'count', 'strand', 'read_lengths', 'read_ids']

    # Process read lengths to calculate IQR
    df['read_lengths'] = df['read_lengths'].str.split(',')
    df['iqr'] = df['read_lengths'].apply(calculate_iqr)

    return df


def bin_and_calculate_rolling_10th_percentile(df):
    """Bins counts into quantiles, calculates the 10th percentile IQR per bin,
    and performs a rolling window smoothing."""
    # Define quantile-based bins
    num_bins = 100  # Number of quantile bins
    df['count_bins'], bin_edges = pd.qcut(df['count'], q=num_bins, retbins=True, duplicates='drop')

    # Calculate the 10th percentile IQR for each bin
    quantiles = df.groupby('count_bins')['iqr'].quantile(0.1).reset_index(name='10th_percentile')

    # Perform rolling window smoothing
    rolling_window = 3
    quantiles['smoothed_10th_percentile'] = quantiles['10th_percentile'].rolling(rolling_window, center=True, min_periods=1).mean()

    # Map smoothed 10th percentile back to the dataframe
    bin_to_smoothed = dict(zip(quantiles['count_bins'], quantiles['smoothed_10th_percentile']))
    df['smoothed_10th_percentile'] = df['count_bins'].map(bin_to_smoothed)

    # Ensure smoothed_10th_percentile is numeric
    df['smoothed_10th_percentile'] = pd.to_numeric(df['smoothed_10th_percentile'], errors='coerce')

    return df


def filter_below_rolling_window(df):
    """Filters entries strictly below the rolling window 10th percentile line."""
    filtered_df = df[df['iqr'] >= df['smoothed_10th_percentile']]
    return filtered_df


def plot_binned_iqr_with_quantile_bins_and_smoothing(df, sample_name):
    """Plots IQRs of read lengths binned by quantiles of `count` and overlays a rolling lower 10% quantile of IQR."""
    num_bins = 100
    df['count_bins'], bin_edges = pd.qcut(df['count'], q=num_bins, retbins=True, duplicates='drop')
    quantiles = df.groupby('count_bins', observed=False)['iqr'].quantile(0.1).reset_index(name='10th_percentile')
    quantiles['smoothed_10th_percentile'] = quantiles['10th_percentile'].rolling(3, center=True, min_periods=1).mean()

    plt.figure(figsize=(12, 6))
    sns.boxplot(x='count_bins', y='iqr', data=df, order=quantiles['count_bins'])
    plt.scatter(
        x=range(len(quantiles)),
        y=quantiles['10th_percentile'],
        color='red',
        label='10th Percentile (per bin)'
    )
    plt.plot(
        range(len(quantiles)),
        quantiles['smoothed_10th_percentile'],
        color='blue',
        label='10th Percentile (Smoothed with Rolling Window)'
    )
    plt.xlabel('Count Bins (Quantile-Based)')
    plt.ylabel('IQR of read lengths')
    plt.title(f"IQR by Quantile Bins with 10th Percentile - {sample_name}")
    plt.xticks(range(len(quantiles)), labels=quantiles['count_bins'].astype(str), rotation=45)
    plt.grid(axis='y')
    plt.legend()
    plt.savefig(f"{sample_name}_IQR_10_quantile.png", dpi=100, bbox_inches='tight')


def main(input_bed, output_bed, sample_name):
    # Step 1: Read and process the BED file
    df = process_bed_file(input_bed)

    # Step 2: Bin data and calculate rolling 10th percentile
    df = bin_and_calculate_rolling_10th_percentile(df)
    print(df.head())  
    # Step 3: Filter entries below the rolling window line
    filtered_df = filter_below_rolling_window(df)
    print(filtered_df.head())
    # Step 4: Save filtered BED file
    filtered_df = filtered_df.drop(columns=['smoothed_10th_percentile', 'count_bins', 'read_lengths'])
    filtered_df.to_csv(output_bed, sep='\t', index=False, header=False)

    print(f"Filtered BED file saved to: {output_bed}")

    # Optional: Generate the plot for visualization
    plot_binned_iqr_with_quantile_bins_and_smoothing(df, sample_name)


def get_args():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Filter a BED file based on IQR and rolling window lower 10% quantile'
    )
    parser.add_argument(
        '-i', '--input-bed', required=True, help='Input BED file with cleavage site data'
    )
    parser.add_argument(
        '-o', '--output-bed', required=True, help='Output filtered BED file'
    )
    parser.add_argument(
        '-s', '--sample-name', required=True, help='Sample name for plotting'
    )

    args = parser.parse_args()
    return args.input_bed, args.output_bed, args.sample_name

if __name__ == "__main__":
    input_bed, output_bed, sample_name = get_args()
    main(input_bed, output_bed, sample_name)
