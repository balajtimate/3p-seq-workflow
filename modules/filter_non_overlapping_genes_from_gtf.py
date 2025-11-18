#!/usr/bin/env python3
"""
Filter non-overlapping, protein-coding genes from a GTF file,
adjust gene-level boundaries based on transcript/exon features,
apply strand-aware downstream extension,
and produce two outputs:
  1. Extended coordinates for all protein-coding genes
  2. Filtered, non-overlapping subset
"""
import argparse
import pandas as pd  # type: ignore
from gtfparse import read_gtf  # type: ignore


def load_and_filter_gtf(gtf_path: str) -> pd.DataFrame:
    print(f"[1/5] Loading GTF from: {gtf_path}")
    df = read_gtf(gtf_path)
    print(f"[1/5] {len(df)} total features read")

    df = df[df["feature"].notna()]

    if "gene_type" in df.columns:
        type_col = "gene_type"
    elif "gene_biotype" in df.columns:
        type_col = "gene_biotype"
    else:
        raise KeyError(
            "Neither 'gene_type' nor 'gene_biotype' found in GTF attributes"
        )

    df = df[df[type_col] == "protein_coding"]
    pc_count = df["gene_id"].nunique()
    print(f"[1/5] {pc_count} protein-coding genes retained")
    return df


def adjust_gene_coordinates(df: pd.DataFrame, extension: int) -> pd.DataFrame:
    print(
        f"[2/5] Computing gene boundaries with {extension} bp downstream extension"
    )
    subset = df[df["feature"] != "gene"]
    coords = (
        subset.groupby("gene_id")
        .agg(
            chrom=("seqname", "first"),
            start=("start", "min"),
            end=("end", "max"),
            strand=("strand", "first"),
        )
        .reset_index()
    )

    if extension:

        def _extend(r):
            if r.strand == "+":
                r.end += extension
            elif r.strand == "-":
                r.start = max(1, r.start - extension)
            return r

        coords = coords.apply(_extend, axis=1)

    print(f"[2/5] {len(coords)} gene boundaries computed")
    return coords


def detect_overlaps(coords: pd.DataFrame) -> set:
    print("[3/5] Detecting overlapping genes")
    sorted_coords = coords.sort_values(["chrom", "start", "end"]).reset_index(
        drop=True
    )
    to_remove = set()

    for i in range(len(sorted_coords)):
        gi, s1, e1, chrom_i = sorted_coords.loc[
            i, ["gene_id", "start", "end", "chrom"]
        ]
        for j in range(i + 1, len(sorted_coords)):
            gj, s2, e2, chrom_j = sorted_coords.loc[
                j, ["gene_id", "start", "end", "chrom"]
            ]
            if chrom_j != chrom_i:
                break
            if s2 > e1:
                break
            if not (e1 < s2 or e2 < s1):
                to_remove.add(gi)
                to_remove.add(gj)
    print(f"[3/5] {len(to_remove)} overlapping genes detected")
    return to_remove


def write_gtf(df: pd.DataFrame, gene_ids: set, out_path: str):
    fixed_cols = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
    ]
    attr_cols = [c for c in df.columns if c not in fixed_cols]
    with open(out_path, "w") as out:
        for _, row in df[df["gene_id"].isin(gene_ids)].iterrows():
            attrs = []
            for c in attr_cols:
                v = row[c]
                if pd.notna(v):
                    attrs.append(f'{c} "{v}";')
            attr_field = " ".join(attrs)

            fields = [
                row.seqname,
                row.source,
                row.feature,
                int(row.start),
                int(row.end),
                row.score if pd.notna(row.score) else ".",
                row.strand,
                row.frame if pd.notna(row.frame) else ".",
                attr_field,
            ]
            out.write("\t".join(map(str, fields)) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Filter non-overlapping, protein-coding genes from a GTF"
    )
    parser.add_argument("--gtf", required=True, help="Input GTF file")
    parser.add_argument(
        "--extended-out",
        required=True,
        help="Output GTF after adjustment and extension",
    )
    parser.add_argument(
        "--filtered-out",
        required=True,
        help="Output filtered GTF (non-overlapping)",
    )
    parser.add_argument(
        "--extend",
        type=int,
        default=0,
        help="Downstream extension in bp (strand-aware)",
    )
    args = parser.parse_args()

    df = load_and_filter_gtf(args.gtf)
    coords = adjust_gene_coordinates(df, extension=args.extend)

    # Overwrite gene-feature rows with adjusted coords
    df = df.merge(
        coords[["gene_id", "start", "end"]],
        on="gene_id",
        how="left",
        suffixes=("", "_new"),
    )
    mask_gene = df["feature"] == "gene"
    df.loc[mask_gene, "start"] = df.loc[mask_gene, "start_new"]
    df.loc[mask_gene, "end"] = df.loc[mask_gene, "end_new"]
    df.drop(columns=["start_new", "end_new"], inplace=True)

    # Step 4/5: write extended-only GTF
    print(f"[4/5] Writing extended GTF to: {args.extended_out}")
    write_gtf(df, set(coords["gene_id"]), args.extended_out)
    print("[4/5] Done.")

    overlaps = detect_overlaps(coords)
    keep_ids = set(coords["gene_id"]) - overlaps
    print(
        f"[3/5] Keeping {len(keep_ids)} of {len(coords)} genes after overlap removal"
    )

    # Step 5/5: write filtered GTF
    print(f"[5/5] Writing filtered GTF to: {args.filtered_out}")
    write_gtf(df, keep_ids, args.filtered_out)
    print("[5/5] Done.")


if __name__ == "__main__":
    main()
