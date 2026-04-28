
# 2026-04-27 JCF

# Script to batch files in a directory by size limit


############################################################################
from pathlib import Path
import re
from collections import defaultdict
import argparse
import sys
###########################################################################

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Batch files in a directory by size limit"
    )
    parser.add_argument(
        "input_directory",
        help="Directory containing files to batch"
    )
    parser.add_argument(
        "--block-size-gb",
        type=float,
        default=60,
        help="Maximum size per batch in GB (default: 60)"
    )
    parser.add_argument(
        "--file-type",
        default=".fastq.gz",
        help="File type pattern to match (default: .fastq.gz)"
    )
    parser.add_argument(
        "--output-file",
        default="",
        help="Output file path. If empty, print to stdout (default: empty)"
    )
    return parser.parse_args()

###########################################################################

def gb_to_bytes(gb):
    """convrt Gb to bytes"""
    return gb * 10**9

def get_file_sizes(root, recursive=True, pattern="*"):
    """
    Return dict of {filepath: size_in_bytes}.

    Parameters
    ----------
    root : str or Path
        Directory to scan
    recursive : bool
        If True, use rglob (walk subdirs). If False, only top level.
    pattern : str
        Glob pattern (e.g. "*.fastq.gz")
    """
    root = Path(root)

    it = root.rglob(pattern) if recursive else root.glob(pattern)

    return {
        str(p): p.stat().st_size
        for p in it
        if p.is_file()
    }

def sum_paired_fastq_sizes(size_dict):
    """
    size_dict: dict {filepath: size_in_bytes}
    returns: dict {sample_key: {"total_size": int, "r1_files": [paths], "r2_files": [paths]}}
    """
    pair_sizes = defaultdict(lambda: {"total_size": 0, "r1_files": [], "r2_files": []})

    for path, size in size_dict.items():
        fname = Path(path).name

        if not fname.endswith(".fastq.gz"):
            continue

        # take everything before _R1_ or _R2_
        m = re.split(r"_R[12]_", fname, maxsplit=1)
        if len(m) < 2:
            continue  # skip if pattern not found

        sample = m[0]  # <- clean sample key
        pair_sizes[sample]["total_size"] += size

        # Track R1 and R2 files
        if "_R1_" in fname:
            pair_sizes[sample]["r1_files"].append(path)
        elif "_R2_" in fname:
            pair_sizes[sample]["r2_files"].append(path)

    return dict(pair_sizes)

def group_samples_by_size(sample_sizes, max_bin_size):
    """
    Greedy bin-packing (First-Fit Decreasing).

    Parameters
    ----------
    sample_sizes : dict {sample: {"total_size": int, "r1_files": [...], "r2_files": [...]}}
        Output of sum_paired_fastq_sizes
    max_bin_size : int
        Maximum total size per group (bytes)

    Returns
    -------
    list of dicts:
        [
          {"samples": [sample1, sample2, ...], "total_size": int},
          ...
        ]
    """
    # sort largest → smallest
    items = sorted(sample_sizes.items(), key=lambda x: x[1]["total_size"], reverse=True)

    bins = []

    for sample, sample_data in items:
        size = sample_data["total_size"]
        if size > max_bin_size:
            raise ValueError(f"{sample} ({size}) exceeds max_bin_size")

        placed = False

        # try to fit into existing bins
        for b in bins:
            if b["total_size"] + size <= max_bin_size:
                b["samples"].append(sample)
                b["total_size"] += size
                placed = True
                break

        # otherwise create new bin
        if not placed:
            bins.append({
                "samples": [sample],
                "total_size": size
            })

    return bins

def format_batches_to_tsv(batches, sample_sizes, out_path=None):
    """
    Convert grouped batches into TSV format:
    batch_id \t sample_id \t size_Gb \t total_batch_size_Gb \t r1_file \t r2_file

    Parameters
    ----------
    batches : list of dicts
        Output of group_samples_by_size()
        e.g. [{"samples": [...], "total_size": int}, ...]
    sample_sizes : dict
        {sample_id: {"total_size": int, "r1_files": [...], "r2_files": [...]}}
    out_path : str or None
        If provided, write TSV to file. Otherwise return string.

    Returns
    -------
    str (TSV content) if out_path is None
    """

    lines = []
    lines.append("batch_id\tsample_id\tsize_Gb\ttotal_batch_size_Gb\tr1_file\tr2_file")

    for i, batch in enumerate(batches, start=1):
        batch_id = f"batch_{i}"
        total_batch_size_gb = batch["total_size"] / 1e9

        for sample in batch["samples"]:
            size_gb = sample_sizes[sample]["total_size"] / 1e9
            r1_file = sample_sizes[sample]["r1_files"][0] if sample_sizes[sample]["r1_files"] else ""
            r2_file = sample_sizes[sample]["r2_files"][0] if sample_sizes[sample]["r2_files"] else ""
            lines.append(f"{batch_id}\t{sample}\t{size_gb:.3f}\t{total_batch_size_gb:.3f}\t{r1_file}\t{r2_file}")

    tsv = "\n".join(lines)

    if out_path:
        with open(out_path, "w") as f:
            f.write(tsv)
    else:
        return tsv

###########################################################################

def main():
    """Main entry point for the script."""
    args = parse_arguments()
    
    sizes = get_file_sizes(args.input_directory, recursive=True, pattern=f"*{args.file_type}")
    size_pairs = sum_paired_fastq_sizes(sizes)
    groups = group_samples_by_size(size_pairs, max_bin_size=gb_to_bytes(args.block_size_gb))
    out = format_batches_to_tsv(groups, size_pairs, out_path=args.output_file if args.output_file else None)
    
    # Print to stdout if no output file specified
    if not args.output_file:
        print(out)


###########################################################################
if __name__ == "__main__":
    main()

