#!/usr/bin/env python3
"""
h5_to_readable.py — Extract data from HDF5 (.h5 / .hdf5) files to CSV, TSV, or TXT.

Usage:
    python h5_to_readable.py input.h5
    python h5_to_readable.py input.h5 --format tsv
    python h5_to_readable.py input.h5 --format txt --output my_output
    python h5_to_readable.py input.h5 --dataset /path/to/dataset
    python h5_to_readable.py input.h5 --list

Requirements:
    pip install h5py pandas numpy
"""

import argparse
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"Missing dependency: {e}")
    print("Install with: pip install h5py pandas numpy")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def list_datasets(filepath: str) -> list[tuple[str, tuple, str]]:
    """Return all dataset paths, shapes, and dtypes in the file."""
    results = []

    def _visitor(name, obj):
        if isinstance(obj, h5py.Dataset):
            results.append((name, obj.shape, str(obj.dtype)))

    with h5py.File(filepath, "r") as f:
        f.visititems(_visitor)

    return results


def print_structure(filepath: str) -> None:
    """Pretty-print the file structure."""
    datasets = list_datasets(filepath)
    if not datasets:
        print("No datasets found.")
        return

    print(f"\nDatasets in '{filepath}':")
    print(f"  {'Path':<50} {'Shape':<20} {'Dtype'}")
    print("  " + "-" * 80)
    for path, shape, dtype in datasets:
        print(f"  /{path:<49} {str(shape):<20} {dtype}")
    print()


def dataset_to_dataframe(dataset: h5py.Dataset) -> pd.DataFrame:
    """Convert an h5py Dataset to a pandas DataFrame."""
    data = dataset[()]  # read all data into memory

    # Scalar
    if data.ndim == 0:
        return pd.DataFrame({"value": [data.item()]})

    # 1-D array → single column
    if data.ndim == 1:
        return pd.DataFrame(data, columns=["value"])

    # 2-D array → rows × columns
    if data.ndim == 2:
        cols = [f"col_{i}" for i in range(data.shape[1])]
        return pd.DataFrame(data, columns=cols)

    # 3-D or higher → flatten to 2-D (rows × all other dims)
    original_shape = data.shape
    flat = data.reshape(data.shape[0], -1)
    cols = [f"dim_{i}" for i in range(flat.shape[1])]
    df = pd.DataFrame(flat, columns=cols)
    df.attrs["original_shape"] = str(original_shape)
    return df

    # Structured / compound dtype → treat fields as columns
    if data.dtype.names:
        return pd.DataFrame({name: data[name].ravel() for name in data.dtype.names})

    raise ValueError(f"Cannot convert dataset with shape {data.shape} and dtype {data.dtype}")


def write_output(df: pd.DataFrame, out_path: Path, fmt: str) -> None:
    """Write a DataFrame to the requested format."""
    if fmt == "csv":
        df.to_csv(out_path, index=False)
    elif fmt == "tsv":
        df.to_csv(out_path, index=False, sep="\t")
    elif fmt == "txt":
        # Human-readable fixed-width table
        with open(out_path, "w") as fh:
            fh.write(df.to_string(index=False))
            fh.write("\n")
    else:
        raise ValueError(f"Unknown format: {fmt}")


def safe_filename(path: str) -> str:
    """Turn a dataset path like /group/sub/data into group_sub_data."""
    return path.strip("/").replace("/", "_").replace(" ", "_")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract HDF5 datasets to CSV / TSV / TXT",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input", help="Path to .h5 / .hdf5 file")
    parser.add_argument(
        "--format", "-f",
        choices=["csv", "tsv", "txt"],
        default="csv",
        help="Output format (default: csv)",
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output file base name (no extension). Default: <input_stem>_<dataset>",
    )
    parser.add_argument(
        "--dataset", "-d",
        default=None,
        help="Extract a specific dataset path (e.g. /group/mydata). "
             "Omit to extract ALL datasets.",
    )
    parser.add_argument(
        "--list", "-l",
        action="store_true",
        help="List all datasets and exit (no extraction)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: file not found: {input_path}")
        sys.exit(1)

    # --list mode
    if args.list:
        print_structure(str(input_path))
        sys.exit(0)

    ext = f".{args.format}"

    with h5py.File(str(input_path), "r") as hf:

        # Decide which datasets to extract
        if args.dataset:
            ds_key = args.dataset.lstrip("/")
            if ds_key not in hf:
                print(f"Error: dataset '/{ds_key}' not found in file.")
                print("Use --list to see available datasets.")
                sys.exit(1)
            targets = [ds_key]
        else:
            targets = [name for name, shape, dtype in list_datasets(str(input_path))]
            if not targets:
                print("No datasets found in file.")
                sys.exit(1)
            print(f"Found {len(targets)} dataset(s). Extracting all…\n")

        # Extract each target
        for ds_key in targets:
            obj = hf[ds_key]
            if not isinstance(obj, h5py.Dataset):
                continue  # skip groups

            print(f"  Processing: /{ds_key}  shape={obj.shape}  dtype={obj.dtype}")

            try:
                df = dataset_to_dataframe(obj)
            except Exception as exc:
                print(f"    ⚠ Skipped (could not convert): {exc}")
                continue

            # Determine output path
            if args.output:
                # If multiple datasets, append dataset name to avoid overwriting
                if len(targets) > 1:
                    out_path = Path(f"{args.output}_{safe_filename(ds_key)}{ext}")
                else:
                    out_path = Path(f"{args.output}{ext}")
            else:
                stem = f"{input_path.stem}_{safe_filename(ds_key)}"
                out_path = input_path.parent / f"{stem}{ext}"

            write_output(df, out_path, args.format)
            print(f"    ✓ Saved → {out_path}  ({len(df):,} rows × {len(df.columns)} cols)")

    print("\nDone.")


if __name__ == "__main__":
    main()