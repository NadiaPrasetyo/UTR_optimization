import argparse
import re
from Bio import SeqIO
import logging
import os
import sys
import subprocess

def setup_logging(verbose=False, output_dir='data/infernal_results'):
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = os.path.join(output_dir, "alifilter.log")
    os.makedirs(output_dir, exist_ok=True)  # <-- add this line
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w') if verbose else logging.NullHandler(),
            logging.StreamHandler(sys.stdout)
        ]
    )
    if verbose:
        logging.info("Logging initialized. Log file: %s", log_file)

def reformat_fasta_to_sto(tool_root, input_file, output, verbose):
    reformat_path = os.path.join(tool_root, "easel", "miniapps", "esl-reformat")
    if not os.path.exists(reformat_path):
        raise FileNotFoundError(f"reformat_path not found: {reformat_path}")
        
    os.makedirs(os.path.join(output, "temp"), exist_ok=True)

    tmp_sto = os.path.join(output, "temp", "tmp.sto")

    if verbose:
        logging.info(f"Converting {input_file} to Stockholm...")

    with open(tmp_sto, "w") as out:
        subprocess.run(
            [reformat_path, "--informat", "afa", "stockholm", input_file],
            stdout=out, check=True
        )

    if verbose:
        logging.info("Done.")

    return tmp_sto


def main(input_file, output, tool_root, genome, verbose):

    if not input_file.endswith(".sto"):
        tmp_sto = reformat_fasta_to_sto(tool_root, input_file, output, verbose)
        input_file = tmp_sto
    # Derive a clean stem from just the filename, not the full path
    input_stem = os.path.splitext(os.path.basename(input_file))[0]

    temp_dir = os.path.join(output, "temp")
    os.makedirs(temp_dir, exist_ok=True)  # ensure temp/ exists for .fasta input too

    cmbuild_path = os.path.join(tool_root, "infernal-1.1.5", "src", "cmbuild")
    cm_path = os.path.join(temp_dir, f"{input_stem}.cm")

    # check that the CM path has not already been built
    if os.path.exists(cm_path):
        if verbose:
            logging.info(f"CM already built for {input_file}. Skipping...")
    else:
        if verbose:
            logging.info(f"Building covariance model for {input_file}...")

        subprocess.run(
            [cmbuild_path, cm_path, os.path.abspath(input_file)],
            check=True
        )

        cmcalibrate_path = os.path.join(tool_root, "infernal-1.1.5", "src", "cmcalibrate")

        if verbose:
            logging.info(f"Calibrating covariance model for {input_file}...")

        subprocess.run(
            [cmcalibrate_path, cm_path],
            check=True
        )

    cmsearch_path = os.path.join(tool_root, "infernal-1.1.5", "src", "cmsearch")

    if verbose:
        logging.info(f"Searching for {input_file} in genome {genome}...")
    
    # run cmscan to output a table and a text file
    subprocess.run(
        [cmsearch_path,
         "--tblout", os.path.join(output, f"{input_stem}_cm_scan_table"),
         "-o",       os.path.join(output, f"{input_stem}_cm_scan"),
         cm_path, genome],
        check=True
    )
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter a MSA using Infernal",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Input fasta or stockholm file")
    parser.add_argument("-o", "--output", type=str, default="data/infernal_results", help="Output file")
    parser.add_argument("--tool_root", type=str, required=True, help="Path to Infernal tools")
    parser.add_argument("--genome", type=str, help="Path to genome")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()
    setup_logging(verbose=args.verbose, output_dir=args.output)
    main(args.input_file, args.output, args.tool_root, args.genome, args.verbose)
