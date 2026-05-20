import pandas as pd
import numpy as np
import sys
import os
import argparse

def main(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    #eTIS strength = 100 − (predicted leaky scanning / maximum predicted leaky scanning) × 100
    max_leaky_scanning = df['predictions_GFP'].max()
    df['eTIS_strength'] = 100 - (df['predictions_GFP'] / max_leaky_scanning) * 100
    df.to_csv(output_file, sep='\t', index=False)

    print(f"Done! Calculated the eTIS strength on {len(df)} rows based on maximum predicted leaky scanning: {max_leaky_scanning}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate eTIS strength')
    parser.add_argument('-i', '--input', required=True, type=str, help='Input TSV file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output TSV file')
    args = parser.parse_args()

    if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    main(args.input, args.output)