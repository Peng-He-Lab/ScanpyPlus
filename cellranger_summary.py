#!/usr/bin/env python3
import os
import argparse
import pandas as pd

def main(base_dir, output):

    sample_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    
    summary_data = []
    
    for sample_dir in sample_dirs:
        summary_path = os.path.join(base_dir, sample_dir, "outs", "summary.csv")
        if os.path.exists(summary_path):
            df = pd.read_csv(summary_path, index_col=0).transpose()
            df.columns = [sample_dir]
            summary_data.append(df)
        else:
            print(f"Warning: {summary_path} not found")
    
    if summary_data:
        combined_df = pd.concat(summary_data, axis=1)
        combined_df.to_csv(output, sep="\t")
        print(f"Output written to {output}")
    else:
        print("No summary data found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge cellranger summary CSV files from multiple samples into one TSV.")
    parser.add_argument("--base_dir", required=True, help="Base directory containing sample subdirectories.")
    parser.add_argument("--output", default="summary.tsv", help="Output TSV file name.")
    
    args = parser.parse_args()
    main(args.base_dir, args.output)

