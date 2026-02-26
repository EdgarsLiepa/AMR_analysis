#!/usr/bin/env python3
import sys
import pandas as pd
from pathlib import Path

def main():
    # 1. Safety check: Ensure we have at least one input and one output
    if len(sys.argv) < 3:
        print("Usage: python mergeCounts.py <input_files...> <output_file>")
        print("Example: python mergeCounts.py in_dir/*_ARG_counts out_dir/Matrix.tsv")
        sys.exit(1)

    # 2. Parse arguments
    # Everything from index 1 up to the second-to-last item are input files
    input_files = sys.argv[1:-1]
    # The very last item in the list is the output file
    output_file = sys.argv[-1]

    print(f"Found {len(input_files)} files to merge...")
    
    dataframes = []

    # 3. Process each file
    for f in input_files:
        try:
            # Extract the sample name from the very first line
            with open(f, 'r') as file:
                sample_name = file.readline().strip()
                
            # Read the rest of the file into a Pandas DataFrame
            # sep=r'\s+' handles any combination of tabs/spaces automatically
            df = pd.read_csv(f, sep=r'\s+', skiprows=1, header=None, names=['Gene', sample_name])
            
            # Set 'Gene' as the index so Pandas aligns rows perfectly
            df.set_index('Gene', inplace=True)
            dataframes.append(df)
            
        except Exception as e:
            print(f"Warning: Could not process {f}. Error: {e}")

    if not dataframes:
        print("Error: No valid dataframes were created. Exiting.")
        sys.exit(1)

    print("Merging data...")
    # 4. Stitch them all together side-by-side
    master_matrix = pd.concat(dataframes, axis=1)

    # Fill any missing genes with 0 and ensure they are whole integers
    master_matrix = master_matrix.fillna(0).astype(int)

    # 5. Ensure the output directory exists before saving
    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # 6. Save the final matrix
    master_matrix.to_csv(out_path, sep='\t')
    print(f"Success! Master matrix saved to: {out_path}")

if __name__ == "__main__":
    main()
