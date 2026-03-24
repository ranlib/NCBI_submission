#!/usr/bin/env python
import argparse
import pandas as pd
import yaml
import os


def load_config(config_path):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def parse_sample_sheet(sample_sheet_path):
    data_section = False
    headers = []
    records = []

    with open(sample_sheet_path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("["):
                data_section = (line == "[BCLConvert_Data]")
                continue

            if data_section:
                if not headers:
                    headers = line.split(",")
                else:
                    values = line.split(",")
                    record = dict(zip(headers, values))
                    records.append(record)

    return records


def generate_fastq_map(records, fastq_dir=None):
    fastq_map = {}

    for i, record in enumerate(records, start=1):
        sample_id = record["Sample_ID"]

        r1 = f"{sample_id}_S{i}_R1_001.fastq.gz"
        r2 = f"{sample_id}_S{i}_R2_001.fastq.gz"

        # Prepend directory if provided
        if fastq_dir:
            r1 = os.path.join(fastq_dir, r1)
            r2 = os.path.join(fastq_dir, r2)

        fastq_map[sample_id] = {
            "R1": r1,
            "R2": r2
        }

    return fastq_map


def main():
    parser = argparse.ArgumentParser(
        description="Augment metadata with FASTQ file names from sample sheet"
    )
    parser.add_argument("-i", "--input_csv", required=True)
    parser.add_argument("-o", "--output_excel", required=True)
    parser.add_argument("-y", "--yaml_config", required=True)
    parser.add_argument("-s", "--sample_sheet", required=True)
    parser.add_argument(
        "-d", "--fastq_dir",
        required=False,
        help="Directory path where FASTQ files are located"
    )

    args = parser.parse_args()

    # Load config
    config = load_config(args.yaml_config)
    columns_to_add = config.get("columns", [])

    # Read metadata CSV
    df = pd.read_csv(args.input_csv)

    # Add columns from config
    if isinstance(columns_to_add, dict):
        for col, value in columns_to_add.items():
            df[col] = value
    else:
        for col in columns_to_add:
            df[col] = None

    # Ensure required columns exist
    df["illumina_sra_file_path_1"] = None
    df["illumina_sra_file_path_2"] = None

    # Set file_location
    df["file_location"] = "local"

    # Set library layout
    df["illumina_library_layout"] = "paired"

    # Parse sample sheet
    records = parse_sample_sheet(args.sample_sheet)
    fastq_map = generate_fastq_map(records, fastq_dir=args.fastq_dir)

    # Map FASTQ names into dataframe
    for idx, row in df.iterrows():
        sample_id = str(row.get("ncbi-spuid-sra"))

        if sample_id in fastq_map:
            df.at[idx, "illumina_sra_file_path_1"] = fastq_map[sample_id]["R1"]
            df.at[idx, "illumina_sra_file_path_2"] = fastq_map[sample_id]["R2"]

    # Write Excel
    df.to_excel(args.output_excel, index=False)


if __name__ == "__main__":
    main()
