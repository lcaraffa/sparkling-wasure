"""Merge occurences counts from ./count_occurences_for_attribute.py
This is intended to be used after running count_occurences_for_attribute.py in parallel
on several files"""

import argparse
import json
import logging
import os
from collections import Counter
from typing import List

from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser("Count points with each value of an attribute.")
    parser.add_argument(
        "--input_files",
        nargs="+",
        type=str,
        help="List of json input files separated by spaces, or directory " + "containing json files",
    )
    parser.add_argument("--output_file", type=str, help="Output json file containing the counts")

    return parser.parse_args()


def merge_counts(input_files: List[str], output_file=""):
    all_counts = Counter()
    # refresh status bar at most every 1/100 iter cf. https://github.com/tqdm/tqdm/issues/1429
    for input_f in tqdm(input_files, miniters=int(len(input_files) / 100), maxinterval=float("inf")):
        with open(input_f, "r") as f:
            count = Counter(json.load(f))

        all_counts += count

    text = ["Number of point per class:"] + [f"Class {k} :: {v:,d}" for k, v in all_counts.items()]
    logging.info("\n".join(text))

    if output_file:
        with open(output_file, "w") as f:
            json.dump(all_counts, f, indent=2)

    return all_counts


def main():
    args = parse_args()
    if len(args.input_files) == 1 and os.path.isdir(args.input_files[0]):
        input_dir = args.input_files[0]
        input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.lower().endswith("json")]
        logging.info(f"Input_files is a directory. Run on {len(input_files)} from this directory.")
    else:
        logging.info(f"Input_files is a list of files. Run on {len(args.input_files)} files.")
        input_files = args.input_files

    merge_counts(input_files, args.output_file)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
