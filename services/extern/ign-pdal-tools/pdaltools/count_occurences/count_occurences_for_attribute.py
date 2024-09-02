"""Count occurences of each value of a given attribute in a set of pointclouds.
Eg. to count points of each class in classified point clouds """

import argparse
import json
import logging
import os
from collections import Counter
from typing import List

import pdal
from tqdm import tqdm

from pdaltools.unlock_file import copy_and_hack_decorator


def parse_args():
    parser = argparse.ArgumentParser("Count points with each value of an attribute.")
    parser.add_argument(
        "--input_files",
        nargs="+",
        type=str,
        help="List of laz input files separated by spaces, or directory " + "containing las/laz files",
    )
    parser.add_argument("--attribute", type=str, default="Classification", help="Attribute on which to count values")
    parser.add_argument("--output_file", type=str, help="Output json file containing the counts")

    return parser.parse_args()


@copy_and_hack_decorator
def compute_count_one_file(filepath: str, attribute: str = "Classification") -> Counter:
    pipeline = pdal.Reader.las(filepath)
    pipeline |= pdal.Filter.stats(dimensions=attribute, count=attribute)
    pipeline.execute()
    # List of "class/count" on the only dimension that is counted
    raw_counts = pipeline.metadata["metadata"]["filters.stats"]["statistic"][0]["counts"]
    split_counts = [c.split("/") for c in raw_counts]
    try:
        # Try to prettify the value by converting it to an integer (eg. for Classification that
        # returns values such as 1.0000 instead of 1 or 1.)
        counts = Counter({str(int(float(value))): int(count) for value, count in split_counts})
    except ValueError:
        # in case value is not a number, float(value) returns a ValueError
        # fallback: use the raw value
        counts = Counter({value: int(count) for value, count in split_counts})

    return counts


def compute_count(input_files: List[str], attribute: str = "Classification", output_file=""):
    all_counts = Counter()
    # refresh status bar at most every 1/100 iter cf. https://github.com/tqdm/tqdm/issues/1429
    for f in tqdm(input_files, miniters=int(len(input_files) / 100), maxinterval=float("inf")):
        logging.debug(f"Counting values of {attribute} for {os.path.basename(f)}")
        all_counts += compute_count_one_file(f, attribute)

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
        input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.lower().endswith(("las", "laz"))]
    else:
        input_files = args.input_files

    compute_count(input_files, args.attribute, args.output_file)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
