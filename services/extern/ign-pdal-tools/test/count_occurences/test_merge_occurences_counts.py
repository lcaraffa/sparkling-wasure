import logging
import os
import shutil
from collections import Counter

from pdaltools.count_occurences.merge_occurences_counts import merge_counts

test_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # File is in subdirectory
tmp_path = os.path.join(test_path, "tmp")
input_dir = os.path.join(test_path, "data/counts")
input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith("json")]
output_file = os.path.join(tmp_path, "merged_counts.json")

attribute = "Classification"
expected_counts = Counter(
    {
        "1": 6830,
        "2": 54740,
        "3": 605,
        "4": 2160,
        "5": 42546,
        "6": 33595,
        "0": 83,
    }
)


def setup_module(module):
    try:
        shutil.rmtree(tmp_path)
    except FileNotFoundError:
        pass
    os.mkdir(tmp_path)


def test_merge_counts():
    count = merge_counts(input_files)
    assert count == expected_counts


def test_merge_counts_with_json():
    count = merge_counts(input_files, output_file)
    assert count == expected_counts
    assert os.path.isfile(output_file)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_merge_counts()
    test_merge_counts_with_json()
