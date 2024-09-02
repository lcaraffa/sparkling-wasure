import json
import logging
import os
import shutil
from collections import Counter

from pdaltools.count_occurences.count_occurences_for_attribute import compute_count

test_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # File is in subdirectory
tmp_path = os.path.join(test_path, "tmp")
input_dir = os.path.join(test_path, "data/classified_laz")
input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(("las", "laz"))]
output_file = os.path.join(tmp_path, "count.json")
single_input_file = os.path.join(input_dir, "test_data_77050_627755_LA93_IGN69.laz")
counts_single_json = os.path.join(test_path, "data", "counts", "count_test_data_77050_627755_LA93_IGN69.json")

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


def test_count_by_attribute_values():
    count = compute_count(input_files, attribute)
    assert count == expected_counts


def test_count_by_attribute_values_with_json():
    count = compute_count(input_files, attribute, output_file)
    assert count == expected_counts
    assert os.path.isfile(output_file)


def test_count_by_attribute_values_one_file():
    compute_count([single_input_file], attribute, output_file)
    with open(counts_single_json, "r") as f:
        expected = Counter(json.load(f))
    with open(output_file, "r") as f:
        count_file = Counter(json.load(f))
    assert count_file == expected


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_count_by_attribute_values()
    test_count_by_attribute_values_with_json()
