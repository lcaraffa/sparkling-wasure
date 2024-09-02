"""Replace values of a given attribute in a las/laz file"""

import argparse
import json
import logging
import os
import tempfile
from collections import Counter
from typing import Dict, List

import pdal

from pdaltools.standardize_format import exec_las2las, get_writer_parameters
from pdaltools.unlock_file import copy_and_hack_decorator


def parse_args():
    parser = argparse.ArgumentParser("Replace values of a given attribute in a las/laz file.")
    parser.add_argument("--input_file", type=str, help="Laz input file")
    parser.add_argument("--output_file", type=str, help="Laz output file.")
    parser.add_argument("--attribute", type=str, default="Classification", help="Attribute on which to count values")
    parser.add_argument(
        "--replacement_map",
        type=str,
        help="Path to a json file that contains the values that we want to "
        + "replace, or string that contains the content of such a file."
        + "It should contain a dict like "
        + "{new_value1: [value_to_replace1, value_to_replace2], "
        + "new_value2: [value_to_replace3, ...]}",
    )
    parser.add_argument(
        "--record_format",
        choices=[6, 8],
        type=int,
        required=True,
        help="Record format: 6 (no color) or 8 (4 color channels)",
    )
    parser.add_argument("--projection", default="EPSG:2154", type=str, help="Projection, eg. EPSG:2154")

    return parser.parse_args()


def check_duplicate_values(d: Dict) -> None:
    """Check that a dict that contains lists of values, eg d = {k1: [value1, value2], k2: [value3]}
    has no duplicate values
    """
    all_values = [elt for value in d.values() for elt in value]
    occurences = Counter(all_values)
    for val, count in occurences.items():
        if count > 1:
            raise ValueError(f"Duplicate value {val} provided more than once (count={count})")


def dict_to_pdal_assign_list(d: Dict, output_attribute: str = "Classification", input_attribute: str = "tmp") -> List:
    """Create an assignment list (to be passed to pdal) from a dictionary of type
    d = {
        output_val1: [input_val1, input_val2],
        output_val2: [input_val3],
    }
    that maps values of input_attribute to the values to assign to output_attribute
    """
    check_duplicate_values(d)
    assignment_list = []
    for output_val, input_values in d.items():
        for input_val in input_values:
            assignment_list.append(f"{output_attribute} = {output_val} WHERE {input_attribute} == {input_val}")

    return assignment_list


@copy_and_hack_decorator
def replace_values(
    input_file: str,
    output_file: str,
    replacement_map: Dict,
    attribute: str = "Classification",
    writer_parameters: Dict = {},
) -> None:
    """
    Replace values of attribute {attribute} using a replacement map
    """
    temp_attribute = "tmp"
    assignment_list = dict_to_pdal_assign_list(replacement_map, attribute, temp_attribute)
    pipeline = pdal.Reader.las(input_file)
    pipeline |= pdal.Filter.ferry(dimensions=f"{attribute} => {temp_attribute}")
    pipeline |= pdal.Filter.assign(value=assignment_list)
    # the temp_attribute dimension should not be written as long as the writer has no "extra_dims"
    # parameter
    pipeline |= pdal.Writer(filename=output_file, forward="all", **writer_parameters)

    pipeline.execute()


def parse_replacement_map_from_path_or_json_string(replacement_map):
    if os.path.isfile(replacement_map):
        with open(replacement_map, "r") as f:
            parsed_map = json.load(f)
    else:
        try:
            parsed_map = json.loads(replacement_map)
        except json.decoder.JSONDecodeError as e:
            print("Invalid json string or file path, please check args.replacement_map input:")
            print(replacement_map)
            raise e

    return parsed_map


def replace_values_clean(
    input_file: str,
    output_file: str,
    replacement_map: Dict,
    attribute: str = "Classification",
    writer_parameters: Dict = {},
):
    filename = os.path.basename(output_file)
    with tempfile.NamedTemporaryFile(suffix=filename) as tmp:
        replace_values(input_file, tmp.name, replacement_map, attribute, writer_parameters)
        exec_las2las(tmp.name, output_file)


def main():
    args = parse_args()
    writer_params_from_parser = dict(dataformat_id=args.record_format, a_srs=args.projection)
    writer_parameters = get_writer_parameters(writer_params_from_parser)
    replacement_map = parse_replacement_map_from_path_or_json_string(args.replacement_map)

    replace_values_clean(args.input_file, args.output_file, replacement_map, args.attribute, writer_parameters)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
