"""Re-write las file with expected format:
    - laz version
    - [TODO] nomenclature ???
    - record format
    - global encoding
    - projection
    - precision
    - no extra-dims
"""

import argparse
import os
import subprocess as sp
import tempfile
from typing import Dict

import pdal

from pdaltools.unlock_file import copy_and_hack_decorator

STANDARD_PARAMETERS = dict(
    major_version="1",
    minor_version="4",  # Laz format version (pdal always write in 1.x format)
    global_encoding=17,  # store WKT projection in file
    compression="true",  # Save to compressed laz format
    extra_dims=[],  # Save no extra_dims
    scale_x=0.01,  # Precision of the stored data
    scale_y=0.01,
    scale_z=0.01,
    offset_x=0,  # No offset
    offset_y=0,
    offset_z=0,
    dataformat_id=6,  # No color by default
    a_srs="EPSG:2154",
)


def parse_args():
    parser = argparse.ArgumentParser("Rewrite laz file with standard format.")
    parser.add_argument("--input_file", type=str, help="Laz input file.")
    parser.add_argument("--output_file", type=str, help="Laz output file")
    parser.add_argument(
        "--record_format", choices=[6, 8], type=int, help="Record format: 6 (no color) or 8 (4 color channels)"
    )
    parser.add_argument("--projection", default="EPSG:2154", type=str, help="Projection, eg. EPSG:2154")
    parser.add_argument(
        "--extra_dims",
        default=[],
        nargs="*",
        type=str,
        help="List of extra dims to keep in the output (default=[], use 'all' to keep all extra dims), "
        "extra_dims must be specified with their type (see pdal.writers.las documentation, eg 'dim1=double')",
    )

    return parser.parse_args()


def get_writer_parameters(new_parameters: Dict) -> Dict:
    """
    Get writer parameters from a set of standard parameters + a new set of parameters that can
    override the standard ones
    """
    params = STANDARD_PARAMETERS | new_parameters

    return params


def rewrite_with_pdal(input_file: str, output_file: str, params_from_parser: Dict) -> None:
    # Update parameters with command line values
    params = get_writer_parameters(params_from_parser)
    pipeline = pdal.Reader.las(input_file)
    pipeline |= pdal.Writer(filename=output_file, forward="all", **params)
    pipeline.execute()


def exec_las2las(input_file: str, output_file: str):
    r = sp.run(["las2las", "-i", input_file, "-o", output_file], stderr=sp.PIPE, stdout=sp.PIPE)
    if r.returncode == 1:
        msg = r.stderr.decode()
        print(msg)
        raise RuntimeError(msg)

    output = r.stdout.decode()
    for line in output.splitlines():
        print(line)


@copy_and_hack_decorator
def standardize(input_file: str, output_file: str, params_from_parser: Dict) -> None:
    filename = os.path.basename(output_file)
    with tempfile.NamedTemporaryFile(suffix=filename) as tmp:
        rewrite_with_pdal(input_file, tmp.name, params_from_parser)
        exec_las2las(tmp.name, output_file)


if __name__ == "__main__":
    args = parse_args()
    params_from_parser = dict(dataformat_id=args.record_format, a_srs=args.projection, extra_dims=args.extra_dims)
    standardize(args.input_file, args.output_file, params_from_parser)
