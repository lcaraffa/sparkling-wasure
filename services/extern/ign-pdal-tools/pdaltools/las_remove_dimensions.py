import argparse
import os

import pdal
from pdaltools.las_info import get_writer_parameters_from_reader_metadata

def remove_dimensions_from_las(input_las: str, dimensions: [str], output_las: str):
    """
    export new las without some dimensions
    """
    pipeline = pdal.Pipeline() | pdal.Reader.las(input_las)
    pipeline.execute()
    points = pipeline.arrays[0]
    input_dimensions = list(points.dtype.fields.keys())
    output_dimensions = [dim for dim in input_dimensions if dim not in dimensions]
    points_pruned = points[output_dimensions]
    params = get_writer_parameters_from_reader_metadata(pipeline.metadata)
    pipeline_end = pdal.Pipeline(arrays=[points_pruned])
    pipeline_end |= pdal.Writer.las(output_las, forward="all", **params)
    pipeline_end.execute()


def parse_args():
    parser = argparse.ArgumentParser("Remove dimensions from las")
    parser.add_argument(
        "--input_las",
        "-i",
        type=str,
        required=True,
        help="Path to the the las for which the dimensions will be removed",
    )
    parser.add_argument(
        "--output_las",
        "-o",
        type=str,
        required=False,
        help="Path to the the output las ; if none, we replace the input las",
    )
    parser.add_argument(
        "--dimensions",
        "-d",
        type=str,
        required=True,
        nargs="+",
        help="The dimension we would like to remove from the point cloud file ; be aware to not remove mandatory "
             "dimensions of las"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    remove_dimensions_from_las(
        input_las=args.input_las,
        dimensions=args.dimensions,
        output_las=args.input_las if args.output_las is None else args.output_las,
    )
