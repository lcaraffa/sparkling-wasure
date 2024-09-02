import argparse
import logging
import os
import tempfile
from functools import wraps
from pathlib import Path
from typing import Callable, List

import pdal

from pdaltools.las_info import (
    get_buffered_bounds_from_filename,
    get_writer_parameters_from_reader_metadata,
)
from pdaltools.las_merge import create_list
from pdaltools.las_remove_dimensions import remove_dimensions_from_las

ORIGINAL_TILE_TAG = "is_in_original"


def create_las_with_buffer(
    input_dir: str,
    tile_filename: str,
    output_filename: str,
    buffer_width: int = 100,
    spatial_ref: str = "EPSG:2154",
    tile_width: int = 1000,
    tile_coord_scale: int = 1000,
    tag_original_tile: bool = False,
):
    """Merge lidar tiles around the queried tile and crop them in order to add a buffer
    to the tile (usually 100m).

    Args:
        input_dir (str): directory of pointclouds (where you look for neighbors)
        tile_filename (str):  full path to the queried LIDAR tile
        output_filename (str): full path to the saved cropped tile
        buffer_width (int, optional): width of the border to add to the tile (in meters).
        Defaults to 100.
        spatial_ref (_type_, optional): Spatial reference to use to override the one from input las.
        Defaults to "EPSG:2154".
        tile_width (int, optional): width of tiles in meters. Defaults to 1000.
        tile_coord_scale (int, optional): scale used in the filename to describe coordinates
        in meters. Defaults to 1000.
        tag_original_tile (bool, optional): if true, add a new "is_in_original" dimension
        to the output las, equal to 1 on points that belong to the original tile, 0 on points
        that belong to the added buffer. Defaults to False.
    """

    bounds = get_buffered_bounds_from_filename(
        tile_filename, buffer_width=buffer_width, tile_width=tile_width, tile_coord_scale=tile_coord_scale
    )

    logging.debug(f"Add buffer of size {buffer_width} to tile.")
    las_merge_and_crop(
        input_dir,
        tile_filename,
        bounds,
        output_filename,
        spatial_ref,
        tile_width=tile_width,
        tile_coord_scale=tile_coord_scale,
        tag_original_tile=tag_original_tile,
    )


def las_merge_and_crop(
    input_dir: str,
    tile_filename: str,
    bounds: List,
    output_filename: str,
    spatial_ref: str = "EPSG:2154",
    tile_width=1000,
    tile_coord_scale=1000,
    tag_original_tile: bool = False,
):
    """Merge and crop las in a single pipeline (for buffer addition)

    For performance reasons, instead of using a pipeline that reads all files, merge them and
    then crop to the desired bbox, what is done is:
    - For each file:
        - read it
        - crop it according to the bounds
        - optionally add a dimension to differentiate points from the central pointscloud
        from those added as a buffer
        - keep the crop in memory
        - delete the pipeline object to release the memory taken by the las reader
    - Merge the already cropped data

    Args:
        input_dir (str): directory of pointclouds (where you look for neighbors)
        tile_filename (str): full path to the queried LIDAR tile
        bounds (List): 2D bounding box to crop to : provided as ([xmin, xmax], [ymin, ymax])
        output_filename (str): full path to the saved cropped tile
        spatial_ref (str, optional): spatial reference for the writer. Defaults to "EPSG:2154".
        tile_width (int, optional): width of tiles in meters (usually 1000m). Defaults to 1000.
        tile_coord_scale (int, optional): scale used in the filename to describe coordinates in meters.
        Defaults to 1000.
        tag_original_tile (bool, optional):  if true, add a new "is_in_original" dimension
        to the output las, equal to 1 on points that belong to the original tile, 0 on points
        that belong to the added buffer. Defaults to False.
    Raises:
        ValueError: if the list of tiles to merge is empty
    """

    # List files to merge
    files_to_merge = create_list(input_dir, tile_filename, tile_width, tile_coord_scale)
    central_file = files_to_merge[-1]
    if len(files_to_merge) > 0:
        # Read and crop each file
        crops = []
        for f in files_to_merge:
            pipeline = pdal.Pipeline()
            pipeline |= pdal.Reader.las(filename=f, override_srs=spatial_ref)
            if tag_original_tile:
                pipeline |= pdal.Filter.ferry(dimensions=f"=>{ORIGINAL_TILE_TAG}")
                pipeline |= pdal.Filter.assign(value=f"{ORIGINAL_TILE_TAG}={int(f == central_file)}")
            pipeline |= pdal.Filter.crop(bounds=str(bounds))
            pipeline.execute()
            if len(pipeline.arrays[0]) == 0:
                logging.warning(f"File {f} ignored in merge/crop: No points in crop bounding box")
            else:
                crops.append(pipeline.arrays[0])

            if f == central_file:
                # Retrieve metadata before the pipeline is deleted
                metadata = pipeline.metadata
            del pipeline

        params = get_writer_parameters_from_reader_metadata(metadata, a_srs=spatial_ref)

        # Merge
        pipeline = pdal.Filter.merge().pipeline(*crops)

        # Write
        pipeline |= pdal.Writer(filename=output_filename, forward="all", **params)

        logging.info(pipeline.toJSON())
        pipeline.execute()
    else:
        raise ValueError("List of valid tiles is empty : stop processing")
    pass


def remove_points_from_buffer(input_file: str, output_file: str):
    """Remove the points that were added as a buffer to a las file using the "is_in_original"
    dimension that has been added by create_las_with_buffer

    Limitation: if any point has been added to the point cloud after adding the buffer, it
    won't be preserved by this operation (only points from the original file are kept)

    Args:
        input_file (str): path to the input file containing the "is_in_original" dimension
        output_file (str): path to the output_file
    """
    with tempfile.NamedTemporaryFile(suffix="_with_additional_dim.las") as tmp_las:
        pipeline = pdal.Pipeline() | pdal.Reader.las(input_file)
        pipeline |= pdal.Filter.range(limits=f"{ORIGINAL_TILE_TAG}[1:1]")
        pipeline |= pdal.Writer.las(filename=tmp_las.name, forward="all", extra_dims="all")
        pipeline.execute()

        remove_dimensions_from_las(tmp_las.name, dimensions=[ORIGINAL_TILE_TAG], output_las=output_file)


def run_on_buffered_las(
    buffer_width: int, spatial_ref: str, tile_width: int = 1000, tile_coord_scale: int = 1000
) -> Callable:
    """Decorator to apply a function that takes a las/laz as input and returns a las/laz output
    on an input with an additional buffer, then remove the buffer points from the output

    The first argument of the decorated function must be an input path
    The second argument of the decorated function must be an output path

    The buffer is added by merging lidar tiles around the queried tile and crop them based
    on their filenames.

    Limitation: if any point has been added to the point cloud by the decorated function, it
    won't be preserved by this operation (only points from the original file are kept)


    Args:
        buffer_width (int): width of the border to add to the tile (in meters)
        spatial_ref (str): spatial reference for the writer. Example: "EPSG:2154".
        tile_width (int, optional): width of tiles in meters (usually 1000m). Defaults to 1000.
        tile_coord_scale (int, optional): scale used in the filename to describe coordinates in meters.
        Defaults to 1000.

    Raises:
        FileNotFoundError: when the first argument of the decorated function is not an existing
            file
        FileNotFoundError:  when the second argument of the decorated function is not a path
            with an existing parent folder

    Returns:
        Callable: decorated function
    """
    """Decorator to run a function that takes a las as input and returns a las output
    on a las with an additional buffer, then remove the buffer points from the buffer points

    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            input_file = args[0]
            output_file = args[1]
            if not Path(input_file).is_file():
                raise FileNotFoundError(
                    f"File {args[0]} not found. The first argument of a function decorated by "
                    "'run_on_buffered_las' is expected to be the path to an existing input file."
                )

            if not Path(output_file).parent.is_dir():
                raise FileNotFoundError(
                    f"Parent folder for file {args[1]} not found. The second argument of a function "
                    "decorated by 'run_on_buffered_las' is expected to be the path to an output "
                    "file in an existing folder"
                )

            with (
                tempfile.NamedTemporaryFile(suffix="_buffered_input.laz", dir=".") as buf_in,
                tempfile.NamedTemporaryFile(suffix="_buffered_output.laz", dir=".") as buf_out,
            ):
                create_las_with_buffer(
                    Path(input_file).parent,
                    input_file,
                    buf_in.name,
                    buffer_width=buffer_width,
                    spatial_ref=spatial_ref,
                    tile_width=tile_width,
                    tile_coord_scale=tile_coord_scale,
                    tag_original_tile=True,
                )
                func(buf_in.name, buf_out.name, *args[2:], **kwargs)

                remove_points_from_buffer(buf_out.name, output_file)

            return

        return wrapper

    return decorator


def parse_args():
    parser = argparse.ArgumentParser("Add a buffer to a las tile by stitching with its neighbors")
    parser.add_argument(
        "--input_dir",
        "-i",
        type=str,
        required=True,
        help="Path to the the folder containing the tile to which you want to add buffer"
        + "as well as its neighbors tiles",
    )
    parser.add_argument(
        "--tile_filename", "-f", type=str, required=True, help="Filename of the input tile (basename only)"
    )
    parser.add_argument("--output_dir", "-o", type=str, required=True, help="Directory folder for saving the outputs")
    parser.add_argument(
        "--buffer_width",
        "-b",
        default=100,
        type=int,
        help="Width (in meter) for the buffer that is added to the tile before interpolation "
        + "(to prevent artefacts)",
    )
    # Optional parameters
    parser.add_argument(
        "--spatial_reference", default="EPSG:2154", help="Spatial reference to use to override the one from input las."
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    create_las_with_buffer(
        input_dir=args.input_dir,
        tile_filename=os.path.join(args.input_dir, args.tile_filename),
        output_filename=os.path.join(args.output_dir, args.tile_filename),
        buffer_width=args.buffer_width,
        spatial_ref=args.spatial_reference,
    )
