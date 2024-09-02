import json
import logging
import os
from typing import Dict, Tuple

import osgeo.osr as osr
import pdal

osr.UseExceptions()


def las_info_metadata(filename: str):
    r = pdal.Reader.las(filename=filename)
    p = r.pipeline()
    metadata = p.quickinfo["readers.las"]

    return metadata


def get_bounds_from_header_info(metadata):
    bounds = metadata["bounds"]
    minx, maxx, miny, maxy = bounds["minx"], bounds["maxx"], bounds["miny"], bounds["maxy"]

    return minx, maxx, miny, maxy


def get_epsg_from_header_info(metadata):
    if "srs" not in metadata.keys():
        raise RuntimeError("EPSG could not be inferred from metadata: No 'srs' key in metadata.")

    else:
        proj = metadata["srs"]
        # use compoundwkt as recommended in https://github.com/PDAL/python/issues/112
        wkt = proj["compoundwkt"]
        osr_crs = osr.SpatialReference()
        osr_crs.ImportFromWkt(wkt)
        authority = osr_crs.GetAttrValue("AUTHORITY", 0)
        if authority == "EPSG":
            proj = osr_crs.GetAttrValue("AUTHORITY", 1)
        else:
            raise RuntimeError("EPSG could not be inferred from metadata: no attribute 'EPSG' found in metadata srs")

    return proj


def las_info_pipeline(filename: str, spatial_ref: str = "EPSG:2154"):
    """Get las info from pdal pipeline with filter.info
    Args:
        filename: input las
        spatial_ref: spatial reference to pass as 'override_srs' argument in las reader
    """
    information = {
        "pipeline": [
            {"type": "readers.las", "filename": filename, "override_srs": spatial_ref, "nosrs": True},
            {"type": "filters.info"},
        ]
    }

    # Create json
    json_info = json.dumps(information, sort_keys=True, indent=4)
    logging.info(json_info)
    pipeline = pdal.Pipeline(json_info)
    pipeline.execute()
    pipeline.arrays
    # Extract metadata
    metadata = pipeline.metadata

    if isinstance(metadata, str):
        metadata = json.loads(metadata)

    return metadata["metadata"]["filters.info"]


def las_get_xy_bounds(filename: str, buffer_width: int = 0, spatial_ref: str = "EPSG:2154"):
    """Get tile bounds (xy only) from las metadata.
    Try getting bounds using las_info_metadata
    As pdal reader does not seem to read spatial reference properly on some data (TerraSolid output for ex),
    fallback to las_info_pipeline with a default spatial reference

    Args:
        filename (str): full path of file for which to get the bounding box
        buffer_width (str): number of pixel to add to the bounding box on each side (buffer size)
        spatial_ref: spatial reference to pass as 'override_srs' argument in las reader
        (Used if using the spatial reference from the las failed)

    Returns:
        bounds(tuple) : Tuple of bounding box from the LIDAR tile with potential buffer
    """
    # Parameters
    _x = []
    _y = []
    bounds = []
    try:
        metadata = las_info_metadata(filename)
        bounds_dict = metadata["bounds"]

    except RuntimeError:
        metadata = las_info_pipeline(filename, spatial_ref)
        bounds_dict = metadata["bbox"]
    if isinstance(metadata, str):
        metadata = json.loads(metadata)
    # Export bound (maxy, maxy, minx and miny), then creating a buffer with 100 m
    _x.append(float((bounds_dict["minx"]) - buffer_width))  # coordinate minX
    _x.append(float((bounds_dict["maxx"]) + buffer_width))  # coordinate maxX
    _y.append(float((bounds_dict["miny"]) - buffer_width))  # coordinate minY
    _y.append(float((bounds_dict["maxy"]) + buffer_width))  # coordinate maxY
    bounds.append(_x)  # [xmin, xmax]
    bounds.append(_y)  # insert [ymin, ymax]

    return tuple(i for i in bounds)


def parse_filename(file: str):
    """Parse filename and return prefix, suffix and coordinates.
    It expects coordinates to be formatted as {prefix1}_{prefix2}_{coordx}_{coordy}_{suffix}
    For example Semis_2021_0000_1111_LA93_IGN69.las"""
    basename = os.path.basename(file)  # Make sure that we work on the base name and not the full path

    prefix1, prefix2, coordx, coordy, suffix = basename.split("_", 4)
    prefix = f"{prefix1}_{prefix2}"

    return prefix, int(coordx), int(coordy), suffix


def get_buffered_bounds_from_filename(
    filename: str, buffer_width: int = 0, tile_width: int = 1000, tile_coord_scale: int = 1000
) -> Tuple:
    """Get tile bounds (xy only) from las metadata.
    Try getting bounds using las_info_metadata
    As command "pdal_info --metadata" does not seem to work properly on some data
    (TerraSolid output for ex), fallback to las_info_pipeline

    Args:
        filename (str): full path of file for which to get the bounding box
        buffer_width (str): number of pixel to add to the bounding box on each side (buffer size)
        tile_width (int): Width of a tile(in the reference unit: 1m)
        tile_coord_scale (int): Scale used in filename to describe coordinates (usually kilometers)
            1000 * 1m (with 1m being the reference)

    Returns:
        bounds(tuple) : Tuple of bounding box from the LIDAR tile with potential buffer
    """
    _, coordX, coordY, _ = parse_filename(filename)
    # Coordinates in the filenames are x_min and y_max
    minX = coordX * tile_coord_scale
    maxX = coordX * tile_coord_scale + tile_width
    maxY = coordY * tile_coord_scale
    minY = coordY * tile_coord_scale - tile_width

    xs = [minX - buffer_width, maxX + buffer_width]
    ys = [minY - buffer_width, maxY + buffer_width]

    return (xs, ys)


def get_writer_parameters_from_reader_metadata(metadata: Dict, a_srs=None) -> Dict:
    """As pdal las writers does not permit to pass easily metadata from one file as
    parameters for a writer, use a trick to generate writer parameters from the
    reader metadata of a previous pipeline:
    This function uses the metadata from the reader of a pipeline to provide parameters
    to pass to the writer of another pipeline

    To be removed once https://github.com/PDAL/python/issues/147 is solved

    Args:
        metadata (Dict): metadata of an executed pipeline (that can be accessed using pipeline.metadata)
    Returns:
        Dict: parameters to pass to a pdal writer
    """

    reader_metadata = metadata["metadata"]["readers.las"]

    params = {
        "major_version": reader_metadata["major_version"],
        "minor_version": reader_metadata["minor_version"],
        "global_encoding": reader_metadata["global_encoding"],
        "extra_dims": "all",
        "scale_x": reader_metadata["scale_x"],
        "scale_y": reader_metadata["scale_y"],
        "scale_z": reader_metadata["scale_z"],
        "offset_x": reader_metadata["offset_x"],
        "offset_y": reader_metadata["offset_y"],
        "offset_z": reader_metadata["offset_z"],
        "dataformat_id": reader_metadata["dataformat_id"],
        "a_srs": a_srs if a_srs else reader_metadata["comp_spatialreference"],
    }
    return params
