import logging
import os

import pdal

from pdaltools.las_info import parse_filename


def create_filenames(file: str, tile_width: int = 1000, tile_coord_scale: int = 1000):
    """Generate the name of the tiles around the input LIDAR tile
    It supposes that the file names are formatted as {prefix1}_{prefix2}_{coordx}_{coordy}_{suffix}
    with coordx and coordy having at least 4 digits

    For example Semis_2021_0000_1111_LA93_IGN69.las

    Args:
        file(str): name of LIDAR file
        tile width (int): width of tiles in meters (usually 1000m)
        tile_coord_scale (int) : scale used in the filename to describe coordinates in meters
                (usually 1000m)
    Returns:
        list_input(list): List of LIDAR's name
    """

    # Create name of LIDAR tiles who cercle the tile
    # # Parameters
    _prefix, coord_x, coord_y, _suffix = parse_filename(file)
    offset = int(tile_width / tile_coord_scale)
    # On left
    _tile_hl = f"{_prefix}_{(coord_x - offset):04d}_{(coord_y + offset):04d}_{_suffix}"
    _tile_ml = f"{_prefix}_{(coord_x - offset):04d}_{coord_y:04d}_{_suffix}"
    _tile_bl = f"{_prefix}_{(coord_x - offset):04d}_{(coord_y - offset):04d}_{_suffix}"
    # On Right
    _tile_hr = f"{_prefix}_{(coord_x + offset):04d}_{(coord_y + offset):04d}_{_suffix}"
    _tile_mr = f"{_prefix}_{(coord_x + offset):04d}_{coord_y:04d}_{_suffix}"
    _tile_br = f"{_prefix}_{(coord_x + offset):04d}_{(coord_y - offset):04d}_{_suffix}"
    # Above
    _tile_a = f"{_prefix}_{coord_x:04d}_{(coord_y + offset):04d}_{_suffix}"
    # Below
    _tile_b = f"{_prefix}_{coord_x:04d}_{(coord_y - offset):04d}_{_suffix}"
    # Return the severals tile's names
    return _tile_hl, _tile_ml, _tile_bl, _tile_a, _tile_b, _tile_hr, _tile_mr, _tile_br


def check_tiles_exist(list_las: list):
    """Check if pointclouds exist
    Args:
        list_las (list): Filenames of the tiles around the LIDAR tile

    Returns:
        li(List): Pruned list of filenames with only existing files
    """
    li = []
    for i in list_las:
        if not os.path.exists(i):
            logging.info(f"NOK : {i}")
            pass
        else:
            li.append(i)
    return li


def create_list(las_dir, input_file, tile_width=1000, tile_coord_scale=1000):
    """Return the paths of 8 tiles around the tile + the input tile
    Args:
        las_dir (str): directory of pointclouds
        input_file (str): path to queried LIDAR tile
        tile_width (int): Width of a tile(in the reference unit: 1m)
        tile_coord_scale (int): Scale used in filename to describe coordinates (usually kilometers)
        1000 * 1m (with 1m being the reference)

    Returns:
        list_files(li): list of tiles
    """

    # Return list 8 tiles around the tile
    list_input = create_filenames(os.path.basename(input_file), tile_width, tile_coord_scale)
    # List pointclouds
    li = [os.path.join(las_dir, e) for e in list_input]
    # Keep only existing files
    li = check_tiles_exist(li)
    # Appending queried tile to list
    li.append(input_file)

    return li


def las_merge(las_dir, input_file, merge_file, tile_width=1000, tile_coord_scale=1000):
    """Merge LIDAR tiles around input_file tile
    Args:
        las_dir (str): directory of pointclouds (to look for neigboprs)
        input_file (str): name of query LIDAR file (with extension)
        output_file (str): path to output
        tile_width (int): Width of a tile(in the reference unit: 1m)
        tile_coord_scale (int): Scale used in filename to describe coordinates (usually kilometers)
        1000 * 1m (with 1m being the reference)
    """
    # List files to merge
    files = create_list(las_dir, input_file, tile_width, tile_coord_scale)
    if len(files) > 0:
        # Merge
        pipeline = pdal.Pipeline()
        for f in files:
            pipeline |= pdal.Reader.las(filename=f)
        pipeline |= pdal.Filter.merge()
        pipeline |= pdal.Writer.las(filename=merge_file, forward="all")
        pipeline.execute()
    else:
        raise ValueError("List of valid tiles is empty : stop processing")
