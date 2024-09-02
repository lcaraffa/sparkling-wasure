import logging
import os
from test.utils import assert_header_info_are_similar

import laspy
import numpy as np

from pdaltools.las_merge import las_merge

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
TMP_PATH = os.path.join(TEST_PATH, "tmp")
INPUT_DIR = os.path.join(TEST_PATH, "data")


# def setup_module(module):
#     try:
#         shutil.rmtree(tmp_path)

#     except (FileNotFoundError):
#         pass
#     os.mkdir(tmp_path)


# Utils functions
def get_nb_points(path):
    """Get number of points in a las file"""
    with laspy.open(path) as f:
        nb_points = f.header.point_count

    return nb_points


def get_2d_bounding_box(path):
    """Get bbox for a las file (x, y only)"""
    with laspy.open(path) as f:
        mins = f.header.mins
        maxs = f.header.maxs

    return mins[:2], maxs[:2]


# Tests
def test_las_merge():
    coord_x = 77055
    coord_y = 627760
    input_file = os.path.join(INPUT_DIR, f"test_data_{coord_x}_{coord_y}_LA93_IGN69.laz")
    output_file = os.path.join(TMP_PATH, "merged.las")
    tile_width = 50
    tile_coord_scale = 10
    input_mins = [coord_x * tile_coord_scale, coord_y * tile_coord_scale - tile_width]
    input_maxs = [coord_x * tile_coord_scale + tile_width, coord_y * tile_coord_scale]
    expected_out_mins = [input_mins[0] - tile_width, input_mins[1] - tile_width]
    expected_out_maxs = [input_maxs[0] + tile_width, input_maxs[1]]  # There is no tile above the tile to merge

    expected_output_nb_points = 405937

    las_merge(INPUT_DIR, input_file, output_file, tile_width=tile_width, tile_coord_scale=tile_coord_scale)

    # check file exists
    assert os.path.isfile(output_file)

    # check difference in bbox
    in_mins, in_maxs = get_2d_bounding_box(input_file)
    out_mins, out_maxs = get_2d_bounding_box(output_file)

    # check number of points
    assert get_nb_points(output_file) == expected_output_nb_points

    # check bounds are the expected ones
    assert np.all(out_mins == expected_out_mins)
    assert np.all(out_maxs == expected_out_maxs)

    # check that the las format is preserved
    assert_header_info_are_similar(output_file, input_file)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_las_merge()
