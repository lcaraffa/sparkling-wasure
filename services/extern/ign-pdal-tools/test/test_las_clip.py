import logging
import os
import shutil
from test.utils import assert_header_info_are_similar

import laspy
import numpy as np

from pdaltools.las_clip import las_crop

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
TMP_PATH = os.path.join(TEST_PATH, "tmp")
INPUT_DIR = os.path.join(TEST_PATH, "data")


def setup_module(module):
    try:
        shutil.rmtree(TMP_PATH)

    except FileNotFoundError:
        pass
    os.mkdir(TMP_PATH)


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
def test_las_crop():
    output_file = os.path.join(TMP_PATH, "cropped.las")

    input_file = os.path.join(INPUT_DIR, "test_data_77055_627760_LA93_IGN69.laz")
    # FYI input min/max are: input_mins = [770550.0, 6277550.0] and input_maxs = [770600.0, 6277600.0]

    expected_output_nb_points = 20862
    expected_out_mins = [770560.0, 6277560.0]
    expected_out_maxs = [770590.0, 6277590.0]
    bounds = ([expected_out_mins[0], expected_out_maxs[0]], [expected_out_mins[1], expected_out_maxs[1]])

    las_crop(input_file, output_file, bounds)

    # check file exists
    assert os.path.isfile(output_file)

    # check difference in bbox
    out_mins, out_maxs = get_2d_bounding_box(output_file)

    assert np.all(out_mins == expected_out_mins)
    assert np.all(out_maxs == expected_out_maxs)

    # check number of points
    assert get_nb_points(output_file) == expected_output_nb_points

    # check that the las format is preserved
    assert_header_info_are_similar(output_file, input_file)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_las_crop()
