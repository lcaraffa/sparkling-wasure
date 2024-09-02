import os
import shutil

import laspy
import pytest

from pdaltools.color import color
from pdaltools.las_info import las_info_metadata
from pdaltools.unlock_file import copy_and_hack_decorator, unlock_file

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
TMPDIR = os.path.join(TEST_PATH, "tmp")
LAZ_FILE = os.path.join(TEST_PATH, "data/test_pdalfail_0643_6319_LA93_IGN69.laz")


def setup_module(module):
    try:
        shutil.rmtree(TMPDIR)
    except FileNotFoundError:
        pass
    os.mkdir(TMPDIR)


@copy_and_hack_decorator
def decorated_las_info_metadata(input_file: str):
    return las_info_metadata(input_file)


def test_copy_and_hack_decorator_simple():
    TMP_FILE = os.path.join(TMPDIR, "copy_and_hack_simple.laz")
    shutil.copy(LAZ_FILE, TMP_FILE)
    ret = decorated_las_info_metadata(TMP_FILE)
    assert ret  # Check that the return value of the decorated function is returned


@pytest.mark.geopf
def test_copy_and_hack_decorator_color():
    # bug during laz opening in pdal (solved with copy_and_hack_decorator)
    LAS_FILE = os.path.join(TMPDIR, "test_pdalfail_0643_6319_LA93_IGN69.las")

    # Color works only when an epsg is present in the header or as a parameter
    color(LAZ_FILE, LAS_FILE, "2154", 1)

    las = laspy.read(LAS_FILE)
    print(las.header)
    print(list(las.point_format.dimension_names))
    print(las.red)
    print(las.green)
    print(las.blue)
    print(las.nir)

    assert os.path.isfile(LAS_FILE)


def test_unlock_file():
    TMP_FILE = os.path.join(TMPDIR, "unlock_file.laz")
    shutil.copy(LAZ_FILE, TMP_FILE)
    with pytest.raises(RuntimeError):
        las_info_metadata(TMP_FILE)
    unlock_file(TMP_FILE)
    las_info_metadata(TMP_FILE)
