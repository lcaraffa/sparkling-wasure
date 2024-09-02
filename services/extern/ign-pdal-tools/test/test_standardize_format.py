import logging
import os
import shutil
import subprocess as sp
from test.utils import EXPECTED_DIMS_BY_DATAFORMAT, get_pdal_infos_summary

import pdal
import pytest

from pdaltools.standardize_format import exec_las2las, rewrite_with_pdal, standardize

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
TMP_PATH = os.path.join(TEST_PATH, "tmp")
INPUT_DIR = os.path.join(TEST_PATH, "data")

MUTLIPLE_PARAMS = [
    {"dataformat_id": 6, "a_srs": "EPSG:2154", "extra_dims": []},
    {"dataformat_id": 8, "a_srs": "EPSG:4326", "extra_dims": []},
    {"dataformat_id": 8, "a_srs": "EPSG:2154", "extra_dims": ["dtm_marker=double", "dsm_marker=double"]},
    {"dataformat_id": 8, "a_srs": "EPSG:2154", "extra_dims": "all"},
]


def setup_module(module):
    try:
        shutil.rmtree(TMP_PATH)

    except FileNotFoundError:
        pass
    os.mkdir(TMP_PATH)


def _test_standardize_format_one_params_set(input_file, output_file, params):
    rewrite_with_pdal(input_file, output_file, params)
    # check file exists
    assert os.path.isfile(output_file)
    # check values from metadata
    json_info = get_pdal_infos_summary(output_file)
    if pdal.info.version < "2.5":
        raise NotImplementedError("This test is not implemented for pdal < 2.5")
    elif pdal.info.version <= "2.5.2":
        metadata = json_info["summary"]["metadata"][1]
    else:
        metadata = json_info["summary"]["metadata"]
    assert metadata["compressed"] is True
    assert metadata["minor_version"] == 4
    assert metadata["global_encoding"] == 17
    assert metadata["dataformat_id"] == params["dataformat_id"]
    # Check that there is no extra dim
    dimensions = set([d.strip() for d in json_info["summary"]["dimensions"].split(",")])
    if params["extra_dims"] == "all":
        assert EXPECTED_DIMS_BY_DATAFORMAT[params["dataformat_id"]].issubset(dimensions)
    else:
        extra_dims_names = [dim.split("=")[0] for dim in params["extra_dims"]]
        assert dimensions == EXPECTED_DIMS_BY_DATAFORMAT[params["dataformat_id"]].union(extra_dims_names)

    # TODO: Check srs
    # TODO: check precision


def test_standardize_format():
    input_file = os.path.join(INPUT_DIR, "test_data_77055_627755_LA93_IGN69_extra_dims.laz")
    output_file = os.path.join(TMP_PATH, "formatted.laz")
    for params in MUTLIPLE_PARAMS:
        _test_standardize_format_one_params_set(input_file, output_file, params)


def exec_lasinfo(input_file: str):
    r = sp.run(["lasinfo", "-stdout", input_file], stderr=sp.PIPE, stdout=sp.PIPE)
    if r.returncode == 1:
        msg = r.stderr.decode()
        print(msg)
        raise RuntimeError(msg)

    output = r.stdout.decode()
    return output


def assert_lasinfo_no_warning(input_file: str):
    errors = [line for line in exec_lasinfo(input_file).splitlines() if "WARNING" in line]

    for line in errors:
        print(line)

    assert errors == [], errors


def test_exec_las2las_error():
    with pytest.raises(RuntimeError):
        exec_las2las("not_existing_input_file", "output_file")


def test_standardize_does_NOT_produce_any_warning_with_Lasinfo():
    # bad file on the store (44 Mo)
    # input_file = (
    #     "/var/data/store-lidarhd/developpement/standaLAS/demo_standardization/Semis_2022_0584_6880_LA93_IGN69.laz"
    # )

    input_file = os.path.join(TEST_PATH, "data/classified_laz/test_data_77050_627755_LA93_IGN69.laz")
    output_file = os.path.join(TMP_PATH, "test_standardize_produce_no_warning_with_lasinfo.las")

    # if you want to see input_file warnings
    # assert_lasinfo_no_warning(input_file)

    standardize(input_file, output_file, MUTLIPLE_PARAMS[0])
    assert_lasinfo_no_warning(output_file)


def test_standardize_malformed_laz():
    input_file = os.path.join(TEST_PATH, "data/test_pdalfail_0643_6319_LA93_IGN69.laz")
    output_file = os.path.join(TMP_PATH, "standardize_pdalfail_0643_6319_LA93_IGN69.laz")
    standardize(input_file, output_file, MUTLIPLE_PARAMS[0])
    assert os.path.isfile(output_file)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_standardize_format()
