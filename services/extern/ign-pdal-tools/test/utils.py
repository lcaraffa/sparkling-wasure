""" Misc tools used in different tests
"""

import json
import subprocess as sp

import numpy as np


def allclose_absolute(a, b, tol: float = 1e-3) -> bool:
    """Check that values are similar with a given precision
    It checks for an absolute tolerance only (unlike np.allclose)
    This enables using it on values with a big offset
    (eg. check millimeter precision on data containing kilometer values)

    Args:
        a (Tuple, or np.array-like): first value to compare
        b (_type_): second value to compare
        tol (float, optional): Tolerance of the comparison.
        Defaults to 1e-3.

    Returns:
        bool: Result of the comparison (Is b-a < tol)
    """
    if isinstance(a, tuple) or isinstance(b, tuple):
        a = np.array(a)
        b = np.array(b)
    return np.all(np.less(np.abs(b - a), 1e-3))


def get_pdal_infos_summary(f: str):
    r = sp.run(["pdal", "info", "--summary", f], stderr=sp.PIPE, stdout=sp.PIPE)
    json_info = json.loads(r.stdout.decode())
    return json_info


def assert_header_info_are_similar(file_to_compare, expected):
    expected_metadata = get_pdal_infos_summary(expected)["summary"]["metadata"]
    file_to_compare_metadata = get_pdal_infos_summary(file_to_compare)["summary"]["metadata"]
    assert file_to_compare_metadata["minor_version"] == expected_metadata["minor_version"]
    assert file_to_compare_metadata["global_encoding"] == expected_metadata["global_encoding"]
    assert file_to_compare_metadata["dataformat_id"] == expected_metadata["dataformat_id"]


EXPECTED_DIMS_BY_DATAFORMAT = {
    6: set(
        [
            "X",
            "Y",
            "Z",
            "Intensity",
            "ReturnNumber",
            "NumberOfReturns",
            "ScanChannel",
            "ScanDirectionFlag",
            "EdgeOfFlightLine",
            "Classification",
            "UserData",
            "ScanAngleRank",
            "PointSourceId",
            "GpsTime",
            "KeyPoint",
            "Overlap",
            "Synthetic",
            "Withheld",
        ]
    ),
    8: set(
        [
            "X",
            "Y",
            "Z",
            "Intensity",
            "ReturnNumber",
            "NumberOfReturns",
            "ScanChannel",
            "ScanDirectionFlag",
            "EdgeOfFlightLine",
            "Classification",
            "UserData",
            "ScanAngleRank",
            "PointSourceId",
            "GpsTime",
            "Red",
            "Green",
            "Blue",
            "Infrared",
            "KeyPoint",
            "Overlap",
            "Synthetic",
            "Withheld",
        ]
    ),
}
