import numpy as np

import convertcloud

ref = np.array([[0.123, 4.567, 8.901],
                [2.345, 6.789, 0.123]], dtype=np.float32)

converter = convertcloud.Converter()


def test_load_xyz():
    converter.load_points("filetest.xyz")
    assert np.all(converter.points == ref)

def test_load_a3d():
    converter.load_points("filetest.a3d")
    assert np.all(converter.points == ref)

def test_load_ply():
    converter.load_points("filetest.ply")
    assert np.all(converter.points == ref)

def test_load_pcd():
    converter.load_points("filetest.pcd")
    assert np.all(converter.points == ref)

