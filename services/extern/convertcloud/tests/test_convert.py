import decimal
import numpy as np

import convertcloud

ref = np.array([[0.123, 4.567, 8.901],
                [2.345, 6.789, 0.123]], dtype=np.dtype(decimal.Decimal))

converter = convertcloud.Converter()
converter.points = ref

def test_convert_xyz(tmpdir):
    output = tmpdir.join("outputtest.xyz")
    converter.convert(str(output))
    with open("filetest.xyz", "r") as ref_file:
        ref_string = ref_file.read()
    assert output.read() == ref_string

def test_convert_a3d(tmpdir):
    output = tmpdir.join("outputtest.a3d")
    converter.convert(str(output))
    with open("filetest.a3d", "r") as ref_file:
        ref_string = ref_file.read()
    assert output.read() == ref_string

def test_convert_ply(tmpdir):
    output = tmpdir.join("outputtest.ply")
    converter.convert(str(output))
    with open("filetest.ply", "r") as ref_file:
        ref_string = ref_file.read()
    assert output.read() == ref_string

def test_convert_pcd(tmpdir):
    output = tmpdir.join("outputtest.pcd")
    converter.convert(str(output))
    with open("filetest.pcd", "r") as ref_file:
        ref_string = ref_file.read()
    assert output.read() == ref_string

def test_pcd2ply(tmpdir):
    converter.load_points("filetest.pcd")
    output = tmpdir.join("outputtest.ply")
    converter.convert(str(output))
    with open("filetest.ply", "r") as ref_file:
        ref_string = ref_file.read()
    assert output.read() == ref_string

def test_ply2pcd(tmpdir):
    converter.load_points("filetest.ply")
    output = tmpdir.join("outputtest.pcd")
    converter.convert(str(output))
    with open("filetest.pcd", "r") as ref_file:
        ref_string = ref_file.read()
    assert output.read() == ref_string

def test_xyz2ply(tmpdir):
    converter.load_points("filetest.xyz")
    output = tmpdir.join("outputtest.ply")
    converter.convert(str(output))
    with open("filetest.ply", "r") as ref_file:
        ref_string = ref_file.read()
    assert output.read() == ref_string

def test_xyz2pcd(tmpdir):
    converter.load_points("filetest.xyz")
    output = tmpdir.join("outputtest.pcd")
    converter.convert(str(output))
    with open("filetest.pcd", "r") as ref_file:
        ref_string = ref_file.read()
    assert output.read() == ref_string


