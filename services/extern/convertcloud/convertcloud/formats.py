import numpy as np
import io
import struct
class Field:
    def __init__(self, name):
        self.name = name
        self.size = None
        self.type = None

    def __str__(self):
        return "Field(name={}, size={}, type={})".format(self.name, self.size, self.type)
    __repr__ = __str__

class Load(object):

    def pcd(self, path):
        """
        Parse a file in .pcd format

        Args:
          path (string): Path to file to be loaded
        Returns:
          points (numpy array): Loaded points
          fields (Field object): Pointcloud field information
        """
        fields = []
        points = []

        _ascii = False
        _binary = False
        _bcompressed = False

        nr_pts = 0
        with open(path, "rb") as f:
            while True:
                line = f.readline()
                if line.startswith(b"#"):
                    pass
                elif line.startswith(b"VERSION"):
                    pass
                elif line.startswith(b"FIELDS"):
                    line = line.strip()
                    line = line.split(b" ")
                    for idx, field in enumerate(line):
                        if idx != 0:
                            fields.append(Field(field))
                elif line.startswith(b"SIZE"):
                    line = line.strip()
                    line = line.split(b" ")
                    for idx, size in enumerate(line):
                        if idx != 0:
                            fields[idx-1].size = int(size)
                elif line.startswith(b"TYPE"):
                    line = line.strip()
                    line = line.split(b" ")
                    for idx, tmp in enumerate(line):
                        if idx != 0:
                            fields[idx-1].type = tmp 
                elif line.startswith(b"COUNT"):
                    pass
                elif line.startswith(b"WIDTH"):
                    pass
                elif line.startswith(b"HEIGHT"):
                    pass
                elif line.startswith(b"VIEWPOINT"):
                    pass
                elif line.startswith(b"POINTS"):
                    line = line.split(b" ")
                    nr_pts = int(line[1])
                elif line.startswith(b"DATA"):
                    line = line.strip()
                    line = line.split(b" ")
                    if line[1] == b"ascii":
                        _ascii = True
                    elif line[1] == b"binary":
                        _binary = True
                    elif line[1] == b"binary_compressed":
                        _bcompressed = True
                        print("Error: binary_compressed format not supported")
                        sys.exit(1)
                    else:
                        print("Error: unknown pcd file format")
                        print(line[1])
                        sys.exit(1)
                    break

            if _ascii:
                for line in f:
                    pt = line.split()
                    points.append(pt)

            if _binary:
                data = f.read()
            if _bcompressed:
                import lzf
                compressed_data = f.read()
                data = lzf.decompress(compressed_data, len(compressed_data))

            #print("Length of data: ", nr_pts)
            #print("Fields: ", fields)

        if _binary or _bcompressed:
            buf = io.BytesIO(data)
            fmt = ""
            size = 0
            for f in fields:
                if f.type == b"F" and f.size == 4:
                    fmt += "f" 
                elif f.type == b"F" and f.size == 8:
                    fmt += "d" 
                elif f.type == b"I" and f.size == 1:
                    fmt += "b" 
                elif f.type == b"I" and f.size == 2:
                    fmt += "h" 
                elif f.type == b"I" and f.size == 4:
                    fmt += "i" 
                elif f.type == b"U" and f.size == 1:
                    fmt += "B" 
                elif f.type == b"U" and f.size == 2:
                    fmt += "H" 
                elif f.type == b"U" and f.size == 4:
                    fmt += "I" 
                else:
                    print("Uknown type: ", f.type)
                size += f.size

            if len(fields) > 3 and fields[3].name == "rgb":
                self._rgb = True
            for _ in range(nr_pts):
                pt = struct.unpack(fmt, buf.read(size))
                break
                points.append(pt)

        return np.array(points), fields

    def ply(self, path):
        """
        Parse a file in .ply format

        Args:
          path (string): Path to file to be loaded
        Returns:
          points (numpy array): Loaded points
          fields (Field object): Pointcloud field information
        """
        points = []
        fields = []

        _ascii = False
        _binary = False

        nr_pts = 0
        with open(path, "rb") as f:
            while True:
                line = f.readline()
                if line.startswith(b"ply"):
                    pass
                elif line.startswith(b"format"):
                    line = line.strip()
                    line = line.split(b" ")
                    ft = line[1]

                    if ft == b"ascii":
                        _ascii = True
                    elif ft == b'binary_little_endian':
                        _binary = True
                        _endianchar = '<'
                    elif ft == b'binary_big_endian':
                        _binary = True
                        _endianchar = '>'

                elif line.startswith(b"comment"):
                    pass
                elif line.startswith(b"element"):
                    line = line.split(b" ")
                    if line[1] == b"vertex":
                        nr_pts = int(line[-1])
                elif line.startswith(b"property"):
                    line = line.strip()
                    line = line.split(b" ")

                    fields.append(Field(line[2].decode()))
                    fields[-1].type = line[1].decode()
                    fields[-1].size = 4

                elif line.startswith(b"end_header"):
                    pass
                    break

            if _ascii:
                for line in f:
                    pt = line.split()
                    points.append(pt)

                    # Do not append faces
                    if len(points) == nr_pts:
                        break

            if _binary:
                data = f.read()

        if _binary:
            buf = io.BytesIO(data)
            fmt = _endianchar
            size = 0
            for f in fields:
                if f.type == b"float" and f.size == 4:
                    fmt += "f"
                else:
                    print("Uknown type: ", f.type)
                size += f.size
            for _ in range(nr_pts):
                pt = struct.unpack(fmt, buf.read(size))
                points.append(pt)

        return np.array(points), fields

    def zdf(self, path):
        """
        Parse a file in ZIVID .zdf format

        Args:
          path (string): Path to file to be loaded
        Returns:
          points (numpy array): Loaded points
          fields (Field object): Pointcloud field information
        """
        from netCDF4 import Dataset

        points = []
        fields = [] 

        f = Dataset(path,'r')
        xyz = f['data']['pointcloud'][:,:,:]
        img = f['data']['rgba_image'][:,:,:]
        f.close()

        pc = np.dstack([xyz, img])
        pc_reshaped = pc.reshape(pc.shape[0]*pc.shape[1], pc.shape[2])

        self._rgba = True

        for val in ["x", "y", "z", "rgb"]:
            field = Field(val)
            field.size = 4
            field.type = "float"
            fields.append(field)

        for pt in pc_reshaped:
            if not np.isnan(pt[0]):
                points.append(pt)
            else:
                points.append([0,0,0,0,0,0,255])

        return np.array(points), fields

    def xyz(self, path):
        """
        Parse a file in .xyz format

        Args:
          path (string): Path to file to be loaded
        Returns:
          points (numpy array): Loaded points
        """
        points = []
        fields = []

        with open(path, 'rb') as f:
            for line in f:
                xyz = line.split()
                if xyz[0] != b'nan':
                    points.append([float(val) for val in xyz])
                else:
                    points.append(len(xyz)*[0.0])

        for val in ["x", "y", "z"]:
            field = Field(val)
            field.size = 4
            field.type = "float"
            fields.append(field)

        if len(xyz) == 6:
            self._rgb = True
            field = Field("rgb")
            field.size = 4
            field.type = "float"
            fields.append(field)

        elif len(xyz) == 7:
            self._rgba = True

        return np.array(points), fields

    def npts(self, path):
        """
        Parse a file in .xyz format

        Args:
          path (string): Path to file to be loaded
        Returns:
          points (numpy array): Loaded points
        """
        points = []
        fields = []

        with open(path, 'rb') as f:
            for line in f:
                xyz = line.split()
                if xyz[0] != b'nan':
                    points.append([float(val) for val in xyz])
                else:
                    points.append(len(xyz)*[0.0])

        for val in ["x", "y", "z"]:
            field = Field(val)
            field.size = 4
            field.type = "float"
            fields.append(field)

        if len(xyz) == 6:
            self._norm = True
            field = Field("norm")
            field.size = 4
            field.type = "float"
            fields.append(field)

        return np.array(points), fields

    def stl(self, path):
        """
        Parse a file in .stl format

        Args:
          path (string): Path to file to be loaded
        Returns:
          points (numpy array): Loaded points
          fields (Field object): Pointcloud field information
        """
        from stl import mesh

        points = []
        fields = []

        stlmesh = mesh.Mesh.from_file(path)
        vects = stlmesh.data["vectors"]
        points = vects.reshape(vects.shape[0]*vects.shape[1], vects.shape[2])

        for val in ["x", "y", "z"]:
            field = Field(val)
            field.size = 4
            field.type = "float"
            fields.append(field)

        return np.array(points), fields

    def a3d(self, path):
        """
        Parse a file in Tordivel Scorpion .a3d format

        Args:
          path (string): Path to file to be loaded
        Returns:
          points (numpy array): Loaded points
        """
        import ast

        points = []
        fields = []

        with open(path, 'r') as f:
            for line in f:
                line = line.replace("\n", "")
                xyzraw = line.split(",")
                xyz = "[{}.{}, {}.{}, {}.{}]".format(xyzraw[0], xyzraw[1], \
                                                     xyzraw[2], xyzraw[3], \
                                                     xyzraw[4], xyzraw[5])
                points.append(ast.literal_eval(xyz))

        for val in ["x", "y", "z"]:
            field = Field(val)
            field.size = 4
            field.type = "float"
            fields.append(field)

        return np.array(points), fields


class Header(object):

    def __init__(self, nr_pts, fields, rgb, rgba,norm,origin):
        """
        Generates header of file to where pointcloud will be saved

        Args:
          nr_pts (int): Number of points of pointcloud
          fields (Fields object): Field information relevant for header
          rgb (bool): True if pointcloud contains RGB information
          rgba (bool): True if pointcloud contains RGBA information
        """

        self._nr_pts = nr_pts
        self._fields = fields
        self._rgb = rgb
        self._rgba = rgba
        self._norm = norm
        self._origin = origin

    def ply(self):

        properties = "property float x\n" \
                   + "property float y\n" \
                   + "property float z\n"

        if self._norm:
            properties += "property float nx\n" \
                        + "property float ny\n" \
                        + "property float nz\n"
        if self._origin:
            properties += "property float x_origin\n" \
                        + "property float y_origin\n" \
                        + "property float z_origin\n"
        if self._rgb:
            properties += "property uchar red\n" \
                        + "property uchar green\n" \
                        + "property uchar blue\n"
        elif self._rgba:
            properties += "property uchar red\n" \
                        + "property uchar green\n" \
                        + "property uchar blue\n" \
                        + "property uchar alpha\n"

        header = 'ply\n' \
               + "format ascii 1.0\n" \
               + "comment https://github.com/SintefRaufossManufacturing/convertcloud\n" \
               + "element vertex {}\n".format(self._nr_pts) \
               + properties \
               + "end_header\n"

        return header

    def pcd(self):
        if self._fields != None:
            fields = ""
            size = ""
            typ = ""
            types = {"float":"F", "int":"I", "uint":"U"}

            for field in self._fields:
                fields += field.name + " "
                size += str(field.size) + " "
                typ += types[field.type] + " "

            # Remove last space
            fields = fields[:-1]
            size = size[:-1]
            typ = typ[:-1]

        else:
            fields = "x y z"
            size = "4 4 4"
            typ = "F F F"
        if self._rgb or self._rgba:
            #TODO calculate rgb value from three R G B values (bitshift)
            fields += " r g b"
            size += " 4 4 4"
            typ += " 4 4 4"
        elif self._rgba:
            self._fields += " r g b a"
            size += " 4 4 4 4"
            typ += " 4 4 4 4"

        header = "# .PCD v0.7 - PointCloud Data file format\n" \
               + "VERSION 0.7\n" \
               + "FIELDS {}\n".format(fields) \
               + "SIZE {}\n".format(size) \
               + "TYPE {}\n".format(typ) \
               + "WIDTH {}\n".format(self._nr_pts) \
               + "HEIGHT 1\n" \
               + "VIEWPOINT 0 0 0 1 0 0 0\n" \
               + "POINTS {}\n".format(self._nr_pts) \
               + "DATA ascii\n"

        return header
