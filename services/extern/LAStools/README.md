_This is a fork of
[the official LAStools project](https://github.com/LAStools/LAStools)
that only keeps the library section of LAStools (that is to say LASlib
and LASzip). This fork is design to be lightweight and self-contained
for people who are only interested in the libraries LASlib and
LASzip._

_If you are interested in the full LAStools software, please
refer to [the official LAStools
project](https://github.com/LAStools/LAStools)._

# LASlib/LASzip

LASlib (with LASzip) is a C++ programming API for reading / writing
LIDAR data stored in standard LAS or in compressed LAZ format (1.0 -
1.3). Both libraries - LASlib with LASzip - are released together
under the terms of the GNU Lesser General Public Licence also known as
LGPL. See [LICENSE.txt](LICENSE.txt) for details.

Internally LASlib includes the stand-alone LASzip compression library that
provides lossless compression for LAS data. This allows you to read or
write *.laz (LASzip compressed *.las files) directly like you would an
uncompressed .las file using the same LASreader and LASwriter classes.

LASlib (with LASzip) is easy-to-use, light-weight, yet extremely fast.

# Compilation

This fork uses CMake to simplify the compilation process. To
compile LASlib/LASzip, use the following set of commands:

```sh
$ mkdir build
$ cd build
$ cmake ..
$ make
```

If you wish to install the library to your operating system's
dedicated directories, do:

```sh
# make install
```

# Example

Here is how simple the programming API is:

```c++
#include "lasreader.hpp"
#include "laswriter.hpp"

int main(int argc, char *argv[])
{
  LASreadOpener lasreadopener;
  lasreadopener.set_file_name("original.las");
  LASreader* lasreader = lasreadopener.open();

  LASwriteOpener laswriteopener;
  laswriteopener.set_file_name("compressed.laz");
  LASwriter* laswriter = laswriteopener.open(&lasreader->header);

  while (lasreader->read_point()) laswriter->write_point(&lasreader->point);

  laswriter->close();
  delete laswriter;

  lasreader->close();
  delete lasreader;

  return 0;
}
```

# Latest update

This fork was updated from [the original
repository](https://github.com/LAStools/LAStools) for the last
time on the 10th August 2017.


--------------------

(c) 2007-2015 martin.isenburg@rapidlasso.com
