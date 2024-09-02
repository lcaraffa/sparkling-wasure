import argparse
import tempfile
import time
from math import ceil

import numpy as np
import pdal
import requests
from osgeo import gdal_array

import pdaltools.las_info as las_info
from pdaltools.unlock_file import copy_and_hack_decorator


def pretty_time_delta(seconds):
    sign_string = "-" if seconds < 0 else ""
    seconds = abs(int(seconds))
    days, seconds = divmod(seconds, 86400)
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    if days > 0:
        return "%s%dd%dh%dm%ds" % (sign_string, days, hours, minutes, seconds)
    elif hours > 0:
        return "%s%dh%dm%ds" % (sign_string, hours, minutes, seconds)
    elif minutes > 0:
        return "%s%dm%ds" % (sign_string, minutes, seconds)
    else:
        return "%s%ds" % (sign_string, seconds)


def retry(times, delay, factor=2, debug=False):
    def decorator(func):
        def newfn(*args, **kwargs):
            attempt = 1
            new_delay = delay
            while attempt <= times:
                need_retry = False
                try:
                    return func(*args, **kwargs)
                except requests.exceptions.ConnectionError as err:
                    print("Connection Error:", err)
                    need_retry = True
                except requests.exceptions.HTTPError as err:
                    if "Server Error" in str(err):
                        print("HTTP Error:", err)
                        need_retry = True
                    else:
                        raise err
                if need_retry:
                    print(f"{attempt}/{times} Nouvel essai après une pause de {pretty_time_delta(new_delay)} .. ")
                    if not debug:
                        time.sleep(new_delay)
                    new_delay = new_delay * factor
                    attempt += 1

            return func(*args, **kwargs)

        return newfn

    return decorator


def is_image_white(filename: str):
    raster_array = gdal_array.LoadFile(filename)
    band_is_white = [np.all(band == 255) for band in raster_array]
    return np.all(band_is_white)


def download_image_from_geoplateforme(
    proj, layer, minx, miny, maxx, maxy, pixel_per_meter, outfile, timeout, check_images
):
    # Give single-point clouds a width/height of at least one pixel to have valid BBOX and SIZE
    if minx == maxx:
        maxx = minx + 1 / pixel_per_meter
    if miny == maxy:
        maxy = miny + 1 / pixel_per_meter

    # for layer in layers:
    URL_GPP = "https://data.geopf.fr/wms-r/wms?"
    URL_FORMAT = "&EXCEPTIONS=text/xml&FORMAT=image/geotiff&SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&STYLES="
    URL_EPSG = "&CRS=EPSG:" + str(proj)
    URL_BBOX = "&BBOX=" + str(minx) + "," + str(miny) + "," + str(maxx) + "," + str(maxy)
    URL_SIZE = (
        "&WIDTH="
        + str(ceil((maxx - minx) * pixel_per_meter))
        + "&HEIGHT="
        + str(ceil((maxy - miny) * pixel_per_meter))
    )

    URL = URL_GPP + "LAYERS=" + layer + URL_FORMAT + URL_EPSG + URL_BBOX + URL_SIZE

    print(URL)
    if timeout < 10:
        print(f"Mode debug avec un timeout à {timeout} secondes")

    req = requests.get(URL, allow_redirects=True, timeout=timeout)
    req.raise_for_status()
    print(f"Ecriture du fichier: {outfile}")
    open(outfile, "wb").write(req.content)

    if check_images and is_image_white(outfile):
        raise ValueError(f"Downloaded image is white, with stream: {layer}")


@copy_and_hack_decorator
def color(
    input_file: str,
    output_file: str,
    proj="",
    pixel_per_meter=5,
    timeout_second=300,
    color_rvb_enabled=True,
    color_ir_enabled=True,
    veget_index_file="",
    check_images=False,
    stream_RGB="ORTHOIMAGERY.ORTHOPHOTOS",
    stream_IRC="ORTHOIMAGERY.ORTHOPHOTOS.IRC",
):
    metadata = las_info.las_info_metadata(input_file)
    minx, maxx, miny, maxy = las_info.get_bounds_from_header_info(metadata)

    if proj == "":
        proj = las_info.get_epsg_from_header_info(metadata)

    pipeline = pdal.Reader.las(filename=input_file)

    writer_extra_dims = "all"

    # apply decorator to retry 3 times, and wait 30 seconds each times
    download_image_from_geoplateforme_retrying = retry(7, 15, 2)(download_image_from_geoplateforme)

    if veget_index_file and veget_index_file != "":
        print(f"Remplissage du champ Deviation à partir du fichier {veget_index_file}")
        pipeline |= pdal.Filter.colorization(raster=veget_index_file, dimensions="Deviation:1:256.0")
        writer_extra_dims = ["Deviation=ushort"]

    tmp_ortho = None
    if color_rvb_enabled:
        tmp_ortho = tempfile.NamedTemporaryFile()
        download_image_from_geoplateforme_retrying(
            proj, stream_RGB, minx, miny, maxx, maxy, pixel_per_meter, tmp_ortho.name, timeout_second, check_images
        )

        pipeline |= pdal.Filter.colorization(
            raster=tmp_ortho.name, dimensions="Red:1:256.0, Green:2:256.0, Blue:3:256.0"
        )

    tmp_ortho_irc = None
    if color_ir_enabled:
        tmp_ortho_irc = tempfile.NamedTemporaryFile()
        download_image_from_geoplateforme_retrying(
            proj, stream_IRC, minx, miny, maxx, maxy, pixel_per_meter, tmp_ortho_irc.name, timeout_second, check_images
        )

        pipeline |= pdal.Filter.colorization(raster=tmp_ortho_irc.name, dimensions="Infrared:1:256.0")

    pipeline |= pdal.Writer.las(
        filename=output_file, extra_dims=writer_extra_dims, minor_version="4", dataformat_id="8", forward="all"
    )

    print("Traitement du nuage de point")
    pipeline.execute()

    # The orthoimages files will be deleted only when their reference are lost.
    # To keep them, make a copy (with e.g. shutil.copy(...))
    # See: https://docs.python.org/2/library/tempfile.html#tempfile.TemporaryFile
    return tmp_ortho, tmp_ortho_irc


def parse_args():
    parser = argparse.ArgumentParser("Colorize tool", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--input", "-i", type=str, required=True, help="Input file")
    parser.add_argument("--output", "-o", type=str, default="", help="Output file")
    parser.add_argument(
        "--proj", "-p", type=str, default="", help="Projection, default will use projection from metadata input"
    )
    parser.add_argument("--resolution", "-r", type=float, default=5, help="Resolution, in pixel per meter")
    parser.add_argument("--timeout", "-t", type=int, default=300, help="Timeout, in seconds")
    parser.add_argument("--rvb", action="store_true", help="Colorize RVB")
    parser.add_argument("--ir", action="store_true", help="Colorize IR")
    parser.add_argument(
        "--vegetation", type=str, default="", help="Vegetation file, value will be stored in Deviation field"
    )
    parser.add_argument("--check-images", "-c", action="store_true", help="Check that downloaded image is not white")
    parser.add_argument(
        "--stream-RGB",
        type=str,
        default="ORTHOIMAGERY.ORTHOPHOTOS",
        help="""WMS raster stream for RGB colorization:
default stream (ORTHOIMAGERY.ORTHOPHOTOS) let the server choose the resolution
for 20cm resolution rasters, use HR.ORTHOIMAGERY.ORTHOPHOTOS
for 50 cm resolution rasters, use ORTHOIMAGERY.ORTHOPHOTOS.BDORTHO""",
    )
    parser.add_argument(
        "--stream-IRC",
        type=str,
        default="ORTHOIMAGERY.ORTHOPHOTOS.IRC",
        help="""WMS raster stream for IRC colorization. Default to ORTHOIMAGERY.ORTHOPHOTOS.IRC
Documentation about possible stream : https://geoservices.ign.fr/services-web-experts-ortho""",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    color(
        input_file=args.input,
        output_file=args.output,
        proj=args.proj,
        pixel_per_meter=args.resolution,
        timeout_second=args.timeout,
        color_rvb_enabled=args.rvb,
        color_ir_enabled=args.ir,
        veget_index_file=args.vegetation,
        check_images=args.check_images,
        stream_RGB=args.stream_RGB,
        stream_IRC=args.stream_IRC,
    )
