# -*- coding: utf-8 -*-
# maintener : MDupays
# version : v.1 06/12/2022
# Extract info from the tile
import json
import logging

import pdal


def las_crop(input_file: str, output_file: str, bounds, spatial_ref: str = "EPSG:2154"):
    """Crop filter removes points that fall inside a cropping bounding box (2D)
    Args:
        input_dir (str): input point cloud file
        output_dir (str): output point cloud file
        bounds : 2D bounding box to crop to : provided as ([xmin, xmax], [ymin, ymax])
    """
    # Parameters
    information = {
        "pipeline": [
            {"type": "readers.las", "filename": input_file, "override_srs": spatial_ref, "nosrs": True},
            {"type": "filters.crop", "bounds": str(bounds)},
            {"type": "writers.las", "a_srs": spatial_ref, "filename": output_file, "forward": "all"},
        ]
    }
    # Create json
    json_crop = json.dumps(information, sort_keys=True, indent=4)
    logging.info(json_crop)
    pipeline = pdal.Pipeline(json_crop)
    pipeline.execute()
