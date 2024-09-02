mkdir tmp

# test flux sur la reunion, flux 20cm:  HR.ORTHOIMAGERY.ORTHOPHOTOS => ok
curl -o tmp/reunion1_20cm_ok.tif "https://data.geopf.fr/wms-r/wms?LAYERS=HR.ORTHOIMAGERY.ORTHOPHOTOS&EXCEPTIONS=text/xml&FORMAT=image/geotiff&SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&STYLES=&CRS=EPSG:2975&BBOX=377000,7654000,377100,7654100&WIDTH=500&HEIGHT=500"

# test flux sur la reunion, flux 50cm: ORTHOIMAGERY.ORTHOPHOTOS.BDORTHO => ok
curl -o tmp/reunion1_50cm_ok.tif "https://data.geopf.fr/wms-r/wms?LAYERS=ORTHOIMAGERY.ORTHOPHOTOS.BDORTHO&EXCEPTIONS=text/xml&FORMAT=image/geotiff&SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&STYLES=&CRS=EPSG:2975&BBOX=377000,7654000,377100,7654100&WIDTH=500&HEIGHT=500"

# test flux sur la reunion, flux qui choisit: ORTHOIMAGERY.ORTHOPHOTOS => timeout
curl -o tmp/reunion1_choix_timeout.txt "https://data.geopf.fr/wms-r/wms?LAYERS=ORTHOIMAGERY.ORTHOPHOTOS&EXCEPTIONS=text/xml&FORMAT=image/geotiff&SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&STYLES=&CRS=EPSG:2975&BBOX=377000,7654000,377100,7654100&WIDTH=500&HEIGHT=500"

# test flux sur la reunion, flux qui choisit: ORTHOIMAGERY.ORTHOPHOTOS => blanc
curl -o tmp/reunion2_choix_blanc.tif "https://data.geopf.fr/wms-r/wms?LAYERS=ORTHOIMAGERY.ORTHOPHOTOS&EXCEPTIONS=text/xml&FORMAT=image/geotiff&SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&STYLES=&CRS=EPSG:2975&BBOX=377000,7655950,377050,7655999.99&WIDTH=250&HEIGHT=250"

# test flux sur la reunion, flux 20cm: HR.ORTHOIMAGERY.ORTHOPHOTOS => ok
curl -o tmp/reunion2_20cm_ok.tif "https://data.geopf.fr/wms-r/wms?LAYERS=HR.ORTHOIMAGERY.ORTHOPHOTOS&EXCEPTIONS=text/xml&FORMAT=image/geotiff&SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&STYLES=&CRS=EPSG:2975&BBOX=377000,7655950,377050,7655999.99&WIDTH=250&HEIGHT=250"

# test flux sur la reunion, flux 50cm: ORTHOIMAGERY.ORTHOPHOTOS.BDORTHO => ok
curl -o tmp/reunion2_50cm_ok.tif "https://data.geopf.fr/wms-r/wms?LAYERS=ORTHOIMAGERY.ORTHOPHOTOS.BDORTHO&EXCEPTIONS=text/xml&FORMAT=image/geotiff&SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&STYLES=&CRS=EPSG:2975&BBOX=377000,7655950,377050,7655999.99&WIDTH=250&HEIGHT=250"