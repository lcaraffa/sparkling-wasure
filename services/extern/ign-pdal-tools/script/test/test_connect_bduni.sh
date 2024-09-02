# Aller dans l'image lidar_prod
# sudo docker run -it --network host --rm  lidar_hd/lidar_prod bash

# installer postgresql-client
# sudo apt install postgresql-client

# Executer ce code
pgsql2shp -f 923000_6308000.shp -h serveurbdudiff.ign.fr -u invite -P 28de# bduni_france_consultation 'SELECT st_setsrid(batiment.geometrie,2154) AS geometry, 1 as presence FROM batiment WHERE batiment.geometrie && ST_MakeEnvelope(922950, 6306950, 924050, 6308050, 2154) and not gcms_detruit'