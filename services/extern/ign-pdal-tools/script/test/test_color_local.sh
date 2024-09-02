# TODO: ajouter le fichier LAS en local
python -u -m pdaltools.color -i ../lidarExpress/data/one_micro_laz/Semis_2021_0785_6378_LA93_IGN69_light.laz -o ./tmp/output.las \
-r 0.1 -t 15 --rvb --ir
