# !/usr/bin/bash
WORK_DIR=/home/christian/git/witx-heatmap/data/routing/valhalla_data
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

rm -f *.osm.pbf
rm -rf valhalla_tiles*

wget http://download.geofabrik.de/europe/denmark-latest.osm.pbf

osmium merge *.osm.pbf -o merged.osm.pbf
find . -name "*.pbf" -not -name "merged.osm.pbf" -delete

/home/christian/git/witx-heatmap/valhalla/install/bin/valhalla_build_tiles -c valhalla.json merged.osm.pbf

find valhalla_tiles | sort -n | tar cf valhalla_tiles.tar --no-recursion -T -