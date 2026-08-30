# !/usr/bin/bash
set -e
set -x
WORK_DIR=/home/christian/git/witx-heatmap/data/routing/
NEW_DIR=${WORK_DIR}/data_new
PROD_DIR=${WORK_DIR}/data

rm -rf ${NEW_DIR}
rm -rf ${PROD_DIR} 
mkdir -p ${NEW_DIR}
cd ${NEW_DIR}
ln -s /home/christian/git/witx-heatmap/osrm-backend/profiles/lib

wget http://download.geofabrik.de/europe/denmark-latest.osm.pbf
#wget http://download.geofabrik.de/europe/switzerland-latest.osm.pbf
#wget https://download.geofabrik.de/europe/austria-latest.osm.pbf
#wget http://download.geofabrik.de/europe/france-latest.osm.pbf
#wget http://download.geofabrik.de/europe/italy-latest.osm.pbf
#wget http://download.geofabrik.de/europe/norway-latest.osm.pbf -O norway-latest.osm.pbf

# Sweden
# wget http://download.geofabrik.de/europe/sweden-latest.osm.pbf

# Germany
#wget http://download.geofabrik.de/europe/germany-latest.osm.pbf

#wget http://download.geofabrik.de/europe/belgium-latest.osm.pbf

# Italy
#wget https://download.geofabrik.de/europe/italy/nord-est-latest.osm.pbf
#wget https://download.geofabrik.de/europe/italy/nord-ovest-latest.osm.pbf
#wget https://download.geofabrik.de/europe/italy/centro-latest.osm.pbf

# France
# wget https://download.geofabrik.de/europe/france/languedoc-roussillon-latest.osm.pbf
#wget https://download.geofabrik.de/europe/france/provence-alpes-cote-d-azur-latest.osm.pbf
#wget https://download.geofabrik.de/europe/france/rhone-alpes-latest.osm.pbf
#wget https://download.geofabrik.de/europe/france/alsace-latest.osm.pbf
#wget https://download.geofabrik.de/europe/france/lorraine-latest.osm.pbf
#wget https://download.geofabrik.de/europe/france/franche-comte-latest.osm.pbf

# Austria
#wget https://download.geofabrik.de/europe/austria-latest.osm.pbf

# Poland
#wget https://download.geofabrik.de/europe/poland-latest.osm.pbf

osmium merge *.osm.pbf -o merged.osm.pbf

find . -name "*.pbf" -not -name "merged.osm.pbf" -delete

cp ${WORK_DIR}/witx.lua ${NEW_DIR}
/home/christian/git/witx-heatmap/build/osrm-backend/osrm-extract -t6 -p ${NEW_DIR}/witx.lua ${NEW_DIR}/merged.osm.pbf 

/home/christian/git/witx-heatmap/build/osrm-backend/osrm-partition -t6 ${NEW_DIR}/merged.osm.pbf
/home/christian/git/witx-heatmap/build/osrm-backend/osrm-customize -t6 ${NEW_DIR}/merged.osm.pbf

#systemctl stop osrm-routing.service
mv ${NEW_DIR} ${PROD_DIR}
