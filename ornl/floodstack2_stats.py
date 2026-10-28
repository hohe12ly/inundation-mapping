# use fiona and shapely to:
# 1. open a geopackage MultiPolygon layer
# 2. count the number of MultiPolygons, polygons, and the total number of coordinates

import fiona
import shapely.geometry
import sys, os
if __name__ == "__main__":
    f_gis = sys.argv[1]
    layer = sys.argv[2]
    with fiona.open(f_gis, "r", layer=layer) as src:
        total_mp = len(src)
        total_polygons = 0
        total_coords = 0
        for feature in src:
            geom = shapely.geometry.shape(feature['geometry'])
            if isinstance(geom, shapely.geometry.MultiPolygon):
                total_polygons += len(geom.geoms)
                for poly in geom.geoms:
                    total_coords += len(poly.exterior.coords)
                    for interior in poly.interiors:
                        total_coords += len(interior.coords)
            else: # error: report and abort
                print(f'ERROR: feature {feature["id"]} is not a MultiPolygon, but a {type(geom)}')
                sys.exit(1)
        print(f'GPKG {f_gis} Layer {layer}: MultiPolygons: {total_mp} polygons: {total_polygons} coordinates: {total_coords}')