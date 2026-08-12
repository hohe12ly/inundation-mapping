# find fathom3m DEM tiles that cover a HUC8
# this script constructs a r-tree from the coverage shapes of fathom3m DEM 
# tiles; takes in a set of huc8 codes; for each huc8 code, find the polygon 
# and then find tiles covering that polygon; finally, it outputs tile file 
# names for each huc8.
import os, sys, time
import numpy as np
import pandas as pd
import shapely
import shapely.geometry
import shapely.errors
import fiona
import fiona.crs
import rasterio
import rasterio.mask
import rasterio.features
import rtree
import math
import csv

# construct r-tree for tiles
def construct_rtree(features, fnfield):
    def generate_items():
        i = 0
        for feature in features:
            box = shapely.geometry.shape(feature["geometry"]).bounds
            yield (i, box, (i, str(feature['properties'][fnfield]).lstrip(), feature))
            i += 1
    return rtree.index.Index(generate_items())

# find huc shapes from wbd
def huc_shapes(l_huc, wbd, wbdlyr, hucfield):
    h_huc = {} # hash table, indexed by huc code

    v_wbd = fiona.open(wbd, mode='r', layer=wbdlyr)
    count = 0
    num = 0
    for f in v_wbd:
        num += 1
        huc = f['properties'][hucfield]
        if not (huc is None or huc == "") and not huc in h_huc:
            h_huc[huc] = f
            count += 1
    print('WBD: loaded', count, 'shapes from ', num, 'features')

    o_hucs = []
    for huc in l_huc:
        if huc in h_huc:
            o_hucs.append((huc, h_huc[huc]))

    print('WBD: found shapes for', len(o_hucs), 'of', len(l_huc), 'HUC codes')
    v_wbd.close()
    return o_hucs

def find_tiles_by_huc(l_h, coverageshp, fnfield, bufsize):
    v_tilecov = fiona.open(coverageshp, mode='r')
    tileitree = construct_rtree(v_tilecov, fnfield)

    o = []
    for huc, feature_h in l_h:
        hucshape = shapely.geometry.shape(feature_h["geometry"])
        hucshape_bufferred = hucshape.buffer(bufsize)
        left, bottom, right, top = hucshape.bounds
        left -= bufsize
        bottom -= bufsize
        right += bufsize
        top += bufsize
        hucbbox = (left, bottom, right, top) # xmin, ymin, xmax, ymax
        hits = tileitree.intersection(hucbbox, objects='raw')
        tiles = []
        for index, fn, feature_tile in hits:
            tileshape = shapely.geometry.shape(feature_tile["geometry"])
            if tileshape.intersects(hucshape_bufferred) or tileshape.within(hucshape_bufferred): # or tileshape.contains(hucshape_bufferred) or tileshape.touches(hucshape_bufferred):
                tiles.append(fn)
        o.append((huc, tiles))
    v_tilecov.close()

    return o

# command:
# python find_tiles_by_huc.py /srv/test/fathom3m/fathom3m_coverage_tacc.shp filename /srv/data/wbd/WBD/WBD.gdb WBDHU8 HUC8 /tmp 0.1 "12010005 12020003 12020006 12020007 12030203 12040201 12040202"
if __name__ == '__main__':
    coverageshp = sys.argv[1] # DEM tile coverage shp
    fnfield = sys.argv[2] # file name field in coverage shp
    wbd = sys.argv[3] # wbd path
    wbdlyr = sys.argv[4] # wbd layer name, e.g., WBDHU8
    hucfield = sys.argv[5] # huc code field in wbd, e.g., 'HUC8'
    odir = sys.argv[6] # output dir
    bufsize = float(sys.argv[7]) # buffer size in the same projection
    huc8s = sys.argv[8] # a string of hucs "huccode1 huccode2..."

    l_h = huc_shapes(huc8s.split(' '), wbd, wbdlyr, hucfield)
    o = find_tiles_by_huc(l_h, coverageshp, fnfield, bufsize)
    for huc, tiles in o:
        print(huc, ':', ' '.join(tiles))
        of = odir + '/' + huc + '_tilelist.txt'
        with open(of, 'w') as f:
            for tile in tiles:
                f.write("%s\n" % tile)

