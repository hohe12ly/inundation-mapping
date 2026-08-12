import os, sys, time
import numpy as np

import pickle
import pandas as pd
import shapely
import shapely.geometry
import shapely.errors
import shapely.ops
import fiona
import fiona.crs

# method 1: 
# further process subsetted streams and catchment files of a HUC
# if a stream/catch's ID (COMID) indicates it is a minor stream,
# we remove it so that later the crosswalk step would not map
# hydroid to a minor stream COMID whose rating curve parameters
# are not representative of the main channel.
def remove_minor_streams_and_catch(f_lookup, f_fabric_streams, lyr_streams, f_fabric_catch, lyr_catch, of_fabric_streams, of_fabric_catch):
    with open(f_lookup, 'rb') as f:
        lookup = pickle.load(f)
        if len(lookup) == 0:
            print('ERROR remove_minor_streams_and_catch:', 'lookup table has not data -', f_lookup)
            return

        streams = fiona.open(f_fabric_streams, mode='r', layer=lyr_streams)
        n_streams = len(streams)
        o_streams = []
        for feat in streams:
            if not feat['properties']['ID'] in lookup: # keep if not a minor braided stream
                o_streams.append(feat)
        with fiona.open(of_fabric_streams,
                       mode = 'w',
                       layer=lyr_streams,
                       driver = streams.driver,
                       crs = streams.crs,
                       schema = streams.schema) as of:
            for feat in o_streams:
                of.write(feat)
        streams.close()
        print('remove_minor_streams_and_catch: read', n_streams, 'streams', ', found', n_streams - len(o_streams), 'minor braided strems', ', wrote', len(o_streams), 'streams into output', of_fabric_streams)

        catch = fiona.open(f_fabric_catch, mode='r', layer=lyr_catch)
        n_catch = len(catch)
        o_catch = []
        for feat in catch:
            if not feat['properties']['ID'] in lookup: # keep if not of a minor braided stream
                o_catch.append(feat)
        with fiona.open(of_fabric_catch,
                       mode = 'w',
                       layer=lyr_catch,
                       driver = catch.driver,
                       crs = catch.crs,
                       schema = catch.schema) as of:
            for feat in o_catch:
                of.write(feat)
        catch.close()
        print('remove_minor_streams_and_catch: read', n_catch, 'catchments', ', found', n_catch - len(o_catch), 'minor braided strems', ', wrote', len(o_catch), 'catchments into output', of_fabric_catch)

# method 2:
# instead of removing minor catchments, we merge minors' catchments 
# with the main channel's catchment. since cross-walrk does the majority
# zone stat, this will increase the change for the main channel to be
# mapped to nearby HydroIDs
def merge_minor_catch(f_lookup, f_fabric_catch, lyr_catch, of_fabric_catch):
    with open(f_lookup, 'rb') as f:
        lookup = pickle.load(f)
        if len(lookup) == 0:
            print('ERROR remove_minor_streams_and_catch:', 'lookup table has not data -', f_lookup)
            return

        # index catchment features using comids
        catch = fiona.open(f_fabric_catch, mode='r', layer=lyr_catch)
        n_catch = len(catch)
        h_flows = {}
        for feat in catch:
            comid = feat['properties']['ID']
            if not comid in h_flows:
                h_flows[comid] = feat
            else:
                print('ERROR', 'comid', comid, 'has two features in catchment data')

        # index minor braided streams using the main channel's comid
        main2minors = {}
        for minor, main in lookup.items():
            if not (minor in h_flows and main in h_flows):
                continue
            if not main in main2minors:
                main2minors[main] = [minor]
            else:
                main2minors[main].append(minor)

        def extract_polygons(geom):
            if geom.geom_type == 'Polygon':
                return [geom]
            elif geom.geom_type == 'MultiPolygon':
                return [ g for g in geom.geoms ]
            else:
                return []

        # foreach main braided stream, merge its catch with its minor streams' catch
        o_catch = []
        for comid, feat in h_flows.items(): # scan all catchments
            if comid in lookup: # a minor braided stream: will be merged with main, so skip
                continue
            elif comid in main2minors: # a main channel, merge with minors
                b_geoms = extract_polygons(shapely.geometry.shape(h_flows[comid]['geometry']))
                for minor in main2minors[comid]:
                    b_geoms += extract_polygons(shapely.geometry.shape(h_flows[minor]['geometry']))
                merged = shapely.geometry.MultiPolygon(b_geoms).buffer(0) # output either polygon or multipolygon
                if merged.geom_type == 'Polygon':
                    merged = shapely.geometry.MultiPolygon([ merged ])
                # create new feature
                new_feat = {
                    'geometry': shapely.geometry.mapping(merged),
                    'properties': feat['properties']
                }
                o_catch.append(new_feat)
            else: # irrelevant, keep
                o_catch.append(feat)

        # save merged catchments
        with fiona.open(of_fabric_catch,
                       mode = 'w',
                       layer=lyr_catch,
                       driver = catch.driver,
                       crs = catch.crs,
                       schema = catch.schema) as of:
            for feat in o_catch:
                of.write(feat)
        catch.close()
        print('merge_minor_catch: read', n_catch, 'catchments', ', found', n_catch - len(o_catch), 'minor braided strems', ', wrote', len(o_catch), 'catchments into output', of_fabric_catch)

if __name__ == '__main__':
    merge_minor_catch(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
