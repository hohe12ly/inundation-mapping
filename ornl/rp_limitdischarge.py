# HUC8 RP and limiting discharge workflow
import dask
from dask.distributed import Client, progress, LocalCluster, performance_report
from dask_jobqueue import PBSCluster
import numcodecs
import numpy as np
import pathlib
import pandas as pd
from rechunker import rechunk
import shutil
import xarray as xr
import zarr as zr
import matplotlib.pyplot as plt
## uda.py: upper drainage area computation from NOAA FIM-based catchment dataset
# Author: Yan Y. Liu <yanliu@ornl.gov>
# Date: October 04, 2022
# This code derives each reach's upper drainage area using the catchment vector
# dataset produced by FIM's HAND process. We construct the river network graph
# using NetworkX and update UDA from headwater reaches and iteratively remove 
# head water reaches from the graph until each reach's UDA is finalized.
# The output is a CSV of columns: HydroID, feature_id (COMID), areasqkm, uda
import os, sys, time
import shapely
import shapely.geometry
import shapely.errors
import fiona
import fiona.crs
import rasterio
import rasterio.mask
import rasterio.features
import affine

from collections import OrderedDict # for output schema

import networkx as nx
import scipy
import scipy.io

import rtree
import math
import csv

cfs2cms = 0.028316846592
cms2cfs = 35.3147 # 1 cms = 35.3147 cfs
m2f = 3.280839895 # 1 meter = 3.280839895 feet
sqmi2sqkm = 2.58999 # 1 square mile 
sqkm2sqmi = 0.386102 # 1 square km

bankfull_threshold = 67 # 67% is bankfull

# func to be applied to each group groupby hydroid
def hydroTable_calc_rp_bf(g, hydroid2fandh, comid2l, l_sf, l_rp, l_cpe, l_epp, l_bankfull, hydroid2lh, lh_H, lh_rp, lh_cpe, lh_epp, lh_bankfull):
    hydroid = g['HydroID'].to_numpy()[0]
    if hydroid in hydroid2fandh:
        comid, huc = hydroid2fandh[hydroid]
        if comid in comid2l: # its comid is not in NOAA NWM Reanalysis
            data_sf = g['discharge_cms'].to_numpy() # for all 84 stages
            i_comid = comid2l[comid] # feature_id level interp
            g['RP_Q_years'] = np.interp(data_sf, l_sf[i_comid], l_rp[i_comid])
            g['CPE_Q'] = np.interp(data_sf, l_sf[i_comid], l_cpe[i_comid])
            g['EPE_Q'] = np.interp(data_sf, l_sf[i_comid], l_epp[i_comid])
            g['BF_Q_cms_67'] = np.full(data_sf.shape[0], l_bankfull[i_comid], dtype=float)

            data_h = g['stage'].to_numpy()
            i_hydroid = hydroid2lh[hydroid]
            g['RP_H_years'] = np.interp(data_h, lh_H[i_hydroid], lh_rp[i_hydroid])
            g['CPE_H'] = np.interp(data_h, lh_H[i_hydroid], lh_cpe[i_hydroid])
            g['EPE_H'] = np.interp(data_h, lh_H[i_hydroid], lh_epp[i_hydroid])
            g['BF_H_m_67'] = np.full(data_h.shape[0], lh_bankfull[i_hydroid], dtype=float)

            return g
        else:
            print("WARN comid", comid, 'not found')
            return g
    else:
        print('WARN hydroid', hydroid, 'not found')
        return g

# derive river network tree structure from NextDownID
# return the directed graph and head water reaches
# i_reachnodes: the indexes of reaches to considered as graph nodes
def graph_by_downstream_id(ii_reachnodes, hydroids, downstream_hydroids):
    h_hydroid = {}
    for ii in ii_reachnodes:
        h_hydroid[hydroids[ii]] = ii
    G = nx.DiGraph()
    for i in range(len(ii_reachnodes)):
        ii = ii_reachnodes[i] # the index of hydroid
        if not downstream_hydroids[ii] in h_hydroid: # departing reach, ignore
            continue
        G.add_edge(ii, h_hydroid[downstream_hydroids[ii]])
    ii_hw = []
    ii_next = []
    for node in G.nodes():
        if len(list(G.predecessors(node))) == 0: # no parent, head water reach
            ii_hw.append(node)
        else: # if not head water reach, it will be considered in the next round
            ii_next.append(node)
    return G, ii_hw, ii_next


def compute_cda(fn_catch, fn_uda, df_uda_table):
    # fn_catch = sys.argv[1] # e.g. '/Users/yxl/test/limiting-discharge/12090204/gw_catchments_reaches_filtered_addedAttributes_crosswalked.gpkg'
    # fn_uda = sys.argv[2] # e.g., '/Users/yxl/test/limiting-discharge/12090204/uda.csv'
    # fn_nhdplusmr_uda = '/home/yxl/test/nwm-reanalyais/pin2flood/limiting-discharge/uda_nwm.csv'

    vec_catch = fiona.open(fn_catch, 'r')
    hydroids = []
    comids = []
    from_nodes = []
    to_nodes = []
    downstream_hydroids = []
    areasqkms = []
    areas = []
    for feature in vec_catch:
        catcharea = shapely.geometry.shape(feature['geometry']).area
        areas.append(catcharea)
        f = feature['properties']
        hydroids.append(int(f['HydroID']))
        comids.append(int(f['feature_id']))
        from_nodes.append(int(f['From_Node']))
        to_nodes.append(int(f['To_Node']))
        downstream_hydroids.append(int(f['NextDownID']))
        areasqkms.append(float(f['areasqkm']))

    num_hydroids = len(hydroids)
    if num_hydroids != np.unique(np.array(hydroids, dtype=int)).size:
        print('ERROR: duplicate HydroIDs in input')
        # return
    print('read', num_hydroids, 'HydroIDs')

    areas = np.array(areas, dtype=float)
    areasqkms = np.array(areasqkms, dtype=float)
    areadiff = np.abs(areas - areasqkms * 1000 * 1000)
    print('diff min:', areadiff.min(), 'diff max:',areadiff.max(), 'diff avg:', areadiff.mean())
    thres1 = 1
    oob1 = areadiff > thres1
    print('num of diff >', thres1, 'meter^2:', oob1.sum())
    areadiff = None

    G, ii_hw, ii_next = graph_by_downstream_id(list(range(len(hydroids))), hydroids, downstream_hydroids)

    # compute external cda
    cda_init = [ (df_uda_table.loc[comid]['TotDASqKM'] - df_uda_table.loc[comid]['AreaSqKM']) if comid in df_uda_table.index and len(df_uda_table.loc[comid])>0 else 0.0 for comid in comids ] # in sqkm

    # now update uda
    uda = areas.copy() / 1000000.0 # sqkm conversion
    # update uda init for flowlines that flow into the huc, i.e., those in ii_hw
    count = 0
    for ii in ii_hw:
        if cda_init[ii] > 0.000001: # ext cda > 1 square meter
            uda[ii] += cda_init[ii]
            print('cda_init:', 'hydroid', hydroids[ii], 'comid', comids[ii], 'ext cda:', cda_init[ii], ', new cda:', uda[ii])
            count += 1
    print("cda_init:", count, "HydroIDs flow into the huc have non-zero cumulative drainage area")

    iterations = 0
    count_braided = 0
    while len(ii_next)>0: # repeat until we cover every hydroid reach
        #print('iteration', iterations, ':', len(ii_hw) + len(ii_next), 'reaches on graph,', len(ii_hw), 'headwater reaches,', len(ii_next), 'non-headwater reaches')
        for ii in ii_hw: # foreach headwater reach
            ii_downstream_reaches = G.successors(ii)
            count = 0
            for ii_downstream_reach in ii_downstream_reaches: # should have just one, but we consider braided cases
                uda[ii_downstream_reach] += uda[ii] # update downstrea reach's uda
                count += 0
            if count>1:
                print('braided from reach', hydroids[ii])
                count_braided += 1
        G, ii_hw, ii_next = graph_by_downstream_id(ii_next, hydroids, downstream_hydroids)

        iterations += 1
    print('uda computation done in', iterations, 'iterations. found', count_braided, 'braided cases')

    # output
    d_o = pd.DataFrame({
        'HydroID': hydroids,
        'feature_id': comids,
        'areasqkm': areas / 1000000.0,
        'uda': uda
    })
    d_o.to_csv(fn_uda, index=False)

    return d_o

def myrun(huc, rdir):
    
    f_annualmax = '/home/yxl/test/nwm-reanalyais/annual_max_sf/annual_max_streamflow_conus_1979-2020.nc'
    f_uda2discharge_table = '/home/yxl/test/nwm-reanalyais/pin2flood/limiting-discharge/peak-sf/env_curve.csv'
    f_uda_table = '/home/yxl/test/nwm-reanalyais/pin2flood/limiting-discharge/uda_nwm.csv'
    
    # March 31 exercise
    # huclist = ['12090201', '12090205', '12090206', '12100201', '12100203', '12090301', '12100202', '12090204']
    # April 21 exercise
    #huclist = ['12070204', '12070103', '12070101', '12030201', '12030108', '12030105']
    # May 11 exercise
    #huclist = ['12070101', '12070102', '12070201', '12070202', '12070203', '12070204', '12070205']
    #huclist = ['12070102', '12070201', '12070202', '12070203', '12070205', '12060203', '12060204']
    #huclist = ['12100201', '12100202', '12100203', '12100304', '12100302', '12110106', '12110107'] # July 01
    #huclist = ['12100101', '12100102', '12100202', '12100203'] # test
    #huclist = ['12100101', '12100102', '12100202', '12100203', '12100204', '12100303', '12100401', '12100402', '12100403', '12100404', '12100405', '12100406'] # August 11
    #huclist = ['12090204']
    #huc = sys.argv[1]
    huclist = [huc]
    
    #rdir = '/home/yxl/test/nwm-reanalyais/pin2flood/2022-July/' # rdir +'/'+ huc +'/hydroTable.csv'
    #rdir = '/home/yxl/test/nwm-reanalyais/pin2flood/limiting-discharge/' 
    #rdir = sys.argv[2] # this is the directory above huc/
    odir = rdir
    
   
    huc_hydrotables = {}
    huc_uniq_comids = {}
    huc_uniq_hydroids = {}
    hydroid2fandh = {} # map HydroID to feature_id and HUC
    for h in huclist:
        df_huc = pd.read_csv(rdir + '/' + str(h) + '/' + 'hydroTable.csv', dtype={'HydroID':np.int32, 'feature_id':np.int32, 'stage':np.float32, 'discharge_cms':np.float32, 'HUC':str, 'LakeID':np.int32})
        huc_hydrotables[h] = df_huc
        huc_uniq_comids[h] = df_huc['feature_id'].unique()
        huc_uniq_hydroids[h] = df_huc['HydroID'].unique()
        for i, row in df_huc.iterrows():
            if not row['HydroID'] in hydroid2fandh:
                hydroid2fandh[row['HydroID']] = (row['feature_id'], h)
    
    comids = np.unique(np.concatenate([fids for h, fids in huc_uniq_comids.items()], axis=0))
    print(len(comids), 'unique feature_ids from', len(huclist), 'HUCs')
    
    count = 0
    for h in huclist:
        huc_count = len(huc_uniq_comids[h])
        count += huc_count
        print('HUC', h, ':', huc_count)
    print(count-len(comids), 'feature_ids are shared among HUCs')
    
    hydroids = np.unique(np.concatenate([hids for h, hids in huc_uniq_hydroids.items()], axis=0))
    
    count = 0
    for h in huclist:
        huc_count = len(huc_uniq_hydroids[h])
        count += huc_count
        print('HUC', h, ':', huc_count)
    print(count-len(hydroids), 'HydroIDs are shared among HUCs')
    
    # set up stage-sf rating curve table for each hydroid
    interp_rc = {}
    for h in huclist:
        rc_by_hydroid = huc_hydrotables[h].groupby('HydroID')
        for hydroid, data in rc_by_hydroid:
            interp_rc[hydroid] = (data['discharge_cms'].to_numpy(), data['stage'].to_numpy())

    ## load NOAA NWM Reanalysis annaul max streamflow data
    annualmax = xr.open_dataset(f_annualmax).streamflow.sel(feature_id=comids)

    ## compute COMID level return period and bankfull data
    l_comid = []
    l_bf_comid = []
    l_sf = []
    l_rp = []
    l_cpe = []
    l_epp = []
    l_bankfull = []
    comid2l = {}
    ii = 0
    nodata_comids = []
    for comid in comids:
        fd = annualmax.sel(feature_id=comid).to_numpy()
        if (np.isnan(fd).sum() == len(fd)):
            nodata_comids.append(comid)
            continue # skip if couldn't find streamflow history
        ob_ranked = np.sort(fd)
        ob_max = ob_ranked[-1]
        maxsf_threshold = ob_max * 2 # twice as max observed in prev 42 years
        n = ob_ranked.shape[0]
        rank = np.array(range(1, n+1))
        T = (n + 1) / (n + 1 - rank)
        xp = ob_ranked
        yp = T # return period
        cpe = rank / (n+1) # cumulative probability estimate
        epp = (n + 1 - rank) / (n + 1) # exceedance probability estimate

        l_comid.append(np.full(xp.shape[0], comid, dtype=np.int32))
        l_sf.append(xp)   # sf sorted
        l_rp.append(yp)   # return period in years (max is 42)
        l_cpe.append(cpe) # cummulative probability estimate
        l_epp.append(epp) # exceed probability 

        l_bf_comid.append(comid)
        l_bankfull.append(np.percentile(xp, bankfull_threshold) )

        comid2l[comid] = ii # map comid to the index of lookup tables
        ii += 1

    print('There are', len(nodata_comids), 'COMIDs missing NOAA NWM Reanalysis data')
    print(nodata_comids)

    ## compute HydroID level return period and bankfull data
    lh_hydroid = []
    lh_bf_hydroid = []
    lh_sf = []
    lh_H = []
    lh_rp = []
    lh_cpe = []
    lh_epp = []
    lh_bankfull = []
    hydroid2lh = {}
    ii = 0
    for hydroid in hydroids:
        comid, h = hydroid2fandh[hydroid] # get its feature_id and huc
        sf_list, stage_list = interp_rc[hydroid] # rating curve
        Q_hist = annualmax.sel(feature_id=comid).to_numpy() # get 42yr annumal max points
        if (np.isnan(Q_hist).sum() == len(Q_hist)):
            continue # skip if couldn't find streamflow history
        H_hist = np.interp(Q_hist, sf_list, stage_list) # convert sf to water depth, in meters
        sorti = np.argsort(H_hist)
        ob_ranked = np.take(H_hist, sorti)
        ob_max = ob_ranked[-1]
        maxsf_threshold = ob_max * 2 # twice as max observed in prev 42 years
        n = ob_ranked.shape[0]
        rank = np.array(range(1, n+1))
        T = (n + 1) / (n + 1 - rank)
        xp = ob_ranked
        yp = T # return period
        cpe = rank / (n+1) # cumulative probability estimate
        epp = (n + 1 - rank) / (n + 1) # exceedance probability estimate

        lh_hydroid.append(np.full(xp.shape[0], hydroid, dtype=np.int32))
        lh_sf.append(np.take(Q_hist, sorti))
        lh_H.append(xp)   # H sorted
        lh_rp.append(yp)   # return period in years (max is 42) by H
        lh_cpe.append(cpe) # cummulative probability estimate by H
        lh_epp.append(epp) # exceed probability by H

        lh_bf_hydroid.append(hydroid)
        lh_bankfull.append(np.percentile(xp, bankfull_threshold) ) # bankfull by H

        hydroid2lh[hydroid] = ii # map comid to the index of lookup tables
        ii += 1

    ## update each HUC's hydroTable
    for h, df_rc in huc_hydrotables.items():
        print('Processing', h, '...')
        df_rc_update = df_rc.groupby('HydroID').apply(hydroTable_calc_rp_bf, 
                                                      hydroid2fandh=hydroid2fandh, 
                                                      comid2l=comid2l, 
                                                      l_sf=l_sf, 
                                                      l_rp=l_rp, 
                                                      l_cpe=l_cpe, 
                                                      l_epp=l_epp, 
                                                      l_bankfull=l_bankfull, 
                                                      hydroid2lh=hydroid2lh, 
                                                      lh_H=lh_H, 
                                                      lh_rp=lh_rp, 
                                                      lh_cpe=lh_cpe, 
                                                      lh_epp=lh_epp, 
                                                      lh_bankfull=lh_bankfull
                                                     ).reset_index(drop=True)
        print('Saving rp bf properties to file.')
        df_rc_update.to_csv(odir + '/' + h + '/' + 'hydroTable_rp_bf.csv', index=False)

    ## computing limiting discharge attributes using hydroid-based upper drainage area (UDA)
    
    # USGS uda-discharge lookup table
    print('Limiting discharge: loading USGS lookup table')
    df_uda2discharge_table = pd.read_csv(f_uda2discharge_table, dtype={'CDA': np.float64, 'PLAUSIBLE_EVELOPMENT': np.float64})
    x = df_uda2discharge_table['CDA'].to_numpy()
    x = x * sqmi2sqkm # converting to sqkm
    y = df_uda2discharge_table['PLAUSIBLE_EVELOPMENT'].to_numpy()
    y = y * cfs2cms # converting to cms
    
    # load comid-uda lookup tabe: COMID, REACHCODE, AreaSqKM, TotDASqKM, DivDASqKM
    print('Limiting discharge: loading NHDPlusV21 MR COMID-TotDASqKM table')
    d_nhdplusmr_cda = pd.read_csv(f_uda_table, index_col = 'COMID', dtype = {'COMID': int, 'REACHCODE': str})
    #d_nhdplusmr_cda_huc = d_nhdplusmr_cda.loc[ d_nhdplusmr_cda['REACHCODE'].str[0:(len(huc) - 4)] == huc[0:(len(huc) - 4)] ] # TODO: for now, only for HUC8
    d_nhdplusmr_cda_huc = d_nhdplusmr_cda # load all to include flowlines in neighboring HUCs

    for h in huclist:
        print('----limiting discharge on UDA for HUC', h, '----')
        df_uda = compute_cda(rdir + '/' + h + '/gw_catchments_reaches_filtered_addedAttributes_crosswalked.gpkg', rdir + '/' + h + '/uda.csv', d_nhdplusmr_cda_huc)
        df_uda.set_index('HydroID', inplace=True)
        f_hydro_table = odir + '/' + h + '/' + 'hydroTable_rp_bf.csv'
        df_hydro_table = pd.read_csv(f_hydro_table, dtype={'discharge_cms': np.float64})
        data_hydroid = df_hydro_table['HydroID'].to_numpy()
        data_discharge = df_hydro_table['discharge_cms'].to_numpy()
        data_uda = np.array([ df_uda.loc[hydroid]['uda'] for hydroid in data_hydroid ], dtype=np.float64) # in sqkm
        data_max_charge = np.interp(data_uda, x, y)
        data_limit_charge_ratio = data_discharge / data_max_charge * 100 # in percentage 

        num_columns = len(df_hydro_table.columns)
        df_hydro_table.insert(num_columns, 'cda_sqkm', data_uda)
        df_hydro_table.insert(num_columns + 1, 'cda_LimDisch_cms', data_max_charge)
        df_hydro_table.insert(num_columns + 2, 'cda_sqmi', data_uda * sqkm2sqmi)
        df_hydro_table.insert(num_columns + 3, 'cda_LimDisch_cfs', data_max_charge * cms2cfs)
        df_hydro_table.insert(num_columns + 4, 'cda_PctLimDisch', data_limit_charge_ratio)

        # rename HydroID to FATSGTID
        df_hydro_table.rename(columns={'HydroID': 'FATSGTID'}, inplace=True)

        df_hydro_table.to_csv(odir + '/' + h + '/' + 'hydroTable_rp_bf_lmtdischarge_cda.csv', index=False)
        print('computed limiting discharge attributes:', h)

# myrun('12090204', '/dev/shm/fim3')
if __name__ == '__main__':
    huc = sys.argv[1]
    rdir = sys.argv[2]
    myrun(huc, rdir)
