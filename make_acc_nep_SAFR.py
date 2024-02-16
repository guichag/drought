"""Functions to make matrices of precipitation accumulation and of non-exceedance probabilities"""

import sys
import os
import dill
import numpy as np
import pandas as pd
import rpy2.robjects as ro

from config import *
from functions import *


### CST ###

r = ro.r
r['source']('kde1d.R')
make_kde1d = ro.globalenv['make_kde1d']    # Load R function for KDE fit
make_pkde1d = ro.globalenv['make_pkde1d']  # Load R function for KDE CDF

mm = np.arange(1, 12+1,1)


### FUNCTIONS ###

    ##################################
    #                                #
    #       LOAD STATIONS DATA       #
    #          FROM R CODE           #
    #                                #
    ##################################


def load_mat_stations(ymin, ymax, nbweeksmin, nbweeksmax):
    """Load matrix of stations precipitation accumulation"""
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    datapath = datadir_stations + '/mat'
    datafile = datapath + '/mat_' + ts

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)


def load_nep_stations(ymin, ymax, typ, nbweeksmin, nbweeksmax):
    """Load matrix of stations NEPs"""
    per = str(ymin) + '-' + str(ymax)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    datapath = DATADIR + '/Fr/stations/' + per + '/nep'
    datafile = datapath + '/nep_' + ts + '_' + typ

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)



    ##################################
    #                                #
    #    PRECIP. ACC. COMPUTATION    #
    #     ON SPATIAL SUB-DOMAINS     #
    #                                #
    ##################################


def make_mat_SAFR_sub(ymin, ymax, nbweeksmin, nbweeksmax, xs=[1116000, 1188000], ys=[1665000, 1737000]):
    """Get matrix of precipitation accumulation on sub-domain"""
    p = get_data_all_SAFR_mean(ymin, ymax, xs, ys)
    p = p.to_dataframe()['product']

    nbweeks = np.arange(nbweeksmin, nbweeksmax+1,1)
    mat = pd.DataFrame(columns=nbweeks)
    for nbweek in nbweeks:
        print(nbweek, end=' : ', flush=True)
        mat[nbweek] = p.rolling(nbweek*7, center=False).sum()
        
    return(mat)


def load_mat_SAFR_sub(reg, ymin, ymax, nbweeksmin, nbweeksmax, xs=[1116000, 1188000], ys=[1665000, 1737000], size=30):
    """Load matrix of SAFRAN precipitation accumulation for sub-domain"""
    per = str(ymin) + '-' + str(ymax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))   # sub-domain size (km x km)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    coords_min_ = str(xs[0]) + '-' + str(ys[0])
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/' + size_ + '/' + coords_min_ + '/mat'
    datafile = datapath + '/mat_' + ts

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)


def load_mat_SAFR(reg, ymin, ymax, nbweeksmin, nbweeksmax):
    """Load matrix of SAFRAN precipitation accumulation for France"""
    per = str(ymin) + '-' + str(ymax)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/mat'
    datafile = datapath + '/mat_' + ts

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)


    ##################################
    #                                #
    #        NEP COMPUTATION         #
    #       ON WHOLE Fr DOMAIN       #
    #                                #
    ##################################


def make_nep_SAFR(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax):
    """Get matrix of Non-Exceedance Probability for France"""
    print('Get matrix of precip. accumulation...')
    mat = load_mat_SAFR(reg, ymin, ymax, nbweeksmin, nbweeksmax)
    nbweeks = np.arange(nbweeksmin, nbweeksmax+1, 1)
    matout = pd.DataFrame(columns=nbweeks, index=mat.index)

    if typ == 'kernel':  # Without deseasonalizing: KDE on each month (Jan, Feb, ...) separately
        print('Computing NEP (type={0})...'.format(typ))
        for m in mm:
            mat_mm = mat[mat.index.month == m]
            print(m)

            for nbweek in nbweeks:
                mat_mm_ts = mat_mm[nbweek]
                values = ro.FloatVector(mat_mm_ts.values.astype(float))
                tmp = make_kde1d(values)
                matout[nbweek][idx_mm_ts] = make_pkde1d(values, tmp)

                
    elif typ == 'kernall':  # With deseasonalizing: KDE on all data with monthly mean removed
        mat2 = pd.DataFrame(columns=nbweeks, index=mat.index)  # to be refilled with deseasonalized values
        print('Deseasonalizing...')
        for m in mm:
            mat_mm = mat[mat.index.month == m]
            idx_mm_ts = mat_mm.index
            print(m, len(idx_mm_ts))

            for nbweek in nbweeks:
                mat2[nbweek][idx_mm_ts] = mat_mm[nbweek] - mat_mm[nbweek].mean(axis=0)

        print('Computing NEP (type={0})...'.format(typ))
        for nbweek in nbweeks:
            print(nbweek, end=' : ', flush=True)
            mat2_ts = mat2[nbweek]
            values = ro.FloatVector(mat2_ts.values.astype(float))
            tmp = make_kde1d(values)
            matout[nbweek] = make_pkde1d(values, tmp)
            
    elif typ == 'empirical':
        print('Not yet implemented')
    
    else:
        print("Incorrect type")
    
    return(matout)


def load_nep_SAFR(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax):
    """Load matrix of SAFRAN NEPs"""
    per = str(ymin) + '-' + str(ymax)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/nep'
    datafile = datapath + '/nep_' + ts + '_' + typ

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)



    ##################################
    #                                #
    #         NEP COMPUTATION        #
    #     ON SPATIAL SUB-DOMAINS     #
    #                                #
    ##################################


#~ Method 1: compute NEPs for each sub-domain and take min across all sub-domains

def make_nep_SAFR_sub(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax, xs=[1116000, 1188000], ys=[1665000, 1737000], size=30):
    """Get matrix of Non-Exceedance Probability on sub-domains with method 1"""
    print('Get matrix of precip. accumulation...')
    mat = load_mat_SAFR_sub(reg, ymin, ymax, nbweeksmin, nbweeksmax, xs, ys, size)
    nbweeks = np.arange(nbweeksmin, nbweeksmax+1, 1)
    matout = pd.DataFrame(columns=nbweeks, index=mat.index)

    if typ == 'kernel':  # Without deseasonalizing: KDE on each month (Jan, Feb, ...) separately
        print('Computing NEP (type={0})...'.format(typ))
        for m in mm:
            mat_mm = mat[mat.index.month == m]
            print(m)

            for nbweek in nbweeks:
                mat_mm_ts = mat_mm[nbweek]
                values = ro.FloatVector(mat_mm_ts.values.astype(float))
                tmp = make_kde1d(values)
                matout[nbweek][idx_mm_ts] = make_pkde1d(values, tmp)

    elif typ == 'kernall':  # With deseasonalizing: KDE on all data with monthly mean removed
        mat2 = pd.DataFrame(columns=nbweeks, index=mat.index)  # to be refilled with deseasonalized values
        print('Deseasonalizing...')
        for m in mm:
            mat_mm = mat[mat.index.month == m]
            idx_mm_ts = mat_mm.index
            print(m, len(idx_mm_ts))

            for nbweek in nbweeks:
                mat2[nbweek][idx_mm_ts] = mat_mm[nbweek] - mat_mm[nbweek].mean(axis=0)

        print('Computing NEP (type={0})...'.format(typ))
        for nbweek in nbweeks:
            print(nbweek, end=' : ', flush=True)
            mat2_ts = mat2[nbweek]
            values = ro.FloatVector(mat2_ts.values.astype(float))
            tmp = make_kde1d(values)
            matout[nbweek] = make_pkde1d(values, tmp)
            
    elif typ == 'empirical':
        print('Not yet implemented')
    
    else:
        print("Incorrect type")

    return(matout)


def load_nep_SAFR_sub(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax, xs=[1116000, 1188000], ys=[1665000, 1737000], size=10):
    """Load matrix of SAFRAN NEPs"""
    per = str(ymin) + '-' + str(ymax)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))   # sub-domain size (km x km)
    coords_min_ = str(xs[0]) + '-' + str(ys[0])
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/' + size_ + '/' + coords_min_ + '/nep'
    datafile = datapath + '/nep_' + ts + '_' + typ

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)


def load_nep_min_SAFR(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax, size=30):
    """Load matrix of SAFRAN minimum NEPs (used for rarity matrix)"""
    per = str(ymin) + '-' + str(ymax)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))   # sub-domain size (km x km)
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/nep_multi-domain'
    datafile = datapath + '/' + ts + '_' + size_

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)



#~ Method 2: take min precip. acc. across all sub-domains and fit 1 kernel / timescale

def make_nep_SAFR_sub2(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax, size=30):
    """Get matrix of Non-Exceedance Probability on sub-domain with method 2"""
    nbweeks = np.arange(nbweeksmin, nbweeksmax+1, 1)
    per = str(ymin) + '-' + str(ymax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))   # sub-domain size (km x km)
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/' + size_
    coords_min = os.listdir(datapath)[:-2]

    if size <= 10:
        fac = 10
    elif 10 < size <= 20:
        fac = 2
    elif size > 20:
        fac = 1 

    mat_all = []
    print('Take min across {0} sub-domains'.format(len(coords_min[::fac])))

    for coord_min in coords_min[::fac]:  # sub-sample across sub-domains to avoid oversized mat_all -> AVOID MEMORY CRASH
        print(coord_min, end=" : ", flush=True)
        xmin_ = int(coord_min[:coord_min.find('-')])
        ymin_ = int(coord_min[coord_min.find('-')+1:])
        xmax_ = xmin_ + size * res - res
        ymax_ = ymin_ + size * res - res
        xs = [xmin_, xmax_]
        ys = [ymin_, ymax_]
        mat = load_mat_SAFR_sub(reg, ymin, ymax, nbweeksmin, nbweeksmax, xs, ys, size)
        mat_all.append(mat.to_numpy())

    mat_all = np.array(mat_all)
    mat_min = np.nanmin(mat_all, axis=0)  # min precip. acc. across sub-domains for each day and timescale
    loc_min = np.argmin(mat_all, axis=0)  # -> dates x nb of weeks with indices of sub-domain with min
    mat_min = pd.DataFrame(mat_min, columns=nbweeks, index=mat.index)

    matout = pd.DataFrame(columns=nbweeks, index=mat.index)

    if typ == 'kernel':  # Without deseasonalizing: KDE on each month (Jan, Feb, ...) separately
        print('Computing NEP (type={0})...'.format(typ))
        for m in mm:
            mat_mm = mat_min[mat_min.index.month == m]
            print(m)

            for nbweek in nbweeks:
                mat_mm_ts = mat_mm[nbweek]
                values = ro.FloatVector(mat_mm_ts.values.astype(float))
                tmp = make_kde1d(values)
                matout[nbweek][idx_mm_ts] = make_pkde1d(values, tmp)

    elif typ == 'kernall':  # With deseasonalizing: KDE on all data with monthly mean removed
        mat2 = pd.DataFrame(columns=nbweeks, index=mat_min.index)  # to be refilled with deseasonalized values
        print('Deseasonalizing...')
        for m in mm:
            mat_mm = mat_min[mat_min.index.month == m]
            idx_mm_ts = mat_mm.index
            print(m, len(idx_mm_ts))

            for nbweek in nbweeks:
                mat2[nbweek][idx_mm_ts] = mat_mm[nbweek] - mat_mm[nbweek].mean(axis=0)

        print('Computing NEP (type={0})...'.format(typ))
        for nbweek in nbweeks:
            print(nbweek, end=' : ', flush=True)
            mat2_ts = mat2[nbweek]
            values = ro.FloatVector(mat2_ts.values.astype(float))
            tmp = make_kde1d(values)
            matout[nbweek] = make_pkde1d(values, tmp)
            
    elif typ == 'empirical':
        print('Not yet implemented')
    
    else:
        print("Incorrect type")

    return(matout, loc_min)


def load_nep_SAFR_sub2(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax, size=30):
    """Load matrix of SAFRAN NEPs (method 2)"""
    per = str(ymin) + '-' + str(ymax)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))   # sub-domain size (km x km)
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/' + size_ + '/method2/nep'
    datafile = datapath + '/nep_' + ts + '_' + typ

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)


def load_nep_min_SAFR2(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax, size=30):
    """Load matrix of SAFRAN minimum NEPs (used for rarity matrix)"""
    per = str(ymin) + '-' + str(ymax)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))   # sub-domain size (km x km)
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/nep_multi-domain_meth2'
    datafile = datapath + '/' + ts + '_' + size_

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)



#~ Method 3: fit 1 kernel / timescale for all sub-domains and take the min NEP

def make_nep_SAFR_sub3(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax, size=30):
    """Get matrix of Non-Exceedance Probability on sub-domain with method 3"""
    nbweeks = np.arange(nbweeksmin, nbweeksmax+1, 1)
    per = str(ymin) + '-' + str(ymax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))   # sub-domain size (km x km)
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/' + size_
    coords_min = os.listdir(datapath)[:-2]

    if size <= 10:
        fac = 20
    elif 10 < size <= 20:
        fac = 5
    elif 20 < size <= 30:
        fac = 3
    elif size > 30:
        fac = 1 

    mat_all = []
    print('Get global distribution from {0} sub-domains'.format(len(coords_min[::fac])))

    for coord_min in coords_min[::fac]:  # sub-sample across sub-domains to avoid oversized mat_all -> AVOID MEMORY CRASH
        print(coord_min, end=" : ", flush=True)
        xmin_ = int(coord_min[:coord_min.find('-')])
        ymin_ = int(coord_min[coord_min.find('-')+1:])
        xmax_ = xmin_ + size * res - res
        ymax_ = ymin_ + size * res - res
        xs = [xmin_, xmax_]
        ys = [ymin_, ymax_]
        mat = load_mat_SAFR_sub(reg, ymin, ymax, nbweeksmin, nbweeksmax, xs, ys, size)
        mat_all.append(mat)

    mat_all = pd.concat(mat_all, axis=0)

    matout = pd.DataFrame(columns=nbweeks, index=mat_all.index)

    if typ == 'kernel':  # Without deseasonalizing: KDE on each month (Jan, Feb, ...) separately
        print('Computing NEP (type={0})...'.format(typ))
        for m in mm:
            mat_mm = mat_all[mat_all.index.month == m]
            print(m)

            for nbweek in nbweeks:
                mat_mm_ts = mat_mm[nbweek]
                values = ro.FloatVector(mat_mm_ts.values.astype(float))
                tmp = make_kde1d(values)
                matout[nbweek][idx_mm_ts] = make_pkde1d(values, tmp)

    elif typ == 'kernall':  # With deseasonalizing: KDE on all data with monthly mean removed
        mat2 = pd.DataFrame(columns=nbweeks, index=mat_all.index)  # to be refilled with deseasonalized values
        print('Deseasonalizing...')
        for m in mm:
            mat_mm = mat_all[mat_all.index.month == m]
            dates_mm_ts = mat_mm.index
            idx_mm_ts = mat_all.index.get_indexer_for(dates_mm_ts)
            idx_mm_ts = np.unique(idx_mm_ts)
            print(m, len(idx_mm_ts))

            for nbweek in nbweeks:
                mat2[nbweek].iloc[idx_mm_ts] = mat_mm[nbweek] - mat_mm[nbweek].mean(axis=0)

        print('Computing NEP (type={0})...'.format(typ))
        for nbweek in nbweeks:
            print(nbweek, end=' : ', flush=True)
            mat2_ts = mat2[nbweek]
            values = ro.FloatVector(mat2_ts.values.astype(float))
            tmp = make_kde1d(values)
            matout[nbweek] = make_pkde1d(values, tmp)

    elif typ == 'empirical':
        print('Not yet implemented')
    
    else:
        print("Incorrect type")

    return(matout)


def load_nep_SAFR_sub3(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax, size=30):
    """Load matrix of SAFRAN NEPs (method 2)"""
    per = str(ymin) + '-' + str(ymax)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))   # sub-domain size (km x km)
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/' + size_ + '/method3/nep'
    datafile = datapath + '/nep_' + ts + '_' + typ

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)


def load_nep_min_SAFR3(reg, ymin, ymax, typ, nbweeksmin, nbweeksmax, size=30):
    """Load matrix of SAFRAN minimum NEPs (used for rarity matrix)"""
    per = str(ymin) + '-' + str(ymax)
    ts = str(nbweeksmin) + '-' + str(nbweeksmax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))   # sub-domain size (km x km)
    datapath = DATADIR + '/' + reg + '/SAFRAN/' + per + '/nep_multi-domain_meth3'
    datafile = datapath + '/' + ts + '_' + size_

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return(out)

