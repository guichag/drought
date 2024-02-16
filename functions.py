"""Functions"""

import xarray as xr
import dill
import numpy as np

from config import *
from utils import transform_latlon_to_lcc


### CST ###


### FUNCTIONS ###

def get_datename_yr_SAFR(year):
    """Get Safran data file dates"""
    datename_yr = str(year) + '080107' + '_' + str(year+1) + '080106'

    return(datename_yr)


def get_datafile_yr_SAFR(year):
    """Get Safran data file name"""
    datename = get_datename_yr_SAFR(year)
    datafile = DATAPATH + '/SAFRAN/ForcPRCP_france_SAFRAN_8Km_1hour_' + datename + '_V2.nc'

    return(datafile)


def get_data_yr_SAFR(year, xs=[None, None], ys=[None, None]):
    """Get year data"""
    datafile = get_datafile_yr_SAFR(year)
    ds = xr.open_dataset(datafile)
    if (xs[0] != None) and (xs[1] != None):
        assert xs[0] < xs[1], 'Incorrect x coordinates'

    if (ys[0] != None) and (ys[1] != None):
        assert ys[0] < ys[1], 'Incorrect y coordinates'

    ds = ds.sel(x=slice(xs[0], xs[1]), y=slice(ys[0], ys[1]))

    return(ds)


def get_data_all_SAFR(ymin, ymax, xs=[None, None], ys=[None, None]):
    """Get SAFRAN data for the period"""
    years = np.arange(ymin, ymax+1, 1)      # consider hydrological years (from 1/08 of year y to 31/07 of year y+1)
    data_all = [get_data_yr_SAFR(y, xs, ys) for y in years]
    data_all = xr.concat(data_all, 'time')
    #data_all = data_all.sel(time=data_all.time.dt.year.isin(range(ymin,ymax+1)))   # not needed for hydrological years
    data_all = data_all.resample(time='1D').mean()  # hourly mean -> daily mean precip. flux
    data_all = data_all * 86400  # kg/m²/s -> mm/d
    data_all = data_all.drop_sel(time=data_all.indexes['time'][-1])  # remove last element (1/08 of last year, averaged over 6 hours)

    return(data_all)


def get_data_yr_SAFR_mean(year, xs=[None, None], ys=[None, None]):
    """Get spatial mean of year data"""
    ds = get_data_yr_SAFR(year, xs, ys)
    ds = ds.mean(dim='x').mean(dim='y')  # spatial mean

    return(ds)


def get_data_all_SAFR_mean(ymin, ymax, xs=[None, None], ys=[None, None]):
    """Get spatial mean of SAFRAN data for the period"""
    years = np.arange(ymin, ymax+1, 1)      # consider hydrological years (from 1/08 of year y to 31/07 of year y+1)
    data_all = [get_data_yr_SAFR_mean(y, xs, ys) for y in years]
    data_all = xr.concat(data_all, 'time')
    data_all = data_all.resample(time='1D').mean()  # hourly mean -> daily mean precip. flux
    data_all = data_all * 86400  # kg/m²/s -> mm/d
    data_all = data_all.drop_sel(time=data_all.indexes['time'][-1])  # remove last element (1/08 of last year, averaged over 6 hours)

    return(data_all)


def load_cum_year(ymin, ymax, year, size=30):
    """Load yearly cumulative rainfall data"""
    per = str(ymin) + '-' + str(ymax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))
    datadir = DATADIR + '/Fr/SAFRAN/' + per + '/cum/' + size_
    datafile = datadir + '/' + str(year)

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return out


def load_pcum(ymin, ymax, size=30):
    """Load cumulative rainfall data"""
    per = str(ymin) + '-' + str(ymax)
    size_ = str(int(size*res/1000)) + 'x' + str(int(size*res/1000))
    datadir = DATADIR + '/Fr/SAFRAN/' + per + '/' + size_
    datafile = datadir + '/pcum'

    with open(datafile, 'rb') as pics:
        out = dill.load(pics)

    return out


