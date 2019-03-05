
# coding: utf-8

# In[5]:


#from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import os
import time
import datetime as dt
import xarray as xr
from datetime import datetime
import pandas
import matplotlib.pyplot as plt
import numpy as np
import math
import cartopy.crs as ccrs
from scipy import stats

####################you will need to change some paths here!#####################
#list of input directories
dir_ccmp='F:/data/sat_data/ccmp/v02.0/Y'
dir_out = 'F:/data/sat_data/ccmp/pdf/'
##where to get the data through opendap, use these directories instead
#dir_cmc = 'https://podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/GDS2/L4/GLOB/CMC/CMC0.1deg/v3/'
#dir_flux = 'http://apdrc.soest.hawaii.edu:80/dods/public_data/WHOI_OAFlux/version3/daily/lh_oaflux/'
#the latest ccmp is from www.remss.com but they do not have an opendap server so you can use this instead:
#dir_ccmp='https://podaac-opendap.jpl.nasa.gov/opendap/allData/ccmp/L3.0/flk/'

#################################################################################
from math import sin, pi
from scipy import interpolate

#functions for running storm data
import sys


# In[6]:



def get_ccmp_month(lyr):
    dir_ccmp='F:/data/sat_data/ccmp/v02.0/monthly/'
    syr=str(lyr)
    fname_tem='/CCMP_Wind_Analysis_' + syr + '_V02.0_L3.0_RSS.nc'
    ccmp_filename = dir_ccmp + fname_tem      
    ds=xr.open_dataset(ccmp_filename,drop_variables=['nobs'])
    ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon').sortby('lat')
    ds['spd']=np.sqrt(ds.uwnd**2+ds.vwnd**2)
    return ds

def get_ccmp(lyr,idyjl):
    dir_ccmp='F:/data/sat_data/ccmp/v02.0/'
    d = dt.datetime(2010,1,1) + dt.timedelta(days=idyjl-1)
    day_of_year, imon, idym = d.timetuple().tm_yday, d.month, d.day
    syr, smon, sdym, sjdy=str(lyr),str(imon),str(idym),str(idyjl)
    fname_tem='/CCMP_Wind_Analysis_' + syr + smon.zfill(2) + sdym.zfill(2) + '_V02.0_L3.0_RSS.nc'
    ccmp_filename = dir_ccmp + '/Y' + syr + '/M' + smon.zfill(2) + fname_tem      
    ds=xr.open_dataset(ccmp_filename,drop_variables=['nobs'])
    ds = ds.rename({'longitude':'lon','latitude':'lat'}) #, inplace = True)  
    ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon').sortby('lat')
    ds['spd']=np.sqrt(ds.uwnd**2+ds.vwnd**2)
    return ds

def get_monthly_sst(lyr):
    dir_data_sst = 'F:/data/sst/cmc/CMC0.2deg/v2/monthly/' 
    filename = dir_data_sst + str(lyr) + 'monthly_average_' + '120000-CMC-L4_GHRSST-SSTfnd-CMC0.2deg-GLOB-v02.0-fv02.0.nc'
    print(filename)
    ds=xr.open_dataset(filename)
    ds.close()
    return ds

lyr=2015
ds = get_monthly_sst(lyr)
ds=ds.sel(time=slice('2015-01-01','2015-05-01'))
land_mask = ds.copy(deep=True)        
ds = get_ccmp_month(lyr)
ds=ds.sel(time=slice('2015-01-01','2015-05-01'))
land_mask2 = land_mask.interp_like(ds,method='nearest')

lons=ds.lon.data
lats=ds.lat.data
#create 2d grid from lats and lons
[lon2d,lat2d]=np.meshgrid(lons,lats)
lon3d = np.stack((lon2d,lon2d,lon2d,lon2d),axis=0)
lat3d = np.stack((lat2d,lat2d,lat2d,lat2d),axis=0)
lon3d.shape

init_data=0
for lyr in range(2000,2019):
    for idyjl in range(1,366):
        ds = get_ccmp(lyr,idyjl)
        land_mask2['time']=ds.time
        ds = ds.where(np.isfinite(land_mask2.mask))

        xdim,ydim,tdim = ds.lon.shape[0],ds.lat.shape[0],ds.time.shape[0]
        pdim=xdim*ydim*tdim

        cbin1 = np.arange(0.001, 30,.1)  #cold wake bins
        bins=cbin1
        data = ds.spd.data
        hist1,mids = np.histogram(data,bins)[0],0.5*(bins[1:]+bins[:-1])

        x1= np.reshape(lat3d.data,(pdim))
        x2= np.reshape(lon3d.data,(pdim))
        x3= np.reshape(ds.spd.data,(pdim))
        print(x1.shape,x2.shape,x3.shape)
        x=np.vstack((x1,x2,x3))
        b1= np.arange(-90,90,.25)
        b2= np.arange(-180,180,.5)
        b3= np.arange(0.001,72,.1)
        #print(b1.shape,b2.shape,b3.shape)
        dbins=np.vstack((b1,b2,b3)).T
        v = np.reshape(ds.spd.data,(pdim))
        hist2=stats.binned_statistic_dd(x.T,v,'count', bins=dbins.T)[0]
        #sum2=stats.binned_statistic_dd(x.T,v, 'sum', bins=dbins.T)[0]    

        if init_data == 0:
            sv_hist1 = hist1
            sv_hist2 = hist2
            init_data=1
        else:
            sv_hist1+= hist1
            sv_hist2+=hist2

    ds=xr.Dataset(data_vars={'hist1': (('spd'),sv_hist1),'hist2': (('lat','lon','spd2'),sv_hist2)},
                  coords={'spd':cbin1[0:-1],'lat':b1[0:-1],'lon':b2[0:-1],'spd2':b3[0:-1]})

    filename=dir_out + str(lyr)+'pdf.nc'
    ds.to_netcdf(filename)



