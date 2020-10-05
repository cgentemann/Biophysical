#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import warnings
# filter some warning messages
warnings.filterwarnings("ignore") 

####################you will need to change some paths here!#####################
#list of input files
filename_origin='f:/data/project_data/NASA_biophysical/collocated_data/NPPSD_GOA_allseabird_full_eddy_info.nc'
# output files
filename_origin_out='f:/data/project_data/NASA_biophysical/collocated_data/NPPSD_GOA_allseabird_full_eddy_info_envdata'
filename_origin_out_nc='f:/data/project_data/NASA_biophysical/collocated_data/NPPSD_GOA_allseabird_full_eddy_info_envdata.nc'


input_data = input("Enter data: aviso,wnd,sst")

ds_bird = xr.open_dataset(filename_origin)
ds_bird.close()
ds_bird['time']=ds_bird.time64
ds_bird['lon'] = (ds_bird['lon'] + 180) % 360 - 180

def get_data():
    
    #climatology years
    cyr1,cyr2='1993-01-01','2018-12-31'
    
    # CCMP test
    dir_pattern_zarr = 'F:/data/sat_data/ccmp/zarr/'
    ds= xr.open_zarr(dir_pattern_zarr)
    ds = ds.rename({'latitude':'lat','longitude':'lon'})
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    ds_ccmp = ds.sortby(ds.lon)
    ds_ccmp = ds_ccmp.drop('nobs')
    for var in ds_ccmp:
        tem = ds_ccmp[var].attrs
        tem['var_name']='ccmp_'+str(var)
        ds_ccmp[var].attrs=tem
    ds_ccmp_clim = ds_ccmp.sel(time=slice(cyr1,cyr2))
    ds_ccmp_clim = ds_ccmp_clim.groupby('time.dayofyear').mean('time',keep_attrs=True,skipna=False)
    
    # AVISO test
    dir_pattern_zarr = 'F:/data/sat_data/aviso/zarr/'
    ds= xr.open_zarr(dir_pattern_zarr)
    ds = ds.rename({'latitude':'lat','longitude':'lon'})
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    ds_aviso = ds.sortby(ds.lon).drop({'lat_bnds','lon_bnds','crs','err'})
    for var in ds_aviso:
        tem = ds_aviso[var].attrs
        tem['var_name']='aviso_'+str(var)
        ds_aviso[var].attrs=tem
    ds_aviso_clim = ds_aviso.sel(time=slice(cyr1,cyr2))
    ds_aviso_clim = ds_aviso_clim.groupby('time.dayofyear').mean('time',keep_attrs=True,skipna=False)    

    #sst
    dir_pattern_zarr = 'F:/data/sat_data/sst/cmc/zarr/'
    ds_sst= xr.open_zarr(dir_pattern_zarr)
    ds_sst = ds_sst.drop({'analysis_error','mask','sea_ice_fraction'})
    tem = ds_sst.analysed_sst.attrs
    tem['var_name']='cmc_sst'
    ds_sst.analysed_sst.attrs=tem
    ds_sst_clim = ds_sst.sel(time=slice(cyr1,cyr2))
    ds_sst_clim = ds_sst_clim.groupby('time.dayofyear').mean('time',keep_attrs=True,skipna=False)
    
    #get bathymetry from ETOPO1
    fname_topo = 'F:/data/topo/ETOPO1_Ice_g_gmt4.grd'
    ds = xr.open_dataset(fname_topo)
    ds_topo = ds.rename_dims({'x':'lon','y':'lat'}).rename({'x':'lon','y':'lat'})
    tem = ds_topo.z.attrs
    tem['var_name']='etopo_depth'
    ds_topo.z.attrs=tem

    #put data into a dictionary
    data_dict={'aviso':ds_aviso,
               'wnd':ds_ccmp,
               'sst':ds_sst,
              'topo':ds_topo}
    clim_dict={'aviso_clim':ds_aviso_clim,
               'wnd_clim':ds_ccmp_clim,
               'sst_clim':ds_sst_clim}
  
    return data_dict,clim_dict

data,clim = get_data()

ilen_bird1 = len(ds_bird.track)
ilen_bird2 = len(ds_bird.observation_number)
for name in data:
    ds_data=data[name]
    if name=='topo':
        continue
    if not name==input_data:
        continue
    print('name',name)   
    for var in ds_data:
        var_tem=ds_data[var].attrs['var_name']
        ds_bird[var_tem]=xr.DataArray(np.nan*np.empty((ilen_bird1,ilen_bird2), 
                                                      dtype=str(ds_data[var].dtype)), 
                                      coords={'track': ds_bird.track,'observation_number':ds_bird.observation_number},
                                      dims=('track','observation_number'))
        ds_bird[var_tem].attrs=ds_data[var].attrs
    print('var',var_tem)
    for i in range(ilen_bird1):
        for j in range(ilen_bird2):
            if np.isnan(ds_bird.lat[i,j]):
                continue
            if ds_bird.time[i,j]<ds_data.time.min():
                continue
            if ds_bird.time[i,j]>ds_data.time.max():
                continue
            t1,t2 = ds_bird.time[i,j]-np.timedelta64(24,'h'), ds_bird.time[i,j]+np.timedelta64(24,'h')
            lat1,lat2=ds_bird.lat[i,j]-.5,ds_bird.lat[i,j]+.5
            lon1,lon2=ds_bird.lon[i,j]-.5,ds_bird.lon[i,j]+.5
            tem = ds_data.sel(time=slice(t1,t2),lat=slice(lat1,lat2),lon=slice(lon1,lon2)).load()
            tem = tem.interp(time=ds_bird.time[i,j],lat=ds_bird.lat[i,j],lon=ds_bird.lon[i,j])
            for var in ds_data:
                var_tem=ds_data[var].attrs['var_name']
                ds_bird[var_tem][i,j]=tem[var].data
        print(i,ilen_bird1)
    df_bird = ds_bird.to_dataframe()
    ds_bird.to_netcdf(filename_origin_out+name+'.nc')


