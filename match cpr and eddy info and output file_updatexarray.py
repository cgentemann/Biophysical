#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import cos, radians
import xarray as xr

####################you will need to change some paths here!#####################
#list of input files
adir = 'f:/data/project_data/NASA_biophysical/'
filename_cpr=adir + 'CPR_data/All CPR Sample catalogue.xlsx'
filename_northpac_eddies=adir + '/aviso/eddy_trajectory_19930101_20170106_north_pacific_2020_10_06a.nc'
filename_cpr_eddy=adir + '/collocated_data/CPR/eddy_cpr_data_north_pacific.nc'
filename_eddy=adir + '/collocated_data/CPR/eddy_ranking_data_north_pacific.nc'
#output files
filename_cpr_expanded=aidr + '/collocated_data/CPR/All CPR Sample catalogue with env info_2020_10_05a'
fname_topo = 'F:/data/topo/ETOPO1_Ice_g_gmt4.grd'
#################################################################################


input_data = input("Enter data: aviso,wnd,sst")


#read in CPR data excell file using pandas library
df = pd.read_excel(filename_cpr)
df = df.rename(columns={'Sample ID':'cpr_sample_id','day':'cpr_sample_day',
                'month':'cpr_sample_month','year':'cpr_sample_year',
                'lat':'cpr_sample_lat','Long':'cpr_sample_lon','Already processed?':'cpr_sample_proc'})
ds_cpr = df.to_xarray()
ds_cpr['index']=ds_cpr.index.astype('int')
ds_cpr
ilen = ds_cpr.cpr_sample_lat.size
tt=np.empty(ilen,dtype='datetime64[ns]') 
for i in range(ilen):
    tstr=str(ds_cpr.cpr_sample_year[i].data)+'-'+str(ds_cpr.cpr_sample_month[i].data).zfill(2)+'-'+str(ds_cpr.cpr_sample_day[i].data).zfill(2)
    tem=np.datetime64(tstr)
    tt[i]=tem
ds_cpr['cpr_sample_time']=xr.DataArray(tt,dims=['index'])

ds_cpr.cpr_sample_lon.min().data,ds_cpr.cpr_sample_lon.max().data,ds_cpr.cpr_sample_lat.min().data,ds_cpr.cpr_sample_lat.max().data

ds_eddy = xr.open_dataset(filename_northpac_eddies).rename({'Longitude':'lon','Latitude':'lat'})

ds_eddy_cpr = xr.open_dataset(filename_cpr_eddy)


#get bathymetry from ETOPO1
ds = xr.open_dataset(fname_topo)
ds_topo = ds.rename_dims({'x':'lon','y':'lat'}).rename({'x':'lon','y':'lat'})
#tem = ds_topo.isel(lat=slice(7000,9500),lon=slice(0,4500))
#tem.z.plot()
tt = ds_topo.z.interp(lat=ds_cpr.cpr_sample_lat,lon=ds_cpr.cpr_sample_lon,method='nearest').data
ds_cpr['ETOPO_depth']= xr.DataArray(tt, coords={'index':ds_cpr.index}, dims=["index"])
ds_cpr['cpr_sample_lon2'] = np.mod(ds_cpr['cpr_sample_lon'],360)

#plt.scatter(ds_cpr.cpr_sample_lon2,ds_cpr.cpr_sample_lat,c=ds_cpr.ETOPO_depth,cmap='coolwarm',vmin=-8000,vmax=8000)

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


# In[ ]:


#ds_cpr = xr.open_dataset(filename_bird_out_eddy_netcdf)
ilen_bird1 = len(ds_cpr.cpr_sample_lon)

clonmin,clonmax = ds_cpr.cpr_sample_lon.min().data,ds_cpr.cpr_sample_lon.max().data
clatmin,clatmax = ds_cpr.cpr_sample_lat.min().data,ds_cpr.cpr_sample_lat.max().data
t1save=0

for name in data:
    ds_data=data[name]
    if name=='topo':
        continue
    if not name==input_data:
        continue
    print('name',name)   
    for var in ds_data:
        var_tem=ds_data[var].attrs['var_name']
        ds_cpr[var_tem]=xr.DataArray(np.nan*np.empty((ilen_bird1), 
                                                      dtype=str(ds_data[var].dtype)), 
                                      coords={'index': ds_cpr.index},
                                      dims=('index'))
        ds_cpr[var_tem].attrs=ds_data[var].attrs
    print('var',var_tem)
    for i in range(ilen_bird1):
        if np.isnan(ds_cpr.cpr_sample_lat[i]):
            continue
        if ds_cpr.cpr_sample_time[i]<ds_data.time.min():
            continue
        if ds_cpr.cpr_sample_time[i]>ds_data.time.max():
            continue
        t1,t2 = ds_cpr.cpr_sample_time[i]-np.timedelta64(24,'h'), ds_cpr.cpr_sample_time[i]+np.timedelta64(24,'h')
        if not t1==t1save:
            tem2 = ds_data.sel(time=slice(t1,t2),lat=slice(clatmin-.5,clatmax+.5),lon=slice(clonmin-.5,clonmax+.5)).load()               
            t1save=t1
            print(i,ilen_bird1)
            #            lat1,lat2=ds_cpr.cpr_sample_lat[i]-.5,ds_cpr.cpr_sample_lat[i,j]+.5
            #            lon1,lon2=ds_cpr.cpr_sample_lon[i]-.5,ds_cpr.cpr_sample_lon[i,j]+.5
            #            tem = ds_data.sel(time=slice(t1,t2),lat=slice(lat1,lat2),lon=slice(lon1,lon2)).load()
        tem = tem2.interp(time=ds_cpr.cpr_sample_time[i],lat=ds_cpr.cpr_sample_lat[i],lon=ds_cpr.cpr_sample_lon[i])
        for var in ds_data:
            var_tem=ds_data[var].attrs['var_name']
            ds_cpr[var_tem][i]=tem[var].data
    #output data
    df_bird = ds_cpr.to_dataframe()
    df_bird.to_csv(filename_cpr_expanded+name+'.csv')
    ds_cpr.to_netcdf(filename_cpr_expanded+name+'.nc')

