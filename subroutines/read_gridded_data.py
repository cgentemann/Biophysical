#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def get_monthly_data(stype):
    from pathlib import Path
    import xarray as xr
    import numpy as np
    
    if np.char.lower(stype)=='oscar':
        dir_data = 'F:/data/sat_data/oscar/L4/oscar_third_deg/oscar_vel*.nc'
        ds=xr.open_mfdataset(dir_data,combine='nested',concat_dim='time').isel(depth=0).drop({'um','vm'}).rename({'latitude':'lat','longitude':'lon'})
        ds = ds.sel(lon=slice(20.0,379.9))
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon').sortby('lat').drop({'year'})
        ds = ds.resample(time='M',keep_attrs=True).mean(skipna=False,keep_attrs=True)
    if np.char.lower(stype)=='godas_mld':
        dir_data = dir_data_mld='F:/data/model_data/godas/*.nc'
        ds=xr.open_mfdataset(dir_data,combine='nested',concat_dim='time')
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon').sortby('lat')
        ds = ds.resample(time='M',keep_attrs=True).mean(skipna=False,keep_attrs=True)
    if np.char.lower(stype)=='ccmp':
        filelist=[]
        dir_data = 'F:/data/sat_data/ccmp/v02.0/'
        for filename in Path(dir_data).rglob('*L3.0_RSS.nc'):
            filelist.append(filename)
            ds=xr.open_mfdataset(filelist,combine='nested',concat_dim='time').rename({'latitude':'lat','longitude':'lon'})
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon').sortby('lat')
        _, index = np.unique(ds['time'], return_index=True)
        ds=ds.isel(time=index)           
        ds = ds.resample(time='M',keep_attrs=True).mean(skipna=False,keep_attrs=True)
    if np.char.lower(stype)=='mur_sst':
        filelist=[]
        dir_data = 'F:/data/sst/jpl_mur/v4.1/'+str(iyr).zfill(4)+'/'
        for filename in Path(dir_data).rglob('*90000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'):
            filelist.append(filename)
        ds=xr.open_mfdataset(filelist,combine='nested',concat_dim='time').drop({'analysis_error'})
        ds = ds.resample(time='M',keep_attrs=True).mean(skipna=False,keep_attrs=True)
    if np.char.lower(stype)=='cmc_sst':
        #combine 0.2 and 0.1 deg data into 0.2 dataset
        filelist=[]
        dir_data = 'F:/data/sst/cmc/CMC0.2deg/v2/data/'
        for filename in Path(dir_data).rglob('*CMC-L4_GHRSST-SSTfnd-CMC0.2deg-GLOB-v02.0-fv02.0.nc'):
            filelist.append(filename)
        ds=xr.open_mfdataset(filelist,combine='nested',concat_dim='time')
        filelist=[]
        dir_data = 'F:/data/sst/cmc/CMC0.1deg/v3/'
        for filename in Path(dir_data).rglob('*CMC-L4_GHRSST-SSTfnd-CMC0.1deg-GLOB-v02.0-fv03.0.nc'):
            filelist.append(filename)
        ds2=xr.open_mfdataset(filelist,combine='nested',concat_dim='time')
        ds2=ds2.interp(lat=ds.lat,lon=ds.lon)
        ds=xr.concat([ds,ds2],dim='time')       
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon').sortby('lat')
        ds = ds.resample(time='M',keep_attrs=True).mean(skipna=False,keep_attrs=True)
    if np.char.lower(stype)=='cmem_sss':
        filelist=[]
        dir_data = 'F:/data/model_data/CMEM/global-reanalysis-phy-001-030-monthly/'
        for filename in Path(dir_data).rglob('subset_mercatorglorys12v1_gl12_mean_*.nc'):
            filelist.append(filename)
        ds=xr.open_mfdataset(filelist,combine='nested',concat_dim='time').rename({'latitude':'lat','longitude':'lon'})
        ds=ds.drop({'zos','bottomT','sithick','usi','vsi','thetao','uo','vo','mlotst'})
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon').sortby('lat')
        ds = ds.resample(time='M',keep_attrs=True).mean(skipna=False,keep_attrs=True)
    if np.char.lower(stype)=='aviso':
        filelist=[]
        dir_data = 'F:/data/sat_data/aviso/data/'
        from pathlib import Path
        for filename in Path(dir_data).rglob('*.nc'):
            filelist.append(filename)
        ds=xr.open_mfdataset(filelist,combine='nested',concat_dim='time').drop({'ugosa','vgosa','err'}).rename({'latitude':'lat','longitude':'lon'})
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon').sortby('lat')
        ds = ds.resample(time='M',keep_attrs=True).mean(skipna=False,keep_attrs=True)
    return ds

