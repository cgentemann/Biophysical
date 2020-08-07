import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import warnings
# filter some warning messages
warnings.filterwarnings("ignore") 
from geopy.distance import geodesic 

####################you will need to change some paths here!#####################
#list of input files
#list of input files

filename_aviso='F:/data/project_data/NASA_biophysical/aviso/eddy_trajectory_19930101_20170106.nc'   #From AVISO  website
filename_bird='f:/data/project_data/NASA_biophysical/collocated_data/NPPSD_GOA_allseabird_wide.csv'
#output files
filename_bird_out='f:/data/project_data/NASA_biophysical/collocated_data/NPPSD_GOA_allseabird_with_sat_and eddy_data.csv'
filename_bird_out_netcdf='f:/data/project_data/NASA_biophysical/collocated_data/NPPSD_GOA_allseabird_wide_sat_and_eddy_data.nc'
#################################################################################
filename_bird_out_eddy_netcdf='f:/data/project_data/NASA_biophysical/collocated_data/NPPSD_GOA_allseabird_with_eddy.nc'

input_data = input("Enter data: aviso,wnd,sst")

#define function to get all the data at once, use same years for climatology for all data
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


ds_bird = xr.open_dataset(filename_bird_out_eddy_netcdf)
ilen_bird = len(ds_bird.lat)
for name in data:
    ds_data=data[name]
    if name=='topo':
        continue
    if not name==input_data:
        continue
    print('name',name)   
    for var in ds_data:
        var_tem=ds_data[var].attrs['var_name']
        ds_bird[var_tem]=xr.DataArray(np.empty(ilen_bird, dtype=str(ds_data[var].dtype)), coords={'index': ds_bird.index}, dims=('index'))
        ds_bird[var_tem].attrs=ds_data[var].attrs
    print('var',var_tem)
    for i in range(len(ds_bird.lat)):
        if ds_bird.time[i]<ds_data.time.min():
            continue
        if ds_bird.time[i]>ds_data.time.max():
            continue
        t1,t2 = ds_bird.time[i]-np.timedelta64(24,'h'), ds_bird.time[i]+np.timedelta64(24,'h')
        lat1,lat2=ds_bird.lat[i]-.5,ds_bird.lat[i]+.5
        lon1,lon2=ds_bird.lon[i]-.5,ds_bird.lon[i]+.5
        tem = ds_data.sel(time=slice(t1,t2),lat=slice(lat1,lat2),lon=slice(lon1,lon2)).load()
        tem = tem.interp(time=ds_bird.time[i],lat=ds_bird.lat[i],lon=ds_bird.lon[i])
        #tem = tem.load()
        for var in ds_data:
            var_tem=ds_data[var].attrs['var_name']
            ds_bird[var_tem][i]=tem[var].data
        if int(i/10000)*10000==i:
            print(i,len(ds_bird.lat))
#at topo info
#interp will create a new 2D array, to avoid that put the lat/lon into dataarrays
if name=='ccmc':
    ds_topo=data['topo']
    new_lat = xr.DataArray(ds_bird.lat.data, dims='z')
    new_lon = xr.DataArray(ds_bird.lon.data, dims='z')
    tem = ds_topo.z.interp(lat=new_lat, lon=new_lon,method='nearest')
    ds_bird['ETOPO_depth'] = xr.DataArray(tem.data, coords={'index': ds_bird.index}, dims=('index'))
#ds_bird = ds_bird.drop('obs')

#output data
df_bird = ds_bird.to_dataframe()
df_bird.to_csv(filename_bird_out+input_data+'.csv')
ds_bird.to_netcdf(filename_bird_out_netcdf+input_data+'.nc')

