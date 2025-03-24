import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
from scipy.ndimage import zoom
# import not installed utils
import sys
sys.path.append('functions/')
from save_fig_util import savefig


def emis_map( ax,  mapdata, latm, lonm, figlb ):

    plt.rcParams.update({'font.size': 16}) 
    ax.set_extent([-170, 170, -60, 80], crs=ccrs.PlateCarree())  # Set longitude from -180 to 180 and latitude from -60 to 90
    ax.coastlines()
    
    lon_grid, lat_grid = np.meshgrid(lonm, latm)
    mapax = ax.pcolormesh(lon_grid, lat_grid, mapdata, transform=ccrs.PlateCarree(), cmap='YlOrBr', vmin=0, vmax=4e-11)

    # Add a legend 
    cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
    cbar = plt.colorbar(mapax, cax=cax, orientation='vertical')
    cbar.set_label(r'SO$_2$ emission (kg m$^{-2}$ s$^{-1}$)')

    ax.text(10, -48, figlb, ha='center', va='center', transform=ccrs.PlateCarree(),fontsize = 20, fontweight = 'bold')
    
    plt.tight_layout()
    plt.show()
  
def diff_map( ax,  mapdata, latm, lonm, figlb ):

    plt.rcParams.update({'font.size': 16}) 
    ax.set_extent([-170, 170, -60, 80], crs=ccrs.PlateCarree())  # Set longitude from -180 to 180 and latitude from -60 to 90
    ax.coastlines()
    
    lon_grid, lat_grid = np.meshgrid(lonm, latm)
    mapax = ax.pcolormesh(lon_grid, lat_grid, mapdata, transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=-1e-11, vmax=1e-11)

    # Add a legend 
    cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
    cbar = plt.colorbar(mapax, cax=cax, orientation='vertical')
    cbar.set_label(r'SO$_2$ emission diff. (kg m$^{-2}$ s$^{-1}$)')

    ax.text(10, -48, figlb, ha='center', va='center', transform=ccrs.PlateCarree(),fontsize = 20, fontweight = 'bold')
    
    plt.tight_layout()
    plt.show()
    

# functions and paths
rootdir = '/storage1/fs1/rvmartin/Active'
savedir = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/01.Emissions/'
cedsfname  = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/CEDS-2024-Gridded/SO2-em-anthro_CMIP_CEDS_2021.nc' 
htapfname  = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/HTAPv3/scale_2018_to_2021/2021/HTAPv3_SO2_0.1x0.1_2021.nc' 
edgarfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/EDGAR/EDGARv81/v8.1_SO2_2021_scaled_w_v61_0.1x0.1.nc'

fig, axes = plt.subplots(3, 1, figsize=(17,21), subplot_kw={'projection': ccrs.PlateCarree()}) 

# read ceds map data
ds = xr.open_dataset(cedsfname) 
so2_vars = [var for var in ds.data_vars if var.lower().startswith('so2')]
ceds = ds[so2_vars].to_array(dim='variable').sum(dim='variable', skipna=True)
ceds = ceds.mean(dim='time', skipna=True)
latc = ds['lat'].values 
lonc = ds['lon'].values
figlb = 'CEDS'
emis_map( axes[0], ceds, latc, lonc, figlb )


### Figure 1 - GCHPwEDGAR - GCHPwCEDS
# read edgar map data
ds = xr.open_dataset(edgarfname) 
so2_vars = [var for var in ds.data_vars if var.lower().startswith('so2')]
edgar = ds[so2_vars].to_array(dim='variable').sum(dim='variable', skipna=True)
edgar = edgar.mean(dim='time', skipna=True)
late = ds['lat'].values 
lone = ds['lon'].values

# regrid ceds to 0.1:
ceds = zoom(ceds, 5, order=0)
mapdata = edgar-ceds
# make the map
figlb = 'EDGAR - CEDS'
diff_map( axes[1], mapdata, late, lone, figlb )



### Figure 2 - GCHPwHTAP - GCHPwCEDS
# read edgar map data
ds = xr.open_dataset(htapfname) 
so2_vars = [var for var in ds.data_vars if var.lower().startswith('so2')]
htap = ds[so2_vars].to_array(dim='variable').sum(dim='variable', skipna=True)
htap = htap.mean(dim='time', skipna=True)

mapdata = htap-ceds
# make the map
figlb = 'HTAP - CEDS'
diff_map( axes[2], mapdata, late, lone, figlb )

# save figure
savefig(f'{savedir}/SO2_emission_map_diff.png', fig)