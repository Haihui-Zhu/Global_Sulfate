import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
# import not installed utils
import sys
sys.path.append('functions/')
from save_fig_util import savefig

def diff_map( ax,  mapdata, latm, lonm, figlb ):

    plt.rcParams.update({'font.size': 16}) 
    ax.set_extent([-170, 170, -60, 80], crs=ccrs.PlateCarree())  # Set longitude from -180 to 180 and latitude from -60 to 90
    ax.coastlines()
    
    lon_grid, lat_grid = np.meshgrid(lonm, latm)
    mapax = ax.pcolormesh(lon_grid, lat_grid, mapdata, transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=-6, vmax=6)

    # Add a legend 
    cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
    cbar = plt.colorbar(mapax, cax=cax, orientation='vertical')
    cbar.set_label('Sulfate Difference (µg/m³)')

    ax.text(10, -48, figlb, ha='center', va='center', transform=ccrs.PlateCarree(),fontsize = 20, fontweight = 'bold')
    
    plt.tight_layout()
    plt.show()
    

# functions and paths
rootdir = '/storage1/fs1/rvmartin/Active'
savedir = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/01.Emissions/'
cedsfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/4.ceds_2021/GCHP_SO2_SO4_PM25_ceds_2021_annual.nc'
htapfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/2.htap_2018/GCHP_SO2_SO4_BC_PM25_htap_2018_annual.nc'
edgarfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/5.edgar_2021/GCHP_SO2_SO4_BC_PM25_edgar_2021_annual.nc'

# # load mask file
# data_mask = loadmat('ceds_scale_2021to2018/mask_fao_ceds_05.mat')
# mask_01 = data_mask['mask_region']
# mlat_01 = data_mask['xlat']
# mlon_01 = data_mask['xlon']
# region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
# region_name = [item[0] for item in region_name_array.flatten()]

# read ceds map data
ds = xr.open_dataset(cedsfname) 
ceds = ds['so4'].values.T # check variable name
latg = ds['latitude'].values 
long = ds['longitude'].values


fig, axes = plt.subplots(2, 1, figsize=(17,14), subplot_kw={'projection': ccrs.PlateCarree()}) 

### Figure 1 - GCHPwEDGAR - GCHPwCEDS
# read edgar map data
ds = xr.open_dataset(edgarfname) 
edgar = ds['so4'].values.T # check variable name
mapdata = edgar-ceds
# make the map
figlb = 'GCHPwEDGAR - GCHPwCEDS'
diff_map( axes[0], mapdata, latg, long, figlb )



### Figure 2 - GCHPwHTAP - GCHPwCEDS
# read edgar map data
ds = xr.open_dataset(htapfname) 
htap = ds['so4'].values.T # check variable name
mapdata = htap-ceds
# make the map
figlb = 'GCHPwHTAP - GCHPwCEDS'
diff_map( axes[1], mapdata, latg, long, figlb )

# save figure
savefig(f'{savedir}/sulfate_map_diff.png', fig)