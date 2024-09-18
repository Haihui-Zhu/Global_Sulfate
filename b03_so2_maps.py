# This script makes maps of global and regional SO2 and SO4

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec

def plot_setup():
    regions = ['Global', 'North America', 'Europe', 'East Asia', 'South Asia', 'Middle East']

    extents = {
        'Global':[-170, 170, -50, 60],
        'North America': [-140, -50, -15, 75],
        'Europe': [-10, 40, 20, 70],
        'East Asia': [95, 145, 5, 55],
        'South Asia': [60, 100, 0, 40],
        'Middle East': [25, 65, 5, 45]
    }

    long_min, long_max, lat_min, lat_max = extents['Global']
    # Calculate the ranges for longitude and latitude
    lat_range = lat_max - lat_min
    long_range = long_max - long_min
    # Calculate the ratio of longitude to latitude range
    ratio =  long_range / lat_range
    fig = plt.figure(figsize=(5*ratio, 15))
    gs = gridspec.GridSpec(3, 1, wspace=0, hspace=0, height_ratios=[1, 1, 1])  

    return regions, extents, fig, gs

def plot_map(data, lat, lon,extents, label, ax, diff_cmap=False):
    lon_min, lon_max, lat_min, lat_max = extents['Global']

    ax.coastlines()
    ax.set_global() 
    ax.set_extent([lon_min+5, lon_max-5,lat_min+5,lat_max-5]) # avoid white space at edges

    # Flatten the latitude and longitude arrays to ensure they are 1D
    lon = lon.flatten()
    lat = lat.flatten()

    lon_range = (lon >= lon_min) & (lon <= lon_max)  
    lat_range = (lat >= lat_min) & (lat <= lat_max) 

    regional_data = data[np.ix_(lat_range, lon_range)]
    regional_total = np.nansum(regional_data)  # Ensure to sum only non-NaN values
    regional_mean = np.nanmean(regional_data)  # Ensure to sum only non-NaN values

    lon_grid, lat_grid = np.meshgrid(lon[lon_range], lat[lat_range])
    if diff_cmap: 
        mesh = ax.pcolormesh(lon_grid, lat_grid, regional_data, transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=-0.2, vmax=0.2)
    else:
        mesh = ax.pcolormesh(lon_grid, lat_grid, regional_data, transform=ccrs.PlateCarree(), cmap='YlOrBr', vmin=0, vmax=0.5)

    
    cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
    cbar = plt.colorbar(mesh, cax=cax, orientation='vertical')
    cbar.set_label('SO2 VCD (DU)',fontsize = 12)
    
    x, y = 0.5, 0.05  # Normalized coordinates (0, 0) is lower left, (1, 1) is upper right
    if diff_cmap:
        text = f"Mean difference: {regional_mean:.2e}"
    else:
        text = f"Sum: {regional_total:.2e}\nMean: {regional_mean:.2e}"

    plt.text(x, y, text, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
             fontsize =14, fontweight = 'bold',backgroundcolor = 'w')
    # adding label of data source
    x, y = 0.2, 0.05
    plt.text(x, y, label, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
             fontsize =14, fontweight = 'bold',backgroundcolor = 'w')
    
def savefig(fname, fig):
    fig.savefig(fname) 
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: '+ fname)
    return


### main starts here ###
simname = 'ceds_2018'
year = 2018
qcstr =  'CF03-SZA60-QA90'
for mn in range(5,13):
    source1 = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/Tropomi_Regrid_{year}{mn:02d}_{qcstr}.nc'
    source2 = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_{mn:02d}.nc'
    savedir = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/'

    ### Figure 1: comparing global so2 from cobra + gchp ###
    ds = xr.open_dataset(source1) 
    so2_1 = ds['so2'].values.T # check variable name
    so2_1 = so2_1/2.69e16 # convert unit
    lat = ds['lat'].values # presumably all sources share the same lat lon
    lon = ds['lon'].values

    # ds = xr.open_dataset(source1)
    # so2_vol = ds['so2_vol'].values.T# adding volcanoes
    # so2_vol = so2_vol/2.69e16
    # # sum with nan:
    # array1_no_nan = np.nan_to_num(so2_pbl, nan=0.0)
    # array2_no_nan = np.nan_to_num(so2_vol, nan=0.0)
    # sum_array = array1_no_nan + array2_no_nan
    # nan_mask = np.isnan(so2_pbl) & np.isnan(so2_vol)
    # so2_1 = np.where(nan_mask, np.nan, sum_array)


    ds = xr.open_dataset(source2) 
    so2_2 = ds['so2'].values.T # check variable name
    so2_2 = so2_2/2.69e16 # convert unit


    # Plot set up
    regions, extents, fig, gs = plot_setup()
    # Plotting the global map 
    ax_global = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    plot_map(so2_1, lat, lon,extents, "TROPOMI-DOAS-QA90 $SO_2$", ax_global)
    ax_global = fig.add_subplot(gs[1, :], projection=ccrs.PlateCarree())
    plot_map(so2_2, lat, lon,extents, "GCHP $SO_2$ ", ax_global)
    ax_global = fig.add_subplot(gs[2, :], projection=ccrs.PlateCarree())
    plot_map(so2_2-so2_1, lat, lon,extents, 'GCHP-TROPOMI', ax_global, diff_cmap=True)

    plt.suptitle(f"$SO_2$ Maps")
    # save figure
    fname = f"{savedir}map_so2_doas-{qcstr}_vs_{simname}_{mn}.png"
    savefig(fname, fig)
