# This script makes maps of global O3

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import cartopy.feature as cfeature

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
    cbar.set_label('O3 VCD (DU)',fontsize = 12)
    
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
year = 2021
qcstr =  'CF03-SZA75-QA75'
# for mn in range(5,13):
indir = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/'
savedir = f'/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/temp'

for mn in range(1,14):
    if mn ==13:
        fname = f"{indir}/Regrid_tropomi_O3_{year}_{qcstr}.nc"
    else:
        fname = f"{indir}/Regrid_tropomi_O3_{year}{mn:02d}_{qcstr}.nc"

    ds = xr.open_dataset(fname) 
    so2_1 = ds['O3'].values.T # check variable name
    so2_1 = so2_1/2.69e16 # convert unit
    lat = ds['lat'].values # presumably all sources share the same lat lon
    lon = ds['lon'].values

    # regrid to 0.5
    res = 0.5
    scale = int(res/0.05)
    shape1 = len(lat)
    shape2 = len(lon)

    so2_1_r = so2_1.reshape(int(shape1/scale),scale,int(shape2/scale),scale)
    so2_1_r = so2_1_r.mean(axis=(1, 3))

    lat_r = lat.reshape(int(shape1/scale),scale)
    lat_r = lat_r.mean(axis=1)
    lon_r = lon.reshape(int(shape2/scale),scale)
    lon_r = lon_r.mean(axis=1)



    # Create a figure and set the projection to a global map with the Pacific centered (180 degrees)
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Plot the concentration data using the meshgrid
    concentration_plot = ax.pcolormesh(lon_r, lat_r, so2_1_r, cmap='RdBu_r', transform=ccrs.PlateCarree(), shading='auto', vmin=300, vmax=500)

    # Add coastlines and features
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    # Add colorbar
    cax = ax.inset_axes([0.16, 0.15, 0.02, 0.3])  # Adjust these parameters as needed for placement
    plt.colorbar(concentration_plot, cax=cax, orientation='vertical',label='O3 VCD (DU)')
        
    # Set gridlines
    ax.gridlines(draw_labels=True)

    # Display the plot
    plt.title(f"$O_3$ Maps")

    # save figure
    fname = f"{savedir}/map_O3_doas-{qcstr}_coarse{res}_{year}{mn}.png"
    savefig(fname, fig)

