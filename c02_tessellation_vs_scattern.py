# This script examine SO2 from TROPOMI and GCHP with focuses on
# 1. how cloud coverage affects sampling and eventually estimated SO2 VCD (03.NASA_L2_RPRO_VCD vs. 04.NASA_L2_RPRO_VCD_CF03)
# 2. how tessellation differ from scatteredInterpolant (05.NASA_RPRO_Tessellation vs. 04.NASA_L2_RPRO_VCD_CF03), including sampling effect related to fake 0s. 
# 3. how tropomi VCD differ from GCHP (05 or 04 vs. 01.GCHP_13.4.0_SO2)

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.io import loadmat
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata

def plot_setup():
    regions = ['Global', 'North America', 'Europe', 'East Asia', 'South Asia', 'Middle East']

    extents = {
        'Global':[-180, 180, -50, 60],
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
    # Calculate the ratio of latitude to longitude range
    ratio =  long_range / lat_range
    # ratio of global map height to regional map height    
    h1 = 5/ratio
    # ratio of figure height to width
    overallh = (1+h1)/5
    fig = plt.figure(figsize=(15, 15*overallh))
    gs = gridspec.GridSpec(2, 5, wspace=0, hspace=0, height_ratios=[h1, 1], width_ratios=[1,1,1,1,1])  # Changed number of columns to 5 if you have 5 regions

    return regions, extents, fig, gs

def plot_map(data, lat, lon,extents, region, ax, vmin, vmax, add_cbar=False):
    lon_min, lon_max, lat_min, lat_max = extents[region]

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

    if vmin == 0:
        mesh = ax.pcolormesh(lon_grid, lat_grid, regional_data, transform=ccrs.PlateCarree(), cmap='YlOrBr', vmin=vmin, vmax=vmax)
    else:
        mesh = ax.pcolormesh(lon_grid, lat_grid, regional_data, transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=vmin, vmax=vmax)

    if add_cbar:
        cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
        cbar = plt.colorbar(mesh, cax=cax, orientation='vertical')
        if vmax == 1:
            cbar.set_label('SO2 VCD (DU)',fontsize = 12)
        else :
            if vmax == 10: 
               cbar.set_label(f'$SO_4$ ($\mu$g/$m^3$)',fontsize = 12) 
        
    
    x, y = 0.5, 0.05  # Normalized coordinates (0, 0) is lower left, (1, 1) is upper right
    text = f"Sum: {regional_total:.2e}\nMean: {regional_mean:.2e}"
    plt.text(x, y, text, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
             fontsize =12, backgroundcolor = 'w')

def savefig(fname, fig):
    fig.savefig(fname) 
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: '+ fname)
    return

def loadso2(fname):
    data = loadmat(fname)

    # Function to normalize variable names
    def normalize_name(name):
        name = name.lower()  # Convert to lower case
        if name in ['so2', 'tso2']:
            return 'so2'
        elif name in ['lat', 'tlat']:
            return 'lat'
        elif name in ['lon', 'tlon']:
            return 'lon'
        return name  # Return the name unchanged if it doesn't match any known pattern

    # Create a new dictionary with normalized keys
    normalized_data = {normalize_name(key): value for key, value in data.items() if not key.startswith('__')}

    # You can now access the variables directly by their normalized names:
    so2 = normalized_data.get('so2', None)
    lat = normalized_data.get('lat', None)
    lon = normalized_data.get('lon', None)

    lon = lon.flatten()
    lat = lat.flatten()
    print(so2.shape)
    print(lat.shape)
    print(lon.shape)

    # only load data within lat [-70 70] and lon [-180 180]
    lat_min = -70
    lat_max = 70
    lat_range = (lat >= lat_min) & (lat <= lat_max)
    print(lat_range.shape)
    lat = lat[lat_range]
    so2 = so2[lat_range,:]

    return so2, lat, lon


### Script starts here ###
indir   = './03.TROPOMI_SO2_Result/'
savedir = './03.TROPOMI_SO2_Result/06.ScatternIn_VS_Tessellation_Figures/'
datadir = ['01.GCHP_13.4.0_SO2', '03.NASA_L2_RPRO_VCD', '04.NASA_L2_RPRO_VCD_CF03', '05.NASA_RPRO_Tessellation']
datanames = ['gchp_13.4.0','tropomi_rpro_scatterin','tropomi_rpro_scatterin_cf03','tropomi_rpro_tess_cf03']

year = 2019

refidx = 1 # index of dataset as a reference, range from 0 to 3
devidx = 2 # index of dataset to evaluate, range from 0 to 3

### Monthly Maps ###

for mn in range(13,14):
    if mn == 13:
        datestr = f'{year}'
    else:
        datestr = f'{year}{mn:02d}'

    file_name =[f'SO2_{datestr}_0.1x0.1deg.mat',
                f'SO2_{datestr}_0.1x0.1deg.mat',
                f'SO2_{datestr}_0.1x0.1deg.mat',
                f'Tropomi_Regrid_{datestr}.mat']
        
    # color scale 
    global_min = 0
    global_max = 1

    ### reading ref dataset  ###
    so2, lat, lon = loadso2(f"{indir}/{datadir[refidx]}/{file_name[refidx]}")
    so2 = so2/2.69e16

    # Plot set up
    regions, extents, fig, gs = plot_setup()
    # Plotting the global map as the first and largest plot
    ax_global = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    plot_map(so2, lat, lon,extents, 'Global', ax_global, global_min, global_max, add_cbar=True)

    for idx, region in enumerate(regions[1:]):
        ax_region = fig.add_subplot(gs[1, idx], projection=ccrs.PlateCarree())
        plot_map(so2, lat, lon, extents,region, ax_region, global_min, global_max)

    plt.suptitle(f"$SO_2$ Maps for {datanames[refidx]}")
    
    # save figure
    fname = f"{savedir}so2_{datanames[refidx]}_{datestr}.png"
    savefig(fname, fig)


    ### reading dev dataset  ###
    so2t, latt, lont  = loadso2(f"{indir}/{datadir[devidx]}/{file_name[devidx]}")
    so2t = so2t/2.69e16
    
    # Plot set up
    regions, extents, fig, gs = plot_setup()
    # Plotting the global map as the first and largest plot
    ax_global = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    plot_map(so2t, latt, lont,extents, 'Global', ax_global, global_min, global_max, add_cbar=True)

    for idx, region in enumerate(regions[1:]):
        ax_region = fig.add_subplot(gs[1, idx], projection=ccrs.PlateCarree())
        plot_map(so2t, latt, lont, extents,region, ax_region, global_min, global_max)

    plt.suptitle(f"$SO_2$ Maps for {datanames[devidx]}")
    
    # save figure
    fname = f"{savedir}so2_{datanames[devidx]}_{datestr}.png"
    savefig(fname, fig)


    # color scale 
    global_min = -0.2
    global_max = 0.2

    ### Compare two datasets ###

    diff = so2t - so2 

    # Plot set up
    regions, extents, fig, gs = plot_setup()
    # Plotting the global map as the first and largest plot
    ax_global = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    plot_map(diff, latt, lont,extents, 'Global', ax_global, global_min, global_max, add_cbar=True)

    for idx, region in enumerate(regions[1:]):
        ax_region = fig.add_subplot(gs[1, idx], projection=ccrs.PlateCarree())
        plot_map(diff, latt, lont, extents,region, ax_region, global_min, global_max)

    plt.suptitle(f"$SO_2$: {datanames[devidx]}-{datanames[refidx]}")
    
    # save figure
    fname = f"{savedir}so2_{datanames[devidx]}-{datanames[refidx]}_{datestr}.png"
    savefig(fname, fig)

