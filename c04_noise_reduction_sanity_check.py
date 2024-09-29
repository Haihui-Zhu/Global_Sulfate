# Sanity check for noise reduction result.
# Sept. 18, 2024
# Haihui Zhu

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec

def savefig(fname, fig):
    fig.savefig(fname) 
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: '+ fname)
    return

def subplot(fig,ps,lon,lat,so2,label,vmin=-0.5,vmax=0.5,diff=False,):
    ax = fig.add_subplot(gs[ps, :], projection=ccrs.PlateCarree(central_longitude=180))

    # Plot the concentration data using the meshgrid
    if diff: 
        concentration_plot = ax.pcolormesh(lon, lat, so2, cmap='RdBu_r', transform=ccrs.PlateCarree(), shading='auto', vmin=-0.1, vmax=0.1)
    else:
        concentration_plot = ax.pcolormesh(lon, lat, so2, cmap='RdBu_r', transform=ccrs.PlateCarree(), shading='auto', vmin=vmin, vmax=vmax)

    # Add coastlines and features
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    # Add colorbar
    cax = ax.inset_axes([0.16, 0.15, 0.02, 0.3])  # Adjust these parameters as needed for placement
    if diff:
        plt.colorbar(concentration_plot, cax=cax, orientation='vertical',label='Abs Diff (DU)')
    else:
        plt.colorbar(concentration_plot, cax=cax, orientation='vertical',label='SO2 VCD (DU)')

    # adding label of data source
    x, y = 0.5, 0.05
    plt.text(x, y, label, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
             fontsize =14, fontweight = 'bold',backgroundcolor = 'w')
    
    return concentration_plot



for mn in range(5,13):
    # load original tropomi
    simname = 'ceds_2018'
    year = 2018
    qcstr =  'CF03-SZA50-QA75'
    source1 = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/Tropomi_Regrid_{year}{mn:02d}_{qcstr}.nc'
    ds = xr.open_dataset(source1) 
    so2_1 = ds['so2'].values.T # check variable name
    so2_1 = so2_1/2.69e16 # convert unit
    lat = ds['lat'].values # presumably all sources share the same lat lon
    lon = ds['lon'].values

    # load noise-reduced tropomi
    source2 = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_{mn:02d}.nc'
    ds = xr.open_dataset(source2) 
    so2_2 = ds['so2_tro'].values.T # check variable name
    so2_2 = so2_2/2.69e16 # convert unit


    # Initialize figure
    fig = plt.figure(figsize=(10, 15))
    gs = gridspec.GridSpec(3, 1, wspace=0, hspace=0, height_ratios=[1, 1, 1])

    plt.title(f'TROPOMI SO2 {year}')

    # subplot 1:  intital tropomi
    ps = 0
    label = 'original'
    subplot(fig,ps,lon,lat,so2_1,label)
    ps = 1
    label = 'noise reduced'
    subplot(fig,ps,lon,lat,so2_2,label)
    ps = 2
    label = 'middle - top'
    subplot(fig,ps,lon,lat,so2_2-so2_1,label,diff=True)

    # Display the plot
    plt.show()

    # save figure
    fname = f"./temp/noise_reduction_sanity_check_{year}_{mn}.png"
    savefig(fname, fig)