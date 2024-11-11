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



simname = 'ceds_2021'
year = 2021
qcstr =  'CF03-SZA50-QA75'

for mn in range(1,13):
    # load original tropomi
    if mn ==13:
        testfname = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_noisereduced_wO3_annual.nc'
        controlfname = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_noisereduced_annual.nc'
    else:
        testfname = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_noisereduced_wO3_{mn:02d}.nc'
        controlfname = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_noisereduced_{mn:02d}.nc'
        # infname = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_noisereduced_annual.nc'
    
    ds = xr.open_dataset(testfname) 
    # load original tropomi
    so2_t = ds['so2_tro_nr'].values.T # check variable name
    so2_t = so2_t/2.69e16 # convert unit
    # load noise-reduced tropomi
    so2_c = ds['so2_tro_nr'].values.T # check variable name
    so2_c = so2_c/2.69e16 # convert unit

    lat = ds['lat'].values # presumably all sources share the same lat lon
    lon = ds['lon'].values


    # Initialize figure 2 - 0.5 x 0.5 average
    res = 0.5
    scale = int(res/0.05)
    shape1 = len(lat)
    shape2 = len(lon)

    so2_t_r = so2_t.reshape(int(shape1/scale),scale,int(shape2/scale),scale)
    so2_t_r = so2_t_r.mean(axis=(1, 3))
    so2_c_r = so2_c.reshape(int(shape1/scale),scale,int(shape2/scale),scale)
    so2_c_r = so2_c_r.mean(axis=(1, 3))

    lat_r = lat.reshape(int(shape1/scale),scale)
    lat_r = lat_r.mean(axis=1)
    lon_r = lon.reshape(int(shape2/scale),scale)
    lon_r = lon_r.mean(axis=1)

    fig = plt.figure(figsize=(10, 15))
    gs = gridspec.GridSpec(3, 1, wspace=0, hspace=0, height_ratios=[1, 1, 1])

    plt.title(f'TROPOMI SO2 {year}-{mn} w and wo O3 filter')

    # subplot 1:  intital tropomi
    ps = 0
    label = 'control'
    subplot(fig,ps,lon_r,lat_r,so2_c_r,label)
    ps = 1
    label = 'with O3 filter'
    subplot(fig,ps,lon_r,lat_r,so2_t_r,label)
    ps = 2
    label = 'Middle - Top'
    subplot(fig,ps,lon_r,lat_r,so2_t_r-so2_c_r,label,diff=True)

    # Display the plot
    plt.show()

    # save figure
    fname = f"./temp/noise_reduction_coarse{res}_{qcstr}_{year}_{mn:02d}_w_vs_wo_O3.png"
    savefig(fname, fig)