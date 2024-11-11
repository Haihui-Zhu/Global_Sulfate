# This script makes maps of global and regional SO2 and SO4

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.interpolate import griddata

def make_plot(region_emisum,emission_labels):
    target_regions_1 = ['Middle East','Asia-East','Asia-South','Asia-Southeast','Africa','America-Central and South'] # regions that will be shown in the plots   
    # target_regions_2 = ['Asia Pacific', 'America-North','Australasia']
    
    # lines = ['-','--','-','--','--'] 
    colors = ['firebrick','goldenrod','steelblue','dimgrey','darkgrey']  # Colors for each site

    # Plot setup
    bar_width = 0.16
    n_groups = len(target_regions_1)
    index = np.arange(n_groups)

    # Create Figure 1: no pattern
    fig, ax1 = plt.subplots(figsize=(12, 5))

    for r, reg in enumerate(target_regions_1): # for each group
        for j, file in enumerate(emission_labels):  # for each inventories
            temi = region_emisum[reg][file]
            ax1.bar(index[r] + j * bar_width, temi, bar_width, color=colors[j],
                label=[f'{file}' if r == 0 else ""])  # Add labels only for one set to avoid duplication in legend 
        
    # Adding labels, title, and axes ticks
    ax1.set_ylabel(f'$SO_2$ VCD (DU)',fontweight='bold')
    ax1.set_xticks(index + 1.5*bar_width)
    ax1.set_xticklabels(target_regions_1)
    # ax1.set_ylim(0, 1e-10)
    
    ax1.legend(loc='upper left', bbox_to_anchor=(1.02, 0.7),fontsize='large')

    plt.tight_layout()
    plt.show()
    return fig
    

def savefig(fname, fig):
    fig.savefig(fname) 
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: '+ fname)
    return


def  spatial_smoothing(mapdata,lat,lon,resout,resin):

    scale = int(resout/resin)
    shape1 = len(lat)
    shape2 = len(lon)

    so2_1_r = mapdata.reshape(int(shape1/scale),scale,int(shape2/scale),scale)
    so2_1_r = so2_1_r.mean(axis=(1, 3))

    lat_r = lat.reshape(int(shape1/scale),scale)
    lat_r = lat_r.mean(axis=1)
    lon_r = lon.reshape(int(shape2/scale),scale)
    lon_r = lon_r.mean(axis=1)

    # # mask out negative values
    # neg_mask = so2_1_r <0
    # so2_1_r[neg_mask] = 0

    return so2_1_r, lat_r, lon_r



### main starts here ###
data_label = ['gchp-ceds', 'gchp-edgar', 'gchp-htap','DOAS','COBRA']
year = 2021
qcstr =  'CF005-SZA40-QA75-vza40'

save_path = '03.TROPOMI_SO2_Result/'

# load mask file
data_mask = loadmat('ceds_scale_2021to2018/mask_fao_ceds_05.mat')
mask_01 = data_mask['mask_region']
mask_01 = mask_01 - 1 # so that id can be use as index
mlat_01 = data_mask['xlat']
mlon_01 = data_mask['xlon']
region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
region_name = [item[0] for item in region_name_array.flatten()]

# need to cut mask to have the same range as tropomi data
lat_range = (mlat_01[:,0] > -70) & (mlat_01[:,0] < 70)   # if change resolution, change here too
lon_range = (mlon_01[0,:] >= -180) & (mlon_01[0,:] <= 180)  
mask = mask_01[np.ix_(lat_range, lon_range)]

region_emisum = {region: {file: [] for file in data_label} for region in region_name} 


for file in data_label:

    if file == 'gchp-ceds':
        fname = f"./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_ceds_2021_noisereduced_annual.nc"

    elif file == 'gchp-edgar':
        fname = f"./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_edgar_2021_noisereduced_annual.nc"

    elif file == 'gchp-htap':
        fname = f"./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_htap_2018_noisereduced_annual.nc"

    elif file == 'DOAS':
        fname = f"./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_ceds_2021_noisereduced_annual.nc"

    elif file == 'COBRA':
        fname = f'02.TROPOMI_SO2_Ref/COBRA/compiled/gchp_so2_cosampled_tropomi_ceds_2019_annual.nc'

    
    ds = xr.open_dataset(fname) # cobra official monthly mean so
    if file == 'COBRA':
        so2 = ds['so2_pbl'].values.T # check variable name
        lat = ds['latitude'].values # presumably all sources share the same lat lon
        lon = ds['longitude'].values
        so2 = so2/2.69e16 # convert unit

        resout,resin = 0.4, 0.05
        [so2, lat, lon] = spatial_smoothing(so2,lat,lon,resout,resin)

        # Flatten the meshgrid and Array A for interpolation
        mask_grid_f = np.vstack([mlat_01.ravel(), mlon_01.ravel()]).T
        mask_f = mask_01.ravel()

        # Create a new meshgrid for the target resolution (Array B)
        grid_x_new, grid_y_new = np.meshgrid(lon, lat)

        # Perform nearest neighbor interpolation
        mask_cobra = griddata(mask_grid_f, mask_f, (grid_y_new, grid_x_new), method='nearest')

        for rgid, region in enumerate(region_name):
            region_mask = np.where(mask_cobra == rgid, 1, 0)
            masked_data = np.where(region_mask, so2, np.nan)
            
            # calc AWM
            awm = np.nanmean(masked_data)
            region_emisum[region][file].append(awm)

    else:
        if file == 'DOAS':
            so2 = ds['so2_tro'].values.T # check variable name
        else:
            so2 = ds['so2_gchp'].values.T # check variable name

        so2 = so2/2.69e16 # convert unit
        lat = ds['lat'].values # presumably all sources share the same lat lon
        lon = ds['lon'].values
        resout,resin = 0.5, 0.05
        [so2, lat, lon] = spatial_smoothing(so2,lat,lon,resout,resin)

        for rgid, region in enumerate(region_name):
            region_mask = np.where(mask == rgid, 1, 0)
            masked_data = np.where(region_mask, so2, np.nan)
            
            # calc AWM
            awm = np.nanmean(masked_data)
            region_emisum[region][file].append(awm)


# Plotting AWM data
fig = make_plot(region_emisum,data_label)

# save figure
fname = f"{save_path}bar_so2_sim_vs_tropomi.png"
savefig(fname, fig)


