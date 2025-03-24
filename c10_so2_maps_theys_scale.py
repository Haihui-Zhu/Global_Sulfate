# This script makes maps of global and regional SO2 and SO4

import numpy as np
import xarray as xr
import xesmf as xe
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import LinearSegmentedColormap

def plot_map(mapdata, latm, lonm, fname,fig_label,sulfate=False):

    # vmin and vmax
    if sulfate is False:
        vmin = -0.5
        vmax = 1.8    
        colors = [
        [0.20392157, 0.60392157, 1.00000000],
        [0.31372549, 0.96078431, 1.00000000],
        [0.84705882, 1.00000000, 1.00000000],
        [0.99215686, 1.00000000, 1.00000000], # white
        [0.83137255, 1.00000000, 0.92941176],
        [0.63137255, 1.00000000, 0.00000000],
        [0.98039216, 0.76078431, 0.00000000],
        [0.96470588, 0.40000000, 0.00000000],
        [0.96078431, 0.16470588, 0.61568627],
        [0.96078431, 0.00000000, 0.95686275],
        [0.88235294, 0.00000000, 0.97254902],
        [0.79607843, 0.00000000, 1.00000000],
        [0.72156863, 0.00000000, 1.00000000],
        [0.62352941, 0.00000000, 1.00000000],
        [0.54117647, 0.00000000, 1.00000000],
        [0.47450980, 0.00000000, 0.87450980]
    ]
        # Create the custom colormap
        custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    else:
        vmin = 0
        vmax = 10  
        colors = [
            (1,1, 1),       # white 
        (0.8078, 0.8784, 0.9725),
        (0.4549, 0.6784, 0.8196),
        (0.2902, 0.7529, 0.2549),
        (0.9961, 0.8784, 0.5647),
        (0.9922, 0.6824, 0.3804),
        (0.9569, 0.4275, 0.2627),
        (0.8431, 0.1882, 0.1529),
        (0.6471, 0,      0.1490)
    ]
        # Create the custom colormap
        custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

    fig, ax = plt.subplots(figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    plt.rcParams.update({'font.size': 16}) 
    ax.set_extent([-170, 170, -60, 80], crs=ccrs.PlateCarree())  # Set longitude from -180 to 180 and latitude from -60 to 90
    ax.coastlines()
    
    lon_grid, lat_grid = np.meshgrid(lonm, latm)
    mesh = ax.pcolormesh(lon_grid, lat_grid, mapdata, transform=ccrs.PlateCarree(), cmap=custom_cmap, vmin=vmin, vmax=vmax)

    # Add a legend 
    cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
    cbar = plt.colorbar(mesh, cax=cax, orientation='vertical')
    cbar.set_label('SO2 VCD (DU)')

    plt.show()
    savefig(fname, fig)
     
def savefig(fname, fig):
    fig.savefig(fname) 
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: '+ fname)
    return


### main starts here ###
simname = 'ceds_2021'
year = 2021
qcstr =  'CF005-SZA40-QA75-vza40'
# for mn in range(5,13):
source1 = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_noisereduced_annual.nc'
# source2 = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/COBRA/compiled/gchp_so2_cosampled_tropomi_ceds_2019_annual.nc'
savedir = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/'


# ### Figure 2: RadiantEarth 2019 Sep to Nov ###
# for mn in range(9,12):
#     source2 = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/COBRA/download_mannual/s5p-l3grd-so2-cobra-pbl-001-month-2019{mn:02d}01-20240702.nc'
#     ds = xr.open_dataset(source2) 
#     so2_1 = ds['SO2_column_number_density'].values # check variable name
#     # so2_1 = so2_1/2.69e16 # convert unit
#     # negative_mask = so2_1 < 0
#     # so2_1[negative_mask] = 0
#     lat = ds['latitude'].values # presumably all sources share the same lat lon
#     lon = ds['longitude'].values

#     # regrid to 0.5
#     res = 0.4
#     scale = int(res/0.05)
#     shape1 = len(lat)
#     shape2 = len(lon)

#     so2_1_r = so2_1.reshape(int(shape1/scale),scale,int(shape2/scale),scale)
#     so2_1_r = so2_1_r.mean(axis=(1, 3))

#     lat_r = lat.reshape(int(shape1/scale),scale)
#     lat_r = lat_r.mean(axis=1)
#     lon_r = lon.reshape(int(shape2/scale),scale)
#     lon_r = lon_r.mean(axis=1)

#     if 'meanso2' in locals():
#         meanso2 += so2_1_r
#     else:
#         meanso2 = so2_1_r.copy()

# meanso2 = meanso2/3
# # Plotting the global map 
# fname = f"{savedir}map_so2_cobra_Sep-Nov_2019_coarse{res}.png"
# plot_map(meanso2, lat_r, lon_r, fname,"TROPOMI-COBRA $SO_2$")

### Figure 3: BIRA-IASB L2 + tessellation 2019 Sep to Nov ###
data_arrays = []
for mn in range(9,12):
    source2 = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/COBRA/COBRA_Tesellation_CF03-SZA60-QA75/COBRA_RG_2019{mn:02d}_CF03-SZA60-QA75.nc'
    ds = xr.open_dataset(source2) 
    so2_1 = ds['so2'].values.T # check variable name # lat x lon
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
    
    data_arrays.append(so2_1_r)
    
print(len(data_arrays))
combined_array = np.stack(data_arrays, axis=0)
print(combined_array.shape)
meanso2 = np.nanmean(combined_array, axis=0) 

# Plotting the global map 
fname = f"{savedir}map_so2_cobra_Sep-Nov_2019_BIRA-IASB_coarse{res}.png"
plot_map(meanso2, lat_r, lon_r, fname,"TROPOMI-COBRA $SO_2$")


""""
### Figure 3: 2021 annual mean, BIRA-IASB L2 + tessellation
qcstr = "CF03-SZA60-QA75"
source1 = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/COBRA/CORBA_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_annual.nc'
ds = xr.open_dataset(source1) 
so2_1 = ds['so2_tro'].values.T # check variable name
so2_1 = so2_1/2.69e16 # convert unit
print(so2_1.shape)
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

print(so2_1_r.shape)
fname = f"{savedir}map_so2_cobra_2021_coarse{res}_annual_complete.png"
plot_map(so2_1_r, lat_r, lon_r, fname,"COBRA $SO_2$ [BIRA-IASB]")



### Figure 4: 2021 annual mean, RadiantEarth
qcstr = "CF03-SZA60-QA75"
source1 = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/COBRA/download_mannual/s5p-l3grd-so2-cobra-pbl-001-year-20210101-20240702.nc'
ds = xr.open_dataset(source1) 
so2_1 = ds['SO2_column_number_density'].values # check variable name
print(so2_1.shape)
lat = ds['latitude'].values # presumably all sources share the same lat lon
lon = ds['longitude'].values
# regrid to 0.5
res = 0.4
scale = int(res/0.05)
shape1 = len(lat)
shape2 = len(lon)

so2_1_r = so2_1.reshape(int(shape1/scale),scale,int(shape2/scale),scale)
so2_1_r = so2_1_r.mean(axis=(1, 3))
lat_r = lat.reshape(int(shape1/scale),scale)
lat_r = lat_r.mean(axis=1)
lon_r = lon.reshape(int(shape2/scale),scale)
lon_r = lon_r.mean(axis=1)

print(so2_1_r.shape)
fname = f"{savedir}map_so2_cobra_2021_coarse{res}_annual_radiantearth.png"
plot_map(so2_1_r, lat_r, lon_r, fname,"COBRA $SO_2$ [RadiantEarth]")
"""