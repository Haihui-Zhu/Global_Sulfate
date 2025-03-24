# making scatter plots comparing SO2 from various sources
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.io import loadmat
# import not installed utils
import sys
sys.path.append('functions/')
from save_fig_util import savefig
from spatial_smoothing import spatial_smoothing
from get_stats import get_stats

def histo_2d(map1, map2, latc, lonc, name1, name2, sfname):
    fz = 14 # font size 
    colors = [
        (1,1, 1),       # white 
        # (0.1922, 0.2118, 0.5843),
        (0.2706, 0.4588, 0.7059),
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
    
    statstrs, slope, intercept, N, x, y =  get_stats(map1,map2)
    
    fig, ax = plt.subplots(figsize=(5, 5))
    # Define bin size
    bin_size = 0.01 
    x_bins = np.arange(min(x), max(x) + bin_size, bin_size)
    y_bins = np.arange(min(y), max(y) + bin_size, bin_size)

    # Create a 2D histogram with specified bin edges
    h = plt.hist2d(x, y, bins=[x_bins, y_bins], cmap=custom_cmap)
    # Set x and y axis limits
    plt.xlim(-0.2, 0.6)  # Set the x-axis limits
    plt.ylim(-0.2, 0.6)  # Set the y-axis limits
    # Labeling
    plt.xlabel(f'{name1} SO2 VCD (DU)',fontsize=fz)
    plt.ylabel(f'{name2} SO2 VCD (DU)',fontsize=fz)

    # Plot the regression line
    x_values = np.linspace(min(x), max(x), 100)
    y_values = intercept + slope * x_values
    plt.plot(x_values, y_values, 'r-', linewidth=1.5)  # Red line for regression
    plt.plot(x_values, x_values, '--', color='grey', linewidth=1.5)  #  1:1
    
    plt.text(0.1, 0.9, statstrs, ha='left', va='center', transform=plt.gca().transAxes, fontsize=fz)
    
    # Create an inset axis for the colorbar inside the main plot
    cbar_ax = fig.add_axes([0.7, 0.15, 0.03, 0.4])  # Example position and size
    cb = fig.colorbar(h[3], cax=cbar_ax)
    cb.set_label('counts in bin',fontsize=fz)
    
    plt.subplots_adjust(left=0.15) 

    plt.show()
    fname = f'{sfname}.png'
    savefig(fname, fig)



# functions and paths
rootdir = '/storage1/fs1/rvmartin/Active'
savedir = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/'
fname = f'{rootdir}/Shared/haihuizhu/SO2_Pandora//figures_2020-2022/site_mean_summary_24hr.csv'

# load mask file
data_mask = loadmat('ceds_scale_2021to2018/mask_fao_ceds_05.mat')
mask_01 = data_mask['mask_region']
mlat_01 = data_mask['xlat']
mlon_01 = data_mask['xlon']
region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
region_name = [item[0] for item in region_name_array.flatten()]
# need to cut mask to have the same range as tropomi data
lat_range = (mlat_01[:,0] > -70) & (mlat_01[:,0] < 70)   # if change resolution, change here too
lon_range = (mlon_01[0,:] >= -180) & (mlon_01[0,:] <= 180)  
mask = mask_01[np.ix_(lat_range, lon_range)]


# reading COBRA SO2 VCD map data 
cobrafname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/02.TROPOMI_SO2_Ref/COBRA/CORBA_Tesellation_CF03-SZA60-QA75/gchp_so2_cosampled_tropomi_ceds_2021_annual.nc' # use monthly mean so far
ds = xr.open_dataset(cobrafname) 
mapdata = ds['so2_tro'].values.T  
mapdata = mapdata/2.69e16 # convert unit  
latg = ds['lat'].values  
long = ds['lon'].values
mapdata, latc, lonc = spatial_smoothing(mapdata,latg,long,0.5,0.05) 
cobra = np.where(mask>0, mapdata, np.nan)

# reading DOAS SO2 VCD map data 
qcstr =  'CF02-SZA40-QA75'
simname = 'ceds_2021'
doasfname = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_noisereduced_annual.nc'
ds = xr.open_dataset(doasfname) 
mapdata = ds['so2_tro'].values.T  
mapdata = mapdata/2.69e16 # convert unit  
latc = ds['lat'].values  
lonc = ds['lon'].values
mapdata, latc, lonc = spatial_smoothing(mapdata,latc,lonc,0.5,0.05) 
doas = np.where(mask>0, mapdata, np.nan)

# reading SO2 VCD map data 
gchpfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/02.TROPOMI_SO2_Ref/COBRA/CORBA_Tesellation_CF03-SZA60-QA75/gchp_so2_cosampled_tropomi_ceds_2021_annual.nc'
ds = xr.open_dataset(gchpfname) 
mapdata = ds['so2_gchp'].values.T # check variable name
mapdata = mapdata/2.69e16 # convert unit
latg = ds['lat'].values 
long = ds['lon'].values
mapdata, latg, long = spatial_smoothing(mapdata,latg,long,0.5,0.05) 
gchp = np.where(mask>0, mapdata, np.nan)


### Figure 1 - COBRA vs DOAS
name1 = 'COBRA'
name2 = 'DOAS'
sfname = f'{savedir}/{name1}_vs_{name2}'
histo_2d(cobra, doas,  latc, lonc, name1, name2, sfname)

### Figure 2 - GCHP vs DOAS
name1 = 'DOAS'
name2 = 'GCHP'
sfname = f'{savedir}/{name1}_vs_{name2}'
histo_2d(doas, gchp, latc, lonc, name1, name2, sfname)


### Figure 1 - GCHP vs COBRA
name1 = 'COBRA'
name2 = 'GCHP'
sfname = f'{savedir}/{name1}_vs_{name2}'
histo_2d(cobra,  gchp, latc, lonc, name1, name2, sfname)
