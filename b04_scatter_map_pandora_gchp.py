import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import xarray as xr
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.io import loadmat
from scipy.ndimage import zoom
# import not installed utils
import sys
sys.path.append('functions/')
from save_fig_util import savefig
from spatial_smoothing import spatial_smoothing
from pop_weighted_mean import pwm

def concentric_map(inner, outer,lats,lons, mapdata, latm, lonm,iname,oname, site, fname,sulfate=False,so4frac=False):
    # remove nan sites:
    valid = ~np.isnan(inner) & ~np.isnan(outer)  # ~ is the bitwise NOT operator, used here to invert the boolean array
    inner = np.array(inner)
    outer = np.array(outer)
    site = np.array(site)
    inner = inner[valid]
    outer = outer[valid]
    lats = lats[valid]
    lons = lons[valid]
    site = site[valid]

    # vmin and vmax
    if sulfate is False:
        cblabel='SO2 VCD (DU)'
        vmin = -0.2
        vmax = 0.6    
        colors = [
        (0.1922, 0.2118, 0.5843),
        (0.8078, 0.8784, 0.9725),
            (1,1, 1),       # white 
        (0.9608, 0.9725, 0.6824),
        (0.9961, 0.8784, 0.5647),
        (0.9922, 0.6824, 0.3804),
        (0.9569, 0.4275, 0.2627),
        (0.8431, 0.1882, 0.1529),
        (0.6471, 0,      0.1490)
    ]
        # Create the custom colormap
        custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
        # save site mean to excel
        data = pd.DataFrame({'site': site,'lat':lats, 'lon': lons, 'pandora': inner, 'tropomi/gchp': outer})
        data.to_csv(f'{fname}.csv', index=False)  

    else:
        if so4frac is True:
            cblabel='Sulfate (%)'
            vmin = 0
            vmax = 40  
        else: # suflate
            cblabel='Sulfate (µg/m³)'
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
        
        # save site mean to excel
        data = pd.DataFrame({'site': site, 'lat':lats, 'lon': lons, 'spartan': inner, 'gchp': outer})
        data.to_csv(f'{fname}.csv', index=False) 

    fig, ax = plt.subplots(figsize=(17, 8.5), subplot_kw={'projection': ccrs.PlateCarree()})
    plt.rcParams.update({'font.size': 16}) 
    ax.set_extent([-170, 170, -60, 80], crs=ccrs.PlateCarree())  # Set longitude from -180 to 180 and latitude from -60 to 90
    ax.coastlines()
    
    lon_grid, lat_grid = np.meshgrid(lonm, latm)
    ax.pcolormesh(lon_grid, lat_grid, mapdata, transform=ccrs.PlateCarree(), cmap=custom_cmap, vmin=vmin, vmax=vmax)

    # Plot sites with more than 10 days data in red
    scatter = ax.scatter(lons, lats, c=outer, cmap=custom_cmap, marker='o', s=300,edgecolor='k',  transform=ccrs.PlateCarree(),vmin=vmin, vmax=vmax)
    ax.scatter(lons, lats, c=inner, cmap=custom_cmap, marker='o', s=120, edgecolor='k', transform=ccrs.PlateCarree(),vmin=vmin, vmax=vmax)

    # Add a legend 
    cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
    cbar = plt.colorbar(scatter, cax=cax, orientation='vertical')
    cbar.set_label(cblabel)
    # Add a legend with custom location (bottom center)
    # plt.legend(loc="lower center", bbox_to_anchor=(0.18, 0.1))
    indian_ocean_lon, indian_ocean_lat = -120, -40
    ax.scatter(indian_ocean_lon, indian_ocean_lat,c=0,cmap=custom_cmap, s=3000, edgecolor='k',vmin=vmin, vmax=vmax)
    ax.scatter(indian_ocean_lon, indian_ocean_lat,c=0,cmap=custom_cmap, s=1200 , edgecolor='k',vmin=vmin, vmax=vmax)
    ax.text(indian_ocean_lon+11, indian_ocean_lat+6, f"{oname}", ha='left', va='center', transform=ccrs.PlateCarree())
    ax.text(indian_ocean_lon+11, indian_ocean_lat, f"{iname}", ha='left', va='center', transform=ccrs.PlateCarree())
    ax.plot([indian_ocean_lon, indian_ocean_lon+10 ], [indian_ocean_lat+6, indian_ocean_lat+6], color='dimgrey')
    ax.plot([indian_ocean_lon, indian_ocean_lon+10], [indian_ocean_lat, indian_ocean_lat], color='dimgrey')

    
    plt.tight_layout()
    plt.show()
    savefig(f"{fname}.png", fig)
    

def read_spt(fname, varname):
    # save(sfname, 'meanmnso2','prc25so2', 'prc75so2', 'Site_cities','latitudes','longitudes','target_sim')
    datain = loadmat(fname)
    meanmnso2 = datain[varname]

    Site_cities = [item[0][0] for item in datain['Site_cities']]
    latitudes = datain['latitudes'].flatten()
    longitudes = datain['longitudes'].flatten()
    site_info = {'site':Site_cities,
                 'lat': latitudes,
                 'lon': longitudes,
                 }
    
    target_sim = ['spt']
    for item in datain['target_sim'][0]:
        target_sim.append(item[0])

    data = {item:{site: [] for site in Site_cities } for item in target_sim}

    for tg, item in enumerate(target_sim):
        for sid, site in enumerate(Site_cities):
            data[item][site] = meanmnso2[:,tg, sid]

    return site_info, data, target_sim 

# functions and paths
rootdir = '/pierce-scratch/haihui/Global_Sulfate/'
savedir1 = f'{rootdir}/figures/'
fname = f'{rootdir}/figures/site_mean_summary_24hr.csv'

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

extents = {
    'Global':{'lon':[-180, 180], 'lat':[-50, 60]},
    'East Asia': {'lon':[95, 145], 'lat':[ 0, 50]},
    'South Asia': {'lon':[60, 100], 'lat':[-5, 35]},
    'Middle East':{'lon': [25, 65], 'lat':[ 0, 40]}
}
    # 'North America':{'lon':[-140, -50], 'lat':[ -15, 75]},
    # 'Europe': {'lon':[-10, 40], 'lat':[20, 70]},

# ### Figure 1 - COBRA vs Pandora
# # read site mean SO2 VCD 
# df = pd.read_csv(fname)
# inner = df['pandora']
# outer = df['cobra-d']
# lats = df['lat']
# lons = df['lon']
# site = df['site']
# # reading SO2 VCD map data 
# cobrafname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/02.TROPOMI_SO2_Ref/COBRA/CORBA_Tesellation_CF03-SZA60-QA75/gchp_so2_cosampled_tropomi_ceds_2021_annual.nc' # use monthly mean so far
# ds = xr.open_dataset(cobrafname) 
# mapdata = ds['so2_tro'].values.T  
# mapdata = mapdata/2.69e16 # convert unit  
# latg = ds['lat'].values  
# long = ds['lon'].values
# mapdata, latc, lonc = spatial_smoothing(mapdata,latg,long,0.5,0.05) 
# mapdata = np.where(mask>0, mapdata, np.nan)
# # make the map
# sfname = f'{savedir1}/pandora_cobra_mask'
# iname = 'Pandora'
# oname = 'COBRA'
# concentric_map(inner, outer,lats,lons, mapdata, latc, lonc, iname, oname, site, sfname)
# for region in extents:
#     print(f'COBRA {region}:')
#     pwm(mapdata, latc, lonc, target_lat=extents[region]['lat'], target_lon=extents[region]['lon'])





### Figure 1 optional - DOAS vs Pandora
# read site mean SO2 VCD 
df = pd.read_csv(fname)
inner = df['pandora']
outer = df['doas-mn']
lats = df['lat']
lons = df['lon']
site = df['site']
qcstr =  'CF005-SZA40-QA75-vza40'
simname = 'ceds_2021'
# reading SO2 VCD map data 
doasfname = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_noisereduced_annual.nc'
ds = xr.open_dataset(doasfname) 
mapdata = ds['so2_tro'].values.T  
mapdata = mapdata/2.69e16 # convert unit  
latc = ds['lat'].values  
lonc = ds['lon'].values
mapdata, latc, lonc = spatial_smoothing(mapdata,latc,lonc,0.5,0.05) 
# mapdata = np.where(mask>0, mapdata, np.nan)
# make the map
sfname = f'{savedir1}/pandora_doas_{qcstr}_mask'
iname = 'Pandora'
oname = 'DOAS'
concentric_map(inner, outer,lats,lons, mapdata, latc, lonc, iname, oname, site, sfname)
for region in extents:
    print(f'DOAS {region}:')
    pwm(mapdata, latc, lonc, target_lat=extents[region]['lat'], target_lon=extents[region]['lon'])


# ### Figure 2 - GCHP vs Pandora
# # read site mean SO2 VCD 
# outer = df['gchp-ceds-daily']
# # reading SO2 VCD map data 
# gchpfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/02.TROPOMI_SO2_Ref/COBRA/CORBA_Tesellation_CF03-SZA60-QA75/gchp_so2_cosampled_tropomi_ceds_2021_annual.nc'
# ds = xr.open_dataset(gchpfname) 
# mapdata = ds['so2_gchp'].values.T # check variable name
# mapdata = mapdata/2.69e16 # convert unit
# latg = ds['lat'].values 
# long = ds['lon'].values
# mapdata, latg, long = spatial_smoothing(mapdata,latg,long,0.5,0.05) 
# mapdata = np.where(mask>0, mapdata, np.nan)
# # make the map
# sfname = f'{savedir1}/pandora_gchp_mask'
# oname = 'GCHP'
# concentric_map(inner, outer, lats,lons, mapdata, latg, long, iname, oname, site,sfname)
# for region in extents:
#     print(f'GCHP-CEDS {region}:')
#     pwm(mapdata, latg, long, target_lat=extents[region]['lat'], target_lon=extents[region]['lon'])

"""
### SI Figure - GCHP-EDGAR vs Pandora
# read site mean SO2 VCD 
outer = df['gchp-ceds-daily']
# reading SO2 VCD map data 
gchpfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/02.TROPOMI_SO2_Ref/COBRA/CORBA_Tesellation_CF03-SZA60-QA75/gchp_so2_cosampled_tropomi_edgar_2021_annual.nc'
ds = xr.open_dataset(gchpfname) 
mapdata = ds['so2_gchp'].values.T # check variable name
mapdata = mapdata/2.69e16 # convert unit
latg = ds['lat'].values 
long = ds['lon'].values
mapdata, latg, long = spatial_smoothing(mapdata,latg,long,0.5,0.05) 
mapdata = np.where(mask>0, mapdata, np.nan)
# make the map
sfname = f'{savedir1}/pandora_gchp-edgar_mask'
oname = 'GCHP'
concentric_map(inner, outer, lats,lons, mapdata, latg, long, iname, oname, site,sfname)
for region in extents:
    print(f'GCHP-CEDS {region}:')
    pwm(mapdata, latg, long, target_lat=extents[region]['lat'], target_lon=extents[region]['lon'])

### SI Figure - GCHP-HTAP vs Pandora
# read site mean SO2 VCD 
outer = df['gchp-ceds-daily']
# reading SO2 VCD map data 
gchpfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/02.TROPOMI_SO2_Ref/COBRA/CORBA_Tesellation_CF03-SZA60-QA75/gchp_so2_cosampled_tropomi_htap_2018_annual.nc'
ds = xr.open_dataset(gchpfname) 
mapdata = ds['so2_gchp'].values.T # check variable name
mapdata = mapdata/2.69e16 # convert unit
latg = ds['lat'].values 
long = ds['lon'].values
mapdata, latg, long = spatial_smoothing(mapdata,latg,long,0.5,0.05) 
mapdata = np.where(mask>0, mapdata, np.nan)
# make the map
sfname = f'{savedir1}/pandora_gchp-htap_mask'
oname = 'GCHP'
concentric_map(inner, outer, lats,lons, mapdata, latg, long, iname, oname, site,sfname)
for region in extents:
    print(f'GCHP-CEDS {region}:')
    pwm(mapdata, latg, long, target_lat=extents[region]['lat'], target_lon=extents[region]['lon'])
"""

"""
### Figure 3 - GCHP vs SPARTAN
# read site mean SPARTAN sulfate
spec = 'SO4';
gchpfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/4.ceds_2021/GCHP_SO2_SO4_PM25_ceds_2021_annual.nc'
fname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/gchp_vs_spartan_{spec}_bysite.mat'
site_info, data, target_sim = read_spt(fname,'meanmnso2')
inner = []
outer = []
for sid, site in enumerate(site_info['site']):
    spt =  np.nanmean(data['spt'][site])
    sim =  np.nanmean(data['ceds-2021'][site])
    inner.append(spt)
    outer.append(sim) 
lats = site_info['lat']
lons = site_info['lon']
# reading sulfate map data 
ds = xr.open_dataset(gchpfname) 
mapdata = ds['so4'].values.T # check variable name
latg = ds['latitude'].values 
long = ds['longitude'].values
# make the map
sfname = f'{savedir2}/spartan_gchp'
iname = 'SPARTAN'
oname = 'GCHP'
concentric_map(inner, outer,lats,lons, mapdata, latg, long,  iname, oname,site_info['site'], sfname, sulfate=True)
for region in extents:
    print(f'GCHP sulfate {region}:')
    pwm(mapdata, latg, long, target_lat=extents[region]['lat'], target_lon=extents[region]['lon'])


spec = 'PM2.5';
fname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/gchp_vs_spartan_{spec}_bysite.mat'
site_info, data, target_sim = read_spt(fname,'meanmnpm')
inner2 = []
outer2 = []
for sid, site in enumerate(site_info['site']):
    spt =  np.nanmean(data['spt'][site])
    sim =  np.nanmean(data['ceds-2021'][site])
    inner2.append(100*inner[sid]/spt)
    outer2.append(100*outer[sid]/sim) 
# inner = 100*inner/inner2
# outer = 100*outer/outer2
# read map data
pm25 = ds['pm25'].values.T 
mapdata = 100*mapdata/pm25
mapdata = zoom(mapdata, 2, order=0)
latg = zoom(latg, 2, order=0)
long = zoom(long, 2, order=0)
mapdata = np.where(mask_01>0, mapdata, np.nan)
# make the map
sfname = f'{savedir2}/spartan_gchp_so4frac'
iname = 'SPARTAN'
oname = 'GCHP'
concentric_map(inner2, outer2,lats,lons, mapdata, latg, long,  iname, oname,site_info['site'], sfname, sulfate=True, so4frac=True)

"""