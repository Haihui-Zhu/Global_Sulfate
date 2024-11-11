# Bar chart showing both so2 emission and so4 concentration

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.ndimage import zoom
import netCDF4 as nc


def savefig(fname, fig):
    fig.savefig(fname)
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: ' + fname)

def make_plot(region_emisum,emission_labels,savedir,fname_label):
    target_regions_1 = ['Middle East','Asia-East','Asia-South','Asia-Southeast','Africa','America-Central and South'] # regions that will be shown in the plots   
    # target_regions_2 = ['Asia Pacific', 'America-North','Australasia']
    
    # lines = ['-','--','-','--','--'] 
    colors = ['firebrick','goldenrod','steelblue','k']  # Colors for each site

    # Plot setup
    bar_width = 0.2
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
    ax1.set_ylabel(f'$SO_2$ emission (kton/year)',fontweight='bold')
    ax1.set_xticks(index + 1.5*bar_width)
    ax1.set_xticklabels(target_regions_1)
    # ax1.set_ylim(0, 1e-10)
    
    ax1.legend(title='Emissions', loc='upper left', bbox_to_anchor=(1.02, 0.7),fontsize='large')

    plt.tight_layout()
    plt.show()
    
    fname = f"{savedir}/{fname_label}.png"
    savefig(fname, fig) 


def find_nearest_index(grid, value):
    idx = np.abs(grid - value).argmin()
    return idx


save_path = '01.Emissions/'

# Load population 
data_population = loadmat('Population-0.1.mat')
npop01 = data_population['npop'].T
latp01 = data_population['tLAT'].squeeze()
lonp01 = data_population['tLON'].squeeze()

# load mask file
data_mask = loadmat('ceds_scale_2021to2018/mask_fao_ceds_01_emis.mat')
mask_01 = data_mask['mask_region']
mask_01 = mask_01 - 1 # so that id can be use as index
mlat_01 = data_mask['xlat']
mlon_01 = data_mask['xlon']
region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
region_name = [item[0] for item in region_name_array.flatten()]

# Area data
fname = 'utils/pixel_area_0.1.nc'
ds = xr.open_dataset(fname)
lat01  = ds['lat'].values
lon01  = ds['lon'].values
area01 = ds['area'].values


###### SO2 emission inventories ####
data_label = ['ceds-2021', 'edgar-2021', 'htap-2018','top-down']
# species = ['BC','CO','NH3','NMVOC','NOx','OC','SO2']
# sectors = ['ENE','IND','TRA','RCO','WST','SHP','AGR','SOL','AVA','AWB','SLV']
sectors = ['ENE','IND','TRA','RCO','WST','SHP'] # sectors with data
Sectors = ['ENE','IND','Other'] 
# AGR SOL all empty
# AVA, AWB available in EDGAR only, and it is unlikely a big source of discrepancies (not used in GCHP either)
# SHP show differences. But its magnitude is much smaller than other ENE and IND.
# SLV in HTAP and only partially. not a big deal. 
daysinmn = [31, 28 ,31, 30, 31 ,30 ,31, 31 ,30, 31, 30 ,31]

spec = 'SO2'
# making maps for each sector
region_emisum = {region: {file: [] for file in data_label} for region in region_name} 


# Load and process data
for file in data_label:
    if file == 'top-down':
        fname = '01.Emissions/MSAQSO2L4_2005-2023_v02-00_20240420.nc'
        ds = nc.Dataset(fname, 'r')
        Catalogue = ds.groups['Catalogue']
        # Useful variables:
        emission = Catalogue.variables['Emissions'][:]
        lat  = ds['Latitude'][:]
        lon  = ds['Longitude'][:]
        sourcetype = ds['SourceType'][:]
        year = ds['Time'][:]

        for yid, yr in enumerate(year):
            if yr == 2021.5:
                temi = emission[:,yid]
        
        region_list = np.zeros(len(lat),dtype='int')
        # if not volcano, assign to the region
        for tid, tlat in enumerate(lat):
            if sourcetype[tid] != 'Volcano    ':
                tlon = lon[tid]
                lat_idx = find_nearest_index(mlat_01[:,1], tlat)
                lon_idx = find_nearest_index(mlon_01[1,:], tlon)
                if abs(mask_01[lat_idx, lon_idx])>=0 : # not a nan
                    region_list[tid] = int(mask_01[lat_idx, lon_idx])
                    region = region_name[region_list[tid]]
                    region_emisum[region][file].append(temi[tid])
        
        # try again to for oceanic sites
        for tid, tlat in enumerate(lat):
            if sourcetype[tid] != 'Volcano    ' and region_list[tid] == 0: # not a volcano and not assigned region
                tlon = lon[tid]
                dis = (lat-tlat)**2 + (lon-tlon)**2
                dis[dis==0] = 999999
                print(min(dis))
                idx = np.abs(dis).argmin()
                print(dis[idx])
                region_list[tid] = region_list[idx]
                print(f'region idx = {region_list[idx]}')
                print(f'lat={tlat:.2f} lon={tlon:.2f}, source type = {sourcetype[tid]} assigned to {region_name[region_list[idx]]}')

        # sum regional total emission (ktons per year):
        for region in region_name:
            region_emisum[region][file] = np.nansum(region_emisum[region][file])
       
        print("\n".join(f"{region}: {region_emisum[region][file]}" for region in region_emisum))

    else:
        if file == 'ceds-2021':
            fname = f"./CEDS-2024-Gridded/{spec}-em-anthro_CMIP_CEDS_2021.nc"

        elif file == 'ceds-2018':
            fname = f"./CEDS-2021/2018/{spec}-em-anthro_CMIP_CEDS_2018.nc"

        elif file == 'edgar-2021':
            fname = f"./EDGAR/EDGARv81/v8.1_{spec}_2021_scaled_w_v61_0.1x0.1.nc"

        elif file == 'edgar-2018':
            fname = f"./EDGAR/EDGARv61/EDGARv6.1_{spec}_2018.0.1x0.1.nc"

        elif file == 'htap-2018':
            # fname = f"./HTAPv3/2018/HTAPv3_{spec}_0.1x0.1_2018.nc"
            fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/HTAPv3_{spec}_0.1x0.1_2018.nc" # file from pfe

        ds = xr.open_dataset(fname)
        lat  = ds['lat'].values
        lon  = ds['lon'].values
        sec_avail = ds.data_vars

        # collect only so2 emissions from desired sectors that used in simulations
        so2_emi = None 
        so2_emi_yr = None
        for tvar in sec_avail:
            if spec in tvar:
                parts = tvar.split('_')
                tsec = parts[1]
                if tsec in sectors or tsec.upper() in sectors:
                    if so2_emi is None:
                        so2_emi = ds[tvar].values
                    else:
                        so2_emi += ds[tvar].values
                    print(f'{tvar} found in {file} and added')
        
        # sum through the time dimension: 
        for month_index in range(so2_emi.shape[0]):  # Assuming the third dimension is time (month)
            if len(lat)==360:
                so2_emi_mn = zoom(so2_emi[month_index, :, :], 5, order=0)
            else:
                so2_emi_mn = so2_emi[month_index, :, :]
            # emission in unit of kg m-2 s-1
            # convert to kg/m2/year:
            if so2_emi_yr is None:
                so2_emi_yr = so2_emi_mn*daysinmn[month_index]*24*3600
            else:
                so2_emi_yr += so2_emi_mn*daysinmn[month_index]*24*3600

        # now apply pixel areas, convert to kg/pixel/year
        so2_emi_yr = so2_emi_yr*area01*1e6 # 1e6 to convert area km2 to m2

        for rgid, region in enumerate(region_name):
            region_mask = np.where(mask_01 == rgid, 1, 0)
            masked_emissions = np.where(region_mask, so2_emi_yr, np.nan)
            
            # calc AWM
            awm = np.nansum(masked_emissions)/1e6 # convert to kton/year 
            region_emisum[region][file].append(awm)


# load spartan vs gchp so4
spec = 'SO4'
# fname = f'06.spartan_gchp/gchp_vs_spartan_{spec}_byregion.mat'
# # sfname, 'meanmnso2',''prc25so2', 'prc75so2', 'region_name', 'existrg','site_label_so4','target_sim'
# so4 = loadmat(fname)

# Plotting AWM data
fname_label = f'bar_emissions_bu_vs_td';
fig = make_plot(region_emisum,data_label,save_path,fname_label)
