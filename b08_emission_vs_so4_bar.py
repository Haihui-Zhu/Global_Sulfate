# Bar chart showing both so2 emission and so4 concentration

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.io import loadmat
from scipy.ndimage import zoom


def savefig(fname, fig):
    fig.savefig(fname)
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: ' + fname)

def make_plot(awm_data,emission_labels,so4,sectors,savedir,fname_label,sptdots=True,sptlines=False,pattern=False):
    target_regions_1 = ['Middle East','Asia-East','Asia-South','Asia-Southeast','Africa','America-Central and South'] # regions that will be shown in the plots   
    # target_regions_2 = ['Asia Pacific', 'America-North','Australasia']
    
    existrg = so4['existrg'].flatten()
    all_regions = so4['region_name']
    all_regions = all_regions[existrg-1]
    all_regions = [item[0][0] for item in all_regions]
    
    # so4_data 
    so4_data = so4['meanmnso2'] # month x sim x region
    so4_data = np.nanmean(so4_data, axis=0)# sim x region
    
    # lines = ['-','--','-','--','--'] 
    colors = ['firebrick','goldenrod','steelblue','k']  # Colors for each site
    patterns = ['///', '\\/...']  # Patterns for each component

    # Plot setup
    n_groups = len(target_regions_1)
    bar_width = 0.25
    index = np.arange(n_groups)

    # Create Figure 1: no pattern
    fig, ax1 = plt.subplots(figsize=(12, 5))

    for r, reg in enumerate(target_regions_1): # for each group
        for j, file in enumerate(emission_labels):  # for each inventories
            bottom = 0
            for i, sec in enumerate(sectors): # stacked component
                temi = np.nanmean(awm_data[reg][file][sec])
                if i ==2 or pattern == False: # for 'other' sectors, no pattern
                    ax1.bar(index[r] + j * bar_width, temi, bar_width,
                        bottom=bottom, color=colors[j],
                        label=f'{sectors[i]}' if j == 0 and r == 0 else "")  # Add labels only for one set to avoid duplication in legend 
                else:
                    ax1.bar(index[r] + j * bar_width, temi, bar_width,
                        bottom=bottom, color=colors[j],
                        hatch=patterns[i],
                        label=f'{sectors[i]}' if j == 0 and r == 0 else "")  # Add labels only for one set to avoid duplication in legend
                       
                bottom += temi # Increment the bottom for the next stack
    # Create secondary y-axis
    if sptdots is True:
        ax2 = ax1.twinx()
        # Plot heights with markers
        for i in range(len(emission_labels)):
            for r, reg in enumerate(target_regions_1):
                for s in range(len(all_regions)):
                    if all_regions[s] == reg:
                        ax2.plot(index[r] + i * bar_width, so4_data[i+1,s], marker='o', color='forestgreen', linestyle=None)
        
        if sptlines is True:
            for r, reg in enumerate(target_regions_1):
                for s in range(len(all_regions)):
                    if all_regions[s] == reg:
                        x = [index[r]-0.5*bar_width,  index[r] + 2.5*bar_width]
                        y = [so4_data[0,s],           so4_data[0,s]]
                        ax2.plot(x, y, linestyle='--', linewidth=1.5, color='dimgrey')


    # Adding labels, title, and axes ticks
    ax1.set_ylabel(f'Annual mean $SO_2$ emission flux (kg m$^{-2}$ s$^{-1}$)',fontweight='bold')
    ax1.set_xticks(index + bar_width)
    ax1.set_xticklabels(target_regions_1)
    ax1.set_ylim(0, 1e-10)
    if sptdots is True:
        ax2.set_ylabel(r'Sulfate concentration (µg/m³)',color='forestgreen',fontweight='bold')
        ax2.set_ylim(0, 15)

    # Create component legend
    legend_components = ax1.legend(title='Sectors', loc='upper left', bbox_to_anchor=(1.08, 1))

    # Add this legend manually to the axes and create another legend for genders
    ax1.add_artist(legend_components)

    # Create emission legend
    emission_patches = [Patch(facecolor=colors[i], label=emission_labels[i]) for i in range(len(emission_labels))]
    ax1.legend(handles=emission_patches, title='Emissions', loc='upper left', bbox_to_anchor=(1.08, 0.7))

    plt.tight_layout()
    plt.show()
    
    fname = f"{savedir}/{fname_label}.png"
    savefig(fname, fig) 




save_path = '01.Emissions/'

# Load population 
data_population = loadmat('Population-0.1.mat')
npop01 = data_population['npop'].T
latp01 = data_population['tLAT'].squeeze()
lonp01 = data_population['tLON'].squeeze()

# load mask file
data_mask = loadmat('ceds_scale_2021to2018/mask_fao_ceds_01.mat')
mask_01 = data_mask['mask_region']
mlat_01 = data_mask['xlat']
mlon_01 = data_mask['xlon']
region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
region_name = [item[0] for item in region_name_array.flatten()]

###### SO2 emission inventories ####
data_label = ['ceds-2021', 'edgar-2021', 'htap-2018']
# species = ['BC','CO','NH3','NMVOC','NOx','OC','SO2']
# sectors = ['ENE','IND','TRA','RCO','WST','SHP','AGR','SOL','AVA','AWB','SLV']
sectors = ['ENE','IND','TRA','RCO','WST','SHP'] # sectors with data
Sectors = ['ENE','IND','Other'] 
# AGR SOL all empty
# AVA, AWB available in EDGAR only, and it is unlikely a big source of discrepancies
# SHP show large differences. EDGAR and HTAP SHP were not used. Need to avoid in this analysis (replace with CEDS SHP)
# SLV in HTAP and only partially. not a big deal. 


spec = 'SO2'
# making maps for each sector
# pwm_data = {region: {file: [] for file in data_label} for region in region_name} 
awm_data = {region: {file: {sector:[] for sector in Sectors} for file in data_label} for region in region_name} 

# Load and process data
for file in data_label:
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
    so2_emi ={sector: None for sector in Sectors}

    for tvar in sec_avail:
        if spec in tvar:
            parts = tvar.split('_')
            tsec = parts[1]
            if tsec in Sectors or tsec.upper() in Sectors:
                if so2_emi[tsec.upper()] is None:
                    so2_emi[tsec.upper()] = ds[tvar].values
                else:
                    so2_emi[tsec.upper()] += ds[tvar].values
                print(f'{tvar} found in {file} and added to {tsec}')
                   
            else: 
                if so2_emi['Other'] is None:
                    so2_emi['Other'] = ds[tvar].values
                else:
                    so2_emi['Other'] += ds[tvar].values
                print(f'{tvar} found in {file} and added to Other')
       

    for rgid, region in enumerate(region_name):
        region_mask = np.where(mask_01 == rgid+1, 1, 0)
        masked_population = np.where(region_mask, npop01, np.nan)
            
        for sec in Sectors:  # Assuming the third dimension is time (month)
            if len(lat)==360:
                so2_emi_sec = zoom(np.nanmean(so2_emi[sec],axis=0), 5, order=0)
            else:
                so2_emi_sec = np.nanmean(so2_emi[sec],axis=0)

            masked_emissions = np.where(region_mask, so2_emi_sec, np.nan)
            
            # calc AWM
            awm = np.nanmean(masked_emissions)
            awm_data[region][file][sec].append(awm)


# load spartan vs gchp so4
spec = 'SO4'
fname = f'06.spartan_gchp/gchp_vs_spartan_{spec}_byregion.mat'
# sfname, 'meanmnso2',''prc25so2', 'prc75so2', 'region_name', 'existrg','site_label_so4','target_sim'
so4 = loadmat(fname)

# Plotting AWM data
# fname_label = f'bar_emissions_vs_{spec}';
# fig = make_plot(awm_data,data_label,so4,Sectors,save_path,fname_label)

fname_label = f'bar_emissions_vs_{spec}_wSPT';
fig = make_plot(awm_data,data_label,so4,Sectors,save_path,fname_label,sptlines=True)

fname_label = f'bar_emissions_vs_{spec}_wPattern';
fig = make_plot(awm_data,data_label,so4,Sectors,save_path,fname_label,pattern=True,sptdots=False)
     

