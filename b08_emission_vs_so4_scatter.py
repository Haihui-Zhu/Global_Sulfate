# access correlation between so2 emission and so4 concentration

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
from scipy.io import loadmat
from scipy.ndimage import zoom

def savefig(fname, fig):
    fig.savefig(fname)
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: ' + fname)

def make_plot(awm_data,emission_labels,so4):
    # awm_data[region][file]
    # emission_labels --> target_sim

    existrg = so4['existrg'].flatten()
    target_regions = so4['region_name']
    target_regions = target_regions[existrg-1]
    target_regions = [item[0][0] for item in target_regions]
    # target_regions = ['Africa','America-Central and South', 'America-North','Asia Pacific','Asia-East','Asia-South','Asia-Southeast','Australasia','Middle East'] # regions that will be shown in the plots   
    
    # so4_data 
    so4_data = so4['meanmnso2'] # month x sim x region
    
    # lines = ['-','--','-','--','--']
    # colors = ['firebrick','firebrick','goldenrod','goldenrod','steelblue']    
    colors = ['firebrick','goldenrod','steelblue','k']  # Colors for each site
    markers = ['o', '^', 's', 'p', '*', '+', 'x', 'D', 'v', '<', 'h']  # Markers for each month

    # Plotting
    fig = plt.figure(figsize=(8.5, 5))
    for ridx, region in enumerate(target_regions):
        for fidx, file in enumerate(emission_labels):
            tso4 = np.nanmean(so4_data[:,fidx+1,ridx])
            eso2 = np.nanmean(awm_data[region][file])
            plt.plot(eso2, tso4,  color=colors[fidx], marker=markers[ridx])
        
        # mark SPARTAN measurement level:
        tso4 = np.nanmean(so4_data[:,0,ridx])
        plt.plot(-1e-11, tso4, color=colors[fidx+1], marker=markers[ridx])
        

    plt.xlabel(f'Annual mean $SO_2$ emission (kg m$^{-2}$ s$^{-1}$)')
    plt.ylabel(r'Sulfate concentration (µg/m³)')

    # Creating a custom legend for sites and another for months
    emission_labels.append('SPARTAN')
    from matplotlib.lines import Line2D
    legend_elements_sites = [Line2D([0], [0], color=color, marker='o', linestyle='', label=file.upper())
                            for color, file in zip(colors,emission_labels)]
    legend_elements_regions = [Line2D([0], [0], color='dimgrey', marker=marker, linestyle='', label=region)
                            for marker, region in zip(markers, target_regions)]

    legend_sites = plt.legend(handles=legend_elements_sites , bbox_to_anchor=(1.05, 1), loc='upper left',title="Data Souce")
    legend_regions = plt.legend(handles=legend_elements_regions, bbox_to_anchor=(1.05, 0.7), loc='upper left',title="Region")
    plt.gca().add_artist(legend_sites)

    plt.grid(True)
    plt.tight_layout()
    plt.show()
    return fig


save_path = '04.Seasonality_Figures/'

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
# AGR SOL all empty
# AVA, AWB available in EDGAR only, and it is unlikely a big source of discrepancies
# SHP show large differences. EDGAR and HTAP SHP were not used. Need to avoid in this analysis (replace with CEDS SHP)
# SLV in HTAP and only partially. not a big deal. 


spec = 'SO2'
# making maps for each sector
# pwm_data = {region: {file: [] for file in data_label} for region in region_name} 
awm_data = {region: {file: [] for file in data_label} for region in region_name} 

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
    so2_emi = None 
    for tvar in sec_avail:
        if spec in tvar:
            for sec in sectors:
                if sec in tvar or sec.lower() in tvar:
                    if so2_emi is None:
                        so2_emi = ds[tvar].values
                    else:
                        so2_emi += ds[tvar].values
                    print(f'{tvar} found in {file} and added')
                              

    for rgid, region in enumerate(region_name):
        region_mask = np.where(mask_01 == rgid+1, 1, 0)
        masked_population = np.where(region_mask, npop01, np.nan)
            
        for month_index in range(so2_emi.shape[0]):  # Assuming the third dimension is time (month)
            if len(lat)==360:
                so2_emi_mn = zoom(so2_emi[month_index, :, :], 5, order=0)
            else:
                so2_emi_mn = so2_emi[month_index, :, :]

            masked_emissions = np.where(region_mask, so2_emi_mn, np.nan)
            
            # # calc PWM
            # total_population = np.nansum(masked_population) # Not all regions guarantee population
            # if total_population > 0:
            #     pwm = np.nansum(masked_emissions * masked_population) / total_population
            # else:
            #     print(f'total_population is 0 for {region}')
            #     pwm = np.nan
            # pwm_data[region][file].append(pwm)

            # calc AWM
            awm = np.nanmean(masked_emissions)
            awm_data[region][file].append(awm)


# # Plotting PWM data
# Ylabel = r'$SO_2$ emission (kg m$^{-2}$ s$^{-1}$)'
# fig = make_plot(pwm_data,Ylabel,data_label)
# fname = f"{save_path}{spec}_emissions_seasonality_pwm.png"
# savefig(fname, fig)

# load spartan vs gchp so4
spec = 'SO4'
fname = f'06.spartan_gchp/gchp_vs_spartan_{spec}_byregion.mat'
# save(sfname, 'meanmnso2',''prc25so2', 'prc75so2', 'region_name', 'existrg','site_label_so4','target_sim')
so4 = loadmat(fname)

# Plotting AWM data
fig = make_plot(awm_data,data_label,so4)
fname = f"{save_path}scatter_emissions_vs_{spec}.png"
savefig(fname, fig)         

