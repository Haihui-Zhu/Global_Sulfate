# compares sectoral emissions of SO2

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
from scipy.io import loadmat
from scipy.ndimage import zoom


def read_emission_data(file_path):
    data = loadmat(file_path)
    return data['emi_so2'], data['lat'], data['lon']

def savefig(fname, fig):
    fig.savefig(fname)
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: ' + fname)

def make_plot(awm_data,Ylabel,emission_labels):
    target_regions = ['Africa','America-Central and South', 'America-North','Asia Pacific','Asia-East','Asia-South','Asia-Southeast','Australasia','Middle East'] # regions that will be shown in the plots   
    markers = ['b-','b--','r-','r--','g--'];
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(3, 3, figure=fig)
    for i, region in enumerate(target_regions):
        ax = fig.add_subplot(gs[i // 3, i % 3])
        for ff, file in enumerate(emission_labels):
            ax.plot(list(range(1, 13)), awm_data[region][file], markers[ff], label=f'{file}')
        ax.set_title(region)
        ax.set_xticks(range(1, 13, 2))
        ax.set_xticklabels(['Jan', 'Mar','May', 'Jul', 'Sep', 'Nov'])
        ax.set_ylim(bottom=0)
        if i%3 == 0:
            ax.set_ylabel(Ylabel)
        ax.grid(True)
        if i == 0:
            ax.legend()

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
data_label = ['CEDS-2021', 'CEDS-2018','EDGAR-2021','EDGAR-2018', 'HTAP-2018']
species = ['BC','CO','NH3','NMVOC','NOx','OC','SO2']
sectors = ['ENE','IND','TRA','RCO','SLV','AGR','AVA','SHP','WST','AWB','SOL']
# sectors = ['ENE','TRA','IND']


spec = 'SO2'
# Load and process data
for sec in sectors:
    # making maps for each sector
    pwm_data = {region: {file: [] for file in data_label} for region in region_name} 
    awm_data = {region: {file: [] for file in data_label} for region in region_name} 

    for file in data_label:
        if file == 'CEDS-2021':
            fname = f"./CEDS-2024-Gridded/{spec}-em-anthro_CMIP_CEDS_2021.nc"
            tvar = f"{spec}_{sec.lower()}"

        elif file == 'CEDS-2018':
            fname = f"./CEDS-2021/2018/{spec}-em-anthro_CMIP_CEDS_2018.nc"
            tvar =  f"{spec}_{sec.lower()}"

        elif file == 'EDGAR-2021':
            fname = f"./EDGAR/EDGARv81/v8.1_{spec}_2021_scaled_w_v61_0.1x0.1.nc"
            tvar = f"{spec}_{sec}"

        elif file == 'EDGAR-2018':
            fname = f"./EDGAR/EDGARv61/EDGARv6.1_{spec}_2018.0.1x0.1.nc"
            tvar = f"{spec}_{sec}"

        elif file == 'HTAP-2018':
            fname = f"./HTAPv3/2018/HTAPv3_{spec}_0.1x0.1_2018.nc"
            tvar = f"{spec}_{sec}"

        ds = xr.open_dataset(fname)
        lat  = ds['lat'].values
        lon  = ds['lon'].values
        sec_avail = ds.data_vars
        
        if tvar in sec_avail:
            so2_emi = ds[tvar].values

            for rgid, region in enumerate(region_name):
                region_mask = np.where(mask_01 == rgid+1, 1, 0)
                masked_population = np.where(region_mask, npop01, np.nan)
                 
                for month_index in range(so2_emi.shape[0]):  # Assuming the third dimension is time (month)
                    if len(lat)==360:
                        so2_emi_mn = zoom(so2_emi[month_index, :, :], 5, order=0)
                    else:
                        so2_emi_mn = so2_emi[month_index, :, :]
                        
                    masked_emissions = np.where(region_mask, so2_emi_mn, np.nan)
                   
                    # calc PWM
                    total_population = np.nansum(masked_population) # Not all regions guarantee population
                    if total_population > 0:
                        pwm = np.nansum(masked_emissions * masked_population) / total_population
                    else:
                        print(f'total_population is 0 for {region}')
                        pwm = np.nan
                    pwm_data[region][file].append(pwm)

                    # calc AWM
                    awm = np.nanmean(masked_emissions)
                    awm_data[region][file].append(awm)

        else:
            print(f"{spec} {sec} not found for {file}")
            for rgid, region in enumerate(region_name):
                for month_index in range(1,13):
                    pwm_data[region][file].append(np.nan)
                    awm_data[region][file].append(np.nan)
    
    # Plotting PWM data
    Ylabel = r'$SO_2$ emission (kg m$^{-2}$ s$^{-1}$)'
    fig = make_plot(pwm_data,Ylabel,data_label)
    fname = f"{save_path}{spec}_emissions_seasonality_pwm_{sec}.png"
    savefig(fname, fig)

    # Plotting AWM data
    fig = make_plot(awm_data,Ylabel,data_label)
    fname = f"{save_path}{spec}_emissions_seasonality_{sec}.png"
    savefig(fname, fig)         

