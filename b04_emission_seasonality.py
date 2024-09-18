import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
from scipy.io import loadmat

def read_emission_data(file_path):
    data = loadmat(file_path)
    return data['emi_so2'], data['lat'], data['lon']

def savefig(fname, fig):
    fig.savefig(fname)
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: ' + fname)

def make_plot(awm_data,Ylabel,emission_labels,regions):
    markers = ['b-','r-','g-'];
    fig = plt.figure(figsize=(10, 15))
    gs = gridspec.GridSpec(5, 3, figure=fig)
    for i, region in enumerate(regions):
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
npop = data_population['npop'].T
latp = data_population['tLAT'].squeeze()
lonp = data_population['tLON'].squeeze()

# load mask file
data_mask = loadmat('ceds_scale_2021to2018/mask_fao_ceds_01.mat')
mask = data_mask['mask_region']
mlat = data_mask['xlat']
mlon = data_mask['xlon']
region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
region_name = [item[0] for item in region_name_array.flatten()]

###### SO2 emission inventories ####
emission_path = '01.Emissions/'
data_label = ['CEDS-2021', 'EDGARv6.1', 'HTAPv3']

# Collect PWM data for each region
pwm_data = {region: {file: [] for file in data_label} for region in region_name}
awm_data = {region: {file: [] for file in data_label} for region in region_name}

# Load and process data
for file in data_label:
    so2_emi, lat, lon = read_emission_data(f"{emission_path}{file}_seasonal_allsectors.mat")
    
    for rgid, region in enumerate(region_name):
        region_mask = np.where(mask == rgid+1, 1, 0)
        for month_index in range(so2_emi.shape[2]):  # Assuming the third dimension is time (month)
            masked_emissions = np.where(region_mask, so2_emi[:, :, month_index], np.nan)
            masked_population = np.where(region_mask, npop, np.nan)
            
            # calc PWM
            total_population = np.nansum(masked_population) # Not all regions guarantee population
            if total_population > 0:
                pwm = np.nansum(masked_emissions * masked_population) / total_population
            else:
                print('total_population is 0')
                pwm = np.nan
            pwm_data[region][file].append(pwm)

            # calc AWM
            awm = np.nanmean(masked_emissions)
            awm_data[region][file].append(awm)

# Plotting PWM data
Ylabel = r'$SO_2$ emission (kg m$^{-2}$ s$^{-1}$)'
fig = make_plot(pwm_data,Ylabel,data_label,region_name)
fname = f"{save_path}so2_emissions_seasonality_pwm.png"
savefig(fname, fig)

# Plotting AWM data
fig = make_plot(awm_data,Ylabel,data_label,region_name)
fname = f"{save_path}so2_emissions_seasonality.png"
savefig(fname, fig)

