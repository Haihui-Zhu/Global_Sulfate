import os
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

def make_plot(awm_data,Ylabel,emission_labels):
    target_regions = ['Africa','America-Central and South', 'America-North','Asia Pacific','Asia-East','Asia-South','Asia-Southeast','Australasia','Middle East'] # regions that will be shown in the plots   
    markers = ['-','--','--','--'];
    colors = ['k','firebrick','goldenrod','steelblue'];
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(3, 3, figure=fig)
    for i, region in enumerate(target_regions):
        ax = fig.add_subplot(gs[i // 3, i % 3])
        print(region)
        for ff, file in enumerate(emission_labels):
            ax.plot(list(range(1, 13)), awm_data[region][file], markers[ff], label=f'{file}', color=colors[ff])
            # # Calculate mean and add as text
            # mean_value = np.mean(awm_data[region][file])
            # if '2021' in file:
            #     ax.text(0.95, 0.31-0.06*ff, f'{mean_value:.2e}', transform=ax.transAxes, 
            #             horizontalalignment='right', verticalalignment='bottom', 
            #             fontsize=10,weight='bold', color=colors[ff])
            # else:
            #     ax.text(0.95, 0.31-0.06*ff, f'{mean_value:.2e}', transform=ax.transAxes, 
            #             horizontalalignment='right', verticalalignment='bottom', 
            #             fontsize=10,color=colors[ff])
            # if ff == 0:
            #     ref_value = mean_value

            # print(f'{mean_value/ref_value:.3f}')

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



def  spatial_smoothing(mapdata,lat,lon,res):

    scale = int(res/0.05)
    shape1 = len(lat)
    shape2 = len(lon)

    so2_1_r = mapdata.reshape(int(shape1/scale),scale,int(shape2/scale),scale)
    so2_1_r = so2_1_r.mean(axis=(1, 3))

    lat_r = lat.reshape(int(shape1/scale),scale)
    lat_r = lat_r.mean(axis=1)
    lon_r = lon.reshape(int(shape2/scale),scale)
    lon_r = lon_r.mean(axis=1)

    # mask out negative values
    neg_mask = so2_1_r <0
    so2_1_r[neg_mask] = 0

    return so2_1_r, lat_r, lon_r

###### SO2 emission inventories ####
data_label  = ['TROPOMI 2021', 'GCHP-CEDS 2021','GCHP-EDGAR 2021', 'GCHP-HTAP 2018'] # options: 'GCHP-CEDS 2019',
data_folder = ['CF03-SZA40-QA75','ceds_2021','edgar_2021','htap_2018'] # options:'CF03-SZA60-QA75' etc; 'ceds_2019' ,'htap_2018','edgar_2018'
flabel = 'doas-SZA40_NR_gchp-all_2021' # change here to update figure name
year = 2021 # year of tropomi data

indir   = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{data_folder[0]}/' # need to update this name ('Tesselation')

save_path = '04.Seasonality_Figures/'

# ## Load population 
# data_population = loadmat('Population-0.05.mat') # loadmat('Population-0.1.mat')
# npop = data_population['npop'].T
# latp = data_population['tLAT'].squeeze()
# lonp = data_population['tLON'].squeeze()
# # need to cut pop to have the same range as tropomi data
# lat_range = (latp  > -70) & (latp  < 70)  # if change resolution, change here too
# lon_range = (lonp  >= -180)   & (lonp <= 180)  
# npop = npop[np.ix_(lat_range, lon_range)]

res = 0.5

## load mask file
data_mask = loadmat('ceds_scale_2021to2018/mask_fao_ceds_05.mat')
mask = data_mask['mask_region']
mlat = data_mask['xlat']
mlon = data_mask['xlon']
region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
region_name = [item[0] for item in region_name_array.flatten()]
# need to cut mask to have the same range as tropomi data
lat_range = (mlat[:,0] > -70) & (mlat[:,0] < 70)   # if change resolution, change here too
lon_range = (mlon[0,:] >= -180) & (mlon[0,:] <= 180)  
mask = mask[np.ix_(lat_range, lon_range)]


## Collect PWM data for each region
pwm_data = {region: {file: [] for file in data_label} for region in region_name}
awm_data = {region: {file: [] for file in data_label} for region in region_name}


## Load and process data
# loop through dataset
for idx, file in enumerate(data_label):
    # loop through months
    for month_index in range(1,13):  
        if file == 'TROPOMI 2021' or file == 'TROPOMI 2018':
            # fname = f"{indir}Tropomi_Regrid_{year}{month_index:02d}_{data_folder[idx]}.nc"
            fname = f"{indir}/gchp_so2_cosampled_tropomi_{data_folder[1]}_noisereduced_{month_index:02d}.nc"
            if os.path.exists(fname):
                print(f'reading {fname}')
                ds = xr.open_dataset(fname)
                mapdata = ds['so2_tro_nr'].values.T
                lat  = ds['lat'].values
                lon  = ds['lon'].values
                mapdata = mapdata/2.69e16
                [mapdata, lat, lon] = spatial_smoothing(mapdata,lat,lon,res)

            else: 
                print(f'{fname} not found')
                for rgid, region in enumerate(region_name):
                    pwm_data[region][file].append(np.nan)
                    awm_data[region][file].append(np.nan)

                continue

        elif file == 'TROPOMI SO2 2019':
            rfname = f"{indir}Tropomi_Regrid_{year}{month_index:02d}.mat"
            print(f'reading {fname}')
            data = loadmat(rfname)
            mapdata, lat, lon = data['so2'], data['tlat'], data['tlon']
            mapdata = mapdata/2.69e16
            [mapdata, lat, lon] = spatial_smoothing(mapdata,lat,lon,res)

        else:
            fname = f"{indir}/gchp_so2_cosampled_tropomi_{data_folder[idx]}_noisereduced_{month_index:02d}.nc"
            if os.path.exists(fname):
                print(f'reading {fname}')
                ds = xr.open_dataset(fname)
                mapdata = ds['so2_gchp'].values.T
                lat  = ds['lat'].values
                lon  = ds['lon'].values
                mapdata = mapdata/2.69e16
                [mapdata, lat, lon] = spatial_smoothing(mapdata,lat,lon,res)
            else: 
                print(f'{fname} not found')
                for rgid, region in enumerate(region_name):
                    pwm_data[region][file].append(np.nan)
                    awm_data[region][file].append(np.nan)

                continue
        # else:
        #     fname = f"{indir}/gchp_so2_cosampled_tropomi_{data_folder[idx]}_{month_index:02d}.nc"
        #     if os.path.exists(fname):
        #         ds = xr.open_dataset(fname)
        #         mapdata = ds['so2'].values.T
        #         lat  = ds['lat'].values
        #         lon  = ds['lon'].values
        #         mapdata = mapdata/2.69e16
        #     else: 
        #         print(f'{fname} not found')
        #         for rgid, region in enumerate(region_name):
        #             pwm_data[region][file].append(np.nan)
        #             awm_data[region][file].append(np.nan)

        #         continue



        # loop through regions    
        for rgid, region in enumerate(region_name):
            region_mask = np.where(mask == rgid+1, 1, 0)

            masked_data = np.where(region_mask, mapdata, np.nan)
            # masked_population = np.where(region_mask, npop, np.nan)
            
            # calc PWM
            # total_population = np.nansum(masked_population) # Not all regions guarantee population
            # if total_population > 0:
            #     pwm = np.nansum(masked_data * masked_population) / total_population
            # else:
            #     print('total_population is 0')
            #     pwm = np.nan
            # pwm_data[region][file].append(pwm)

            # calc AWM
            awm = np.nanmean(masked_data)
            awm_data[region][file].append(awm)


## Plotting PWM data
Ylabel = r'$SO_2$ VCD (DU)'
# fig = make_plot(pwm_data,Ylabel,data_label)
# fname = f"{save_path}so2_vcd_seasonality_pwm_{flabel}_nr.png"
# savefig(fname, fig)

## Plotting AWM data
fig = make_plot(awm_data,Ylabel,data_label)
fname = f"{save_path}so2_vcd_seasonality_awm_{flabel}_nr.png"
savefig(fname, fig)

