import numpy as np
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

def make_plot(awm_data,Ylabel,emission_files,regions):
    markers = ['b-','b-.','r-','g-'];
    fig = plt.figure(figsize=(8,4))
    gs = gridspec.GridSpec(1, 2, figure=fig)
    for i, region in enumerate(regions):
        ax = fig.add_subplot(gs[0,1])
        for ff, file in enumerate(emission_files):
            ax.plot(list(range(1, 13)), awm_data[region][file], markers[ff], label=f'{file}')
        ax.set_title(region)
        ax.set_xticks(range(1, 13, 2))
        ax.set_xticklabels(['Jan', 'Mar','May', 'Jul', 'Sep', 'Nov'])
        ax.set_ylabel(Ylabel)
        ax.grid(True)
        if i == 0:
            ax.legend()

    plt.tight_layout()
    plt.show()
    return fig

def make_plot_2(awm_data,Ylabel1,Ylabel2,emission_files,regions):
    markers = ['b-','b-.','r-.'];
    fig = plt.figure(figsize=(8,4))
    gs = gridspec.GridSpec(1, 2, figure=fig)
    for i, region in enumerate(regions):
        ax = fig.add_subplot(gs[0,i])
        for ff, file in enumerate(emission_files):
            if file == 'GCHP SO4':
                ax2 = ax.twinx()
                ax2.plot(list(range(1, 13)), awm_data[region][file], markers[ff], label=f'{file}')
                ax2.set_ylim(bottom = 0)
                ax2.set_ylabel(Ylabel2)
            else:
                ax.plot(list(range(1, 13)), awm_data[region][file], markers[ff], label=f'{file}')
                ax.set_title(region)
                ax.set_xticks(range(1, 13, 2))
                ax.set_xticklabels(['Jan', 'Mar','May', 'Jul', 'Sep', 'Nov'])
                ax.set_ylim(bottom = 0)
                ax.grid(True)
                ax.set_ylabel(Ylabel1)
        if i == 0:
            # Combining legends
            handles, labels = ax.get_legend_handles_labels()
            handles2, labels2 = ax2.get_legend_handles_labels()
            handles.extend(handles2)
            labels.extend(labels2)
            ax.legend(handles, labels, loc='lower center')

    plt.tight_layout()
    plt.show()
    return fig

def make_plot_3(awm_data,Ylabel1,Ylabel2,emission_files,regions):
    markers = ['b-','b-.','r-.'];
    fig = plt.figure(figsize=(12,4))
    gs = gridspec.GridSpec(1, 3, figure=fig)
    for i, region in enumerate(regions):
        ax = fig.add_subplot(gs[0,i])
        for ff, file in enumerate(emission_files):
            if file == 'GCHP SO4':
                ax2 = ax.twinx()
                ax2.plot(list(range(1, 13)), awm_data[region][file], markers[ff], label=f'{file}')
                ax2.set_ylim(bottom = 0, top = 12)
                ax2.set_ylabel(Ylabel2)
            else:
                ax.plot(list(range(1, 13)), awm_data[region][file], markers[ff], label=f'{file}')
                ax.set_title(region)
                ax.set_xticks(range(1, 13, 2))
                ax.set_xticklabels(['Jan', 'Mar','May', 'Jul', 'Sep', 'Nov'])
                ax.set_ylim(bottom = 0, top = 0.8)
                ax.grid(True)
                ax.set_ylabel(Ylabel1)
        if i == 0:
            # Combining legends
            handles, labels = ax.get_legend_handles_labels()
            handles2, labels2 = ax2.get_legend_handles_labels()
            handles.extend(handles2)
            labels.extend(labels2)
            ax.legend(handles, labels, loc='lower center')

    plt.tight_layout()
    plt.show()
    return fig

def plot_setup_2():
    fig = plt.figure(figsize=(12,5))
    gs = gridspec.GridSpec(1, 3, wspace=0, hspace=0,width_ratios =[1,1,1])  # Changed number of columns to 5 if you have 5 regions
    return fig, gs

def plot_setup():
    fig = plt.figure(figsize=(10,5))
    gs = gridspec.GridSpec(1, 2, wspace=0, hspace=0,width_ratios =[1,1])  # Changed number of columns to 5 if you have 5 regions
    return fig, gs

def plot_map(so2_emi, lat, lon,extents, region, ax, vmin, vmax, add_cbar=0):
    lon_min, lon_max, lat_min, lat_max = extents[region]

    ax.coastlines()
    ax.set_global() 
    ax.set_extent([lon_min, lon_max,lat_min,lat_max]) # avoid white space at edges

    # Flatten the latitude and longitude arrays to ensure they are 1D
    lon = lon.flatten()
    lat = lat.flatten()

    lon_range = (lon >= lon_min) & (lon <= lon_max)  
    lat_range = (lat >= lat_min) & (lat <= lat_max) 

    regional_data = so2_emi[np.ix_(lat_range, lon_range)]
    regional_total = np.nansum(regional_data)  # Ensure to sum only non-NaN values
    regional_mean = np.nanmean(regional_data)  # Ensure to sum only non-NaN values

    lon_grid, lat_grid = np.meshgrid(lon[lon_range], lat[lat_range])
    mesh = ax.pcolormesh(lon_grid, lat_grid, regional_data, transform=ccrs.PlateCarree(), cmap='YlOrBr', vmin=vmin, vmax=vmax)
    
    if add_cbar==1:
        cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
        cbar = plt.colorbar(mesh, cax=cax, orientation='vertical')
        if vmax ==0.5:
            cbar.set_label(f'$SO_2$ VCD (DU)',fontsize = 12)
        if vmax ==2:
            cbar.set_label(f'$SO_4$ ($\mu$g/$m^3$)',fontsize = 12)
    
    x, y = 0.5, 0.05  # Normalized coordinates (0, 0) is lower left, (1, 1) is upper right
    text = f"Sum: {regional_total:.2e}\nMean: {regional_mean:.2e}"
    plt.text(x, y, text, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
             fontsize =14, fontweight = 'bold',backgroundcolor = 'w')
    return regional_mean


def plot_map_2(so2_emi, lat, lon,extents, region, ax, vmin, vmax, add_cbar=0):
    lon_min, lon_max, lat_min, lat_max = extents[region]

    ax.coastlines()
    ax.set_global() 
    ax.set_extent([lon_min, lon_max,lat_min,lat_max]) # avoid white space at edges

    # Flatten the latitude and longitude arrays to ensure they are 1D
    lon = lon.flatten()
    lat = lat.flatten()

    lon_range = (lon >= lon_min-5) & (lon <= lon_max+5)  # make it larger to ensure correct result
    lat_range = (lat >= lat_min-5) & (lat <= lat_max+5) 

    regional_data = so2_emi[np.ix_(lat_range, lon_range)]
    regional_total = np.nansum(regional_data)  # Ensure to sum only non-NaN values
    regional_mean = np.nanmean(regional_data)  # Ensure to sum only non-NaN values

    lon_grid, lat_grid = np.meshgrid(lon[lon_range], lat[lat_range])
    mesh = ax.pcolormesh(lon_grid, lat_grid, regional_data, transform=ccrs.PlateCarree(), cmap='YlOrBr', vmin=vmin, vmax=vmax)
    
    if add_cbar==1:
        cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
        cbar = plt.colorbar(mesh, cax=cax, orientation='vertical')
        if vmax ==0.5:
            cbar.set_label(f'$SO_2$ VCD (DU)',fontsize = 12)
        if vmax ==2:
            cbar.set_label(f'$SO_4$ ($\mu$g/$m^3$)',fontsize = 12)
    
    x, y = 0.5, 0.05  # Normalized coordinates (0, 0) is lower left, (1, 1) is upper right
    text = f"Sum: {regional_total:.2e}\nMean: {regional_mean:.2e}"
    plt.text(x, y, text, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
             fontsize =14, fontweight = 'bold',backgroundcolor = 'w')
    return regional_mean






save_path = '04.Seasonality_Figures/'
regions = ['Middle East','Persian Gulf']
extents = {
    'Middle East': [20, 70, 5, 45],
    'Persian Gulf': [45, 60, 20.5, 32.5]
}

# At SPARTAN sites
target_sites = ['Abu Dhabi', 'Haifa','Rehovot']
extents_sites = {
    'Abu Dhabi': [53.5, 55.5, 23.5, 25.5], 
    'Haifa': [34.5, 35.5, 32.5, 33.5],
    'Rehovot': [ 34.5, 35.5, 31.5, 32.5]
}


data_population = loadmat('Population-0.1.mat')
npop = data_population['npop'].T
latp = data_population['tLAT'].squeeze()
lonp = data_population['tLON'].squeeze()

data_population = loadmat('Population-1x1.mat')
npop1x1 = data_population['npop'].T
latp1x1 = data_population['tLAT'].squeeze()
lonp1x1 = data_population['tLON'].squeeze()

###### SO2 VCD and SO4 conc ####
# Indir_tropomi      = './03.TROPOMI_SO2_Result/02.NASA_L2_OFFLINE_VCD/'
Indir_tropomi      = './03.TROPOMI_SO2_Result/05.NASA_RPRO_Tessellation/'
Indir_gchp         = './03.TROPOMI_SO2_Result/01.GCHP_13.4.0_SO2/'
Indir_gchp_sulfate = './GCHP_13.4.0_Aerosol/'
in_files = ['TROPOMI SO2 tess', 'GCHP SO2', 'GCHP SO4']
year = 2019

# Collect PWM data for each region
pwm_data = {region: {file: [] for file in in_files} for region in regions}
awm_data = {region: {file: [] for file in in_files} for region in regions}
pwm_data_site = {site: {file: [] for file in in_files} for site in target_sites}

# Load and process data
for file in in_files:
    for month_index in range(1,13):
        global_min = 0
        global_max = 0.5  
        if file == 'TROPOMI SO2 tess':
            rfname = f"{Indir_tropomi}Tropomi_Regrid_{year}{month_index:02d}.mat"
            data = loadmat(rfname)
            mapdata, lat, lon = data['so2'], data['tlat'], data['tlon']
            mapdata = mapdata/2.69e16
            # mapdata[mapdata<0] = np.nan
        if file == 'GCHP SO2':
            data = loadmat(f"{Indir_gchp}SO2_GCHP_CEDS_{year}{month_index:02d}.mat");
            mapdata, lat, lon = data['so2'], data['tLAT'], data['tLON']
            mapdata = mapdata/2.69e16
        if file == 'GCHP SO4':
            data = loadmat(f"{Indir_gchp_sulfate}PM25_AOD_{year}{month_index:02d}_ParaReff_C90.mat");
            mapdata, lat, lon = data['tTRAC'], data['tLAT'], data['tLON']
            mapdata = mapdata[:,:,0,0] 
            global_max = 10

        # making maps 
        fig, gs = plot_setup() # Plot set up
        # regional map
        for idx, region in enumerate(regions): 
            ax_region = fig.add_subplot(gs[idx], projection=ccrs.PlateCarree())
            regional_mean = plot_map(mapdata, lat, lon, extents,region, ax_region, global_min, global_max, add_cbar=idx)
            awm_data[region][file].append(regional_mean)
            
        # save figure
        fname = f"{save_path}{file}_{year}{month_index:02d}_ME.png"
        savefig(fname, fig)

        # calculate regional PWM 
        for region in regions:
            lon_min, lon_max, lat_min, lat_max = extents[region]
            lon = lon.flatten()
            lat = lat.flatten()
            lat_range = (lat >= lat_min) & (lat <= lat_max)
            lon_range = (lon >= lon_min) & (lon <= lon_max)
            # region_mask = np.outer(lat_mask, lon_mask)
            if file == 'TROPOMI SO2 tess':
                masked_emissions = mapdata[np.ix_(lat_range, lon_range)]
                masked_population = npop[np.ix_(lat_range, lon_range)]
            else:
                masked_emissions = mapdata[np.ix_(lat_range, lon_range)]
                masked_population = npop1x1[np.ix_(lat_range, lon_range)]

            # calc PWM
            total_population = np.sum(masked_population)
            if total_population > 0:
                pwm = np.nansum(masked_emissions * masked_population) / total_population
            else:
                pwm = np.nan
            pwm_data[region][file].append(pwm)
            

        # making maps at SPARTAN sites
        fig, gs = plot_setup_2() # Plot set up
        # regional map
        for idx, site in enumerate(target_sites): 
            ax_region = fig.add_subplot(gs[idx], projection=ccrs.PlateCarree())
            regional_mean = plot_map_2(mapdata, lat, lon, extents_sites,site, ax_region, global_min, global_max, add_cbar=idx)
            # awm_data[site][file].append(regional_mean)
            
        # save figure
        fname = f"{save_path}{file}_{year}{month_index:02d}_SPARTAN_Sites.png"
        savefig(fname, fig)

        # reading data at SPARTAN sites
        for site in target_sites:
            lon_min, lon_max, lat_min, lat_max = extents_sites[site]
            lon = lon.flatten()
            lat = lat.flatten()
            lat_range = (lat >= lat_min) & (lat <= lat_max)
            lon_range = (lon >= lon_min) & (lon <= lon_max)
            # region_mask = np.outer(lat_mask, lon_mask)
            if file == 'TROPOMI SO2 tess':
                masked_emissions = mapdata[np.ix_(lat_range, lon_range)]
                masked_population = npop[np.ix_(lat_range, lon_range)]
            else:
                masked_emissions = mapdata[np.ix_(lat_range, lon_range)]
                masked_population = npop1x1[np.ix_(lat_range, lon_range)]

            # calc PWM
            total_population = np.sum(masked_population)
            if total_population > 0:
                pwm = np.nansum(masked_emissions * masked_population) / total_population
            else:
                pwm = np.nan
            pwm_data_site[site][file].append(pwm)


# Plotting PWM data
Ylabel1 = r'$SO_2$ VCD (DU)'
Ylabel2 = r'$SO_4$ ($\mu$g/$m^3$)'
fig = make_plot_2(pwm_data,Ylabel1,Ylabel2,in_files,regions)
fname = f"{save_path}so2_vcd_so4_seasonality_pwm_ME.png"
savefig(fname, fig)

# Plotting AWM data
fig = make_plot_2(awm_data,Ylabel1,Ylabel2,in_files,regions)
fname = f"{save_path}so2_vcd_so4_seasonality_ME.png"
savefig(fname, fig)

# Plotting PWM data
fig = make_plot_3(pwm_data_site,Ylabel1,Ylabel2,in_files,target_sites)
fname = f"{save_path}so2_vcd_so4_seasonality_SPARTAN_sites.png"
savefig(fname, fig)
