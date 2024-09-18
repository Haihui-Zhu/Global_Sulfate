import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.io import loadmat
import matplotlib.gridspec as gridspec

def read_emission_data(file_path):
    data = loadmat(file_path)
    return data['emi_so2'], data['lat'], data['lon']

def plot_map(so2_emi, lat, lon,extents, region, ax, vmin, vmax, add_cbar=False):
    lon_min, lon_max, lat_min, lat_max = extents[region]

    ax.coastlines()
    ax.set_global() 
    ax.set_extent([lon_min+5, lon_max-5,lat_min+5,lat_max-5]) # avoid white space at edges

    # calc annual mean emission
    data = np.nanmean(so2_emi, axis=2)
    # Flatten the latitude and longitude arrays to ensure they are 1D
    lon = lon.flatten()
    lat = lat.flatten()

    lon_range = (lon >= lon_min) & (lon <= lon_max)  
    lat_range = (lat >= lat_min) & (lat <= lat_max)  

    regional_data = data[np.ix_(lat_range, lon_range)]
    regional_total = np.nansum(regional_data)  # Ensure to sum only non-NaN values
    regional_mean = np.nanmean(regional_data)  # Ensure to sum only non-NaN values

    lon_grid, lat_grid = np.meshgrid(lon[lon_range], lat[lat_range])
    mesh = ax.pcolormesh(lon_grid, lat_grid, regional_data, transform=ccrs.PlateCarree(), cmap='YlOrBr', vmin=vmin, vmax=vmax)
    
    if add_cbar:
        cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
        cbar = plt.colorbar(mesh, cax=cax, orientation='vertical')
        cbar.set_label('SO2 emission (kg m-2 s-1)',fontsize = 12)
    
    x, y = 0.5, 0.05  # Normalized coordinates (0, 0) is lower left, (1, 1) is upper right
    text = f"Sum: {regional_total:.2e}\nMean: {regional_mean:.2e}"
    plt.text(x, y, text, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
             fontsize =14, fontweight = 'bold',backgroundcolor = 'w')


def plot_map_diff(so2_emi, lat, lon,extents, region, ax, vmin, vmax, add_cbar=False):
    ax.coastlines()
    ax.set_global() 

    lon_min, lon_max, lat_min, lat_max = extents[region]

    ax.set_extent([lon_min+5, lon_max-5,lat_min+5,lat_max-5]) # avoid white space at edges

    # Flatten the latitude and longitude arrays to ensure they are 1D
    lon = lon.flatten()
    lat = lat.flatten()

    lon_range = (lon >= lon_min) & (lon <= lon_max)  
    lat_range = (lat >= lat_min) & (lat <= lat_max)  

    regional_data = so2_emi[np.ix_(lat_range, lon_range)]
    regional_mean_emission = np.nanmean(regional_data)  # Ensure to sum only non-NaN values


    lon_grid, lat_grid = np.meshgrid(lon[lon_range], lat[lat_range])
    mesh = ax.pcolormesh(lon_grid, lat_grid, regional_data, transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=vmin, vmax=vmax)
    
    if add_cbar:
        cax = ax.inset_axes([0.05, 0.1, 0.02, 0.5])  # Adjust these parameters as needed for placement
        cbar = plt.colorbar(mesh, cax=cax, orientation='vertical')
        cbar.set_label('SO2 emission (kg m-2 s-1)',fontsize = 12)
    
    x, y = 0.5, 0.05  # Normalized coordinates (0, 0) is lower left, (1, 1) is upper right
    text = f"Mean: {regional_mean_emission:.2e}"
    plt.text(x, y, text, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
             fontsize =14, fontweight = 'bold',backgroundcolor = 'w')

    # ax.set_title(f"{region} - Tot: {np.nansum(so2_emi):.2e} molec cm-2")

def main():
    emission_path = '01.Emissions/'
    save_path = '01.Emissions/'
    emission_files = ['CEDS-2021', 'EDGARv6.1','HTAPv3']  # Add more files as needed
    regions = ['Global', 'North America', 'Europe', 'East Asia', 'South Asia', 'Middle East']
    

    extents = {
        'Global':[-180, 180, -50, 60],
        'North America': [-140, -50, -15, 75],
        'Europe': [-10, 40, 20, 70],
        'East Asia': [95, 145, 0, 50],
        'South Asia': [60, 100, -5, 35],
        'Middle East': [25, 65, 0, 40]
    }

    # Initialize a list
    ratios = []

    # Loop through each region in the regions list
    for region in regions:
        # Extract longitude and latitude bounds from the dictionary
        long_min, long_max, lat_min, lat_max = extents[region]
        
        # Calculate the ranges for longitude and latitude
        lat_range = lat_max - lat_min
        long_range = long_max - long_min
        
        # Calculate the ratio of latitude to longitude range
        ratio =  long_range / lat_range
        
        ratios.append(ratio)  # Append each ratio to the list

    # ratio of global map height to regional map height    
    h1 = sum(ratios[1:])/ratios[0]

    overallh = (1+h1)/sum(ratios[1:])

    for file in emission_files:
        so2_emi, lat, lon = read_emission_data(f"{emission_path}{file}_seasonal_allsectors.mat")
        # color scale 
        global_min = 1e-12
        global_max = 4e-11

        fig = plt.figure(figsize=(15, 15*overallh))
        gs = gridspec.GridSpec(2, 5, wspace=0, hspace=0, height_ratios=[h1, 1], width_ratios=ratios[1:])  # Changed number of columns to 5 if you have 5 regions

        # Plotting the global map as the first and largest plot
        ax_global = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
        plot_map(so2_emi, lat, lon,extents, 'Global', ax_global, global_min, global_max, add_cbar=True)

        for idx, region in enumerate(regions[1:]):
            ax_region = fig.add_subplot(gs[1, idx], projection=ccrs.PlateCarree())
            plot_map(so2_emi, lat, lon, extents,region, ax_region, global_min, global_max)

        plt.suptitle(f"Emission Maps for {file}")
        # save figure
        fname = f"{save_path}so2_{file}.png"
        fig.savefig(fname) 
        plt.close(fig)  # Close the plot to free memory
        print('Figure saved: '+ fname)


    # comparing inventories
    for file in emission_files[1:]:
        so2_emi_1, lat, lon = read_emission_data(f"{emission_path}{emission_files[0]}_seasonal_allsectors.mat")
        so2_emi_2, lat, lon = read_emission_data(f"{emission_path}{file}_seasonal_allsectors.mat")
        so2_emi_1_ann = np.nanmean(so2_emi_1, axis=2)
        so2_emi_2_ann = np.nanmean(so2_emi_2, axis=2)
        so2_emi = (so2_emi_2_ann-so2_emi_1_ann)
        global_min = -1e-11
        global_max = 1e-11

        fig = plt.figure(figsize=(15, 15*overallh))
        gs = gridspec.GridSpec(2, 5, wspace=0, hspace=0, height_ratios=[h1, 1], width_ratios=ratios[1:])  # Changed number of columns to 5 if you have 5 regions

        # Plotting the global map as the first and largest plot
        ax_global = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
        plot_map_diff(so2_emi, lat, lon,extents, 'Global', ax_global, global_min, global_max, add_cbar=True)

        for idx, region in enumerate(regions[1:]):
            ax_region = fig.add_subplot(gs[1, idx], projection=ccrs.PlateCarree())
            plot_map_diff(so2_emi, lat, lon, extents,region, ax_region, global_min, global_max)

        plt.suptitle(f"{file} - CEDS")
        # save figure
        fname = f"{save_path}so2_{file}-{emission_files[0]}.png"
        fig.savefig(fname) 
        plt.close(fig)  # Close the plot to free memory
        print('Figure saved: '+ fname)



    # seasonal variation
    for file in emission_files:
        so2_emi, lat, lon = read_emission_data(f"{emission_path}{file}_seasonal_allsectors.mat")
        so2_emi_diff = so2_emi[:,:,7]-so2_emi[:,:,1]
        # color scale 
        global_min = -1e-11
        global_max = 1e-11

        fig = plt.figure(figsize=(15, 15*overallh))
        gs = gridspec.GridSpec(2, 5, wspace=0, hspace=0, height_ratios=[h1, 1], width_ratios=ratios[1:])  # Changed number of columns to 5 if you have 5 regions

        # Plotting the global map as the first and largest plot
        ax_global = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
        plot_map_diff(so2_emi_diff, lat, lon,extents, 'Global', ax_global, global_min, global_max, add_cbar=True)

        for idx, region in enumerate(regions[1:]):
            ax_region = fig.add_subplot(gs[1, idx], projection=ccrs.PlateCarree())
            plot_map_diff(so2_emi_diff, lat, lon, extents,region, ax_region, global_min, global_max)

        plt.suptitle(f"{file}: July-Jan")
        # save figure
        fname = f"{save_path}so2_{file}_season.png"
        fig.savefig(fname) 
        plt.close(fig)  # Close the plot to free memory
        print('Figure saved: '+ fname)


if __name__ == "__main__":
    main()
