# this script create a mask file for pacific 
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def savefig(fname, fig):
    fig.savefig(fname) 
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: '+ fname)
    return


# Define latitude and longitude ranges (global, from -90 to 90 for latitude, -180 to 180 for longitude)
lat_min, lat_max = -90, 90
lon_min, lon_max = -180, 180
resolution = 0.05

# Create arrays for latitude and longitude values
latitudes = np.linspace(lat_min+0.5*resolution, lat_max+0.5*resolution, (lat_max-lat_min)/resolution) # max + 0.5xres -> end needs to be outside of the last target point
longitudes = np.linspace(lon_min+0.5*resolution, lon_max+0.5*resolution, (lon_max-lon_min)/resolution)
lon_grid, lat_grid = np.meshgrid(longitudes, latitudes)  # shape (n, m)

# Initialize a 2D array of NaN values with shape (latitudes, longitudes)
array = np.full((len(latitudes), len(longitudes)), np.nan)


print(f"Latitude grid shape: {lat_grid.shape}")
print(f"Longitude grid shape: {lon_grid.shape}")

### SQUARE 1 ###
pacific_lat_min, pacific_lat_max = 4, 45  # Rough estimate for pacific lat range
pacific_lon_min, pacific_lon_max = 150, 180  # Rough estimate for pacific lon range

array = np.where( ( (lat_grid >= pacific_lat_min) & (lat_grid <= pacific_lat_max) & \
                    (lon_grid >= pacific_lon_min) & (lon_grid <= pacific_lon_max) ), \
                  1, array)

### SQUARE 2 ###
pacific_lat_min, pacific_lat_max = 0, 60  # Rough estimate for pacific lat range
pacific_lon_min, pacific_lon_max = -180, -140  # Rough estimate for pacific lon range

array = np.where( ( (lat_grid >= pacific_lat_min) & (lat_grid <= pacific_lat_max) & \
                    (lon_grid >= pacific_lon_min) & (lon_grid <= pacific_lon_max) ), \
                  1, array)

### SQUARE 3 ###
pacific_lat_min, pacific_lat_max = 0, 30  # Rough estimate for pacific lat range
pacific_lon_min, pacific_lon_max = -140, -125  # Rough estimate for pacific lon range

array = np.where( ( (lat_grid >= pacific_lat_min) & (lat_grid <= pacific_lat_max) & \
                    (lon_grid >= pacific_lon_min) & (lon_grid <= pacific_lon_max) ), \
                  1, array)

### SQUARE 4 ###
pacific_lat_min, pacific_lat_max = -60, 0  # Rough estimate for pacific lat range
pacific_lon_min, pacific_lon_max = -170, -90  # Rough estimate for pacific lon range

array = np.where( ( (lat_grid >= pacific_lat_min) & (lat_grid <= pacific_lat_max) & \
                    (lon_grid >= pacific_lon_min) & (lon_grid <= pacific_lon_max) ), \
                  1, array)

### SQUARE 5 ### Exclude Hawaii
pacific_lat_min, pacific_lat_max = 10, 26  # Rough estimate for pacific lat range
pacific_lon_min, pacific_lon_max = -180, -150  # Rough estimate for pacific lon range

array = np.where( ( (lat_grid >= pacific_lat_min) & (lat_grid <= pacific_lat_max) & \
                    (lon_grid >= pacific_lon_min) & (lon_grid <= pacific_lon_max) ), \
                  np.nan, array)

### SQUARE 6 ### Exclude Island Isabela
pacific_lat_min, pacific_lat_max = -10, 0  # Rough estimate for pacific lat range
pacific_lon_min, pacific_lon_max = -120, -90  # Rough estimate for pacific lon range

array = np.where( ( (lat_grid >= pacific_lat_min) & (lat_grid <= pacific_lat_max) & \
                    (lon_grid >= pacific_lon_min) & (lon_grid <= pacific_lon_max) ), \
                  np.nan, array)

## Save mask in a nc file
# Step 1: Create an xarray Dataset
ds = xr.Dataset(
    {
        "mask": (["latitude", "longitude"], array),  # Add the mask data
        "lat_mesh": (["latitude", "longitude"], lat_grid),  # Add the mask data
        "lon_mesh": (["latitude", "longitude"], lon_grid)  # Add the mask data
    },
    
    coords={
        "latitude": latitudes,  # Add the latitude values
        "longitude": longitudes  # Add the longitude values
    }
)
# Step 2: Save the dataset to a NetCDF file
output_nc_file = "pacific_mask.nc"
ds.to_netcdf(output_nc_file)


##  Make a map for validation
# load gchp so2 vcd
simname = 'ceds_2018'
year = 2018
qcstr =  'CF03-SZA60-QA90'
source2 = f'./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_{simname}_annual.nc'
ds = xr.open_dataset(source2) 
so2 = ds['so2'].values.T 
so2 = so2/2.69e16 # convert unit
so2_lat = ds['lat'].values
so2_lon = ds['lon'].values

print(f"so2 grid shape: {so2.shape}")

# Create a figure and set the projection to a global map with the Pacific centered (180 degrees)
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))

# Plot the concentration data using the meshgrid
concentration_plot = ax.pcolormesh(so2_lon, so2_lat, so2, cmap='RdBu_r', transform=ccrs.PlateCarree(), shading='auto', vmin=-0.5, vmax=0.5)

# Plot the mask data, with a color map and handling NaN values
masked_array = np.ma.masked_where(np.isnan(array), array)

# Plot the data
# img = ax.pcolormesh(longitudes, latitudes, masked_array, cmap='Blues', transform=ccrs.PlateCarree(), shading='auto')
mask_plot = ax.pcolormesh(lon_grid, lat_grid, masked_array, cmap='viridis', alpha=0.2, transform=ccrs.PlateCarree(), shading='auto')

# Add coastlines and features
ax.coastlines()
ax.add_feature(cfeature.BORDERS, linestyle=':')

# Add colorbar
cax = ax.inset_axes([0.16, 0.15, 0.02, 0.3])  # Adjust these parameters as needed for placement
plt.colorbar(concentration_plot, cax=cax, orientation='vertical',label='SO2 VCD (DU)')
    
# Set gridlines
ax.gridlines(draw_labels=True)

# Display the plot
plt.title('Global Map of the Mask (Pacific-Centered)')
plt.show()

# save figure
fname = f"./temp/map_pacific_mask.png"
savefig(fname, fig)