
import numpy as np
from pyproj import Geod
import xarray as xr

# Create arrays for latitudes and longitudes
res = 0.1
lats = np.arange(-90+res/2, 90+res/2, res)  # Latitude array from -90 to 89 degrees
lons = np.arange(-180+res/2, 180+res/2, res)  # Longitude array from -180 to 179 degrees

# Initialize the Geod object for WGS84 (World Geodetic System 1984 (WGS84) is a global coordinate system 
# and geodetic datum that defines the Earth's shape, size, gravity, and geomagnetic fields)
geod = Geod(ellps='WGS84')
# geod = Geod('+a=6378137 +f=0.0033528106647475126')
# Initialize an empty 2D array to store the areas of each grid cell
area_grid = np.zeros((len(lats), len(lons)))

# Loop through latitudes and longitudes to calculate area of each grid cell
for i, lat in enumerate(lats):
    for j, lon in enumerate(lons):
        # Define the corners of the grid cell (1-degree x 1-degree box)
        lon_corners = [lon-res/2, lon-res/2, lon + res/2, lon+res/2]
        lat_corners = [lat-res/2, lat+res/2, lat + res/2, lat-res/2]
        
        # Calculate the area of the grid cell (in square meters)
        area, _ = geod.polygon_area_perimeter(lon_corners, lat_corners)
        
        # Store the absolute value of the area (since area can be negative)
        area_grid[i, j] = abs(area)

# Convert area from square meters to square kilometers
area_grid_km2 = area_grid / 1e6



# output emission data
ds_out = xr.Dataset({'lat': lats, 'lon': lons})

data_array = xr.DataArray(area_grid_km2, coords=[lats, lons], dims=['lat', 'lon'], \
    attrs=dict(long_name=f'Pixel Area for Lat Lon grid with Res={res}', \
        units='km2',missing_value = 1.e20))
ds_out['area'] = data_array
ds_out.lon.attrs = dict(units='degrees_east', axis='X', standard_name='Longitude', long_name='Longitude')
ds_out.lat.attrs = dict(units='degrees_north', axis='Y', standard_name='Latitude', long_name='Latitude')

ds_out.to_netcdf(f"utils/pixel_area_{res}.nc", engine='netcdf4', format='NETCDF4', unlimited_dims='time')
print('file saved')
