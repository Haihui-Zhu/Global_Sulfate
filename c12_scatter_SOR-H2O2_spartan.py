import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


P = 101325
R = 8.314
T = 298.15 # will need to use ambient T and P 

# 1. Read site details
sites_df = pd.read_csv('spartan/Site_details.csv', usecols=['Site_Code', 'Latitude', 'Longitude'])

# 2. Read NetCDF file
ds = xr.open_dataset('temp/GCHP_monthly_htap_202101.nc')
# Extract variables and squeeze singleton dimensions
so2 = ds['so2_s'].squeeze()
h2o2 = ds['h2o2_s'].squeeze()
so4 = ds['so4_s'].squeeze()
unitconv = ds['unit_converter_s'].squeeze().values # molec/cm2

# Calculate latitude and longitude from dataset
lats = ds['lat'].values
lons = ds['lon'].values

# 3. Unit conversion for SO2 (convert to µg/m³)
with np.errstate(divide='ignore', invalid='ignore'):
    so2_ugm3 = (so2 / unitconv) * P / (R * T) * 64 * 1e6
    h2o2_ugm3 = (h2o2 / unitconv) * P / (R * T) * 34 * 1e6


# 4. Create interpolation function using xarray's nearest neighbor
def get_nearest_values(lat, lon):
    try:
        return {
            'so2_ugm3': so2_ugm3.sel(lat=lat, lon=lon, method='nearest').item(),
            'h2o2_ugm3': h2o2_ugm3.sel(lat=lat, lon=lon, method='nearest').item(),
            'so4': so4.sel(lat=lat, lon=lon, method='nearest').item()
        }
    except KeyError:
        return {'so2_ugm3': np.nan, 'h2o2_ugm3': np.nan, 'so4': np.nan}

# Apply interpolation and filter invalid points
sites_df[['so2_ugm3', 'h2o2_ugm3', 'so4']] = sites_df.apply(
    lambda row: pd.Series(get_nearest_values(row['Latitude'], row['Longitude'])),
    axis=1
)
sites_df.dropna(subset=['so2_ugm3', 'h2o2_ugm3', 'so4'], inplace=True)

# 5. Create scatter plot
plt.figure(figsize=(10, 8))
sc = plt.scatter(
    sites_df['so2_ugm3'],
    sites_df['so4'],
    c=sites_df['h2o2_ugm3'],
    cmap='viridis',
    s=60,
    edgecolor='w',
    linewidth=0.5
)

# Formatting
plt.colorbar(sc, label='H₂O₂')
plt.xlabel('SO₂', fontsize=12)
plt.ylabel('SO₄', fontsize=12)
plt.title('SPARTAN Sites GCHP-HTAP-Jan', fontsize=14)

# Save figure
plt.savefig('figures/Scatter_spartan_SO2_SO4_H2O2.png', dpi=300, bbox_inches='tight')
plt.close()

# Optional: Save the table to CSV
sites_df.to_csv('spartan/SPARTAN_site_SO4_SO2_H2O2.csv', index=False)



# 5. SOR vs H2O2/SO2 plot
plt.figure(figsize=(10, 8))
sc = plt.scatter(
    sites_df['h2o2_ugm3'] / sites_df['so2_ugm3'],
    sites_df['so4'] / (sites_df['so4'] + sites_df['so2_ugm3']),
    c=sites_df['so4'],
    cmap='viridis',
    s=60,
    edgecolor='w',
    linewidth=0.5
)

# Formatting
plt.colorbar(sc, label='SO₄')
plt.xlabel('H₂O₂/SO₂', fontsize=12)
plt.ylabel('SOR', fontsize=12)
plt.title('SPARTAN Sites GCHP-HTAP-Jan', fontsize=14)

# Save figure
plt.savefig('figures/Scatter_spartan_SOR_H2O2.png', dpi=300, bbox_inches='tight')
plt.close()
