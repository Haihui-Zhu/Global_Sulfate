import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.io import loadmat
import pandas as pd
import numpy as np
import xarray as xr
import os
from tqdm import tqdm  # For progress bars


# Configuration
indir = '/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/'
outdir = '/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/all_data_by_site/'
os.makedirs(outdir, exist_ok=True)

# Satellite data version
qc_str = 'CF03-SZA60-QA75'  # Update this with your QC string
year = 2021

# Load SPARTAN data
data_version = 'pm25_so4_202410'
fname = f'{indir}/SPT_{data_version}.mat'
datain = loadmat(fname)
data = datain['TOT']

# Process metadata
D3_SiteCode = [item[0][0] for item in datain['D3_SiteCode']]
lats = datain['latitudes'].flatten()
lons = datain['longitudes'].flatten()
D1_Dates = datain['D1_Dates']
n_sites = len(D3_SiteCode)

# Process dates
date_str = [f'{row[0]}-{row[1]:02d}-{row[2]:02d}' for row in D1_Dates]
dates = pd.to_datetime(date_str)
n_dates = len(dates)

# Initialize DataFrames
pm25_df = pd.DataFrame(index=dates, columns=D3_SiteCode)
so4_df = pd.DataFrame(index=dates, columns=D3_SiteCode)

cobra_so2 = pd.DataFrame(index=dates, columns=D3_SiteCode)
gchp_so2 = pd.DataFrame(index=dates, columns=D3_SiteCode)
gchp_so4 = pd.DataFrame(index=dates, columns=D3_SiteCode)
gchp_pm25 = pd.DataFrame(index=dates, columns=D3_SiteCode)

# Load species data
for spid, spec in enumerate(datain['Species'][0]):
    spec_name = spec[0]
    if spec_name == 'PM2.5':
        pm25_df.values[:] = data[:, spid, :]
    elif spec_name == 'SO4':
        so4_df.values[:] = data[:, spid, :]

# functions
def find_nearest_index(grid_points, target):
    """Find index of nearest value in 1D array"""
    return np.abs(grid_points - target).argmin()

for date in tqdm(dates, desc='Processing dates'):
    # Process COBRA data
    cobra_path = f"/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/COBRA/\
COBRA_Tesellation_{qc_str}/COBRA_RG_{date.strftime('%Y%m%d')}_{qc_str}.nc"
    if os.path.exists(cobra_path):
        with xr.open_dataset(cobra_path) as ds:
            try:
                if 'idxclat' not in locals():
                    # Get grid and data
                    lat_grid = ds.lat.values
                    lon_grid = ds.lon.values
                    idxclat = np.zeros(n_sites, dtype=int)
                    idxclon = np.zeros(n_sites, dtype=int)
                    for sid, site in enumerate(D3_SiteCode):
                        idxclat[sid] = find_nearest_index(lat_grid, lats[sid])
                        idxclon[sid] = find_nearest_index(lon_grid, lons[sid])
                   
                so2_data = ds.so2.values.T
                cobra_so2.loc[date] = so2_data[idxclat, idxclon]
            except Exception as e:
                print(f"COBRA interpolation error {date}: {str(e)}")
                cobra_so2.loc[date] = np.nan
                exit()
    else:
        cobra_so2.loc[date] = np.nan
        
    # Process GCHP-CEDS data
    gchp_path = f"/storage1/fs1/rvmartin2/Active/haihuizhu/5.GEOS-Chem/\
5.nasa_run_full_in_archive/4.ceds_2021/processed/GCHP_SO2_SO4_PM25_ceds_{date.year}_{date.month:02d}{date.day:02d}.nc"
    if os.path.exists(gchp_path):
        with xr.open_dataset(gchp_path) as ds:
            try:
                if 'idxglat' not in locals():
                    # Get grid and data
                    lat_grid = ds.latitude.values
                    lon_grid = ds.longitude.values
                    idxglat = np.zeros(n_sites, dtype=int)
                    idxglon = np.zeros(n_sites, dtype=int)
                    for sid, site in enumerate(D3_SiteCode):
                        idxglat[sid] = find_nearest_index(lat_grid, lats[sid])
                        idxglon[sid] = find_nearest_index(lon_grid, lons[sid])
                
                so2_data = ds.so2.values.T
                so4_data = ds.so4.values.T
                pm25_data = ds.pm25.values.T
                gchp_so2.loc[date] = so2_data[idxglat, idxglon]
                gchp_so4.loc[date] = so4_data[idxglat, idxglon]
                gchp_pm25.loc[date] = pm25_data[idxglat, idxglon]
                
            except Exception as e:
                print(f"GCHP interpolation error {date}: {str(e)}")
                gchp_so2.loc[date] = np.nan
                gchp_so4.loc[date] = np.nan
                gchp_pm25.loc[date] = np.nan
                exit()
    else:
        gchp_so2.loc[date] = np.nan
        gchp_so4.loc[date] = np.nan
        gchp_pm25.loc[date] = np.nan
    
# write to csv by site
for site in tqdm(D3_SiteCode, desc='Processing sites'):
    site_data = {
        'date': date_str,
        'pm25': pm25_df[site],
        'so4': so4_df[site],
        'cobra_so2': cobra_so2[site],
        'gchp-ceds_so2': gchp_so2[site],
        'gchp-ceds_so4': gchp_so4[site],
        'gchp-ceds_pm25':gchp_pm25[site],
    }
    
    # Create and save DataFrame
    df = pd.DataFrame(site_data)
    df.to_csv(os.path.join(outdir, f"{site}_pm25_so2_so4.csv"), index=False)

    
print("Processing complete! All files saved to:", outdir)