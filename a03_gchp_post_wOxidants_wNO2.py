import os
import xesmf as xe
import numpy as np
import xarray as xr
from scipy.interpolate import griddata
from datetime import datetime

# Constants
R = 8.314  # m3⋅Pa⋅K−1⋅mol−1
Na = 6.02e23
SO4MW = 96
BCMW = 12
RH_growth_sia = 1.10

def cube2latlon2d(var1, lat, lon, xLAT, xLON):
    # print(var1.shape, lat.shape, lon.shape, xLAT.shape, xLON.shape)
    points = np.column_stack((lon.ravel(), lat.ravel()))
    values = var1.ravel()
    # print(points.shape, values.shape, xLON.shape, xLAT.shape)
    var_interp = griddata(points, values, (xLON, xLAT), method='linear', fill_value=0)
    return var_interp

def savenc(sfname, so2, h2o2, o3, oh, no2s, pm25, so4, so2s, h2o2s, o3s, ohs, unit_conv, tLAT, tLON):
    """Save data to netCDF file with all variables"""
    ds = xr.Dataset(
        {
            'so2': (('lat', 'lon'), so2),
            'h2o2': (('lat', 'lon'), h2o2),
            'o3': (('lat', 'lon'), o3),
            'oh': (('lat', 'lon'), oh),
            'no2_s': (('lat', 'lon'), no2s),
            'pm25_s': (('lat', 'lon'), pm25),
            'so4_s': (('lat', 'lon'), so4),
            'so2_s': (('lat', 'lon'), so2s),
            'h2o2_s': (('lat', 'lon'), h2o2s),
            'o3_s': (('lat', 'lon'), o3s),
            'oh_s': (('lat', 'lon'), ohs),
            'unit_converter_s': (('lat', 'lon'), unit_conv),
        },
        coords={'lat': tLAT, 'lon': tLON}
    )
    
    # Add variable attributes
    ds.lat.attrs['units'] = 'degrees_north'
    ds.lon.attrs['units'] = 'degrees_east'
    
    ds.so2.attrs = {'units': 'molec/cm2', 
                   'long_name': 'Sulfur Dioxide Vertical Column Density at S5P overpass time'}
    ds.h2o2.attrs = {'units': 'molec/cm2',
                    'long_name': 'Hydrogen Peroxide Vertical Column Density at S5P overpass time'}
    ds.o3.attrs = {'units': 'molec/cm2',
                  'long_name': 'Ozone Vertical Column Density at S5P overpass time'}
    ds.oh.attrs = {'units': 'molec/cm2',
                  'long_name': 'Hydroxyl Radical Vertical Column Density at S5P overpass time'}
    
    ds.no2_s.attrs = {'units': 'mol/mol dry air',
                     'long_name': 'Surface Nitrogen Dioxide at S5P overpass time'}
    ds.unit_converter_s.attrs = {'units': 'molec/cm2',
                                'long_name': 'Surface unit converter from mol/mol dry air to molec/cm2 at satellite overpass time'}
    
    ds.pm25_s.attrs = {'units': 'ug/m3', 
                      'long_name': '24 hr mean surface PM25 concentration at RH 35%'}
    ds.so4_s.attrs = {'units': 'ug/m3',
                     'long_name': '24 hr mean surface sulfate concentration at RH 35%'}
    ds.so2_s.attrs = {'units': 'molec/cm2',
                     'long_name': '24 hr mean surface Sulfur Dioxide concentration'}
    ds.h2o2_s.attrs = {'units': 'molec/cm2',
                      'long_name': '24 hr mean surface Hydrogen Peroxide concentration'}
    ds.o3_s.attrs = {'units': 'molec/cm2',
                    'long_name': '24 hr mean surface Ozone concentration'}
    ds.oh_s.attrs = {'units': 'molec/cm2',
                    'long_name': '24 hr mean surface Hydroxyl Radical concentration'}
    
    ds.to_netcdf(sfname)
    print(f"Saved {sfname}")

# Configuration
YEARS = [2021]
SimName = 'htap_2021'
InDir = '/storage1/fs1/rvmartin2/Active/haihuizhu/5.GEOS-Chem/1.aws/2.htap2021_met_2021_all/OutputDir/'
OutDir = '/storage1/fs1/rvmartin2/Active/haihuizhu/5.GEOS-Chem/1.aws/2.htap2021_met_2021_all/processed/'
Mns = range(1, 13)
daysinmonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# Prepare grid
tLAT = np.linspace(-89.5, 89.5, 180)
tLON = np.linspace(-179.5, 179.5, 360)
xLON, xLAT = np.meshgrid(tLON, tLAT)

# Get original grid from sample file
sample_file = f"{InDir}/01/GEOSChem.ACAG.{YEARS[0]}0101_0030z.nc4"
with xr.open_dataset(sample_file) as ds:
    lat = ds['lats'].values.astype(float)
    lon = ds['lons'].values.astype(float)
    lon[lon > 180] -= 360  # Convert to -180-180

for Yr in YEARS:
    # Annual accumulators
    so2 = np.zeros((len(tLAT), len(tLON)))
    annual_vars = {
        'so2': np.zeros_like(so2),
        'h2o2': np.zeros_like(so2),
        'o3': np.zeros_like(so2),
        'oh': np.zeros_like(so2),
        'no2s': np.zeros_like(so2),
        
        'pm25': np.zeros_like(so2),
        'so4': np.zeros_like(so2),
        'so2s': np.zeros_like(so2),
        'h2o2s': np.zeros_like(so2),
        'o3s': np.zeros_like(so2),
        'ohs': np.zeros_like(so2),
        
        'unit_conv': np.zeros_like(so2),
    }
    
    for Mn in Mns:
        # Monthly accumulators
        monthly_vars = {k: np.zeros_like(v) for k, v in annual_vars.items()}
        
        for Dy in range(1, daysinmonth[Mn-1]+1):
            sfname = f"{OutDir}/GCHP_daily_{SimName}{Mn:02d}{Dy:02d}.nc"
            if os.path.exists(sfname):
                print(f"{sfname} exists, reading")
                ds = xr.open_dataset(sfname)
                # Accumulate annual
                for var in monthly_vars:
                    monthly_vars[var] += ds[var]
                continue
            
            # Daily overpass time [melec/cm2 or mol/mol dry air]
            daily_overpass_vars = {
                'so2': np.zeros_like(so2),
                'h2o2': np.zeros_like(so2),
                'o3': np.zeros_like(so2),
                'oh': np.zeros_like(so2),
                'no2s': np.zeros_like(so2),
                'unit_conv': np.zeros_like(so2),
            }
            # daily 24 hr mean [ug/m3]
            daily_24hr_vars = {
                'pm25': np.zeros_like(so2),
                'so4': np.zeros_like(so2),
                'so2s': np.zeros_like(so2),
                'h2o2s': np.zeros_like(so2),
                'o3s': np.zeros_like(so2),
                'ohs': np.zeros_like(so2),
            }
            
            for Hr in range(24):
                # Read files
                fname = f"{InDir}/{Mn:02d}/GEOSChem.ACAG.{Yr}{Mn:02d}{Dy:02d}_{Hr:02d}30z.nc4"
                
                # TROPOMI sampling mask
                loc = (13 - Hr) * 15
                if loc > 180: loc -= 360
                lower = loc - 7.5
                upper = loc + 7.5
                mask = ((tLON >= lower) & (tLON < upper)) if lower <= upper else ((tLON >= lower) | (tLON < upper))
                LonInd = np.where(mask)[0]
                
                with xr.open_dataset(fname) as ds:
                    
                    # Load 4D variables
                    hBxH = ds['Met_BXHEIGHT'].values.squeeze()
                    tT = ds['Met_T'].values.squeeze()
                    tP = ds['Met_PMID'].values.squeeze() * 100
                    
                    # Calculate unit conversion factors
                    vol = R * tT / tP
                    mol_conv = hBxH / vol * Na * 1e-4
                    
                    # Process column variables
                    for var in ['so2', 'h2o2', 'o3', 'oh']:
                        data = ds[f'SpeciesConcVV_{var.upper()}'].values.squeeze()
                        data = (data * mol_conv)
                        data_sfc = data[0, ...]
                        data_col = data.sum(axis=0)  # Sum vertical layers
                        
                        data_sfc_regrid = cube2latlon2d(data_sfc, lat, lon, xLAT, xLON)
                        data_col_regrid = cube2latlon2d(data_col, lat, lon, xLAT, xLON)
                        daily_24hr_vars[f'{var}s'] += data_sfc_regrid
                        daily_overpass_vars[var][:, LonInd] += data_col_regrid[:, LonInd]
                        
                    
                    # NO2: surface,  not unit convertion
                    no2_sfc = ds['SpeciesConcVV_NO2'].values.squeeze()[0, ...]
                    no2_regrid = cube2latlon2d(no2_sfc, lat, lon, xLAT, xLON)
                    daily_overpass_vars['no2s'][:, LonInd] += no2_regrid[:, LonInd]
                    # Unit converter surface
                    unit_conv_sfc = cube2latlon2d(mol_conv[0, ...], lat, lon, xLAT, xLON)
                    daily_overpass_vars['unit_conv'][:, LonInd] += unit_conv_sfc[:, LonInd]  # Only keep last hour's value
                    
                    # PM25 surface
                    pm25_sfc = ds['PM25'].values.squeeze()[0, ...]
                    daily_24hr_vars['pm25'] += cube2latlon2d(pm25_sfc, lat, lon, xLAT, xLON)
                    
                    # SO4 surface
                    so4 = ds['SpeciesConcVV_SO4'].values.squeeze()[0, ...]
                    unitconv = (1e6 / R) * (tP / tT)[0, ...]
                    so4_sfc = so4 * unitconv * SO4MW * RH_growth_sia
                    daily_24hr_vars['so4'] += cube2latlon2d(so4_sfc, lat, lon, xLAT, xLON)
                
                
            # Create daily averages
            for var in daily_24hr_vars:
                daily_24hr_vars[var] /= 24
            
            # Save daily data
            savenc(sfname,
                   daily_overpass_vars['so2'], daily_overpass_vars['h2o2'], daily_overpass_vars['o3'], daily_overpass_vars['oh'],
                   daily_overpass_vars['no2s'],
                   daily_24hr_vars['pm25'], daily_24hr_vars['so4'],
                   daily_24hr_vars['so2s'], daily_24hr_vars['h2o2s'], daily_24hr_vars['o3s'], daily_24hr_vars['ohs'], daily_overpass_vars['unit_conv'],
                   tLAT, tLON)
            
            # Accumulate monthly
            for var in monthly_vars:
                if var in daily_overpass_vars:
                    monthly_vars[var] += daily_overpass_vars[var]
                else:
                    monthly_vars[var] += daily_24hr_vars[var]
        
        # Calculate monthly averages
        for var in monthly_vars:
            monthly_vars[var] /= daysinmonth[Mn-1]
        
        # Save monthly data
        sfname_month = f"{OutDir}/GCHP_monthly_{SimName}{Mn:02d}.nc"
        savenc(sfname_month, *monthly_vars.values(), tLAT, tLON)
        
        # Accumulate annual
        for var in annual_vars:
            annual_vars[var] += monthly_vars[var]
    
    # Calculate annual averages
    for var in annual_vars:
        annual_vars[var] /= len(Mns)
    
    # Save annual data
    sfname_annual = f"{OutDir}/GCHP_annual_{SimName}.nc"
    savenc(sfname_annual, *annual_vars.values(), tLAT, tLON)

print("Done")