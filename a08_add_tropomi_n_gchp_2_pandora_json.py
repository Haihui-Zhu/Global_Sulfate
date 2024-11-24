# This script reads the compiled pandora SO2 json file and populate TROPOMI and GCHP SO2
# 
import json
import xarray as xr
import numpy as np
from scipy.interpolate import griddata
# import not installed utils
import sys
sys.path.append('functions/')
from spatial_smoothing import spatial_smoothing


def interp_site(latd,lond,doas,site_info):
    # interpolate
    lon_mesh, lat_mesh = np.meshgrid(lond, latd)
    points = np.vstack([ lon_mesh.ravel(),lat_mesh.ravel()]).T
    values = doas.ravel()
    doas_intp = griddata(points, values, (site_info['lon'],site_info['lat']), method='nearest')
    return doas_intp




indir = '/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/'
qa = 'aq01-21'
ystr = '2020-2024'
infname1 = f'{indir}/compiled_pandora_so2_{ystr}_s5popt_{qa}.json'
outfname1 = f'{indir}/compiled_pandora_so2_{ystr}_s5popt_{qa}_with_cobra-d_gchp-d.json'
infname2 = f'{indir}/compiled_pandora_so2_{ystr}_24hr_{qa}.json'
outfname2 = f'{indir}/compiled_pandora_so2_{ystr}_24hr_{qa}_with_cobra-d_gchp-d.json'

# step 1: read pandora data
with open(infname1, 'r') as file:
    data1 = json.load(file)

with open(infname2, 'r') as file:
    data2 = json.load(file)

# step 2: put lat, lon, site name in one dic for easy call
site_info = {'site':[], 'lat':[], 'lon':[]}
for site, coordinates in data1.items():
    site_info['site'].append(site)
    site_info['lat'].append(float(data1[site]['lat']))
    site_info['lon'].append(float(data1[site]['lon']))

# step 3: loop through days in the year to read tropomi and gchp
# - load nc files
# - read data through interpolation
# - loop through sites to find populate data
yr = 2021
daysinmn  = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; 
qcstr =  'CF005-SZA40-QA75-vza40'
simname = 'ceds_2021'


for sid, site in enumerate(site_info['site']):
    # data1[site]['tropomi-doas'] = [None]*len(data1[site]['date'])
    data1[site]['tropomi-cobra'] = [None]*len(data1[site]['date'])
    data1[site]['gchp-ceds'] = [None]*len(data1[site]['date'])
    data1[site]['gchp-edgar'] = [None]*len(data1[site]['date'])
    data1[site]['gchp-htap'] = [None]*len(data1[site]['date'])


    # data2[site]['tropomi-doas'] = [None]*len(data2[site]['date'])
    data2[site]['tropomi-cobra'] = [None]*len(data2[site]['date'])
    data2[site]['gchp-ceds'] = [None]*len(data2[site]['date'])
    data2[site]['gchp-edgar'] = [None]*len(data2[site]['date'])
    data2[site]['gchp-htap'] = [None]*len(data2[site]['date'])


for mn in range(1,13):
    print(f'process mn-{mn}')

    # doasfname = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/gchp_so2_cosampled_tropomi_ceds_{yr}_noisereduced_{mn:02d}.nc'
    # # read doas:
    # ds = xr.open_dataset(doasfname) 
    # doas = ds['so2_tro'].values.T # check variable name
    # doas = doas/2.69e16 # convert unit
    # latd = ds['lat'].values 
    # lond = ds['lon'].values
    # # spatial smoothing
    # doas, latd, lond = spatial_smoothing(doas,latd,lond,0.5,0.05)
    # # interpolate:
    # doas = interp_site(latd,lond,doas,site_info)

    for dy in range(1,daysinmn[mn-1]+1):
    # for dy in range(1,2):
        # doasfname = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/Tropomi_Regrid_{yr}{mn:02d}{dy:02d}_{qcstr}.nc'
        # # read doas:
        # ds = xr.open_dataset(doasfname) 
        # doas = ds['so2'].values.T # check variable name
        # doas = doas/2.69e16 # convert unit
        # latd = ds['lat'].values 
        # lond = ds['lon'].values
        # # spatial smoothing
        # doas, latd, lond = spatial_smoothing(doas,latd,lond,0.5,0.05)
        # # interpolate:
        # doas = interp_site(latd,lond,doas,site_info)

        cobrafname = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/COBRA/download/s5p-l3grd-so2-cobra-pbl-001-day-{yr}{mn:02d}{dy:02d}-20240619.nc' # use monthly mean so far
        # read COBRA:
        ds = xr.open_dataset(cobrafname) 
        cobra = ds['SO2_column_number_density'].values # check variable name
        latc = ds['latitude'].values # presumably all sources share the same lat lon
        lonc = ds['longitude'].values
        # spatial smoothing
        cobra, latc, lonc = spatial_smoothing(cobra,latc,lonc,0.4,0.05)
        # interpolate:
        cobra = interp_site(latc,lonc,cobra,site_info)

        # GCHP-CEDS
        gchpfname = f'/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/4.ceds_2021/GCHP_SO2_SO4_PM25_ceds_{yr}_{mn:02d}{dy:02d}.nc'
        ds = xr.open_dataset(gchpfname) 
        gchp = ds['so2'].values.T # check variable name
        gchp = gchp/2.69e16 # convert unit
        latg = ds['latitude'].values 
        long = ds['longitude'].values
        # interpolate:
        ceds = interp_site(latg,long,gchp,site_info)

        # GCHP-EDGAR
        gchpfname = f'/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/5.edgar_2021/GCHP_SO2_SO4_BC_PM25_edgar_{yr}_{mn:02d}{dy:02d}.nc'
        ds = xr.open_dataset(gchpfname) 
        gchp = ds['so2'].values.T # check variable name
        gchp = gchp/2.69e16 # convert unit
        latg = ds['latitude'].values 
        long = ds['longitude'].values
        # interpolate:
        edgar = interp_site(latg,long,gchp,site_info)
        
        # GCHP-HTAP
        gchpfname = f'/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/2.htap_2018/GCHP_SO2_SO4_BC_PM25_htap_2018_{mn:02d}{dy:02d}.nc'
        ds = xr.open_dataset(gchpfname) 
        gchp = ds['so2'].values.T # check variable name
        gchp = gchp/2.69e16 # convert unit
        latg = ds['latitude'].values 
        long = ds['longitude'].values
        # interpolate:
        htap = interp_site(latg,long,gchp,site_info)

        # populate to the compiled dictionary
        datestr0 = f'2020-{mn:02d}-{dy:02d}'
        datestr1 = f'{yr}-{mn:02d}-{dy:02d}'
        datestr2 = f'2022-{mn:02d}-{dy:02d}'
        datestr3 = f'2023-{mn:02d}-{dy:02d}'
        datestr4 = f'2024-{mn:02d}-{dy:02d}'
        for sid, site in enumerate(site_info['site']):
            for i, tdate in enumerate(data1[site]['date']):
                if tdate == datestr1 or tdate == datestr2 or tdate == datestr3 or tdate == datestr4 or tdate == datestr0:
                    data1[site]['tropomi-cobra'][i] = float(cobra[sid])
                    data1[site]['gchp-ceds'][i] = float(ceds[sid])
                    data1[site]['gchp-edgar'][i] = float(edgar[sid])
                    data1[site]['gchp-htap'][i] = float(htap[sid])
                
            for i, tdate in enumerate(data2[site]['date']):
                if tdate == datestr1 or tdate == datestr2 or tdate == datestr3 or tdate == datestr4 or tdate == datestr0:
                    data2[site]['tropomi-cobra'][i] = float(cobra[sid])
                    data2[site]['gchp-ceds'][i] = float(ceds[sid])
                    data2[site]['gchp-edgar'][i] = float(edgar[sid])
                    data2[site]['gchp-htap'][i] = float(htap[sid])

# Step 5: save as a JSON file
with open(outfname1, 'w') as json_file:
    # json_file.write(data)
    json.dump(data1, json_file, indent=4)
print(f'{outfname1} saved')

with open(outfname2, 'w') as json_file:
    # json_file.write(data)
    json.dump(data2, json_file, indent=4)
print(f'{outfname2} saved')