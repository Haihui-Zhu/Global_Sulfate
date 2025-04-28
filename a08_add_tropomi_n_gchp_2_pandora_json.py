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


indir = '/pierce-scratch/haihui/Global_Sulfate/'
qa = 'aq01-21'
ystr = '2020-2024'
infname1 = f'{indir}/SO2_Pandora/compiled/compiled_pandora_so2_{ystr}_s5popt_{qa}.json'
outfname1 = f'{indir}/SO2_Pandora/compiled/compiled_pandora_so2_{ystr}_s5popt_{qa}_with_omipca_d.json'
infname2 = f'{indir}/SO2_Pandora/compiled/compiled_pandora_so2_{ystr}_24hr_{qa}.json'
outfname2 = f'{indir}/SO2_Pandora/compiled/compiled_pandora_so2_{ystr}_24hr_{qa}_with_omipca_d.json'

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
qcstr =  'CF03-SZA75-QA0'


for sid, site in enumerate(site_info['site']):
    data1[site]['omi-pca'] = [None]*len(data1[site]['date'])
    # data1[site]['tropomi-doas'] = [None]*len(data1[site]['date'])
    # data1[site]['gchp-ceds'] = [None]*len(data1[site]['date'])
    # data1[site]['gchp-edgar'] = [None]*len(data1[site]['date'])
    # data1[site]['gchp-htap'] = [None]*len(data1[site]['date'])


    data2[site]['omi-pca'] = [None]*len(data2[site]['date'])
    # data2[site]['tropomi-doas'] = [None]*len(data2[site]['date'])
    # data2[site]['gchp-ceds'] = [None]*len(data2[site]['date'])
    # data2[site]['gchp-edgar'] = [None]*len(data2[site]['date'])
    # data2[site]['gchp-htap'] = [None]*len(data2[site]['date'])


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

        omipcafname = f'{indir}/02.TROPOMI_SO2_Ref/OMI_SO2_Tesellation_CF03-SZA75-QA0/Regrid_omi_so2_{yr}{mn:02d}{dy:02d}_{qcstr}.nc' # use monthly mean so far
        # read omipca:
        ds = xr.open_dataset(omipcafname) 
        omipca = ds['SO2'].values # check variable name
        latc = ds['lat'].values # presumably all sources share the same lat lon
        lonc = ds['lon'].values
        # spatial smoothing
        omipca, latc, lonc = spatial_smoothing(omipca,latc,lonc,0.5,0.1)
        # interpolate:
        omipca = interp_site(latc,lonc,omipca,site_info)


        # populate to the compiled dictionary
        datestr0 = f'2020-{mn:02d}-{dy:02d}'
        datestr1 = f'2021-{mn:02d}-{dy:02d}'
        datestr2 = f'2022-{mn:02d}-{dy:02d}'
        datestr3 = f'2023-{mn:02d}-{dy:02d}'
        datestr4 = f'2024-{mn:02d}-{dy:02d}'
        for sid, site in enumerate(site_info['site']):
            for i, tdate in enumerate(data1[site]['date']):
                if tdate == datestr1 or tdate == datestr2 or tdate == datestr3 or tdate == datestr4 or tdate == datestr0:
                    data1[site]['omi-pca'][i] = float(omipca[sid])
                
            for i, tdate in enumerate(data2[site]['date']):
                if tdate == datestr1 or tdate == datestr2 or tdate == datestr3 or tdate == datestr4 or tdate == datestr0:
                    data2[site]['omi-pca'][i] = float(omipca[sid])

# Step 5: save as a JSON file
with open(outfname1, 'w') as json_file:
    # json_file.write(data)
    json.dump(data1, json_file, indent=4)
print(f'{outfname1} saved')

with open(outfname2, 'w') as json_file:
    # json_file.write(data)
    json.dump(data2, json_file, indent=4)
print(f'{outfname2} saved')