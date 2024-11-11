# This script reads the compiled pandora SO2 json file and populate TROPOMI and GCHP SO2
# 
import json
import xarray as xr
import numpy as np
from scipy.interpolate import griddata

def interp_site(latd,lond,doas,site_info):
    # interpolate
    lon_mesh, lat_mesh = np.meshgrid(lond, latd)
    points = np.vstack([ lon_mesh.ravel(),lat_mesh.ravel()]).T
    values = doas.ravel()
    doas_intp = griddata(points, values, (site_info['lon'],site_info['lat']), method='nearest')
    return doas_intp

def spatial_smoothing(mapdata,lat,lon,resout,resin):

    scale = int(resout/resin)
    shape1 = len(lat)
    shape2 = len(lon)

    so2_1_r = mapdata.reshape(int(shape1/scale),scale,int(shape2/scale),scale)
    so2_1_r = so2_1_r.mean(axis=(1, 3))

    lat_r = lat.reshape(int(shape1/scale),scale)
    lat_r = lat_r.mean(axis=1)
    lon_r = lon.reshape(int(shape2/scale),scale)
    lon_r = lon_r.mean(axis=1)

    # # mask out negative values
    # neg_mask = so2_1_r <0
    # so2_1_r[neg_mask] = 0

    return so2_1_r, lat_r, lon_r



indir = '/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/'
infname1 = f'{indir}/compiled_pandora_so2_2021_s5popt.json'
outfname1 = f'{indir}/compiled_pandora_so2_2021_s5popt_with_trop_gchp.json'
infname2 = f'{indir}/compiled_pandora_so2_2021_24hr.json'
outfname2 = f'{indir}/compiled_pandora_so2_2021_24hr_with_trop_gchp.json'

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
    data1[site]['tropomi-doas'] = [None]*len(data1[site]['date'])
    data1[site]['tropomi-cobra'] = [None]*len(data1[site]['date'])
    data1[site]['gchp-ceds'] = [None]*len(data1[site]['date'])

    data2[site]['tropomi-doas'] = [None]*len(data2[site]['date'])
    data2[site]['tropomi-cobra'] = [None]*len(data2[site]['date'])
    data2[site]['gchp-ceds'] = [None]*len(data2[site]['date'])


for mn in range(1,13):
    print(f'process mn-{mn}')
    cobrafname = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/COBRA/compiled/gchp_so2_cosampled_tropomi_ceds_2019_{mn:02d}.nc' # use monthly mean so far
    # read doas:
    ds = xr.open_dataset(cobrafname) 
    cobra = ds['so2_pbl'].values.T # check variable name
    cobra = cobra/2.69e16 # convert unit  
    latc = ds['latitude'].values # presumably all sources share the same lat lon
    lonc = ds['longitude'].values
    # spatial smoothing
    cobra, latc, lonc = spatial_smoothing(cobra,latc,lonc,0.4,0.05)
    # interpolate:
    cobra = interp_site(latc,lonc,cobra,site_info)

    for dy in range(1,daysinmn[mn-1]+1):
    # for dy in range(1,2):
        doasfname = f'/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_{qcstr}/Tropomi_Regrid_{yr}{mn:02d}{dy:02d}_{qcstr}.nc'
        gchpfname = f'/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/4.ceds_2021/GCHP_SO2_SO4_PM25_ceds_{yr}_{mn:02d}{dy:02d}.nc'
        
        # read doas:
        ds = xr.open_dataset(doasfname) 
        doas = ds['so2'].values.T # check variable name
        doas = doas/2.69e16 # convert unit
        latd = ds['lat'].values 
        lond = ds['lon'].values
        # spatial smoothing
        doas, latd, lond = spatial_smoothing(doas,latd,lond,0.5,0.05)
        # interpolate:
        doas = interp_site(latd,lond,doas,site_info)


        ds = xr.open_dataset(gchpfname) 
        gchp = ds['so2'].values.T # check variable name
        gchp = gchp/2.69e16 # convert unit
        latg = ds['latitude'].values 
        long = ds['longitude'].values
        # interpolate:
        gchp = interp_site(latg,long,gchp,site_info)

        # populate to the compiled dictionary
        datestr = f'{yr}-{mn:02d}-{dy:02d}'
        for sid, site in enumerate(site_info['site']):
            for i, tdate in enumerate(data1[site]['date']):
                if tdate == datestr:
                    data1[site]['tropomi-doas'][i] = doas[sid]
                    data1[site]['tropomi-cobra'][i] = cobra[sid]
                    data1[site]['gchp-ceds'][i] = gchp[sid]
                    break

            for i, tdate in enumerate(data2[site]['date']):
                if tdate == datestr:
                    data2[site]['tropomi-doas'][i] = doas[sid]
                    data2[site]['tropomi-cobra'][i] = cobra[sid]
                    data2[site]['gchp-ceds'][i] = gchp[sid]
                    break

# Step 5: save as a JSON file
with open(outfname1, 'w') as json_file:
    # json_file.write(data)
    json.dump(data1, json_file, indent=4)
print(f'{outfname1} saved')

with open(outfname2, 'w') as json_file:
    # json_file.write(data)
    json.dump(data2, json_file, indent=4)
print(f'{outfname2} saved')