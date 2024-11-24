# Long term trend of sulfate at spartan sites
import pandas as pd
import json
import os


# read XRF sulfur and IC sulfate from SPARTAN database 
# compare OS (sulfur - S in sulfate) across sites -> if there is any differences between the global south and the global north

rootdir = '/storage1/fs1/rvmartin/Active'
indir = f'{rootdir}/SPARTAN-shared/Public_Data/Chemical_Filter_Data/PM25'
indir2 = f'{rootdir}/SPARTAN-shared/Public_Data/RCFM'
outfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/sulfur_2020_forward.json'

# site info
site_details = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/Site_details.csv'
df = pd.read_csv(site_details)
Site_cities = df['City'].values.tolist()
code = df['Site_Code'].values.tolist()
latitudes = df['Latitude'].values.tolist()
longitudes = df['Longitude'].values.tolist()

site_info = {'site':Site_cities, 'code':code, 'lat':latitudes, 'lon':longitudes}


# loop through site code and 
data = {scode:{'filter_id':[], 'xrf_sulfur':[], 'ic_so4':[], 'org_frac':[]} for scode in site_info['code']}

for scode in site_info['code']: 
    file = f'{indir}/{scode}_PM25_speciation.csv'
    if os.path.exists(file):
        # Read the file starting from the 4th line (skip first 3 lines)
        df = pd.read_csv(file, skiprows=3)
         
        # find sulfur and assign to data
        df_washu = df[df['Start_Year_local'] > 2019]
        sulfur = df_washu[df_washu['Parameter_Name'] == ' Sulfur PM2.5' ]
        data[scode]['filter_id'] = sulfur['Filter_ID'].values.tolist()
        data[scode]['xrf_sulfur'] = sulfur['Value'].values.tolist() # ng/m3

        # find sulfate and assign to data
        so4 = df[df['Parameter_Name'] == ' Sulfate Ion PM2.5'] 
        data[scode]['ic_so4'] = [None]*len(data[scode]['xrf_sulfur'])
        for fid, tfilter in enumerate(data[scode]['filter_id']):
            for so4id in so4.index:
                so4filter = so4['Filter_ID'][so4id]
                if so4filter == tfilter:
                    data[scode]['ic_so4'][fid] = so4['Value'][so4id]
    
    
    file = f'{indir2}/{scode}_PM25_RCFM.csv'
    if os.path.exists(file):
        # Read the file starting from the 4th line (skip first 3 lines)
        df = pd.read_csv(file, skiprows=3)
         
        # find sulfur and assign to data
        df_washu = df[df['Start_Year_local'] > 2019]
        dfrm = df_washu[df_washu['Parameter_Name'] == 'Residual Matter' ]
        rm = dfrm['Value'].mean()
        dfoc = df_washu[df_washu['Parameter_Name'] == 'OC PM2.5' ] 
        oc = dfoc['Value'].mean()
        dfpm = df_washu[df_washu['Parameter_Name'] == 'Filter PM2.5 mass' ] 
        pm = dfpm['Value'].mean()
        data[scode]['org_frac'] = (rm+oc)/pm


# save as a JSON file
with open(outfname, 'w') as json_file:
    json.dump(data, json_file, indent=4)
print(f'{outfname} saved')

outfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/site_info.json'
with open(outfname, 'w') as json_file:
    json.dump(site_info, json_file, indent=4)
print(f'{outfname} saved')
