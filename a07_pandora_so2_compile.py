# This script compile Pandora SO2 VCD during desired period to a single JSON file 

import glob
import json
import pandas as pd
import numpy as np
from datetime import timedelta

# Step 0: Define paths and func:
# indir = '/Volumes/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/'
indir = '/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/'
outfilename1 = f'{indir}/compiled_pandora_so2_2020-2024_s5popt_aq01-21.json'
outfilename2 = f'{indir}/compiled_pandora_so2_2020-2024_24hr_aq01-21.json'

def find_start_index(file_path):
    dash_count = 0
    with open(file_path, 'r', encoding='ISO-8859-1') as file:
        for i, line in enumerate(file):
            if line.strip() == '-' * len(line.strip()):  # Checks if the line is only dashes
                dash_count += 1
                if dash_count == 2:
                    return i + 1  # Returns the line after the second dash line
    return None

def convert_to_local_solar_time(utc_time_str, longitude):
    time_offset = timedelta(hours = float(longitude) / 15)
    # Remove the 'Z' and convert to datetime
    utc_datetime = pd.to_datetime(utc_time_str.replace('Z', ''), format='%Y%m%dT%H%M%S.%f')
    # Apply the time offset
    local_solar_time = utc_datetime + time_offset
    return local_solar_time




# Step 1: Define lists and dictionaries
site_info_string = [
    'Short location name:',
    'Location latitude [deg]:',
    'Location longitude [deg]:',
    'File name:'
]
site_info = {var: [] for var in site_info_string}

column_num = [1, 3, 4, 5, 16, 36, 39, 43] 
column_num = [item-1 for item in column_num] #  -1 because index starts from 0
column_info = [
    'UT date and time',
    'Effective duration',
    'Solar zenith angle',
    'Solar azimuth',
    'Mean value of measured data',
    'L2 data quality flag for sulfur dioxide',
    'Sulfur dioxide total vertical column amount',
    'Total uncertainty of sulfur dioxide'
]
column_info_short = ['date','duration','sza', 'saa', 'mean', 'L2qa', 'vcd', 'unc_tot']

dic_name = ['lat', 'lon', 'date', 'duration','sza', 'saa', 'mean', 'L2qa', 'vcd', 'unc_tot',\
            'tropomi-doas','tropomi-cobra','gchp-ceds']

# Step 2: Read each file
for filename in glob.glob(f'{indir}*.txt'):
    text_temp = []
    with open(filename, 'r', encoding='ISO-8859-1') as file:
        for _ in range(25):  # Limit to the first 100 lines
            line = file.readline()
            if not line:
                break
            text_temp.append(line)

    # Extract site information from the first 25 lines only
    for line in text_temp[:25]:  # Limit search to the first 25 lines
        for key in site_info_string:
            if line.startswith(key):
                site_info[key].append(line[len(key):].strip())
    """
    # check if key variables are in the expeced column (only need to do once for each version of column_info)
    # (will need to revise the code if columns not match)
    text_temp = []
    with open(filename, 'r', encoding='ISO-8859-1') as file:
        # read column information
        for _ in range(20,100):  # Limit to the lines 20 to about 100
            line = file.readline()
            if not line:
                break
            text_temp.append(line)
    col_number = []            
    for item_id, info_string in enumerate(column_info):
        for line in text_temp:
            if info_string in line:
                # Extract the number after 'Column' and before ':'
                col_number.append(int(line.split(':')[0].split(' ')[-1]))
                if col_number[item_id] != column_num[item_id]:
                    print(f"Warning for {filename}: Expected '{info_string}' in column {column_num[item_id]}, found in {col_number}")
                break
        else:
            print(f"Warning for {filename}: '{info_string}' not found")
    """

# Step 3: create a nested dictionary to store final results
site_name = site_info.get('Short location name:', 'Unknown')
data1 = {site: {var: [] for var in dic_name} for site in site_name}
data2 = {site: {var: [] for var in dic_name} for site in site_name}


# Step 4: read data and compile 
# - find data start line
# - reformate date and time
# - find time zone
# - [otptional] find data for S5P overpass time
# - quatlity checks:
#   - L2qa == 1, 2, 11, 12, 21, 22,
#   - for qualified points, using duration weighted mean for each day. 
# - assign all site info and daily data
for sid, site in enumerate(site_name):
    print(f'processing site-{sid}: {site}')
    # assigning lat, and lon
    data1[site]['lat'] = site_info['Location latitude [deg]:'][sid]
    data1[site]['lon'] = site_info['Location longitude [deg]:'][sid]
    data2[site]['lat'] = site_info['Location latitude [deg]:'][sid]
    data2[site]['lon'] = site_info['Location longitude [deg]:'][sid]
    
    tfile = site_info['File name:'][sid]
    filename = f'{indir}{tfile}'

    # Find the start index
    start_index = find_start_index(filename)

    # Read the data into a DataFrame from the start_index 
    if start_index is not None:
        data_temp = pd.read_csv(filename, skiprows=start_index, delim_whitespace=True, header=None, encoding='ISO-8859-1')

    # Convert UT to local solar time
    data_temp[0] = data_temp[0].apply(lambda x: convert_to_local_solar_time(x, data1[site]['lon']))

    # ====== [optional 1] find data within S5P overpass time (13:30) =======
    mask =  (data_temp[0].dt.hour == 13) & (data_temp[0].dt.year >= 2020) 
    data_temp_1 = data_temp[mask]
    
    # quality checks
    for i, item in enumerate(column_info):
        if item == 'L2 data quality flag for sulfur dioxide':
            tcol = column_num[i]
            mask = (data_temp_1[tcol]%10) < 1
            data_temp_1 = data_temp_1[mask]

    # duration weighted daily mean
    unique_dates = data_temp_1[0].dt.date.unique()  
    for item in unique_dates:
        mask = data_temp_1[0].dt.date==item
        t_date_temp = data_temp_1[mask][column_num]
        t_date_temp.columns = column_info_short

        # calc & populate daily data:
        data1[site][dic_name[2]].append(str(item)) # Date
        data1[site][dic_name[3]].append(np.sum(t_date_temp[dic_name[3]])) # sum of duration 
        #  dic_name = ['lat', 'lon', 'UT', 'sza', 'saa', 'mean', 'L2qa', 'vcd', 'unc_tot',\
        # 'tropomi-doas','tropomi-cobra','gchp-ceds']
        for i in range(4, len(dic_name)-3):
            data1[site][dic_name[i]].append( np.average(t_date_temp[dic_name[i]], weights=t_date_temp['duration']) )
    


    # ====== [optional 2] calculate 24 hr mean data ============================
    mask = (data_temp[0].dt.year >= 2020) 
    data_temp_2 = data_temp[mask]
    # data_temp_2 = data_temp
    
    # quality checks
    for i, item in enumerate(column_info):
        if item == 'L2 data quality flag for sulfur dioxide':
            tcol = column_num[i]
            mask = (data_temp_2[tcol]%10) < 1
            data_temp_2 = data_temp_2[mask]

    # duration weighted daily mean
    unique_dates = data_temp_2[0].dt.date.unique()  
    for item in unique_dates:
        mask = data_temp_2[0].dt.date==item
        t_date_temp = data_temp_2[mask][column_num]
        t_date_temp.columns = column_info_short

        # calc & populate daily data:
        data2[site][dic_name[2]].append(str(item)) # Date
        data2[site][dic_name[3]].append(np.sum(t_date_temp[dic_name[3]])) # sum of duration 
        #  dic_name = ['lat', 'lon', 'UT', 'sza', 'saa', 'mean', 'L2qa', 'vcd', 'unc_tot',\
        # 'tropomi-doas','tropomi-cobra','gchp-ceds']
        for i in range(4, len(dic_name)-3):
            data2[site][dic_name[i]].append( np.average(t_date_temp[dic_name[i]], weights=t_date_temp['duration']) )
    

# Step 5: save as a JSON file
with open(outfilename1, 'w') as json_file:
    json.dump(data1, json_file, indent=4)
    
with open(outfilename2, 'w') as json_file:
    json.dump(data2, json_file, indent=4)