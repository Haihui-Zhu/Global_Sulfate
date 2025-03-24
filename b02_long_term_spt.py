# Long term trend of sulfate at spartan sites
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.io import loadmat
import pandas as pd
import numpy as np
from scipy.stats import linregress
# import not installed utils
import sys
sys.path.append('functions/')
from save_fig_util import savefig

def make_plot(df,site,spec,savedir,single_timeseries=True):
    # cut exceeding indices
    first_valid_index = df[site].first_valid_index()
    last_valid_index = df[site].last_valid_index()
    df = df.loc[first_valid_index:last_valid_index]

    ## Calculate the time span
    start_date = pd.to_datetime(df.index.min())
    end_date = pd.to_datetime(df.index.max())
    time_span = end_date - start_date
    days = time_span.days

    # Check if the time span is longer than 2 years
    two_years = pd.Timedelta(days=2 * 365)  # Roughly 2 years
    if time_span > two_years:
        # make a copy of df to distinguish interpolated points
        df_cp = df.copy()
        # interpolate to fill missing values
        df = df.interpolate(method='linear')
        
        # calculate linear trend
        df['date'] = pd.to_datetime(df.index) # datetinme format
        df['date_ordinal'] = pd.to_datetime(df.index).map(pd.Timestamp.toordinal) # Date number
        # coefficients = np.polyfit(df['date_ordinal'], df[site], 1)
        # slope, intercept = coefficients
        x = df['date_ordinal'].values
        X = x.reshape(-1, 1)
        y = df[site].values
        slope1, intercept1, r_value1, p_value1, std_err1 = linregress(X.flatten(), y)

        if single_timeseries is True:
            linear_trend1 = np.polyval([slope1, intercept1], df['date_ordinal'])

            fig = plt.figure(figsize=(7.5, 5))

            # plt.plot(df.index, df[site],  marker='o',color='lightgrey',markerfacecolor='none')  # Plot a single column
            plt.plot(df.index, df_cp[site],label=site, marker='o',color='grey')  # Plot a single column
            plt.plot(df.index, linear_trend1,'--', label='Linear Trend daily', color='green')
            
            # plt.xlabel('Date')
            if spec == 'SO4_frac':
                plt.ylabel(f'Sulfate Fraction (unitless)',fontsize=14 )
                # Format the coefficient label text
                label_text = f"{slope1*365:.2e}±{std_err1*365:.2e} /yr "   
                
            else:
                plt.ylabel(f'Sulfate Concentration (µg/m³)',fontsize=14)
                # Format the coefficient label text
                label_text = f"{slope1*365:.2f}±{std_err1*365:.2f} µg/m³/yr "   

            plt.text(0.02, 0.9, label_text, ha='left', va='center', transform=plt.gca().transAxes, color='green',fontsize=16,\
                     bbox=dict(facecolor='white', alpha=0.7, edgecolor='lightgrey'))

            plt.title(site)
            # plt.legend(loc='upper left')
            plt.grid()
            if site == 'Beijing':
                plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=6))  
            else:
                plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))  
            # plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))  # Format dates
            
            plt.xticks(rotation=45)  # Rotate for readability
            plt.subplots_adjust(bottom=0.19) 

            plt.tight_layout()
            plt.show()
            
            fname = f'{savedir}/single_{spec}_{site}.png'
            savefig(fname, fig)
    else:
        print(f'Skip {site}, less than two years data. ')
        slope1 = np.nan

    return slope1, days
  

def find_nearest_index(grid, value):
    idx = np.abs(grid - value).argmin()
    return idx

def make_bar_chart(indata,n,site_name,regions,region_name,ylabel, fig_name_lb, savedir):
    global_south = ['Africa','America-Central and South','Asia-East','Asia-South','Asia-Southeast','Middle East'] # regions that will be shown in the plots   
    
    indata = np.array(indata)
    regions = np.array(regions,dtype='int')
    site_name = np.array(site_name)
    n = np.array(n)

    # Remove NaN values
    valid_indices = ~np.isnan(indata)  
    indata = indata[valid_indices]
    regions = regions[valid_indices]
    site_name = site_name[valid_indices]
    n = n[valid_indices]

    # Sort the filtered data in descending order
    sorted_indices = np.argsort(indata)[::-1]
    indata = indata[sorted_indices]
    regions = regions[sorted_indices]
    site_name = site_name[sorted_indices]
    n = n[sorted_indices]

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(indata))
    for rig, tdata in enumerate(indata):
        if region_name[regions[rig]] in global_south:
            bars = ax.bar(x[rig], tdata, color='firebrick')
        else:
            bars = ax.bar(x[rig], tdata, color='steelblue')
        
        for i, bar in enumerate(bars):
            if tdata>0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{int(n[rig])}',
                        ha='center', va='bottom')
            else:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{int(n[rig])}',
                        ha='center', va='top')

    # Set x-axis labels
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(site_name, rotation=45, ha='right')

    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()
    fname = f'{savedir}/bar_{fig_name_lb}_trend.png'
    savefig(fname, fig)

def sum_site_info(so4_data, D1_Dates, site, site_info):
    year1 = 2019
    if site == 'Beijing' or site == 'Rehovot' or site == 'Abu Dhabi' or site == 'Halifax' :
        year1 = 2012
        
    # find the start and end date (cut at 2020 and 2023)
    idx = (D1_Dates[:, 0] > year1) & (D1_Dates[:, 0] < 2024)
    D1_Dates = D1_Dates[idx,:]
    so4_data = so4_data[idx]
    
    non_nan_indices = np.where(~np.isnan(so4_data))[0]
    first_date = D1_Dates[non_nan_indices[0]] if non_nan_indices.size > 0 else None
    last_date = D1_Dates[non_nan_indices[-1]] if non_nan_indices.size > 0 else None
    
    # calculate statistics
    non_nan_values = so4_data.dropna()
    
    for sid,tsite in enumerate(site_info['site']):
        if tsite is site:
            if first_date is not None:
                site_info['Start Date'][sid] = f'{first_date[0]}-{first_date[1]}-{first_date[2]}'
                site_info['End Date'][sid] = f'{last_date[0]}-{last_date[1]}-{last_date[2]}'
                site_info['N'][sid] = f'{len(non_nan_values):.0f}'
                site_info['Mean'][sid] = f'{np.nanmean(non_nan_values):.2f}'
                site_info['Median'][sid] = f'{np.nanmedian(non_nan_values):.2f}'
                site_info['std'][sid] = f'{np.nanstd(non_nan_values):.3f}'
    return site_info


# Read data
indir = '/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/'
outdir =  '/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/long-term_figures'
data_version = 'pm25_so4_202410'
simyear = 2021;
simname = 'ceds';
single_timeseries = False;

fname = f'{indir}/SPT_{data_version}.mat'
# 'TOT','Site_cities','D3_SiteCode','latitudes','longitudes','Species','DatesNum','D1_Dates')
datain = loadmat(fname)
data = datain['TOT']

Site_cities = [item[0][0] for item in datain['Site_cities']]
D3_SiteCode = [item[0][0] for item in datain['D3_SiteCode']] 
latitudes = datain['latitudes'].flatten()
longitudes = datain['longitudes'].flatten()
Species = [item[0] for item in datain['Species'][0]]
DatesNum = datain['DatesNum'].flatten()
D1_Dates = datain['D1_Dates']

site_info = {
'site_code':D3_SiteCode,
'site': Site_cities,
'latitudes': latitudes,
'longitudes':longitudes,
'Start Date': [None] * len(D3_SiteCode),
'End Date': [None] * len(D3_SiteCode),
'N': [None] * len(D3_SiteCode),
'Mean': [None] * len(D3_SiteCode),
'Median':[None] * len(D3_SiteCode),
'std': [None] * len(D3_SiteCode),
}

date_str = []
for i in range(len(D1_Dates)):
    date_str.append(f'{D1_Dates[i,0]}-{D1_Dates[i,1]:02d}-{D1_Dates[i,2]:02d}')

for spid, spec in enumerate(Species):
    if spec == 'PM2.5':
        A = data[:, spid, :]
        pm25 = pd.DataFrame(A, index=date_str, columns=site_info['site'])
        
    elif spec == 'SO4':
        B = data[:, spid, :]
        so4 = pd.DataFrame(B, index=date_str, columns=site_info['site'])

so4_frac = pd.DataFrame(B/A, index=date_str, columns=site_info['site'])

# Loop through sites and make a long-term trend plot
so4slope = []
so4days = []
pm25slope = []
pm25days = []
so4fracslope = []
so4fracdays = []
for sid, site in enumerate(site_info['site']):
    so4_data = so4[[site]]
    slope, days = make_plot(so4_data,site,'SO4',outdir,single_timeseries)
    so4slope.append( slope )
    so4days.append(days)
    
    site_info = sum_site_info(so4_data, D1_Dates, site, site_info)
    
    pm_data = pm25[[site]]
    slope, days = make_plot(pm_data,site,'PM2.5',outdir,single_timeseries)
    pm25slope.append( slope )
    pm25days.append(days)

    so4frac_data = so4_frac[[site]]
    slope, days = make_plot(so4frac_data,site,'SO4_frac',outdir,single_timeseries)
    so4fracslope.append( slope )
    so4fracdays.append(days)

# save site info as a csv
filename = f'{indir}/site_info.csv'
df = pd.DataFrame(site_info)
df.to_csv(filename, index=False)

# assign regions 
data_mask = loadmat('ceds_scale_2021to2018/mask_fao_ceds_01_emis.mat')
mask_01 = data_mask['mask_region']
mask_01 = mask_01 - 1 # so that id can be use as index
mlat_01 = data_mask['xlat']
mlon_01 = data_mask['xlon']
region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
region_name = [item[0] for item in region_name_array.flatten()]

regionid = []
for sid,lat in enumerate(latitudes):
    lon = longitudes[sid]

    lat_idx = find_nearest_index(mlat_01[:,1], lat)
    lon_idx = find_nearest_index(mlon_01[1,:], lon)
    if abs(mask_01[lat_idx, lon_idx])>=0 : # not a nan
        regionid.append(int(mask_01[lat_idx, lon_idx]))
    else:
        tsite = site_info['site'][sid]
        if tsite == 'Halifax':
            regionid.append(2)
        elif tsite == 'Singapore':
            regionid.append(7)
        else:
            print(f'region not assigned: {tsite}')


# # make bar for each region showing trends in PM2.5 and SO4
# for rgid, region in enumerate(region_name):
#     tslopes_so4 =[]
#     tslopes_pm  = []
#     tsitename = []
#     for ii in range(len(regionid)):
#         if regionid[ii]==rgid:
#             tslopes_so4.append(so4slope[ii])
#             tslopes_pm.append(pm25slope[ii])
#             tsitename.append(site_info['site'][ii])

#     # make a bar chart
#     if len(tslopes_so4)>0:
#         make_bar_chart(tslopes_so4,tslopes_pm,tsitename,region,outdir)

# make bar chart for each species
make_bar_chart(so4slope, so4days, site_info['site'],regionid,region_name, 'sulfate (ug/m3)','so4', outdir)
make_bar_chart(pm25slope,pm25days, site_info['site'],regionid,region_name, 'PM2.5 (ug/m3)','PM', outdir)
make_bar_chart(so4fracslope,so4fracdays, site_info['site'],regionid,region_name, 'sulfate fraction (unitless)','so4frac', outdir)
