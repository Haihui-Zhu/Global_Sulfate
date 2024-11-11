# Long term trend of sulfate at spartan sites
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.io import loadmat
import pandas as pd
import numpy as np

def make_plot(df,site,spec,savedir,single_timeseries=True ):
    # cut exceeding indices
    first_valid_index = df[site].first_valid_index()
    last_valid_index = df[site].last_valid_index()
    df = df.loc[first_valid_index:last_valid_index]
    # make a copy of df to distinguish interpolated points
    df_cp = df.copy()
    # interpolate to fill missing values
    df = df.interpolate(method='linear')
    
    # calculate moving average
    df['MA'] = df[site].rolling(window=90).mean()  # Adjust window size as needed
    # calculate linear trend
    df['date_ordinal'] = pd.to_datetime(df.index).map(pd.Timestamp.toordinal)
    coefficients = np.polyfit(df['date_ordinal'], df[site], 1)
    slope, intercept = coefficients

    if single_timeseries is True:
        linear_trend = np.polyval(coefficients, df['date_ordinal'])

        fig = plt.figure(figsize=(10, 8))

        plt.plot(df.index, df[site],  marker='o',color='lightgrey',markerfacecolor='none')  # Plot a single column
        plt.plot(df.index, df_cp[site],label=site, marker='o',color='grey')  # Plot a single column
        plt.plot(df.index, df['MA'], label='90-Day Moving Average', color='orange')
        plt.plot(df.index, linear_trend,'--', label='Linear Trend', color='green')
        
        # Format the coefficient label text
        label_text = f"y = {slope:.2e}x + {intercept:.2e}"   
        # Add the label text to the center-top area of the plot
        plt.text(0.5, 0.9, label_text, ha='center', va='center', transform=plt.gca().transAxes, color='green',fontsize=14)
        
        plt.xlabel('Date')
        plt.ylabel(f'{spec} concentration (µg/m³)')
        plt.title(site)
        plt.legend(loc='upper left')
        plt.grid()
        plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))  
        # plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))  # Format dates
        
        plt.xticks(rotation=45)  # Rotate for readability

        plt.show()
        
        fname = f'{savedir}/single_{spec}_{site}.png'
        savefig(fname, fig)
    return slope


def savefig(fname, fig):
    fig.savefig(fname) 
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: '+ fname)
    

def find_nearest_index(grid, value):
    idx = np.abs(grid - value).argmin()
    return idx

def make_bar_chart(so4,pm,sitename,region,outdir):

    fig = plt.figure(figsize=(10, 8))

    # Step 1: Sort the data by heights in descending order
    sorted_indices = sorted(range(len(pm)), key=lambda i: pm[i], reverse=True)
    sorted_names = [sitename[i] for i in sorted_indices]
    sorted_pm = [pm[i] for i in sorted_indices]
    sorted_so4 = [so4[i] for i in sorted_indices]

    # Step 2: Create a bar chart
    x = np.arange(len(sorted_names))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    bars1 = ax.bar(x - width/2, sorted_pm, width, label='PM2.5', color='skyblue')
    bars2 = ax.bar(x + width/2, sorted_so4, width, label='Sulfate', color='salmon')

    # Step 3: Add labels, title, and custom x-axis tick labels
    ax.set_xlabel('Site')
    ax.set_ylabel('Slope')
    ax.set_title(region)
    ax.set_xticks(x)
    ax.set_xticklabels(sorted_names)
    ax.legend()
    plt.subplots_adjust(bottom=0.15) 

    # Display the plot
    plt.xticks(rotation=45)
    plt.show()
    
    fname = f'{outdir}/bar_regional_{region}.png'
    savefig(fname, fig)

    

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
}

date_str = []
for i in range(len(D1_Dates)):
    date_str.append(f'{D1_Dates[i,0]}-{D1_Dates[i,1]:02d}-{D1_Dates[i,2]:02d}')

# Step 1: Create a MultiIndex from date and site
index = pd.MultiIndex.from_product([date_str, site_info['site_code']], names=['Date', 'Site'])
for spid, spec in enumerate(Species):
    if spec == 'PM2.5':
        A = data[:, spid, :]
        pm25 = pd.DataFrame(A, index=date_str, columns=site_info['site_code'])
        
    elif spec == 'SO4':
        A = data[:, spid, :]
        so4 = pd.DataFrame(A, index=date_str, columns=site_info['site_code'])

# Loop through sites and make a long-term trend plot
so4slope = []
pm25slope = []
for sid, site in enumerate(site_info['site_code']):
    so4_data = so4[[site]]
    so4slope.append( make_plot(so4_data,site,'SO4',outdir,single_timeseries) )
    
    pm_data = pm25[[site]]
    pm25slope.append( make_plot(pm_data,site,'PM2.5',outdir,single_timeseries) )

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


# make bar for each region showing trends in PM2.5 and SO4
for rgid, region in enumerate(region_name):
    tslopes_so4 =[]
    tslopes_pm  = []
    tsitename = []
    for ii in range(len(regionid)):
        if regionid[ii]==rgid:
            tslopes_so4.append(so4slope[ii])
            tslopes_pm.append(pm25slope[ii])
            tsitename.append(site_info['site'][ii])

    # make a bar chart
    if len(tslopes_so4)>0:
        make_bar_chart(tslopes_so4,tslopes_pm,tsitename,region,outdir)

