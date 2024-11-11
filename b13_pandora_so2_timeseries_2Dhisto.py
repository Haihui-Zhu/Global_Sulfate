# Making plots to compare SO2 VCD from folowing sources
# - Pandora
# - TROPOMI
# - GCHP
# Plots include: 
# - time series for each site when there are more than 10 ground measurements
# - 2d histogram to visualize overall agreemene among sources
#   - daily
#   - site mean

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.stats import linregress

# user switches
s_timeseries = 1
s_histo = 1

# functions and paths
indir = '/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/'
savedir = '/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/figures/'
files =  ['s5popt','24hr']
mpsm2du = 6.22e23/1e4/2.69e16 # [moles per square meter] to DU


def timeseries_plot(date, pandora,  doas, cobra,gchp,site, pdr_label, savedir ):

    fig = plt.figure(figsize=(10, 8))

    plt.plot(date, pandora,  label=pdr_label, color='lightgrey')  
    plt.plot(date, doas, label='doas', color='steelblue')  
    plt.plot(date, cobra, label='cobra', color='orange')
    plt.plot(date, gchp, label='gchp-ceds', color='orchid')
    
    # calcualte correlation between sources
    label_text = get_r_str(pandora,doas)
    plt.text(0.5, 0.9, label_text, ha='center', va='center', transform=plt.gca().transAxes, color='steelblue',fontsize=14)
    # Format the coefficient label text
    label_text = get_r_str(pandora,cobra)
    plt.text(0.5, 0.85, label_text, ha='center', va='center', transform=plt.gca().transAxes, color='orange',fontsize=14)
    # Format the coefficient label text
    label_text = get_r_str(pandora,gchp)
    plt.text(0.5, 0.8, label_text, ha='center', va='center', transform=plt.gca().transAxes, color='orchid',fontsize=14)
    
    plt.xlabel('Date')
    plt.ylabel(f'SO2 VCD (DU)')
    plt.title(site)
    plt.legend(loc='upper left')
    plt.grid()
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))   
    plt.xticks(rotation=45)  # Rotate for readability

    plt.show()
    
    fname = f'{savedir}/timesereis_{pdr_label}_{site}.png'
    savefig(fname, fig)


def scatter_plots_site(x1,y1,ylabel,freq_label,savedir):
    valid = ~np.isnan(x1) & ~np.isnan(y1)  # ~ is the bitwise NOT operator, used here to invert the boolean array
    x1 = np.array(x1)
    y1 = np.array(y1)
    x = x1[valid]
    y = y1[valid]
    # Calculate the correlation coefficient
    correlation = np.corrcoef(x, y)[0, 1]

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    text_str = f"Line equation: y = {slope:.2f}x + {intercept:.2f}\nr = {correlation:.2f}"

    # Plotting the 2D histogram
    if text_str is None:
        print(f'skip {site}, less than 10 data point')
    else:
        fig = plt.figure(figsize=(10, 8))
        plt.scatter(x, y)
        
        # Plot the regression line
        x_values = np.linspace(min(x), max(x), 100)
        y_values = intercept + slope * x_values
        plt.plot(x_values, y_values, 'r-', linewidth=1.5)  # Red line for regression
        plt.plot(x_values, x_values, '--', color='grey', linewidth=1.5)  # 1:1
        
        plt.text(0.08, 0.9, text_str, ha='left', va='center', transform=plt.gca().transAxes,fontsize=14)
        
        # Labeling
        plt.xlabel('Pandora SO2 (DU)')
        plt.ylabel(ylabel)

        plt.show()
        fname = f'{savedir}/scatter_{freq_label}_{ylabel}.png'
        savefig(fname, fig)

def scatter_plots(x1,y1,ylabel,freq_label,savedir):
    valid = ~np.isnan(x1) & ~np.isnan(y1)  # ~ is the bitwise NOT operator, used here to invert the boolean array
    x1 = np.array(x1)
    y1 = np.array(y1)
    x = x1[valid]
    y = y1[valid]
    # Calculate the correlation coefficient
    correlation = np.corrcoef(x, y)[0, 1]

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    text_str = f"Line equation: y = {slope:.2f}x + {intercept:.2f}\nr = {correlation:.2f}"

    # Plotting the 2D histogram
    if text_str is None:
        print(f'skip {site}, less than 10 data point')
    else:
        fig = plt.figure(figsize=(10, 8))
        plt.hist2d(x, y, bins=30, cmap='PuRd')
        cb = plt.colorbar()  # Adds a colorbar to show the density scale
        cb.set_label('counts in bin')

        # Plot the regression line
        x_values = np.linspace(min(x), max(x), 100)
        y_values = intercept + slope * x_values
        plt.plot(x_values, y_values, 'r-', linewidth=1.5)  # Red line for regression
        plt.plot(x_values, x_values, '--', color='grey', linewidth=1.5)  #  1:1
        
        plt.text(0.08, 0.9, text_str, ha='left', va='center', transform=plt.gca().transAxes,fontsize=14)
        
        # Labeling
        plt.xlabel('Pandora SO2 (DU)')
        plt.ylabel(ylabel)

        plt.show()
        fname = f'{savedir}/histo2D_{freq_label}_{ylabel}.png'
        savefig(fname, fig)



def get_r_str(x,y):
    x = np.array([np.nan if v is None else v for v in x])
    y = np.array([np.nan if v is None else v for v in y])
    valid = ~np.isnan(x) & ~np.isnan(y)  # ~ is the bitwise NOT operator, used here to invert the boolean array
    x_clean = x[valid]
    y_clean = y[valid]
    if len(x_clean)>5:
        r = np.corrcoef(x_clean, y_clean)[0, 1]
        label_text = f"r = {r:.2f}"  
    else:
        label_text = None
    return label_text

def savefig(fname, fig):
    fig.savefig(fname) 
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: '+ fname)
    
def add_items(y1doas,y2doas, data):
    doas = np.array([np.nan if v is None else v for v in data])
    for item in doas:
        y1doas.append(item)

    y2doas.append(np.nanmean(doas))

    return y1doas, y2doas

# start processing data
for file in files:
    fullname = f'{indir}/compiled_pandora_so2_2021_{file}_with_trop_gchp.json'
    # load data
    with open(fullname, 'r') as json_file:
        data = json.load(json_file)


    site_info = {'site':[], 'lat':[], 'lon':[]}
    for site, coordinates in data.items():
        site_info['site'].append(site)
        site_info['lat'].append(float(data[site]['lat']))
        site_info['lon'].append(float(data[site]['lon']))

    # initialize histogram data
    x1 = [] # all data point
    y1doas = [] 
    y1cobra = []
    y1gchp = []
    x2 = [] # site mean
    y2doas = [] 
    y2cobra = []
    y2gchp = []

    for site in site_info['site']:
        # making time series plots
        label = f'{file}-vcd'
        savedir1 = f'{savedir}/vcd'
        pandora = [item*mpsm2du for item in data[site]['vcd'] ]
        if len(pandora)<=10:
            print(f'skip {site}, less than 10 data point')
        else:
            if s_timeseries == 1:
                timeseries_plot(data[site]['date'], pandora, \
                                data[site]['tropomi-doas'], data[site]['tropomi-cobra'], data[site]['gchp-ceds'], \
                                site,label, savedir1)
            # collecting data for histogram
            x1, x2           = add_items(x1,x2, pandora)
            y1doas, y2doas   = add_items(y1doas,y2doas, data[site]['tropomi-doas'])
            y1cobra, y2cobra = add_items(y1cobra,y2cobra, data[site]['tropomi-cobra'])
            y1gchp, y2gchp   = add_items(y1gchp,y2gchp, data[site]['gchp-ceds'])

            """
            # for debugging purpose: 
            # gchp = data[site]['gchp-ceds']
            # for v in gchp:
            #     if v is None:
            #         print(f'None in gchp data for {site}')
            # y1gchp.append(data[site]['gchp-ceds'])
            # y2gchp.append(np.nanmean(data[site]['gchp-ceds']))
            
        # label = f'{file}-mean'
        # savedir1 = f'{savedir}/mean'
        # if len(pandora)>5:
        #     pandora = [item*mpsm2du for item in data[site]['mean'] ]
        #     timeseries_plot(data[site]['date'], pandora, \
        #                     data[site]['tropomi-doas'], data[site]['tropomi-cobra'], data[site]['gchp-ceds'], \
        #                     site,label, savedir1)
          """  

    # making 2d histogram through out all data 
    if s_histo == 1:
        # daily 
        scatter_plots(x1,y1doas,'doas',f'{file}-daily',savedir)
        scatter_plots(x1,y1cobra,'cobra',f'{file}-daily',savedir)
        scatter_plots(x1,y1gchp,'gchp',f'{file}-daily',savedir)

        # site mean
        scatter_plots_site(x2,y2doas,'doas',f'{file}-site',savedir)
        scatter_plots_site(x2,y2cobra,'cobra',f'{file}-site',savedir)
        scatter_plots_site(x2,y2gchp,'gchp',f'{file}-site',savedir)
