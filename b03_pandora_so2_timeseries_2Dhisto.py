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
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.stats import linregress
# import not installed utils
import sys
sys.path.append('functions/')
from save_fig_util import savefig
from interp_site_util import interp_site

# user switches
s_timeseries = 0
s_histo = 1

# functions and paths
rootdir = '/storage1/fs1/rvmartin/Active'
indir = f'{rootdir}/Shared/haihuizhu/SO2_Pandora/'
# files =  ['24hr', 's5popt']
files = ['s5popt']
qa = 'aq01-21'
ystr = '2020-2024'
savedir = f'{rootdir}/Shared/haihuizhu/SO2_Pandora/figures_{ystr}/'
mpsm2du = 6.22e23/1e4/2.69e16 # [moles per square meter] to DU
maskfname = f'{rootdir}/haihuizhu/4.SPARTAN_SO4/ceds_scale_2021to2018/mask_fao_ceds_05.mat'

def timeseries_plot(date, pandora, cobra, doas, ceds, edgar, htap, l2qa, site, pdr_label, savedir ):


    fig, axes = plt.subplots(2, 1, figsize=(14, 16))
    # Set the global font size
    plt.rcParams.update({'font.size': 18})  # Adjust the value to set the desired font size

    
    axes[0].plot(date, pandora,  label=pdr_label, color='grey')
    axes[1].plot(date, pandora,  label=pdr_label, color='grey')
    # for pandora, use a darker color to show quality assured data:
    pandora2 = pandora.copy()
    for i, item in enumerate(l2qa):
        if item > 5: # unassured data points
            pandora2[i] = np.nan

    axes[0].plot(date, pandora2, color='k', linewidth=3) 
    axes[0].plot(date, cobra, label='cobra', color='firebrick')
    axes[0].plot(date, doas, label='doas', color='steelblue')
    
    axes[1].plot(date, pandora2, color='k', linewidth=3) 
    axes[1].plot(date, ceds, label='gchp-ceds', color='firebrick')
    axes[1].plot(date, edgar, label='gchp-edgar', color='goldenrod')
    axes[1].plot(date, htap, label='gchp-htap', color='steelblue')
    
    m = np.mean(pandora)
    pstd = np.std(pandora)
    # Calculate correlation between sources
    # For the first subplot (axes[0])
    axes[0].text(0.5, 0.95, f'Pandora M = {m:.3f} ({pstd:.4f})', ha='center', va='top', transform=axes[0].transAxes, color='k')
    label_text = get_r_str(pandora, cobra)
    axes[0].text(0.5, 0.9, label_text, ha='center', va='top', transform=axes[0].transAxes, color='firebrick')
    label_text = get_r_str(pandora, doas)
    axes[0].text(0.5, 0.85, label_text, ha='center', va='top', transform=axes[0].transAxes, color='steelblue')

    # For the second subplot (axes[1])
    axes[1].text(0.5, 0.95, f'Pandora M = {m:.3f} ({pstd:.4f})', ha='center', va='top', transform=axes[1].transAxes, color='k')
    label_text = get_r_str(pandora, ceds)
    axes[1].text(0.5, 0.9, label_text, ha='center', va='top', transform=axes[1].transAxes, color='firebrick')
    label_text = get_r_str(pandora, edgar)
    axes[1].text(0.5, 0.85, label_text, ha='center', va='top', transform=axes[1].transAxes, color='goldenrod')
    label_text = get_r_str(pandora, htap)
    axes[1].text(0.5, 0.8, label_text, ha='center', va='top', transform=axes[1].transAxes, color='steelblue')
    
    # axes[0].set_xlabel('Date')
    axes[0].set_ylabel(f'SO2 VCD (DU)')
    axes[0].legend(loc='upper left')
    axes[0].grid()
    axes[0].xaxis.set_major_locator(mdates.MonthLocator(interval=3)) 
    axes[0].tick_params(axis='x', rotation=45)  # Rotate for readability
    
    # axes[1].set_xlabel('Date')
    axes[1].set_ylabel(f'SO2 VCD (DU)')
    axes[1].legend(loc='upper left')
    axes[1].grid()
    axes[1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))  
    axes[1].tick_params(axis='x', rotation=45)  # Rotate for readability 
    
    # plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=10))   
    fig.subplots_adjust(bottom=0.15) 

    plt.show()
    
    fname = f'{savedir}/timesereis_{pdr_label}_{site}.png'
    savefig(fname, fig)

def subscatter(axe, x, y, datalb, ylabel,colors):

    # Calculate the correlation coefficient
    correlation = np.corrcoef(x, y)[0, 1]
    meanx = np.mean(x)
    meany = np.mean(y)

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    it = '-+'
    itv = abs(intercept)
    text_str = f"y = {slope:.2f}x {it[int(intercept>0)]} {itv:.2f}\nr = {correlation:.2f}\nN = {len(x)}"
    text_2 = f"{datalb}\nPandora M = {meanx:.3f}\n{ylabel} M = {meany:.3f}"
    
    # make the scatter
    axe.scatter(x, y,  s=20, color=colors[2])
    

    # Plot the regression line
    x_values = np.linspace(-0.5,2.1, 10)
    y_values = intercept + slope * x_values
    axe.plot(x_values, y_values, 'r-', linewidth=1.5)  # Red line for regression
    axe.plot(x_values, x_values, '--', color='grey', linewidth=1.5)  # 1:1

    # axe.text(0.05, 0.95, text_str, ha='left', va='top', transform=plt.gca().transAxes)
    # axe.text(0.9, 0.05, text_2, ha='right', va='bottom', transform=plt.gca().transAxes) 
    # Add text annotations, compare to the two line above for difference between plot and subplot
    axe.text(0.05, 0.95, text_str, ha='left', va='top', transform=axe.transAxes, fontsize=14)
    axe.text(0.9, 0.05, text_2, ha='right', va='bottom', transform=axe.transAxes, fontsize=14)

    # Set x and y axis limits
    if max(y)>1:
        axe.set_xlim(-0.3, 2.1)  # Set the x-axis limits
        axe.set_ylim(-0.3, 2.1)  # Set the y-axis limits
    else:
        axe.set_xlim(-0.1, 1.1)  # Set the x-axis limits
        axe.set_ylim(-0.1, 1.1)  # Set the y-axis limits

    axe.set_xlabel(f'Pandora SO_2 (DU)')
    axe.set_ylabel(f'{ylabel} SO_2 (DU)')

def scatter_plots_site(x1,y1,regions,region_name,ylabel,freq_label,savedir):
    valid = ~np.isnan(x1) & ~np.isnan(y1)  # ~ is the bitwise NOT operator, used here to invert the boolean array
    x1 = np.array(x1)
    y1 = np.array(y1)
    x = x1[valid]
    y = y1[valid]
    regions = regions[valid]
    
    # figure style
    colors = ['firebrick','goldenrod','steelblue','k']  
    
    # subplot 1 NA sites
    nax = []
    nay = []
    for rg, region in enumerate(regions):
        if region_name[region] == 'America-North':
            nax.append(x[rg])
            nay.append(y[rg])

    # subplot 1 EU sites
    eux = []
    euy = []
    for rg, region in enumerate(regions):
        if np.char.find(region_name[region], "Europe") >= 0:
            eux.append(x[rg])
            euy.append(y[rg])

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    # Set the global font size
    plt.rcParams.update({'font.size': 14})  # Adjust the value to set the desired font size
   
    subscatter(axes[0], nax, nay, 'North America',ylabel,colors) 
    subscatter(axes[1], eux, euy, 'Europe',ylabel,colors)

    # start final plot
    target_regions = ['Africa','America-Central and South','Asia-East','Asia-South','Asia-Southeast','Middle East','Asia Pacific'] # regions that will be shown in the plots   
    
    for rid, region in enumerate(region_name):
        tx = x[regions==rid]
        ty = y[regions==rid]
        # plt.scatter(x, y)
        if region in target_regions: # global south
            axes[2].scatter(tx, ty,  s=20,  color=colors[0] )
        else:
            axes[2].scatter(tx, ty,  s=20,  color=colors[2] )

    # Calculate the correlation coefficient
    correlation = np.corrcoef(x, y)[0, 1]
    meanx = np.mean(x)
    meany = np.mean(y)

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    it = '-+'
    itv = abs(intercept)
    text_str = f"y = {slope:.2f}x {it[int(intercept>0)]} {itv:.2f}\nr = {correlation:.2f}\nN = {len(x)}"
    text_2 = f"All sites\nPandora M = {meanx:.3f}\n{ylabel} M = {meany:.3f}"
    

    # Plot the regression line
    x_values = np.linspace(-0.5,2.1, 10)
    y_values = intercept + slope * x_values
    axes[2].plot(x_values, y_values, 'r-', linewidth=1.5)  # Red line for regression
    axes[2].plot(x_values, x_values, '--', color='grey', linewidth=1.5)  # 1:1
    
    axes[2].text(0.05, 0.95, text_str, ha='left', va='top', transform=plt.gca().transAxes )
    axes[2].text(0.9, 0.05, text_2, ha='right', va='bottom', transform=plt.gca().transAxes )
    
    # # Set x and y axis limits
    # if max(y)>1:
    #     axes[2].set_xlim(-0.3, 2.1)  # Set the x-axis limits
    #     axes[2].set_ylim(-0.3, 2.1)  # Set the y-axis limits
    # else:
    #     axes[2].set_xlim(-0.1, 1.1)  # Set the x-axis limits
    #     axes[2].set_ylim(-0.1, 1.1)  # Set the y-axis limits

    # Labeling
    axes[2].set_xlabel(f"Pandora $SO_2$ (DU)")
    axes[2].set_ylabel(f"{ylabel} $SO_2$ (DU)")
    # plt.subplots_adjust(right=0.6) 
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.show()
    fname = f'{savedir}/scatter_{freq_label}_{ylabel}.png'
    savefig(fname, fig)
    
    
def scatter_plots_site2(x1,y1,regions,region_name,ylabel,freq_label,savedir):
    valid = ~np.isnan(x1) & ~np.isnan(y1)  # ~ is the bitwise NOT operator, used here to invert the boolean array
    x1 = np.array(x1)
    y1 = np.array(y1)
    x = x1[valid]
    y = y1[valid]
    regions = regions[valid]
    
    # figure style
    colors = ['firebrick','goldenrod','steelblue','k']  
    
    fig,ax = plt.subplots(figsize=(5, 5))
    # Set the global font size
    plt.rcParams.update({'font.size': 14})  # Adjust the value to set the desired font size
   
    # start final plot
    target_regions = ['Africa','America-Central and South','Asia-East','Asia-South','Asia-Southeast','Middle East','Asia Pacific'] # regions that will be shown in the plots   
    
    for rid, region in enumerate(region_name):
        tx = x[regions==rid]
        ty = y[regions==rid]
        # plt.scatter(x, y)
        if region in target_regions: # global south
            ax.scatter(tx, ty,  s=20,  color=colors[0] )
        else:
            ax.scatter(tx, ty,  s=20,  color=colors[2] )

    # Calculate the correlation coefficient
    correlation = np.corrcoef(x, y)[0, 1]
    correlation = correlation**2
    meanx = np.mean(x)
    meany = np.mean(y)

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    it = '-+'
    itv = abs(intercept)
    text_str = f"y = {slope:.2f}x {it[int(intercept>0)]} {itv:.2f}\n$r^2$ = {correlation:.2f}\nN = {len(x)}"
    text_2 = f"Pandora M = {meanx:.3f}\n{ylabel} M = {meany:.3f}"
    

    # Plot the regression line
    x_values = np.linspace(-0.5,2.1, 10)
    y_values = intercept + slope * x_values
    ax.plot(x_values, y_values, 'r-', linewidth=1.5)  # Red line for regression
    ax.plot(x_values, x_values, '--', color='grey', linewidth=1.5)  # 1:1
    
    ax.text(0.05, 0.95, text_str, ha='left', va='top', transform=plt.gca().transAxes )
    ax.text(0.9, 0.05, text_2, ha='right', va='bottom', transform=plt.gca().transAxes )
    
    # Set x and y axis limits
    ax.set_xlim(-0.4, 2)  # Set the x-axis limits
    ax.set_ylim(-0.4, 2)  # Set the y-axis limits

    # Labeling
    ax.set_xlabel(f"Pandora $SO_2$ (DU)")
    ax.set_ylabel(f"{ylabel} $SO_2$ (DU)")
    plt.subplots_adjust(left=0.2,bottom=0.2) 
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.show()
    fname = f'{savedir}/scatter_{freq_label}_{ylabel}_single.png'
    savefig(fname, fig)

def histo2d(x1,y1,ylabel,freq_label,savedir):
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
        # Define bin size
        bin_size = 0.05 
        x_bins = np.arange(min(x), max(x) + bin_size, bin_size)
        y_bins = np.arange(min(y), max(y) + bin_size, bin_size)

        # Create a 2D histogram with specified bin edges
        plt.hist2d(x, y, bins=[x_bins, y_bins], cmap='PuRd')
        # Set x and y axis limits
        plt.xlim(-0.5, 1.5)  # Set the x-axis limits
        plt.ylim(-0.5, 1.5)  # Set the y-axis limits
        cb = plt.colorbar()  # Adds a colorbar to show the density scale
        cb.set_label('counts in bin')
        # Set the global font size
        plt.rcParams.update({'font.size': 16})  # Adjust the value to set the desired font size


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
    m = np.mean(y_clean)
    tstd = np.std(y_clean)
    if len(x_clean)>5:
        r = np.corrcoef(x_clean, y_clean)[0, 1]
        label_text = f"r = {r:.2f} M = {m:.3f} ({tstd:.4f})"  
    else:
        label_text = None
    return label_text

def add_items(y1doas,y2doas,tstd, data):
    doas = np.array([np.nan if v is None else v for v in data])
    for item in doas:
        y1doas.append(item)

    y2doas.append(np.nanmean(doas))
    tstd.append(np.nanstd(doas))

    return y1doas, y2doas, tstd

# start processing data
for file in files:
    fullname = f'{indir}/compiled_pandora_so2_{ystr}_{file}_{qa}_with_trop-d_gchp-d.json'
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
    y1cobra = []
    y1doas = []
    y1ceds = []
    y1edgar = []
    y1htap = []
    x2 = [] # site mean
    x2std = []
    y2cobra = []
    y2cobrastd = []
    y2doas = []
    y2doasstd = []
    y2ceds = []
    y2cedsstd = []
    y2edgar = []
    y2edgarstd = []
    y2htap = []
    y2htapstd = []
    

    site_out = []
    lat_out = []
    lon_out = []

    for sid, site in enumerate(site_info['site']):
        # making time series plots
        label = f'{file}-vcd'
        savedir1 = f'{savedir}/vcd'
        pandora = [item*mpsm2du for item in data[site]['vcd'] ]
        if len(pandora)<=10:
            print(f'skip {site}, less than 10 data point')
        else:
            if s_timeseries == 1:
                timeseries_plot(data[site]['date'], pandora, \
                                data[site]['tropomi-cobra'],  data[site]['tropomi-doas'],\
                                data[site]['gchp-ceds'], data[site]['gchp-edgar'], data[site]['gchp-htap'],\
                                data[site]['L2qa'], site,label, savedir1)
                
            # collecting data for histogram
            x1, x2, x2std      = add_items(x1,x2,x2std, pandora)
            y1cobra, y2cobra, y2cobrastd = add_items(y1cobra,y2cobra, y2cobrastd, data[site]['tropomi-cobra'])
            y1doas, y2doas,   y2doasstd   = add_items(y1doas, y2doas,y2doasstd, data[site]['tropomi-doas'])
            y1ceds, y2ceds,   y2cedsstd   = add_items(y1ceds, y2ceds,y2cedsstd, data[site]['gchp-ceds'])
            y1edgar, y2edgar, y2edgarstd = add_items(y1edgar, y2edgar, y2edgarstd, data[site]['gchp-edgar'])
            y1htap, y2htap,   y2htapstd   = add_items(y1htap, y2htap,y2htapstd, data[site]['gchp-htap'])

            site_out.append(site)
            lat_out.append(site_info['lat'][sid])
            lon_out.append(site_info['lon'][sid])


    # making 2d histogram through out all data 
    if s_histo == 1:
        # daily 
        histo2d(x1,y1cobra,'cobra',f'{file}-daily',savedir)
        histo2d(x1,y1doas,'doas',f'{file}-daily',savedir)
        histo2d(x1,y1ceds,'gchp-ceds',f'{file}-daily',savedir)
        histo2d(x1,y1edgar,'gchp-edgar',f'{file}-daily',savedir)
        histo2d(x1,y1htap,'gchp-htap',f'{file}-daily',savedir)

        # assign regions
        regions,region_name = interp_site(maskfname,lat_out,lon_out)

        # site mean
        scatter_plots_site(x2,y2cobra,regions,region_name,'cobra',f'{file}-site',savedir)
        scatter_plots_site(x2,y2doas,regions,region_name,'doas',f'{file}-site',savedir)
        scatter_plots_site(x2,y2ceds,regions,region_name,'gchp-ceds',f'{file}-site',savedir)
        scatter_plots_site(x2,y2edgar,regions,region_name,'gchp-edgar',f'{file}-site',savedir)
        scatter_plots_site(x2,y2htap,regions,region_name,'gchp-htap',f'{file}-site',savedir)
        
        scatter_plots_site2(x2,y2cobra,regions,region_name,'cobra',f'{file}-site',savedir)
        scatter_plots_site2(x2,y2doas,regions,region_name,'doas',f'{file}-site',savedir)
        scatter_plots_site2(x2,y2ceds,regions,region_name,'gchp-ceds',f'{file}-site',savedir)
        scatter_plots_site2(x2,y2edgar,regions,region_name,'gchp-edgar',f'{file}-site',savedir)
        scatter_plots_site2(x2,y2htap,regions,region_name,'gchp-htap',f'{file}-site',savedir)

    # save site mean to excel
    data = pd.DataFrame({'site': site_out,'lat':lat_out, 'lon': lon_out, 'pandora': x2, 'pandora-std':x2std, \
                        'cobra-d': y2cobra, 'cobra-std':y2cobrastd, 'doas-d': y2doas, 'doas-std': y2doasstd, \
                        'gchp-ceds': y2ceds, 'ceds-std':y2cedsstd, \
                        'gchp-edgar': y2edgar,'edgar-std':y2edgarstd, \
                        'gchp-htap': y2htap, 'htap-std':y2htapstd})
    data.to_csv(f'{savedir}site_mean_summary_{file}.csv', index=False)  
