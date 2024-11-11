import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import json



# functions and paths
indir = '/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/'
savedir = '/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/figures/'

fullname = f'{indir}/compiled_pandora_so2_2021_24hr.json'
# load data
with open(fullname, 'r') as json_file:
    data = json.load(json_file)


site_info = {'site':[], 'lat':[], 'lon':[],'numofdays':[]}
for site, coordinates in data.items():
    site_info['site'].append(site)
    site_info['lat'].append(float(data[site]['lat']))
    site_info['lon'].append(float(data[site]['lon']))
    numofdays = len(data[site]['vcd'])
    site_info['numofdays'].append(numofdays)

# Create a map with cartopy
fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())  # Set longitude from -180 to 180 and latitude from -60 to 90
ax.coastlines()

# Plot sites with more than 10 days data in red
r = 0
b = 0
for i, site in enumerate(site_info['site']):
    if site_info['numofdays'][i] > 10:
        ax.scatter(site_info['lon'][i], site_info['lat'][i], color='firebrick', marker='o', s=30,  transform=ccrs.PlateCarree())
        r = r+1
    else:
        ax.scatter(site_info['lon'][i], site_info['lat'][i], color='steelblue', marker='o', s=30, transform=ccrs.PlateCarree())
        b = b+1

label = f'Sites with > 10 days data: N = {r}'
x, y = 0.6, 0.10
plt.text(x, y, label, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
            fontsize =12, color= 'firebrick')
label = f'Sites with < 10 days data: N = {b}'
x, y = 0.6, 0.05
plt.text(x, y, label, transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom',\
            fontsize =12, color= 'steelblue')

# Add a legend and gridlines
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

plt.title('Pandora Sites')
plt.show()

# Save figure
savedir = '/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/figures/'
fname = f'{savedir}/map_location.png'
fig.savefig(fname) 
plt.close(fig)  # Close the plot to free memory
print('Figure saved: '+ fname)
    
