
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.ndimage import zoom
import sparselt.esmf
import sparselt.xr

def savefig(fname, fig):
    fig.savefig(fname)
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: ' + fname)

def make_plot(diag_emi,input_emi,emission_labels,ylabel):
    
    # lines = ['-','--','-','--','--']
    # colors = ['firebrick','firebrick','goldenrod','goldenrod','steelblue']    
    colors = ['firebrick','goldenrod','steelblue','k']  # Colors for each site

    # Plotting
    fig = plt.figure(figsize=(7, 5))
    for fidx, file in enumerate(emission_labels):
        eso2 = np.nanmean(diag_emi[file])
        plt.plot(input_emi[file], eso2,  color=colors[fidx], marker='o',label=file)
    
    plt.ylabel('Diag')
    plt.xlabel('Input')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    return fig
        
def cal_region(region_mask,data_in):
    masked_emissions = np.where(region_mask, data_in, np.nan)
    awm = np.nanmean(masked_emissions)
    return awm

save_path = '04.Seasonality_Figures/'

###### SO2 emission inventories ####
data_label = ['ceds-2018', 'edgar-2018', 'htap-2018']
# species = ['BC','CO','NH3','NMVOC','NOx','OC','SO2']
sectors = ['ENE','IND','TRA','RCO','WST','SHP'] # sectors with data 
# sectors = ['ENE','IND','TRA','RCO','WST','SHP','AGR','SOL','AVA','AWB','SLV']
prodvar = ['ProdSO4fromH2O2inCloud','ProdSO4fromO2inCloudMetal','ProdSO4fromO3inCloud',\
           'ProdSO4fromO3inSeaSalt','ProdSO4fromHOBrInCloud','ProdSO4fromSRO3','ProdSO4fromSRHOBr','ProdSO4fromO3s']
spec = 'SO2'

# making maps for each sector
diag_emi_so2 = {file: [] for file in data_label}  
diag_emi_so4 = {file: [] for file in data_label}  

input_emi_so2 = {file: [] for file in data_label}  
input_emi_so4 = {file: [] for file in data_label}   

# Read GCHP diag data
for file in data_label:
    print(file)
    for Dy in range(1,2):
        for Hr in range(0,24):
            print(f"{Dy}-{Hr}")
            if file == 'ceds-2021':
                fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/4.ceds_2021/OutputDir/GEOSChem.ACAG.202105{Dy:02d}_{Hr:02d}30z.nc4"

            elif file == 'ceds-2018':
                fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/1.ceds_2018/OutputDir/GEOSChem.ACAG.201805{Dy:02d}_{Hr:02d}30z.nc4"

            elif file == 'edgar-2021':
                fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/5.edgar_2021/OutputDir/GEOSChem.ACAG.202105{Dy:02d}_{Hr:02d}30z.nc4"

            elif file == 'edgar-2018':
                fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/3.edgar_2018/OutputDir/GEOSChem.ACAG.201805{Dy:02d}_{Hr:02d}30z.nc4"

            elif file == 'htap-2018':
                fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/2.htap_2018/OutputDir/GEOSChem.ACAG.201805{Dy:02d}_{Hr:02d}30z.nc4"

            ds_cube = xr.open_dataset(fname)
            # gridspec: cube to latxlon
            weight_file = '~/my-projects/8.regrid_file/Weight_C90_1x1.nc'
            grid_in = 90
            lat_out = 180
            lon_out = 360
            # Create a linear transformer object from an ESMF weights file
            transform = sparselt.esmf.load_weights(
                    weight_file,                 
                    input_dims=[('nf', 'Ydim', 'Xdim'), (6, grid_in, grid_in)],                        
                    output_dims=[('lat', 'lon'), (lat_out, lon_out)],             
                )
            # Apply the transform to ds
            ds = sparselt.xr.apply(transform, ds_cube)

            lat  = ds['lat'].values
            lon  = ds['lon'].values
            sec_avail = ds.data_vars

            so2_emi = ds['EmisSO2_Total'].values   

            # sum emi and pro to get 2D column total emi and prod
            so2_emi = np.nanmean(so2_emi,axis=(0,2,3)) 
            so2_emi = np.nansum(so2_emi,axis=(0)) # sum through elevation   

            # so2 emi
            diag_emi_so2[file].append(so2_emi)


# Read emission input
print('reading input emission:')
for file in data_label:
    if file == 'ceds-2021':
        fname = f"./CEDS-2024-Gridded/{spec}-em-anthro_CMIP_CEDS_2021.nc"

    elif file == 'ceds-2018':
        fname = f"./CEDS-2021/2018/{spec}-em-anthro_CMIP_CEDS_2018.nc"

    elif file == 'edgar-2021':
        fname = f"./EDGAR/EDGARv81/v8.1_{spec}_2021_scaled_w_v61_0.1x0.1.nc"

    elif file == 'edgar-2018':
        fname = f"./EDGAR/EDGARv61/EDGARv6.1_{spec}_2018.0.1x0.1.nc"

    elif file == 'htap-2018':
        # fname = f"./HTAPv3/2018/HTAPv3_{spec}_0.1x0.1_2018.nc"
        fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/HTAPv3_{spec}_0.1x0.1_2018.nc" # file from pfe

    ds = xr.open_dataset(fname)
    lat  = ds['lat'].values
    lon  = ds['lon'].values
    sec_avail = ds.data_vars

    # collect only so2 emissions from desired sectors that used in simulations
    so2_emi = None 
    for tvar in sec_avail:
        for sec in sectors:
            if sec in tvar or sec.lower() in tvar:
                if so2_emi is None:
                    so2_emi = ds[tvar].values
                else:
                    so2_emi += ds[tvar].values
                print(f'{tvar} found in {file} and added')


    # sum emi and pro to get 2D column total emi and prod
    so2_emi = np.nanmean(so2_emi[4,:,:],axis=(0,1))         

    # so2 emi
    input_emi_so2[file].append(so2_emi)

# Plotting emi data
yaxislabel = f'$SO_2$ emission (kg m$^{-2}$ s$^{-1}$)'
fig = make_plot(diag_emi_so2,input_emi_so2,data_label,yaxislabel)
fname = f"{save_path}bar_{spec}_emission.png"
savefig(fname, fig)  


