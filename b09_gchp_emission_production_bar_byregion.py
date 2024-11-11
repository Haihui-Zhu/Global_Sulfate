# compares sectoral emissions of SO2

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

def make_plot(awm_emi,emission_labels,ylabel):
    
    target_regions = ['Africa','America-Central and South', 'America-North','Asia Pacific','Asia-East','Asia-South','Asia-Southeast','Australasia','Middle East'] # regions that will be shown in the plots   
    
    # lines = ['-','--','-','--','--']
    # colors = ['firebrick','firebrick','goldenrod','goldenrod','steelblue']    
    colors = ['firebrick','goldenrod','steelblue','k']  # Colors for each site
    markers = ['o', '^', 's', 'p', '*', '+', 'x', 'D', 'v', '<', 'h']  # Markers for each month

    # Plotting
    fig = plt.figure(figsize=(7, 5))
    for ridx, region in enumerate(target_regions):
        for fidx, file in enumerate(emission_labels):
            eso2 = np.nanmean(awm_emi[region][file])
            plt.plot(ridx, eso2,  color=colors[fidx], marker=markers[ridx],label=file)

    plt.xticks(ticks=range(len(target_regions)), labels=target_regions, rotation=90)  # `rotation` for better label visibility
    plt.ylabel(ylabel)

    from matplotlib.lines import Line2D
    legend_elements_sites = [Line2D([0], [0], color=color, marker='o', linestyle='', label=file.upper())
                            for color, file in zip(colors,emission_labels)]
    plt.legend(handles=legend_elements_sites , bbox_to_anchor=(1.05, 1), loc='upper left',title="Data Souce")

    plt.grid(True)
    plt.tight_layout()
    plt.show()
    return fig
        
def cal_region(region_mask,data_in):
    masked_emissions = np.where(region_mask, data_in, np.nan)
    awm = np.nanmean(masked_emissions)
    return awm

save_path = '04.Seasonality_Figures/'

# load mask file
data_mask = loadmat('ceds_scale_2021to2018/mask_fao_ceds_05.mat')
mask_01 = data_mask['mask_region']
mlat_01 = data_mask['xlat']
mlon_01 = data_mask['xlon']
region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
region_name = [item[0] for item in region_name_array.flatten()]

###### SO2 emission inventories ####
data_label = ['ceds-2018', 'edgar-2018', 'htap-2018']
species = ['BC','NOx','OC','SO2','NH3']
diagvar = {'BC':['SpeciesConcVV_BCPI','SpeciesConcVV_BCPO'],
           'NOx':['SpeciesConcVV_NO2','SpeciesConcVV_NO'],
           'OC':['SpeciesConcVV_OCPI','SpeciesConcVV_OCPO'],
           'SO2':['EmisSO2_Total'],
           'NH3':['SpeciesConcVV_NH3','SpeciesConcVV_NH4'],
}
# sectors = ['ENE','IND','TRA','RCO','WST','SHP','AGR','SOL','AVA','AWB','SLV']
prodvar = ['ProdSO4fromH2O2inCloud','ProdSO4fromO2inCloudMetal','ProdSO4fromO3inCloud',\
           'ProdSO4fromO3inSeaSalt','ProdSO4fromHOBrInCloud','ProdSO4fromSRO3','ProdSO4fromSRHOBr','ProdSO4fromO3s']


# making maps for each sector
awm_emi_so2 = {spec: {region: {file: [] for file in data_label} for region in region_name} for spec in species}
awm_emi_so4 = {spec: {region: {file: [] for file in data_label} for region in region_name} for spec in species}
awm_prod_so4 = {spec: {region: {file: [] for file in data_label} for region in region_name} for spec in species}

# Load and process data
for file in data_label:
    if file == 'ceds-2021':
        fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/4.ceds_2021/OutputDir/GEOSChem.ACAG.20210501_1230z.nc4"

    elif file == 'ceds-2018':
        fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/1.ceds_2018/OutputDir/GEOSChem.ACAG.20180501_1230z.nc4"

    elif file == 'edgar-2021':
        fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/5.edgar_2021/OutputDir/GEOSChem.ACAG.20210501_1230z.nc4"

    elif file == 'edgar-2018':
        fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/3.edgar_2018/OutputDir/GEOSChem.ACAG.20180501_1230z.nc4"

    elif file == 'htap-2018':
        fname = f"~/my-projects2/5.GEOS-Chem/5.nasa_run/2.htap_2018/OutputDir/GEOSChem.ACAG.20180501_1230z.nc4"

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

    for spec in species:
        for idx, varname in enumerate(diagvar[spec]):
            if idx == 0:
                so2_emi = ds[varname].values
            else:
                so2_emi = so2_emi + ds[varname].values  

        # sum emi and pro to get 2D column total emi and prod
        if spec =='SO2':
            so2_emi = np.nansum(so2_emi[0,:,:,:],axis=0)
        else:
            so2_emi = so2_emi[0,0,:,:] 
        print(so2_emi.shape)       

        for rgid, region in enumerate(region_name):
            region_mask = np.where(mask_01 == rgid+1, 1, 0)
            
            so2_emi_mn = zoom(so2_emi, 2, order=0)

            # so2 emi
            awm = cal_region(region_mask,so2_emi_mn)
            awm_emi_so2[spec][region][file].append(awm)


# Plotting emi data
for spec in species:
    yaxislabel = f'{spec} emission/conc'
    fig = make_plot(awm_emi_so2[spec],data_label,yaxislabel)
    fname = f"{save_path}bar_{spec}_emission_byregion.png"
    savefig(fname, fig)  


