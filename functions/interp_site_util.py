
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import griddata

def interp_site(region_file,lat_out,lon_out):

    # load mask file
    data_mask = loadmat(region_file)
    mask = data_mask['mask_region']
    mask = mask - 1 # so that id can be use as index
    lat_mesh = data_mask['xlat']
    lon_mesh = data_mask['xlon']
    region_name_array = data_mask['region_name'] # a total of 14 regions + 1 international
    region_name = [item[0] for item in region_name_array.flatten()]

    points = np.vstack([ lon_mesh.ravel(),lat_mesh.ravel()]).T
    values = mask.ravel()

    valid  = ~np.isnan(values)
    values = values[valid]
    points = points[valid,:]

    region = griddata(points, values, (lon_out,lat_out), method='nearest')

    region = np.array([int(x) for x in region])
    
    return region, region_name

