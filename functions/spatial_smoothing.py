
import xarray as xr
import numpy as np

def spatial_smoothing(mapdata,lat,lon,resout,resin):

    scale = int(resout/resin)
    shape1 = len(lat)
    shape2 = len(lon)

    so2_1_r = mapdata.reshape(int(shape1/scale),scale,int(shape2/scale),scale)
    so2_1_r = np.nanmean(so2_1_r, axis=(1, 3))

    lat_r = lat.reshape(int(shape1/scale),scale)
    lat_r = lat_r.mean(axis=1)
    lon_r = lon.reshape(int(shape2/scale),scale)
    lon_r = lon_r.mean(axis=1)

    # # mask out negative values
    # neg_mask = so2_1_r <0
    # so2_1_r[neg_mask] = 0

    return so2_1_r, lat_r, lon_r
