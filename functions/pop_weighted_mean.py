
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import griddata

def pwm(mapdata, lat_coarse, lon_coarse, target_lat=None, target_lon=None):
    # Provide default values if None
    if target_lat is None:
        target_lat = [-90, 90]
    if target_lon is None:
        target_lon = [-180, 180]
        
    # Flatten lat/lon in case they're 2D
    lat_coarse = lat_coarse.flatten()
    lon_coarse = lon_coarse.flatten()
    
    # Select appropriate population dataset
    if lat_coarse[3] - lat_coarse[2] >= 0.1:
        popdata = loadmat('/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/Population-0.1.mat')
        population = popdata['npop'].T
        poplat = popdata['tLAT'].flatten()
        poplon = popdata['tLON'].flatten()
        
        # interpolate
        lon_mesh, lat_mesh = np.meshgrid(lon_coarse, lat_coarse)
        mappoints = np.vstack([ lon_mesh.ravel(),lat_mesh.ravel()]).T
        mapvalues = mapdata.ravel()
        
        lon_mesh, lat_mesh = np.meshgrid(poplon, poplat)
        poppoints = np.vstack([ lon_mesh.ravel(),lat_mesh.ravel()]).T
        population = population.ravel()
        
        mapdata = griddata(mappoints, mapvalues, poppoints, method='nearest')
        # print(mapdata.shape)
    
    # elif lat_coarse[3] - lat_coarse[2] == 0.05:
    #     print('loading 0.05 pop data')
    #     popdata = loadmat('/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/Population-0.05.mat')
    #     population = popdata['npop'].T
    #     poplat = popdata['tLAT'].flatten()
    #     poplon = popdata['tLON'].flatten()
        
    else:
        print('Pop data not available, please add')
    
    
    indices = np.where((poppoints[:,1] >= target_lat[0]) & (poppoints[:,1] <= target_lat[1]) \
                    & (poppoints[:,0] >= target_lon[0]) & (poppoints[:,0]<= target_lon[1]) )[0]
    mapdata = mapdata[np.ix_(indices)]
    # print(mapdata.shape)
    population = population[np.ix_(indices)]
    # print(population.shape)
    
    # Calculate weighted mean
    population_weighted = mapdata * population
    weighted_mean = np.nansum(population_weighted) / np.nansum(population)

    print(f"Population-weighted mean: {weighted_mean}")
    return weighted_mean
