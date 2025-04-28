import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import binned_statistic_2d

P = 101325
R = 8.314
T = 298.15 # will need to use ambient T and P 


def plot_histo2d(x_data, y_data, bins, ranges, xlabel, ylabel, title, fname):
    
    fig = plt.figure(figsize=(10, 8))

    # Create a colormap copy and set 'under' color to grey (for zero-count bins)
    cmap = plt.cm.viridis.copy()
    cmap.set_under('grey')

    # Create 2D histogram with 20 bins per axis and grey-masked zeros
    plt.hist2d(
        x_data, 
        y_data, 
        bins=bins, 
        range=ranges,
        cmap=cmap,
        vmin=1  # Values <1 (i.e., zeros) will use 'under' color (grey)
    )

    # Add colorbar and labels
    plt.colorbar(label='Number of Pixels')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    plt.show()

    fig.savefig(fname) 
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: '+ fname)


# Read the NetCDF file
ds = xr.open_dataset('temp/GCHP_monthly_htap_202101.nc')

# Extract variables and remove singleton dimensions (e.g., time, level if present)
so2 = ds['so2_s'].squeeze().values # molec/cm2
h2o2 = ds['h2o2_s'].squeeze().values # molec/cm2
so4 = ds['so4_s'].squeeze().values # ug/m3
unitconv = ds['unit_converter_s'].squeeze().values # molec/cm2
so2_ugm3 = so2/unitconv *P /R/T * 64 * 1e6 # Convert SO2 to ug/m3 for consistency

# Explore SO2 concentration distribution
# The 0 to 0.08e16 molec/cm2 bin takes 64092/64800 pixels.
# Try serveral numbers that is lower than that. 
source_threshold = 0.05e16
'''''
# Flatten and filter NaN/inf values
so2_flat = so2.flatten()
valid_mask = np.isfinite(so2_flat)
so2_valid = so2_flat[valid_mask]
so2_valid = so2_valid * 1e-16  # Convert to ppbv
print(f"\nTotal valid points: {len(so2_valid)}")

# Create histogram
n_bins = 20  # Adjust this number as needed
counts, bin_edges = np.histogram(so2_valid, bins=n_bins)

# Print bin information
print("Bin#\tLower Bound\tUpper Bound\tCount")
for i, (count, left, right) in enumerate(zip(counts, bin_edges[:-1], bin_edges[1:])):
    print(f"{i+1}\t{left:.2f}\t\t{right:.2f}\t\t{count}")

# Create plot
plt.figure(figsize=(10, 6))
plt.hist(so2_valid, bins=bin_edges, edgecolor='k', alpha=0.7)
plt.title('SO₂ Concentration Distribution', fontsize=14)
plt.xlabel('SO₂ Concentration (x10^16 moecl/cm2)', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.grid(True, alpha=0.3)

# Set logarithmic scale if needed
plt.yscale('log')

plt.tight_layout()
plt.savefig('figures/SO2_histogram.png', dpi=300)
plt.close()

print(f"Plot saved as SO2_histogram.png")
'''''

# Calculate the H2O2/SO2 ratio and SO4/(SO4+SO2) ratio, handling division by zero
with np.errstate(divide='ignore', invalid='ignore'):
    h2o2_so2_ratio = h2o2 / so2
    so4_ratio = so4 / (so4 + so2_ugm3)

# FIGIRE 1: SOR vs H2O2/SO2 
valid_mask = ~np.isnan(h2o2_so2_ratio) & ~np.isnan(so4_ratio) & (so2 > source_threshold)
x_data = h2o2_so2_ratio[valid_mask].flatten()
y_data = so4_ratio[valid_mask].flatten()
bins = [20,20]
ranges = [[0.1, 0.6], [0, 1]] 
xlabel = 'H2O2/SO2'
ylabel = 'SOR'
title = 'Global January Mean'
fname = 'figures/histo2d_SOR-H2O2SO2_global_202102_HTAP.png'
plot_histo2d(x_data, y_data, bins, ranges, xlabel, ylabel, title, fname)


# FIGURE 2: SO4 vs SO2
valid_mask = (so4 > 0) & (so2_ugm3 > 0) & (so2 > source_threshold)
x_data = so2_ugm3[valid_mask].flatten()
y_data = so4[valid_mask].flatten()
bins = [20,20]
ranges = [[0, 220], [0, 30]] 
xlabel = 'SO2'
ylabel = 'SO4'
title = 'Global January Mean'
fname = 'figures/histo2d_SO4-SO2_global_202102_HTAP.png'
plot_histo2d(x_data, y_data, bins, ranges, xlabel, ylabel, title, fname)







# %% FIGURE 3: SO4 vs H2O2
h2o2_ugm3 = h2o2/unitconv *P /R/T * 34 * 1e6 # Convert H2O2 to ug/m3 for consistency
mask = (so4 > 0) & (so2_ugm3 > 0) & (h2o2_ugm3 > 0) & (so2 > source_threshold)
so2_flat = so2_ugm3[mask].flatten()
so4_flat = so4[mask].flatten()
h2o2_flat = h2o2_ugm3[mask].flatten()

# Create bins for SO2 and SO4
so2_bins = np.linspace(0, 150, 20)
so4_bins = np.linspace(0, 32, 20)

# Calculate mean H2O2 in each bin
statistic, x_edge, y_edge, _ = binned_statistic_2d(
    so2_flat,
    so4_flat,
    h2o2_flat,
    statistic='mean',
    bins=[so2_bins, so4_bins]
)

fig = plt.figure(figsize=(10, 8))
cmap = plt.cm.viridis.copy()
# cmap.set_bad('grey', 1.0)  # Grey color for NaN/missing values

# Plot using pcolormesh to handle bin edges correctly
pc = plt.pcolormesh(
    x_edge,
    y_edge,
    statistic.T,
    cmap=cmap
)

# Formatting
plt.colorbar(pc, label='Mean H2O2 [mol/mol]')
# plt.xscale('log')
# plt.yscale('log')
plt.xlabel('SO2 [ug/m3]')
plt.ylabel('SO4 [ug/m3]')
plt.title('SO2 vs SO4')

# Adjust grid and layout
plt.grid(True, which='both', alpha=0.3)
plt.tight_layout()
plt.show()

fname = 'figures/histo2d_SO4-SO2-H2O2_global_202102_HTAP.png'
fig.savefig(fname) 
plt.close(fig)  # Close the plot to free memory
print('Figure saved: '+ fname)
