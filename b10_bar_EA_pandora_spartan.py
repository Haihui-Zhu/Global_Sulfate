import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
# import not installed utils
import sys
sys.path.append('functions/')
from save_fig_util import savefig

# File path to the CSV file
csv_file = "/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/figures_2020-2024/site_mean_summary_s5popt.csv"  # Replace with your actual file name
savedir = "/storage1/fs1/rvmartin/Active/Shared/haihuizhu/SO2_Pandora/figures_2020-2024/"

# Parameters
site = ['Beijing-RADI','Seoul-SNU', 'Ulsan','Seosan','Tsukuba','Yokosuka']  # Sites of interest
columns = ['pandora', 'cobra-d', 'doas-d','gchp-ceds','gchp-edgar','gchp-htap']  # Data columns to plot
column_present = ['Pandora', 'COBRA', 'DOAS','GCHP-CEDS','GCHP-EDGAR','GCHP-HTAP']  # Data columns to plot
std_col = ['pandora-std', 'cobra-std', 'doas-std','ceds-std','edgar-std','htap-std']  # Corresponding standard deviation columns

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_file)

# Filter the DataFrame for the selected sites
filtered_df = df[df['site'].isin(site)]

# change order to the same as 'site' list
filtered_df['site'] = pd.Categorical(filtered_df['site'], categories=site, ordered=True)
filtered_df = filtered_df.sort_values('site').reset_index(drop=True)
print(filtered_df)

# Extract data for plotting
plot_data = filtered_df[columns].values  # Transpose for bar plot (columns as groups)
std_data = filtered_df[std_col].values  # Standard deviation (transpose)

# Create the bar chart
x = np.arange(len(site))  # Positions for groups (one per site)
bar_width = 0.15  # Width of each bar
colors = ['Black', 'darkgreen', 'limegreen','firebrick','goldenrod','steelblue']  # Colors for columns

fig, ax = plt.subplots(figsize=(10, 6))
plt.rcParams.update({'font.size': 16}) 

# Add bars for each column
for i, col_name in enumerate(columns):
    ax.bar(x + i * bar_width, plot_data[:, i], width=bar_width, yerr=std_data[:, i],
        capsize=5, label=column_present[i], color=colors[i])


# Add labels, legend, and title
ax.set_xticks(x + bar_width*2.5)  # Center ticks
ax.set_xticklabels(filtered_df['site'])  # Set site names as x-tick labels
ax.set_ylabel('SO2 VCD (DU)')
ax.set_title('SO2 VCD at East Asian Sites')
ax.legend()
plt.tight_layout()

# Show the plot
plt.show()
fname = f'{savedir}/bar_EastAsian_withTROPOMI.png'
savefig(fname, fig)


columns = ['pandora', 'gchp-ceds','gchp-edgar','gchp-htap']  # Data columns to plot
column_present = ['Pandora','GCHP-CEDS','GCHP-EDGAR','GCHP-HTAP']  # Data columns to plot
std_col = ['pandora-std', 'ceds-std','edgar-std','htap-std']  # Corresponding standard deviation columns

# Extract data for plotting
plot_data = filtered_df[columns].values  # Transpose for bar plot (columns as groups)
std_data = filtered_df[std_col].values  # Standard deviation (transpose)
# print(plot_data)
# print(std_data)

# Create the bar chart
x = np.arange(len(site))  # Positions for groups (one per site)
bar_width = 0.15  # Width of each bar
colors = ['Black','firebrick','goldenrod','steelblue']  # Colors for columns

fig, ax = plt.subplots(figsize=(10.5, 6))
plt.rcParams.update({'font.size': 16}) 

# Add bars for each column
for i, col_name in enumerate(columns):
    ax.bar(x + i * bar_width, plot_data[:, i], width=bar_width, yerr=std_data[:, i],
           capsize=5, label=column_present[i], color=colors[i])

# Add labels, legend, and title
ax.set_xticks(x + bar_width*1.5)  # Center ticks
ax.set_xticklabels(filtered_df['site'])  # Set site names as x-tick labels
ax.set_ylabel(r'SO$_2$ VCD (DU)')
ax.set_title('SO2 VCD at East Asian Sites')
ax.legend()
plt.tight_layout()

# Show the plot
plt.show()
fname = f'{savedir}/bar_EastAsian_noTROPOMI.png'
savefig(fname, fig)




# making plots for spartan sulfate
fname = "/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/gchp_vs_spartan_SO4_bysite.mat"
spt = loadmat(fname)
site = spt['Site_cities']
target_sites = ['Beijing','Seoul','Ulsan','Kaohsiung','Taipei']
so4 = np.nanmean(spt['meanmnso2'],axis=0)
# per25 = np.nanmean(spt['prc25so2'],axis=0)
# per75 = np.nanmean(spt['prc75so2'],axis=0)
stds = spt['stdso2']

print(stds.shape)
label = ['SPARTAN','GCHP-CEDS','GCHP-EDGAR','GCHP-HTAP']

# Extract data for plotting
plot_data = filtered_df[columns].values  # Transpose for bar plot (columns as groups)
std_data = filtered_df[std_col].values  # Standard deviation (transpose)
# print(plot_data)
# print(std_data)

# Create the bar chart
x = np.arange(len(target_sites))  # Positions for groups (one per site)
bar_width = 0.15  # Width of each bar
colors = ['Black','firebrick','goldenrod','steelblue']  # Colors for columns

fig, ax = plt.subplots(figsize=(10.5, 6))
plt.rcParams.update({'font.size': 16}) 

# Add bars for each column
n = 0
for sid, site_name in enumerate(site):
    if site_name in target_sites:
        n = n+1
        for i, col_name in enumerate(label):
            if n == 1:
                ax.bar(x[n-1] + i * bar_width, so4[i, sid], width=bar_width, yerr=stds[i, sid],
                    capsize=5, label=label[i], color=colors[i])
            else:
                ax.bar(x[n-1] + i * bar_width, so4[i, sid], width=bar_width, yerr=stds[i, sid],
                    capsize=5, color=colors[i])
                

# Add labels, legend, and title
ax.set_xticks(x + bar_width*1.5)  # Center ticks
ax.set_xticklabels(target_sites)  # Set site names as x-tick labels
ax.set_ylabel('Sulfate (µg/m³)')
ax.set_title('Surface Sulfate East Asian Sites')
ax.legend()
plt.tight_layout()

# Show the plot
plt.show()
fname = f'{savedir}/bar_EastAsian_spartan_sulfate.png'
savefig(fname, fig)

