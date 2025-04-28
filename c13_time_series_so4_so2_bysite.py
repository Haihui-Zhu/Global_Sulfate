import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import os
from pathlib import Path

# Configuration
input_dir = './06.spartan_gchp/all_data_by_site/'
output_dir = './06.spartan_gchp/long-term_figures/spartan_cobra_gchp-ceds/'
os.makedirs(output_dir, exist_ok=True)

# Date ranges
full_range = ('2019-01-01', '2023-01-01')
focus_range = ('2021-01-01', '2021-12-31')

# Initialize data collection
all_stats = []

# Get all site files
site_files = [f for f in os.listdir(input_dir) if f.endswith('_pm25_so2_so4.csv')]

for site_file in site_files:
    site = site_file.split('_')[0]
    df = pd.read_csv(os.path.join(input_dir, site_file), parse_dates=['date'])
    
    # Filter data to date range
    mask = (df['date'] >= full_range[0]) & (df['date'] <= full_range[1])
    dff = df[mask].copy()
    
    # Filter data
    dffso4 = dff[(dff['so4'].notna())]
    
    # Requirement 5: Skip sites with <50 SO4 obs
    if dffso4['so4'].count() < 50:
        continue
    
    # Create figure
    fig = plt.figure(figsize=(7.5, 10))
    ax1 = plt.subplot(2, 1, 1)
    ax2 = plt.subplot(2, 1, 2)
    
    # Plotting function
    def plot_data(ax, data):
        # Common scatter parameters
        obs_params = {
            'so4': {'marker': 'o', 's': 30, 'c': 'black',  'label': 'Obs SO$_4$'},
            'so2': {'marker': 's', 's': 30, 'c': 'darkgrey', 'label': 'Obs SO$_2$'}
        }
        
        sim_params = {
            'so4': {'marker': 'o', 's': 30, 'facecolors': 'firebrick', 'label': 'Sim SO$_4$'},
            'so2': {'marker': 's', 's': 30, 'facecolors': 'lightcoral', 'label': 'Sim SO$_2$'}
        }

        # Create twin axis for SO2
        axr = ax.twinx()
        
        # Plot simulations as markers (not lines)
        valid_sim = data.dropna(subset=['gchp-ceds_so4', 'gchp-ceds_so2'], how='all')
        ax.scatter(valid_sim['date'], valid_sim['gchp-ceds_so4'], **sim_params['so4'])
        axr.scatter(valid_sim['date'], valid_sim['gchp-ceds_so2'], **sim_params['so2'])
        
        # Plot observations
        valid_obs = data.dropna(subset=['so4'])
        ax.scatter(valid_obs['date'], valid_obs['so4'], **obs_params['so4'])
        valid_so2_obs = data.dropna(subset=['cobra_so2'])
        axr.scatter(valid_so2_obs['date'], valid_so2_obs['cobra_so2'], **obs_params['so2'])
        
        # Formatting
        # ax.set_xlim(pd.to_datetime(date_range))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_ylabel('SO$_4$ (µg/m³)', fontsize=14)
        axr.set_ylabel('SO$_2$ (molec/cm²)', fontsize=14)
        
        return ax, axr
    
    # Plot full range
    ax_top, ax_top2 = plot_data(ax1, dff)
    
    # Plot focus range
    mask = (df['date'] >= focus_range[0]) & (df['date'] <= focus_range[1])
    dff2021 = df[mask].copy()
    ax_bot, ax_bot2 = plot_data(ax2, dff2021)
    
    # Calculate statistics
    def calc_stats(tdf, co_sampled=False):
        stats = {}
        
        # SO4
        stats['so4_obs'] = tdf['so4'].mean()
        stats['so4_sim'] = tdf['gchp-ceds_so4'].mean()
        # adding count of so4 obs
        stats['so4_obs_count'] = tdf['so4'].count()
        
        # SO2
        stats['so2_obs'] = tdf['cobra_so2'].mean()
        stats['so2_sim'] = tdf['gchp-ceds_so2'].mean()
        # adding count of so2 obs
        stats['so2_obs_count'] = tdf['cobra_so2'].count()
        
        if co_sampled:
            mask = tdf[['so4', 'cobra_so2']].notna().all(axis=1)
            stats['so4_sim_cs'] = tdf.loc[mask, 'gchp-ceds_so4'].mean()
            stats['so4_count_cs'] = tdf.loc[mask, 'so4'].count()
            
            stats['so2_sim_cs'] = tdf.loc[mask, 'gchp-ceds_so2'].mean()
            stats['so2_count_cs'] = tdf.loc[mask, 'cobra_so2'].count()
            
        return stats
    
    # Add statistics annotations
    stats_full = calc_stats(dff)
    stats_focus = calc_stats(dff2021, co_sampled=True)
    
    # Top plot text
    text = f"Mean SO$_4$ (Obs): {stats_full['so4_obs']:.2f} µg/m³ [N = {stats_full['so4_obs_count']}]\n" \
           f"Mean SO$_4$ (Sim): {stats_full['so4_sim']:.2f} µg/m³\n" \
           f"Mean SO$_2$ (Obs): {stats_full['so2_obs']:.2e} molec/cm² [N = {stats_full['so2_obs_count']}]\n" \
           f"Mean SO$_2$ (Sim): {stats_full['so2_sim']:.2e} molec/cm²"
    ax_top.text(0.02, 0.95, text, transform=ax_top.transAxes, va='top', fontsize=14)
    
    # Bottom plot text
    text = f"Obs SO$_4$: {stats_focus['so4_obs']:.2f} µg/m³ [N = {stats_focus['so4_count_cs']}]\n" \
           f"Obs SO$_2$: {stats_focus['so2_obs']:.2e} molec/cm² [N = {stats_focus['so2_count_cs']}]\n"\
           f"Co-sampled sim SO$_4$: {stats_focus['so4_sim_cs']:.2f} µg/m³\n" \
           f"Co-sampled sim SO$_2$: {stats_focus['so2_sim_cs']:.2e} molec/cm²"
    ax_bot.text(0.02, 0.95, text, transform=ax_bot.transAxes, va='top', fontsize=14)
    
    # Legend
    handles, labels = ax_top.get_legend_handles_labels()
    handles2, labels2 = ax_top2.get_legend_handles_labels()
    ax_top.legend(handles + handles2, labels + labels2, loc='lower center', 
                bbox_to_anchor=(0.5, 1), ncol=2, fontsize=14)
    
    # Final formatting
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{site}_all_long_term.png"), dpi=300)
    plt.close()
    print(f"Saved figure for {site}")
    
    
    # Calculate Oxidation Ratio (OR)
    if stats_focus['so4_count_cs']>0:
        obs_OR = 1e16 * (stats_focus['so4_obs'] / stats_focus['so2_obs'])
        sim_OR = 1e16 * (stats_focus['so4_sim_cs'] / stats_focus['so2_sim_cs'])
        if site != 'KRSE':
            all_stats.append({
                'site': site,
                'obs_OR': obs_OR,
                'sim_OR': sim_OR,
                'N': stats_focus['so4_count_cs']
            })
    
    

# Sort by observed OR
all_stats.sort(key=lambda x: x['obs_OR'], reverse=True)

# Plotting
fig, ax = plt.subplots(figsize=(15, 8))
x = np.arange(len(all_stats))
bar_width = 0.35

# Create bars
bars1 = ax.bar(x - bar_width/2, [s['obs_OR'] for s in all_stats], 
               bar_width, label='Observed OR', color='navy')
bars2 = ax.bar(x + bar_width/2, [s['sim_OR'] for s in all_stats], 
               bar_width, label='Simulated OR', color='firebrick')

# Add labels and formatting
ax.set_ylabel('Oxidation Ratio (1e16 × SO4/SO2)', fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels([s['site'] for s in all_stats], rotation=45, ha='right', fontsize=10)
ax.legend()

# Add N values above bars
for i, (bar1, bar2) in enumerate(zip(bars1, bars2)):
    height = max(bar1.get_height(), bar2.get_height())
    ax.text((bar1.get_x() + bar2.get_x())/2, height * 1.02, 
            f"N={all_stats[i]['N']}", 
            ha='center', va='bottom', fontsize=8)

# Add horizontal grid
ax.yaxis.grid(True, linestyle='--', alpha=0.7)
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'SO4-SO2-ratio_site_rank_2021_cosampled.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"Saved OR comparison plot to {output_dir}")