import pandas as pd
import matplotlib.pyplot as plt

def savefig(fname, fig):
    fig.savefig(fname)
    plt.close(fig)  # Close the plot to free memory
    print('Figure saved: ' + fname)

def find_emission(emissions_data,country_names,ctn):
    # Merge the country names with the emissions data on the country code
    emissions_data_merged = pd.merge(emissions_data, country_names[['iso', 'Country_Name']], left_on='country', right_on='iso', how='left')

    # Filter for United Arab Emirates
    uae_emissions = emissions_data_merged[emissions_data_merged['Country_Name'] == ctn]

    # Extract the relevant years
    uae_emissions_years = uae_emissions[['X2017', 'X2018', 'X2019', 'X2020', 'X2021', 'X2022']].transpose()
    uae_emissions_years.columns = ['Emissions']
    uae_emissions_years.index = [2017, 2018, 2019, 2020, 2021, 2022]
    return uae_emissions_years




in_path = 'CEDS-2024/'
save_path = '01.Emissions/'
# countrynames = ['Saudi Arabia','Islamic Republic of Iran','Iraq','United Arab Emirates','Kuwait','Qatar']
target_countrynames = pd.read_csv(f'{in_path}Master_Country_List.csv')
#### NOTE: the code assumes country list follows the same order of emission countries. need to doubel check!####

fig = plt.figure(figsize=(6, 4))
for idx, ctn in enumerate(target_countrynames): 
    # Load the emissions data
    emissions_data = pd.read_csv(f'{in_path}CEDS_v_2024_04_01_aggregate/BC_CEDS_emissions_by_country_v2024_04_01.csv')

    # Load the country names data
    country_names = pd.read_csv(f'{in_path}Master_Country_List.csv')
    # print(ctn)
    uae_emissions_years = find_emission(emissions_data,country_names,ctn)
    
    if idx == 0:
        total_emission =  uae_emissions_years.copy()
    else:
        total_emission = total_emission + uae_emissions_years.copy()

    # Plot the emissions data for United Arab Emirates from 2017 to 2022
    plt.plot(uae_emissions_years.index, uae_emissions_years['Emissions'], marker='o', linestyle='-',label=ctn)

plt.title(f'BC Emissions China (2017-2022)')
plt.xlabel('Year')
plt.ylabel('Emissions (ktC)')
plt.grid(True)
plt.legend(loc='center left', bbox_to_anchor=(0, 0.5))
plt.show()

# save figure
fname = f"{save_path}BC_CEDS-2024_China.png"
savefig(fname, fig)

# # Plot the emissions data for United Arab Emirates from 2017 to 2022
# fig = plt.figure(figsize=(6, 4))
# plt.plot(uae_emissions_years.index, total_emission['Emissions'], marker='o', linestyle='-', color='b')
# plt.title(f'$SO_2$ Emissions near Persian Gulf (2017-2022)')
# plt.xlabel('Year')
# plt.ylabel('Emissions (ktSO2)')
# plt.grid(True)
# plt.show()

# # save figure
# fname = f"{save_path}SO2_CEDS-2024_Persian_Gulf.png"
# savefig(fname, fig)