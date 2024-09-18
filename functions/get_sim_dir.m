function SimDir = get_sim_dir(SimName, simyear)

switch SimName
    case 'ceds'
        if simyear == 2018
            SimDir = '~/my-projects/4.SPARTAN_SO4/05.GCHP_outputs/1.ceds_2018';
        elseif simyear == 2021
            SimDir = '~/my-projects/4.SPARTAN_SO4/05.GCHP_outputs/4.ceds_2021';
        end
    case 'htap'
        SimDir = '~/my-projects/4.SPARTAN_SO4/05.GCHP_outputs/2.htap_2018';
    case 'edgar'
        SimDir = '~/my-projects/4.SPARTAN_SO4/05.GCHP_outputs/3.edgar_2018';
end

