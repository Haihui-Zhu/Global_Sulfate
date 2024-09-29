function SimDir = get_sim_dir(SimName, simyear, RootDir)

switch SimName
    case 'ceds'
        if simyear == 2018
            SimDir = sprintf('%s/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/1.ceds_2018', RootDir);
        elseif simyear == 2021
            SimDir = sprintf('%s/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/4.ceds_2021', RootDir);
        end
    case 'htap'
        SimDir = sprintf('%s/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/2.htap_2018', RootDir);
    case 'edgar'
        SimDir = sprintf('%s/haihuizhu/4.SPARTAN_SO4/05.GCHP_outputs/3.edgar_2018', RootDir);
end

