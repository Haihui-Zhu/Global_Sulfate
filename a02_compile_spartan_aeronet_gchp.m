% July-18-2024
% Haihui Zhu

% collecting public SPARTAN data (PM2.5 and Sulfate only) and collect a one
% year gchp sim in the same data format. 

clear 
close all

% RootDir = '/Volumes/rvmartin/Active/'; % run locally % doesn't work since soft links can't be found locally
RootDir = '/storage1/fs1/rvmartin/Active/'; % run on compute1
addpath(sprintf('%s/haihuizhu/1.code',RootDir))
addpath(sprintf('%s/haihuizhu/4.SPARTAN_SO4/functions',RootDir))
SaveDir = sfmkdir(sprintf('%s/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/',RootDir)); % save to compute1

masterfname = sprintf('pm25_so4_202410'); % file name to save

% SPARTAN species
Species     = {'PM2.5','SO4',  'NO3',  'NH4',  'BC',   'OM',   'dust', 'SS','SampledHour'}; 
Specie_code = {'28101','28401','28402','28802','28202','28306','28305','28304',''};
% Data Dimensions 
DatesNum = datenum([2013 01 01]):1:today;
D1_Dates = datevec(DatesNum); 
D1_Dates = D1_Dates(:,1:3); % 'YEAR','MONTH','DATE'

% Sections
S_RawFiltered = 0; % Read raw filter-based SPARTAN data

S_gchp = 1;
    simyear = 2021;
    SimName = 'ceds';
    SimDir = get_sim_dir(SimName, simyear,RootDir);

S_tropomi = 0;
    year = 2021;
    qcstr = 'CF03-SZA50-QA75';

site_details = readtable(sprintf('%s/SPARTAN-shared/Site_Sampling/Site_details.xlsx',RootDir));
Site_cities = table2array(site_details(:,3));
latitudes = table2array(site_details(:,5));
longitudes = table2array(site_details(:,6));
timezone_GMT = table2array(site_details(:,12));


%% Read Raw filter-based SPARTAN data from compute1
% will creat a master file with a 'TOT' variable containing daily PM2.5/spec
% data from all sites
if S_RawFiltered == 1

    diary(sprintf('%s/%s_%s_record_SPT.txt',SaveDir,masterfname,datestr(today,'yyyy-mm-dd')))
    fprintf('%s \n', datestr(today))

    D3_SiteCode = table2array(site_details(:,1));

    FltData = sprintf('%s/SPARTAN-shared/Public_Data/Chemical_Filter_Data/PM25',RootDir);
    RCFMData = sprintf('%s/SPARTAN-shared/Public_Data/RCFM',RootDir);

    TOT = NaN.*zeros(length(D1_Dates),length(Species),length(D3_SiteCode));
    
    for st = 1:numel(Site_cities)
        PMNum = 0;
        % PM2.5 SO4 NO3 and NH4 from speciation
        Specfname = sprintf('%s/%s_PM25_speciation.csv',FltData,D3_SiteCode{st});
        if exist(Specfname,'file') ==2
            fprintf('Reading %s_PM25_speciation.csv\n',D3_SiteCode{st})

            ThisTable = readtable(Specfname);
            TableVars = ThisTable.Properties.VariableNames;
            
            % consider PM2.5 first
            sp = 1; spec = 'PM2.5';
            ThisData = GetSpecMat(sp,Specie_code,ThisTable,TableVars); % n x 3 mat: datenum, spec conc, and sampled hours.
            % ReDistribute only operates the first 2 columns
            TOT(:,findcol(spec, Species),st)          = ReDistribute(ThisData,DatesNum); 
            TOT(:,findcol('SampledHour', Species),st) = ReDistribute(ThisData(:,[1 3]),DatesNum);
            
            PMNum = size(ThisData,1);

            % SO4 NO3 NH4 since they are not available in RCFM files
            for sp = 2:4
                spec = sprintf('%s',Species{sp});
                ThisData = GetSpecMat(sp,Specie_code,ThisTable,TableVars);
                TOT(:,findcol(spec, Species),st) = ReDistribute(ThisData,DatesNum);
            end
            
        else
            fprintf('%s_PM25_speciation.csv not found \n',D3_SiteCode{st})
        end
        
        % BC Dust SS and OM from reconstructed
        RCFMfname = sprintf('%s/%s_PM25_RCFM.csv',RCFMData,D3_SiteCode{st});
        if exist(RCFMfname,'file') ==2
            fprintf('Reading %s_PM25_RCFM.csv\n',D3_SiteCode{st})

            ThisTable = readtable(RCFMfname);
            
            for sp = 5:numel(Species)-1
                spec = sprintf('%s',Species{sp});
                ThisData = GetSpecMat(sp,Specie_code,ThisTable,TableVars);
                TOT(:,findcol(spec, Species),st) = ReDistribute(ThisData,DatesNum);
            end
            
        else
            fprintf('%s_PM25_RCFM.csv not found \n',D3_SiteCode{st})
        end
        fprintf('%d PM2.5 measurement added for %s\n\n',PMNum,D3_SiteCode{st})
    end
    save(sprintf('%s/SPT_%s.mat',SaveDir,masterfname),...
        'TOT','Site_cities','D3_SiteCode','latitudes','longitudes','Species','DatesNum','D1_Dates')
end

%% Read simulation from GCHP 13.4.0 C90 with NIT cut to 50% and AOD add 0.04 (Haihui)
if S_gchp == 1

    diary(sprintf('%s/%s_%s_record_sim_%s-%d.txt',SaveDir,masterfname,datestr(today,'yyyy-mm-dd'), SimName,simyear))
    fprintf('%s \n', datestr(today))
    
    load(sprintf('%s/SPT_%s.mat',SaveDir,masterfname)) % load the corresponding spartan data

    SimDatesNum = datenum([simyear 01 01]) : 1 : datenum([simyear 12 31]);
    SimDatesVec = datevec(SimDatesNum);
    
    pm25_sim = NaN.*zeros(length(D1_Dates),length(D3_SiteCode));
    so4_sim = NaN.*zeros(length(D1_Dates),length(D3_SiteCode));
    so2_sim = NaN .* zeros(length(D1_Dates), length(D3_SiteCode));
    
    for d = 1:length(SimDatesNum)
        tic
        % fname = sprintf('%s/GCHP_SO2_SO4_BC_PM25_%s_%d_%.2d%.2d.nc',SimDir,SimName,simyear,SimDatesVec(d,2),SimDatesVec(d,3));
        fname = sprintf('%s/GCHP_SO2_SO4_PM25_%s_%d_%.2d%.2d.nc',SimDir,SimName,simyear,SimDatesVec(d,2),SimDatesVec(d,3));

        dd = find(DatesNum == SimDatesNum(d));
        if exist(fname,'file')
            
            lat=ncread(fname,'latitude');
            lon=ncread(fname,'longitude');

            pm25=ncread(fname,'pm25');
            pm25_sim(dd,:) = interp2(lon',lat,pm25,longitudes,latitudes);

            so4=ncread(fname,'so4');
            so4_sim(dd,:) = interp2(lon',lat,so4,longitudes,latitudes); 

            so2 = ncread(fname, 'so2');
            so2_sim(dd, :) = interp2(lon', lat, so2, longitudes, latitudes);
            
        else
            fprintf('%s not found \n',fname)
        end
        toc
    end
    fprintf('Done reading Sim from %s.\n',SimName)

    save(sprintf('%s/SIM_%s_%s-%d.mat',SaveDir,masterfname,SimName,simyear),...
        'pm25_sim','so4_sim','so2_sim','Site_cities','D3_SiteCode','latitudes','longitudes','Species','DatesNum','D1_Dates')

end

%% Read TROPOMI SO2
if S_tropomi == 1

    diary(sprintf('%s/%s_%s_record_tropomi_%s-%d.txt', SaveDir, masterfname, datestr(today, 'yyyy-mm-dd'), SimName, year))
    fprintf('%s \n', datestr(today))

    load(sprintf('%s/SPT_%s.mat', SaveDir, masterfname)) % load the corresponding spartan data

    SatDatesNum = datenum([year 01 01]):1:datenum([year 12 31]);
    SatDatesVec = datevec(SatDatesNum);

    so2_tro = NaN .* zeros(length(D1_Dates), length(D3_SiteCode));

    for d = 1:length(SatDatesNum)
        tic
        fname = sprintf('./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_%s/Tropomi_Regrid_%d%.2d%.2d_%s.nc', ...
                        qcstr, year, SatDatesVec(d, 2), SatDatesVec(d, 3), qcstr);

        dd = find(DatesNum == SatDatesNum(d));

        if exist(fname, 'file')

            lat = ncread(fname, 'lat');
            lon = ncread(fname, 'lon');

            so2 = ncread(fname, 'so2');
            so2_tro(dd, :) = interp2(lon', lat, so2, longitudes, latitudes);

        else
            fprintf('%s not found \n', fname)
        end

        toc
    end

    fprintf('Done reading Sim from tropomi %d.\n', year)

    save(sprintf('%s/TROPOMI_%s-%d.mat', SaveDir, masterfname, year), ...
        'so2_tro', 'Site_cities', 'D3_SiteCode', 'latitudes', 'longitudes','DatesNum', 'D1_Dates')

end

%%
diary off
clear S_*

%% functions


% re-distribute measured data to D1_Dates dimension 
function out = ReDistribute(ThisMat,DatesNum)
    out = NaN.* DatesNum;
    for i = 1:size(ThisMat,1)
        Ind = find(DatesNum == ThisMat(i,1));
        if isnan( out(Ind) )
            out(Ind) = ThisMat(i,2); 
        else
            out(Ind) = 0.5*( out(Ind) + ThisMat(i,2) );
        end
    end

end


% read dates and spec conc from raw table
function [ ThisData ] = GetSpecMat(sp,Specie_code,ThisTable,TableVars)

ThisPcode = ThisTable.Parameter_Code;

ThisSpec = find( ThisPcode == str2double(Specie_code{sp}) );
SampHr = table2array(ThisTable(ThisSpec,findcol('Hours_sampled', TableVars)));
SampHr_Ft = find(SampHr>3); % could change to 24 to be more strick on filter data

DateCols = findcol('Start_Year_local', TableVars) : findcol('End_hour_local', TableVars);
ThisDates = table2array(ThisTable(ThisSpec(SampHr_Ft),DateCols));
ThisMeas = table2array(ThisTable(ThisSpec(SampHr_Ft),findcol('Value', TableVars)));
SampHr = SampHr(SampHr_Ft);

ThisData = zeros(1,3);
for d = 1:size(ThisDates,1)
    Sdate = datenum(ThisDates(d,1:3));
    Edate = datenum(ThisDates(d,5:7));

    %%%% this is to avoid double counting the same day. Or give two PM2.5 measurement to the same day.
%     if ThisDates(d,4) > 12
%         Sdate = Sdate+1; % start from the next day
%     end
%     if ThisDates(d,8) < 12
%         Edate = Edate-1; % ends from the previous day
%     end

    % 2022-04-29: commented out the above code. That average later if there
    % are more than 1 measurement for the same day. 
     Dates = Sdate:Edate;
     ThisData(end+1:end+length(Dates), :) = [Dates' repmat(ThisMeas(d),[length(Dates),1]) repmat(SampHr(d)/length(Dates),[length(Dates),1])];

     % give a warning if there are filters sampled for only 1-2 day, as
     % said by Randall, those filters should be valued more.
     if length(Dates)<3
        fprintf('%d %d %d to %d %d %d sampled %d hours\n', ThisDates(d,1:3),ThisDates(d,5:7),SampHr(d))
     end
   
end
ThisData = ThisData(2:end,:); % delete the first empty row 

end


    