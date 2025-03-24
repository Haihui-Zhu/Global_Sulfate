% Sept-26-2024
% Haihui Zhu
%
% This script makes the time-series plots that show SPARTAN measurement and
% multiple GCHP sims by regions .
clear
close all

% RootDir = '/Volumes/rvmartin/Active/'; % read from compute1
RootDir = '/storage1/fs1/rvmartin/Active/'; % read from compute1
addpath(sprintf('%s/haihuizhu/1.code',RootDir))
addpath(sprintf('%s/haihuizhu/4.SPARTAN_SO4/functions',RootDir))
SaveDir = sfmkdir(sprintf('%s/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/',RootDir)); % save to compute1

masterfname = sprintf('pm25_so4_202410'); % file name to save

target_sim = {'ceds-2021', 'edgar-2021', 'htap-2018'};

ObsYear = 2019:2023; % years that spartan data can be included (WashU)
simyear = 2021;
SimName = 'ceds';
Species = {'SO4',  'NO3',  'NH4',   'OM', 'BC',  'dust', 'SS'};

% calculate monthly mean: both SPT and AERONET exist at the same day for at
% least 3 days in a month,ignore hours of samples
thres_Mn = 5; % at least X days available for a month;
thres_Hr = 24*6; % if there are filters sampled via 1-day period instead of 9-day, 6 filters might suffice;
thres_Ann = 9; % at least X months available for a year


%% data processing begins

load(sprintf('%s/SPT_%s.mat',SaveDir,masterfname)) % load the corresponding spartan data

% loading regional mask
load('./ceds_scale_2021to2018/mask_fao_ceds_005.mat')
region_id = interp2(xlon(1,:),xlat(:,1),mask_region,longitudes,latitudes,'nearest');
region_id(7) = 3;
region_id(28) = 8;
if sum(find(isnan(region_id)))>0
    error('Site not assigned to a region:%s',Site_cities(isnan(region_id)))
end

existrg = unique(region_id);
numofrg = length(existrg);
site_label_pm = cell(numofrg,1);
site_label_so4 = cell(numofrg,1);


% ==== Data prep ====
% remove unneeded years and species
cols = [findcol('PM2.5', Species)
        findcol('SO4', Species)
        findcol('SampledHour', Species) ];
rows = find(ismember(D1_Dates(:,1),ObsYear)); % select one year of data to plot
dates_selected = D1_Dates(rows,:);
sptall = squeeze( TOT(rows,cols,:) );

for sid = 1:numofrg
    tsites = find(sid);
    sitelabel1 = '';
    sitelabel2 = '';
    for sid = 1:length(tsites)
        nofobs = sum(~isnan(sptall(:,2,tsites(sid))),'all');
        if nofobs == 0
            fprintf('No SO4 measurement for %s\n',Site_cities{tsites(sid)})
            continue
        else
            sitelabel1 = sprintf('%s\n%s (%d) ', sitelabel1,Site_cities{tsites(sid)},sum(~isnan(sptall(:,1,tsites(sid))),'all'));
            sitelabel2 = sprintf('%s\n%s (%d) ', sitelabel2,Site_cities{tsites(sid)},nofobs);
        end
    end
    site_label_pm{sid} = sitelabel1;
    site_label_so4{sid} = sitelabel2;
end

% initiate Mn data
meanmnpm = NaN.*zeros(12,length(target_sim)+1,length(latitudes)); % Mn mean PM2.5 at each site
prc25pm  = NaN.*zeros(12,length(target_sim)+1,length(latitudes)); % 25% percentile
prc75pm  = NaN.*zeros(12,length(target_sim)+1,length(latitudes)); % 75% percentile
stdpm    = NaN.*zeros(length(target_sim)+1,length(latitudes));
meanmnso2= NaN.*zeros(12,length(target_sim)+1,length(latitudes)); % Mn mean SO4 at each site
prc25so2 = NaN.*zeros(12,length(target_sim)+1,length(latitudes)); % 25% percentile
prc75so2 = NaN.*zeros(12,length(target_sim)+1,length(latitudes)); % 75% percentile
stdso2   = NaN.*zeros(length(target_sim)+1,length(latitudes));

for tg = 1:length(target_sim)
    % ==== load measurement + coincident sim ====
    load(sprintf('%s/SIM_%s_%s.mat',SaveDir,masterfname,target_sim{tg}),'so4_sim','pm25_sim')

    % ==== before filter data for dates/species, repeat SIM for other years ===
    if contains(target_sim{tg}, '2018')
        rowsimyear = find(ismember(D1_Dates(:, 1), 2018)); % select one year of data to plot
    else
        rowsimyear = find(ismember(D1_Dates(:, 1), simyear)); % select one year of data to plot
    end

    for yr = ObsYear
        trow = find(ismember(D1_Dates(:, 1), yr));
        if length(trow) == 366
            trow = trow(1:365);
        elseif length(trow)<365 % last year in record
            rowsimyear = rowsimyear(1:length(trow));
        end
        pm25_sim(trow,  :) = pm25_sim(rowsimyear,  :);
        so4_sim(trow,  :)  = so4_sim(rowsimyear, :);
    end
    clear rowsimyear  yr trow


    % ==== Data prep ====
    % remove unneeded years and species [so that sim has the same size as obs]
    pm25_sim = squeeze( pm25_sim(rows,:) );
    so4_sim = squeeze( so4_sim(rows,:) );

    % ==== collect useful data ====

    for sid = 1:length(latitudes)
        % total std
        Ind = find(~isnan(sptall(:,1,sid)));
        stdpm(1,sid) = std(sptall(Ind,1,sid));
        stdpm(tg+1,sid) = std(pm25_sim(Ind,sid));
        Ind = find(~isnan(sptall(:,2,sid)));
        stdso2(1,sid) = std(sptall(Ind,2,sid));
        stdso2(tg+1,sid) = std(so4_sim(Ind,sid));

        % monthly mean, 25 prct and 75 prct
        for Mn = 1:12
            Ind = find(dates_selected(:,2) == Mn); % ignore interannual differences
            tmnPM = sptall(Ind,1,sid);
            tmnSO4 = sptall(Ind,2,sid);
            tmnPMsim = pm25_sim(Ind,sid);
            tmnSO4sim = so4_sim(Ind,sid);

            ind = find(~isnan(tmnPM));
            if length(ind)>thres_Mn
                meanmnpm(Mn,1,sid) = mean(tmnPM(ind),'omitnan');
                prc25pm(Mn,1,sid) = prctile(tmnPM(ind),25);
                prc75pm(Mn,1,sid) = prctile(tmnPM(ind),75);

                meanmnpm(Mn,tg+1,sid) = mean(tmnPMsim(ind),'omitnan');
                prc25pm(Mn,tg+1,sid) = prctile(tmnPMsim(ind),25);
                prc75pm(Mn,tg+1,sid) = prctile(tmnPMsim(ind),75);

            end

            ind = find(~isnan(tmnSO4));
            if length(ind)>thres_Mn
                meanmnso2(Mn,1,sid) = mean(tmnSO4(ind),'omitnan');
                prc25so2(Mn,1,sid) = prctile(tmnSO4(ind),25);
                prc75so2(Mn,1,sid) = prctile(tmnSO4(ind),75);

                meanmnso2(Mn,tg+1,sid) = mean(tmnSO4sim(ind),'omitnan');
                prc25so2(Mn,tg+1,sid) = prctile(tmnSO4sim(ind),25);
                prc75so2(Mn,tg+1,sid) = prctile(tmnSO4sim(ind),75);
            end

        end
    end
end


%% save data

% SO24
spec = 'SO4';
sfname = sprintf('%s/gchp_vs_spartan_%s_bysite.mat', SaveDir,spec);
save(sfname, 'meanmnso2','prc25so2', 'prc75so2', 'Site_cities','latitudes','longitudes','target_sim','stdso2')

% PM2.5
spec = 'PM2.5';
sfname = sprintf('%s/gchp_vs_spartan_%s_bysite.mat', SaveDir, spec);
save(sfname, 'meanmnpm', 'prc25pm', 'prc75pm', 'Site_cities','latitudes','longitudes','target_sim','stdpm')


