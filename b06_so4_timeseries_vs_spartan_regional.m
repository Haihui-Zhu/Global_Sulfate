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

for rgn = 1:numofrg
    tsites = find(region_id==existrg(rgn));
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
    site_label_pm{rgn} = sitelabel1;
    site_label_so4{rgn} = sitelabel2;
end

% initiate Mn data
meanmnpm = NaN.*zeros(12,length(target_sim)+1,numofrg); % Mn mean PM2.5 at each site
prc25pm  = NaN.*zeros(12,length(target_sim)+1,numofrg); % 25% percentile
prc75pm  = NaN.*zeros(12,length(target_sim)+1,numofrg); % 75% percentile
meanmnso2= NaN.*zeros(12,length(target_sim)+1,numofrg); % Mn mean SO4 at each site
prc25so2 = NaN.*zeros(12,length(target_sim)+1,numofrg); % 25% percentile
prc75so2 = NaN.*zeros(12,length(target_sim)+1,numofrg); % 75% percentile

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

    for rgn = 1:numofrg

        if isempty(site_label_pm{rgn}) % if no sites have any data
            continue
        end

        for Mn = 1:12
            Ind = find(dates_selected(:,2) == Mn); % ignore interannual differences
            tmnPM = sptall(Ind,1,region_id==existrg(rgn));
            tmnSO4 = sptall(Ind,2,region_id==existrg(rgn));
            tmnPMsim = pm25_sim(Ind,region_id==existrg(rgn));
            tmnSO4sim = so4_sim(Ind,region_id==existrg(rgn));

            ind = find(~isnan(tmnPM));
            if length(ind)>thres_Mn
                meanmnpm(Mn,1,rgn) = mean(tmnPM(ind),'omitnan');
                prc25pm(Mn,1,rgn) = prctile(tmnPM(ind),25);
                prc75pm(Mn,1,rgn) = prctile(tmnPM(ind),75);

                meanmnpm(Mn,tg+1,rgn) = mean(tmnPMsim(ind),'omitnan');
                prc25pm(Mn,tg+1,rgn) = prctile(tmnPMsim(ind),25);
                prc75pm(Mn,tg+1,rgn) = prctile(tmnPMsim(ind),75);

            end

            ind = find(~isnan(tmnSO4));
            if length(ind)>thres_Mn
                meanmnso2(Mn,1,rgn) = mean(tmnSO4(ind),'omitnan');
                prc25so2(Mn,1,rgn) = prctile(tmnSO4(ind),25);
                prc75so2(Mn,1,rgn) = prctile(tmnSO4(ind),75);

                meanmnso2(Mn,tg+1,rgn) = mean(tmnSO4sim(ind),'omitnan');
                prc25so2(Mn,tg+1,rgn) = prctile(tmnSO4sim(ind),25);
                prc75so2(Mn,tg+1,rgn) = prctile(tmnSO4sim(ind),75);
            end

        end
    end
end


%% making the time series plots

% SO2
spec = 'SO4';
make_timeseries(meanmnso2,prc25so2,prc75so2,region_name,existrg, site_label_so4, target_sim,spec,SaveDir)
sfname = sprintf('%s/gchp_vs_spartan_%s_byregion.mat', SaveDir,spec);
save(sfname, 'meanmnso2','prc25so2', 'prc75so2', 'region_name', 'existrg','site_label_so4','target_sim')

% PM2.5
spec = 'PM2.5';
make_timeseries(meanmnpm,prc25pm,prc75pm,region_name,existrg, site_label_pm, target_sim,spec,SaveDir)
sfname = sprintf('%s/gchp_vs_spartan_%s_byregion.mat', SaveDir, spec);
save(sfname, 'meanmnpm', 'prc25pm', 'prc75pm', 'region_name', 'existrg', 'site_label_pm', 'target_sim')


%% FUNCTIONS
function make_timeseries(meanmnall,prc25all,prc75all,region_name, existrg, site_label_so4, target_sim, spec, SaveDir)
Color =[0 0 0;
    0.6350 0.0780 0.1840;
    0.9290 0.6940 0.1250;
    0 0.4470 0.7410]; 
lines = {'-','--','-.',':'};
fz = 14; % fontsize
lw = 1.2; % line width
fal = 0.25; % confidence face alpha
legends = [];

c1 = 0.05; c2 = 0.35; c3 = 0.65;
l1 = 0.65;
l2 = 0.35;
l3 = 0.05;
fw = 0.27;
fh = 0.25;
% the order is strange because we want to show ME as first and other GS in the first 6, and GN at the end
positions =[c2,l1, fw,fh; % AF
            c3, l1, fw, fh; % CAm
            c2, l3, fw, fh; % NAm
            c1, l3, fw, fh; % As P
            c1, l2, fw, fh; % As E
            c2, l2, fw, fh; % As S
            c3, l2, fw, fh; % As SE
            c3, l3, fw, fh; % Au
            c1, l1, fw, fh; ]; % ME

figure('Position',[0, 0, 1200, 1200])

for rgn = 1:length(existrg)
    sitelabel = site_label_so4{rgn};
    meandata = meanmnall(:,:,rgn);
    prc25 = prc25all(:,:,rgn);
    prc75 = prc75all(:,:,rgn);

    subplot('Position',positions(rgn,:))

    for sid = size(meandata,2):-1:1 % loop through simulations/measurement
        tmean = meandata(:,sid);
        tprc25 = prc25(:,sid);
        tprc75 = prc75(:,sid);
        % to deal with missing data (part 1: using mean for missing data)
        nanind =  find(isnan(tprc25));
        if ~isempty(nanind)
            mnvalid = 1:12;
            % valid only matrices
            mnvalid(nanind)=[]; % index for month with valid data
            tprc25valid = tprc25; tprc25valid(nanind)=[];
            tprc75valid = tprc75; tprc75valid(nanind)=[];
            meanvalid = tmean; meanvalid(nanind)=[];

            % expand to 12 month
            tmean = interp1(mnvalid,meanvalid,1:12,'nearest','extrap')';
            tprc25 = interp1(mnvalid,tprc25valid,1:12,'nearest','extrap')';
            tprc75 = interp1(mnvalid,tprc75valid,1:12,'nearest','extrap')';
        end
        

        % ================= Plot Begins ===========================

        % adding shadow
        xconf = [1:12 12:-1:1] ;
        yconfsim = [tprc75' tprc25(end:-1:1)'];
        p = fill(xconf,yconfsim(1,:),Color(sid,:),'FaceAlpha',fal);
        p.EdgeColor = 'none';
        hold on

        % adding the center line
        if sid ==1
            legends(sid) = plot(1:12,tmean(1:12),lines{sid},'Color',Color(sid,:),'LineWidth',lw,'DisplayName','Measurement');
        else
            legends(sid) = plot(1:12, tmean(1:12), lines{sid}, 'Color', Color(sid, :), 'LineWidth', lw, 'DisplayName', target_sim{sid - 1});
        end
        hold on

        % to deal with missing data (part 2: using white mask)
        mu = max(tprc75); % top of the mask
        md = min(tprc25); % bottom of the mask
        if ~isempty(nanind)
            for ii = 1:length(nanind)
                sind = nanind(ii)-0.5;
                eind = nanind(ii)+0.5;
                if sind<0
                    sind = 0;
                end
                if eind>12
                    eind =12;
                end

                xconf =    [       sind:0.5:eind        eind:-0.5:sind] ;
                yconfsim = [repmat(mu,[1,length(xconf)/2])  repmat(md,[1,length(xconf)/2]) ];

                p = fill(xconf,yconfsim(1,:),[1,1,1]); % white mask
                p.EdgeColor = 'none';
            end
            clear nanind
        end
        hold on

        % add y axis label
        if rgn == 9 || rgn == 4 || rgn == 5 % strange location since plots are placed in a strange order (see 'positions')
            switch spec
                case 'PM2.5'
                    ylabel(sprintf('%s (%s)','PM_{2.5}','\mug/m^{3}'));
                case 'SO4'
                    ylabel(sprintf('%s (%s)','Sulfate','\mug/m^{3}'));
            end
            hold on
        end

    end

    title(region_name{existrg(rgn)})
    xticks(1:12)
    xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
    xlim([1 12])
    xtickangle(45)
    set(gca ,'FontName', 'Arial', 'FontSize',fz);
    if rgn == 1
        legend(legends,'location','northeast')
    end

    % adding site label
    % sitelabel = '';
    % for nn = 1:length(tsites)
    %     sitelabel =  sprintf('%s\n%s',sitelabel,tsites{nn});
    % end
    text(0.05,1,sitelabel,'Units','normalized','VerticalAlignment','top','fontSize',fz+2)


end
% save figure
fname = sprintf('%s/TimeSeries_Region_%s.png',SaveDir,spec);
saveas(gcf,fname)
fprintf('%s saved\n',fname)

close

end