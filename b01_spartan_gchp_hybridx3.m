% July-18-2024
% Haihui Zhu
% 
% This script makes the 3 essential figures comparing SPARTAN data and
% other data source.
clear 
close all

% RootDir = '/Volumes/rvmartin/Active/'; % read from compute1
RootDir = '/storage1/fs1/rvmartin/Active/'; % read from compute1
addpath(sprintf('%s/haihuizhu/1.code',RootDir))
addpath(sprintf('%s/haihuizhu/4.SPARTAN_SO4/functions',RootDir))
SaveDir = sfmkdir(sprintf('%s/haihuizhu/4.SPARTAN_SO4/06.spartan_gchp/',RootDir)); % save to compute1

masterfname = sprintf('pm25_so4_202407'); % file name to save

S_fig_statis = 0; % run to get hybridx3 figure and the statistics
S_bar = 1; % run to get bar including multiple sims
    % if s_bar is true:
    target_sim = {'ceds_2021', 'edgar_2018', 'htap_2018'};

ObsYear = 2019:2023; % years that spartan data can be included (WashU)
simyear = 2018;
SimName = 'edgar' ; 
SimDir = get_sim_dir(SimName, simyear,RootDir);
Species = {'SO4',  'NO3',  'NH4',   'OM', 'BC',  'dust', 'SS'}; 

% calculate monthly mean: both SPT and AERONET exist at the same day for at
% least 3 days in a month,ignore hours of samples
thres_Mn = 5; % at least X days available for a month; 
thres_Hr = 24*6; % if there are filters sampled via 1-day period instead of 9-day, 6 filters might suffice; 
thres_Ann = 9; % at least X months available for a year


% get region id 
rgidfname = sprintf('%s/region_id.mat',SaveDir);
if ~exist(rgidfname,'file')
    load(sprintf('%s/SPT_%s.mat',SaveDir,masterfname)) % load the corresponding spartan data

    % loading regional mask
    load('./ceds_scale_2021to2018/mask_fao_ceds_005.mat')
    region_id = interp2(xlon(1,:),xlat(:,1),mask_region,longitudes,latitudes,'nearest');
    region_id(7) = 3;
    region_id(28) = 8;
    if sum(find(isnan(region_id)))>0
        error('Site not assigned to a region:%s',Site_cities(isnan(region_id)))
    end
    % trim region_name according to spartan data availability
    unirg = unique(region_id);
    region_id2 = nan.*region_id;
    region_name2 = cell(length(unirg),1);
    for ii = 1:length(unirg)
        region_id2(region_id == unirg(ii)) = ii;
        region_name2{ii} = region_name{unirg(ii)};
    end
    region_name = region_name2;
    region_id = region_id2;
    save(rgidfname,'region_id','Site_cities','region_name')
else
    load(rgidfname)
end






%% Make Site Mean Contour + Bar + Monthly Scatter:  ETA, PM, and AOD
if S_fig_statis == 1
% ==== Map background data ====
fname = sprintf('%s/GCHP_SO2_SO4_BC_PM25_%s_%d_annual.nc',SimDir,SimName,simyear);
% fname = sprintf('%s/GCHP_SO2_SO4_PM25_%s_%d_annual.nc',SimDir,SimName,simyear);

simlat=ncread(fname,'latitude');
simlon=ncread(fname,'longitude');

pm25map = ncread(fname,'pm25');
so4map = ncread(fname,'so4');

% ==== load measurement + coincident sim ====
load(sprintf('%s/SPT_%s.mat',SaveDir,masterfname)) % load the corresponding spartan data
load(sprintf('%s/SIM_%s_%s-%d.mat',SaveDir,masterfname,SimName,simyear),'so4_sim','pm25_sim')

D2_Titles = Species;

% ==== before filter data for dates/species, repeat SIM for other years ===
rowsimyear = find(ismember(D1_Dates(:, 1), simyear)); % select one year of data to plot
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
% remove unneeded years and species
cols = [findcol('PM2.5', D2_Titles) 
        findcol('SO4', D2_Titles) 
        findcol('SampledHour', D2_Titles) ];
rows = find(ismember(D1_Dates(:,1),ObsYear)); % select one year of data to plotcl
D1_Dates = D1_Dates(rows,:);
sptall = squeeze( TOT(rows,cols,:) );
pm25_sim = squeeze( pm25_sim(rows,:) );
so4_sim = squeeze( so4_sim(rows,:) );

% initiate Mn data
sptmn = NaN.*zeros(12,4,length(Site_cities)); % Mn mean PM2.5, SO4 at each site
stdsmn = NaN.*zeros(12,4,length(Site_cities)); % standard deviation; not sure if needed
      
for st = 1:length(Site_cities)

    if sum(~isnan(sptall(:,1,st))) == 0
        fprintf('No SPT measurement for %s\n',Site_cities{st})
        continue
    end

    for Mn = 1:12
        Ind = find(D1_Dates(:,2) == Mn); % ignore interannual differences
        tmnPM = sptall(Ind,1,st);
        tmnSO4 = sptall(Ind,2,st);
        tmnPMsim = pm25_sim(Ind,st);
        tmnSO4sim = so4_sim(Ind,st);
        
        ind = find(~isnan(tmnPM));
        if length(ind)>thres_Mn
            sptmn(Mn,1,st) = mean(tmnPM(ind),'omitnan');
            stdsmn(Mn,1,st) = std(tmnPM(ind),'omitnan');

            sptmn(Mn,3,st) = mean(tmnPMsim(ind),'omitnan');
            stdsmn(Mn,3,st) = std(tmnPMsim,'omitnan');
        end

        ind = find(~isnan(tmnSO4));
        if length(ind)>thres_Mn
            sptmn(Mn,2,st) = mean(tmnSO4(ind),'omitnan');
            stdsmn(Mn,2,st) = std(tmnSO4(ind),'omitnan');

            sptmn(Mn,4,st) = mean(tmnSO4sim(ind),'omitnan');
            stdsmn(Mn,3,st) = std(tmnSO4sim,'omitnan');
        end

    end
end


% ==== PM25 ====
Legends = {'SPARTAN',sprintf('gchp-%s',SimName)}; 
spec = 'PM25';
% Scatter 
Mn_obs =  squeeze(sptmn(:,1,:)) ; % PM25_SPT 
Mn_sim =  squeeze(sptmn(:,3,:)) ; % PM25_sim 
% Bar
Mn_obs_num =  sum(~isnan(Mn_obs),1); % num of Mn at each site; 1xlength(site_cities)
Mn_obs_site = mean(Mn_obs,1,'omitnan'); % site mean
Mn_sim_site = mean(Mn_sim,1,'omitnan');
Mn_obs_std = std(Mn_obs,1,'omitnan'); % site std
Mn_sim_std = std(Mn_sim,1,'omitnan');
% making figure
[slope,b,r2,nrmsd,nmb,n,obsm]= mkfigure_hybridx3(pm25map, simlat,simlon, ... % map 
                  Mn_obs_site,Mn_obs_std, Mn_sim_site,Mn_sim_std, Mn_obs_num,... % contour and bar
                  Mn_obs, Mn_sim, Site_cities,latitudes,longitudes,spec,Legends, region_id,region_name); % scatter and labels
% save figure
sfname = sprintf('%s/Hybridx3_%s_%s_%d.png',SaveDir,spec,SimName,simyear);
saveas(gcf, sfname)
fprintf('saving %s\n',sfname)
close 
% save statistics
sfname = sprintf('%s/Statis_%s_%s_%d.mat',SaveDir,spec,SimName,simyear);
save(sfname,'slope','b','r2','nrmsd','nmb','n','obsm','region_name')
fprintf('saving %s\n',sfname)


% ==== SO4 ====
Legends = {'SPARTAN',sprintf('gchp-%s',SimName)}; 
spec = 'SO4';
% Scatter 
Mn_obs =  squeeze(sptmn(:,2,:)) ; % SO4_SPT 
Mn_sim =  squeeze(sptmn(:,4,:)) ; % SO4_sim 
% Bar
Mn_obs_num =  sum(~isnan(Mn_obs),1); % num of Mn at each site; 1xlength(site_cities)
Mn_obs_site = mean(Mn_obs,1,'omitnan'); % site mean
Mn_sim_site = mean(Mn_sim,1,'omitnan');
Mn_obs_std = std(Mn_obs,1,'omitnan'); % site std
Mn_sim_std = std(Mn_sim,1,'omitnan');
% making figure
[slope,b,r2,nrmsd,nmb,n,obsm] = mkfigure_hybridx3(so4map, simlat,simlon, ... % map 
                  Mn_obs_site,Mn_obs_std, Mn_sim_site,Mn_sim_std, Mn_obs_num,... % contour and bar
                  Mn_obs, Mn_sim, Site_cities,latitudes,longitudes,spec,Legends, region_id,region_name); % scatter and labels
% save figure
sfname = sprintf('%s/Hybridx3_%s_%s_%d.png',SaveDir,spec,SimName,simyear);
saveas(gcf, sfname)
fprintf('saving %s\n',sfname)
close 
% save statistics
sfname = sprintf('%s/Statis_%s_%s_%d.mat',SaveDir,spec,SimName,simyear);
save(sfname,'slope','b','r2','nrmsd','nmb','n','obsm','region_name')
fprintf('saving %s\n',sfname)


diary off

end

%% making bar
if S_bar == 1

    note = 'annual'; % any note on the figure explaining the data

    spec = 'PM25';
    load_data_make_bar(SaveDir,spec,target_sim,note);
    sfname = sprintf('%s/Bars_%s.png',SaveDir,spec); % change figure name
    saveas(gcf,sfname)
    fprintf('saving %s\n',sfname)
    close


    spec = 'SO4';
    load_data_make_bar(SaveDir,spec,target_sim,note);
    sfname = sprintf('%s/Bars_%s.png',SaveDir,spec);
    saveas(gcf,sfname)
    fprintf('saving %s\n',sfname)
    close

end

%% Function
function fig = load_data_make_bar(SaveDir,spec,target_sim,note)
    % four subplots
    r2_all = [];
    slo_all = [];
    nrmsd_all = [];
    nmb_all = [];
    % two used in all plots
    obsm_all = [];
    n_all = [];
    
    for sid = 1:length(target_sim)
        sfname = sprintf('%s/Statis_%s_%s.mat',SaveDir,spec,target_sim{sid});
        load(sfname)
        % save(sfname,'slope','b','r2','nrmsd','nmb','n','obsm','region_name')

        r2_all(:,sid) = r2;
        nmb_all(:,sid) = nmb;
        slo_all(:,sid) = slope;
        nrmsd_all(:,sid) = nrmsd;

        n_all(:,sid) = n;
        obsm_all(:,sid) = obsm;
    end
    region_name{end+1} = 'All'; % statistics contain one more column as 'all'

    % making stat# of plots in a figure. each plot compares serveral sims
    % on serveral regions for one statistics
    PaperSize = [0 0 800 2400]*1.2;
    fig = figure('Position',PaperSize);
 
    subplot(4,1,  1)
    stats_lable = 'r2';
    yrng = [-0.18 1.0];
    subbar(r2_all,n_all,obsm_all,target_sim,stats_lable,region_name,yrng,note)

    subplot(4,1,  2)
    stats_lable = 'nmb';
    yrng = [-0.6 2];
    subbar(nmb_all,n_all,obsm_all,target_sim,stats_lable,region_name,yrng,note)

    subplot(4,1,  3)
    stats_lable = 'slope';
    yrng = [0 2.6];
    subbar(slo_all,n_all,obsm_all,target_sim,stats_lable,region_name,yrng,note)

    subplot(4,1,  4)
    stats_lable = 'nrmsd';
    yrng = [0 2.6];
    subbar(nrmsd_all,n_all,obsm_all,target_sim,stats_lable,region_name,yrng,note)
    
end

function subbar(statistics,N,ObsM,target_sims,stats_lable,region_name,yrng,note)
    N = N(:,1); % rg x 1
    ObsM = ObsM(:,1);

    % move 'all' to the beginning
    idx = [size(N,1), 1:size(N,1)-1];
    region_name = region_name(idx); % rg x 1
    statistics = statistics(idx,:); % rg x sims
    ObsM = ObsM(idx,:); % rg x sims
    N = N(idx,:); % rg x sims

    FaceColor = {[0.6350 0.0780 0.1840],[0.9290 0.6940 0.1250],[0 0.4470 0.7410],...
                [0.4660 0.6740 0.1880],[0 0 0],[1 1 0],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330]};
    X = categorical(region_name);
    X = reordercats(X,region_name);
    
    for rg = 1:numel(region_name)
        s = bar(X(rg),statistics(rg,:),1);
     
        % fomatting and style
        for ssi = 1:size(statistics,2) % statistics = regions x sims
            s(ssi).FaceColor = FaceColor{ssi};
            s(ssi).LineWidth = 1.5;
        end


        if ismember(stats_lable,'r2') % only label N and M for R2
            if size(statistics,2)>1
                xtips1 = rg; % center is the index of spec
            else
                xtips1 = s(1).XEndPoints;
            end

            numstr = sprintf('N=%d', N(rg));
            if ~isempty(ObsM)
                numstr = sprintf('%s\nM=%.2f',numstr,ObsM(rg));
            end
            labels1 = string(numstr);

            if statistics(rg)>0
                ytips1 = -0.05;
                text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
                    'VerticalAlignment','top','fontsize',12)
            else
                ytips1 = 0.05;
                text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
                    'VerticalAlignment','bottom','fontsize',12)
            end

            notes = sprintf('N = number of observations\nM = mean of the observations');
            text(3.2,0.95,notes, ...
                'horizontalalignment','left','verticalalignment','top','color','k','fontsize',18);
%             text(3,0.95,note, ...
%                 'horizontalalignment','left','verticalalignment','top','color','k','fontsize',18,'FontWeight','bold');
        end

        hold on
    end

    % adding axis label and legend

    % format legend
    for lg = 1:length(target_sims)
        tlg = target_sims{lg};
        tlg(tlg=='_') = '-';
        target_sims{lg} = tlg;
    end

    switch stats_lable
        case 'r2'
             ylabel('Coefficient of determination');
             legend(target_sims,'Location','northwest') % only give r2 and nrmsd legends
             xticklabels({}) % hide x tick label for manuscript
        case 'nmb'
            ylabel('Normalized Mean Bias')
        case 'slope'
            ylabel('Slope')
            y = ones;
            plot(X,y,':k')
            xticklabels({})% hide x tick label for manuscript
        case 'nrmsd'
            ylabel('Normalized Root Mean Square Deviation')
            legend(target_sims,'Location','northeast') % only give r2 and nrmsd legends
    end

    ylim(yrng)
    set(gca ,'FontName', 'Arial', 'FontSize',18);



end


function [Ylabels, units, maprange,ylimbar,rng]  = getRngNLabels(specie)

switch specie
    case 'PM25'
        Ylabels='PM_{2.5}';
        units = '\mug/m^3';
        maprange = [0 80] ;
        ylimbar = [0 80];
        rng = [0 200]; % scatter range

    case 'AOD'
        Ylabels='AOD';
        units = 'unitless';
        maprange = [0 0.6] ;
        ylimbar = [0 1];
        rng = [0 2.4]; % scatter range

    case 'ETA'
        Ylabels='\eta';
        units = '\mug/m^3';
        maprange = [0 200] ;
        ylimbar = [0 200];
        rng = [0 500]; % scatter range

    case 'SO4'
        Ylabels='Sulfate';
        units = '\mug/m^3';
        maprange = [0 10] ;
        ylimbar = [0 15];
        rng = [0 28]; % scatter range
    case 'NO3'
        Ylabels='Nitrate';
        units = '\mug/m^3';
        maprange = [0 8] ;
        ylimbar = [0 20];
        rng = [0 30]; % scatter range
    case 'NH4'
        Ylabels='Ammonium';
        units = '\mug/m^3';
        maprange = [0 5] ;
        ylimbar = [0 8];
        rng = [0 20]; % scatter range
    case 'BC'
        Ylabels='Black Carbon';
        units = '\mug/m^3';
        maprange = [0 5] ;
        ylimbar = [0 8];
        rng = [0 12]; % scatter range
    case 'OM'
        Ylabels='Organics';
        units = '\mug/m^3';
        maprange = [0 20] ;
        ylimbar = [0 40];
        rng = [0 60]; % scatter range

    case 'dust'
        Ylabels='Dust';
        units = '\mug/m^3';
        maprange = [0 20] ;
        ylimbar = [0 50];
        rng = [0 80]; % scatter range

    case 'SS'
        Ylabels='Sea Salt';
        units = '\mug/m^3';
        maprange = [0 4] ;
        ylimbar = [0 2];
        rng = [0 4]; % scatter range
end

end

function [m,b,r2,nrmsd,nmb,n,obsm] = mkfigure_hybridx3(mapdata, simlat,simlon, bar_obs_site,bar_obs_std, bar_sim_site,bar_sim_std, bar_obs_num,...
                            scatter_obs, scatter_sim, Site_cities,latitudes,longitudes,specie,Legends, region_id,region_name)

[Ylabels, units, maprange,ylimbar,rng]  = getRngNLabels(specie);
% Ylabels = species print name; rng = range for scatter 

PaperSize = [0 0 1000 600]*1.2;

fontsize = 12;
Styles = {'o','>','h','s','^','d','<','p','v','*'};
% Colors = {'#D95319','#4DBEEE','#77AC30','#FFFF00','#0072BD'}; %,
Colors = {[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.9290 0.6940 0.1250],...
            [0.4660 0.6740 0.1880],[1 1 0],[0.4940 0.1840 0.5560],...
            [0.3010 0.7450 0.9330],[ 0.4667 0.6745 0.1882],[ 0.8510 0.3255 0.0980],[0,0,0]};
latmin = -45;
latmax = 60;
lonmin = -155;
lonmax = 165;
fw = 0.65;
fh = fw * ((latmax-latmin)/180) / ((lonmax-lonmin)/360) ; % inner figure height 
Positions = [0.1   0.45  fw   fh ;
             0.1   0.15  0.4  0.3;
             0.5   0.15  0.5  0.3]; % change the width of scatter 

for i = 1:length(latitudes)
    if isnan(bar_obs_site(i)) % no any value available
        region_id(i) = NaN; % so that can be deleted later 
    end
end

% ==== Sort data by region & delete empty site ==== 
[region_id, I] = sort(region_id) ;
bar_obs_site = sortNdelete(bar_obs_site,I,region_id);
bar_obs_std = sortNdelete(bar_obs_std,I,region_id);
bar_sim_site  = sortNdelete(bar_sim_site,I,region_id);
bar_sim_std =  sortNdelete(bar_sim_std,I,region_id);
bar_obs_num = sortNdelete(bar_obs_num,I,region_id);
latitudes = sortNdelete(latitudes,I,region_id);
longitudes = sortNdelete(longitudes,I,region_id);
Site_cities = sortNdelete(Site_cities,I,region_id);

scatter_obs = scatter_obs(:,I);
scatter_sim = scatter_sim(:,I);
scatter_obs(:,isnan(region_id)) = [];
scatter_sim(:,isnan(region_id)) = [];

region_id(isnan(region_id)) = [];

figure('Position',PaperSize);
% ==== make background map ==== 
subplot('Position',Positions(1,:))
h1 = worldmap([latmin  latmax],[lonmin lonmax]);
setm(gca,'Grid','off','MapProjection','miller','parallellabel','off','meridianlabel','off')

% mask out ocean 
load('LandMask-0.1.mat') % wMASK_land
mapdata = interp2(simlon',simlat,mapdata,Lon,Lat');
surfm(Lat, Lon, mapdata.*wMASK_land);
load coastlines
plotm(coastlat,coastlon)
hold on

% colormap jet 
cm = flipud(cbrewer('div','RdYlBu',100,'spline'));
cm(cm>1) = 1;
cm(cm<0) = 0;
colormap(cm);
set(gca,'clim',maprange);

% ==== make contour ==== 
s1 = scatterm(latitudes, longitudes, 150, bar_sim_site,'filled','MarkerEdgeColor','k');
s2 = scatterm(latitudes, longitudes, 60, bar_obs_site,'filled','MarkerEdgeColor','k');
set(h1,'clim', maprange)

cb1=colorbar('vertical', 'fontsize',fontsize, 'fontweight', 'bold');
set(get(cb1,'YLabel'),  'string', sprintf('%s %s',Ylabels,units),...
    'fontsize',14,'fontweight','bold','FontName','Helvetica');

% ==== make bar ====
subplot('Position',Positions(2,:))
bardata = [bar_obs_site' bar_sim_site'];
X = categorical(Site_cities);
X = reordercats(X,Site_cities);
for st = 1:numel(Site_cities)
    s = bar(X(st),bardata(st,:),1);
    s(1).FaceColor = Colors{1};
    s(2).FaceColor = Colors{2};
    s(2).LineStyle = '-.';
    alpha(s(2),0.6)

    xtips1 = s(2).XEndPoints;
    ytips1 = max([s.YEndPoints]);
    numstr = sprintf('N=%d', bar_obs_num(st));
    text(xtips1-0.1,ytips1,numstr,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','fontsize',fontsize-3.5)
    hold on
end
set(gca, 'LineWidth', 1) % box line width
ylim(ylimbar)
ylabel(sprintf('%s (%s)',Ylabels,units),'FontWeight','bold');
set(gca ,'FontName', 'Arial', 'FontSize',fontsize);
legend(Legends,'Location','northwest','FontSize',fontsize+2)

% ==== make scatter ====
unirg = unique(region_id);
subplot('Position',Positions(3,:))
for rg = 1:length(unirg)
    X = reshape(scatter_obs(:,region_id == unirg(rg)),[],1);
    Y = reshape(scatter_sim(:,region_id == unirg(rg)),[],1);
    C = rg.*ones(size(X));
    s = scatter(X, Y, 24,C,Styles{rg},'filled');
    [~,m(rg),b(rg),r2(rg),nrmsd(rg),nmb(rg),obsm(rg),n(rg)] = get_all_Statics(X,Y);
    hold on
end
colormap(gca ,flip(lines(numel(Colors)))) % change color

% Grid and framing
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , 'k', ...
    'YColor'      , 'k', ...
    'Xlim'        , (rng), ...
    'Ylim'        , (rng), ...
    'YTick'       , 0:(rng(2)/4):rng(2), ...
    'XTick'       , 0:(rng(2)/4):rng(2), ...
    'LineWidth'   , 1         );
pbaspect([1 1 1])

% Statistics and labels
X = reshape(scatter_obs,[],1);
Y = reshape(scatter_sim,[],1);
[ta1,m(rg+1),b(rg+1),r2(rg+1),nrmsd(rg+1),nmb(rg+1),obsm(rg+1),n(rg+1)] = get_all_Statics(X,Y);
% label on the left top
xl = 0:0.01*rng(2):rng(2);
yl = m(rg+1)*xl + b(rg+1);

text(rng(2)-0.95*(rng(2)-rng(1)),rng(1)+0.96*(rng(2)-rng(1)),ta1, ...
    'horizontalalignment','left','verticalalignment','top','color','k','fontsize',fontsize);
plot(xl,yl,'-r','Linewidth',1)

% calculate poplulation weighted mean
if ~exist('tpop','Var')
    load('Population-0.05.mat')
    tpop = interp2(tLAT',tLON,npop,latitudes,longitudes);
    clear pop popm LAT LON
end

popWmobs = sum(tpop.*bar_obs_site','omitnan')/(sum(tpop));
popWmsim = sum(tpop.*bar_sim_site','omitnan')/(sum(tpop));
notes = sprintf('PWM_{GM} = %4.1f\nPWM_{SIM} = %4.1f',popWmobs,popWmsim);
% label on the right bottom
text(0.95*(rng(2)-rng(1)),rng(1)+0.04*(rng(2)-rng(1)),notes, ...
    'horizontalalignment','Right','verticalalignment','bottom','color','k','fontsize',fontsize);

% Tick font size
set(gca,'FontName','Arial','FontSize', fontsize-2);

% X and Y label
hYLabel = ylabel(sprintf('%s_{SIM}(%s)',Ylabels,units));
hXLabel = xlabel(sprintf('%s_{GM}(%s)',Ylabels,units));
set([hXLabel, hYLabel], ...
    'FontName'   , 'Arial',...
    'FontSize'   , fontsize,...
    'fontweight', 'bold');

% 1:1 line
y3=xl;
plot(xl,y3,'-k')
% 2:1 line
y3=0.5*xl;
plot(xl,y3,'--k')
% 1:2 line
y3=2*xl;
plot(xl,y3,'--k')
legend_label = region_name;
legend_label(end+1:end+3) = {'Fit Line','1:1 Line','2:1 Line'};
legend(legend_label,'Location','eastoutside')

pause(3)
end

function out = sortNdelete(in,I,rgID)
    out = in(I);
    out(isnan(rgID)) = [];
end