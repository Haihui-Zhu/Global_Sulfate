% This script analyzes SO2 emission from 3 emission inventories, particularly 
% emission from anthropogentic sources. 
% It creates plots for 
% 1. comparison of the annual PWM SO2 (all sectors) among the inventories
% 2. seasonality of global PWM SO2 
% 3. seasonality of PWM SO2 in the Middle East, also the SPARTAN sites

clear
close all
addpath('/storage1/fs1/rvmartin/Active/haihuizhu/1.code')
addpath('/storage1/fs1/rvmartin2/Active/haihuizhu/10.SpatialETA')
RootDir = '/storage1/fs1/rvmartin/Active';

SaveDir = sfmkdir(sprintf('./01.Emissions'));

S_CEDS_2023 = 1;
S_CEDS_2021 = 0;
S_EDGARv61 = 0;
S_HTAPv3 = 0;

S_maps = 0;
    S_global_jan2jul_ratio = 1;
S_seasonal_plots = 0;

SimYr = 2018;

% load SPARTAN site information 
site_details = readtable(sprintf('%s/SPARTAN-shared/Site_Sampling/Site_details.xlsx',RootDir),'VariableNamingRule' , 'preserve');
Site_cities = table2array(site_details(:,3));
latitudes = table2array(site_details(:,5));
longitudes = table2array(site_details(:,6));

tunit = 'kg m-2 s-1';

%% CEDS emission 2023
if S_CEDS_2023 == 1

    EmissionName = 'CEDS-2023';
    CEDSDir = './CEDS-2023';
    Sectors = {'agr','ene','ind','rco',   'shp','slv','tra','wst'};
    
    Spec = {'SO2','NO','N2O'}; % NOx + SO2
    
    for sp = 1%:length(Spec)

        Infname = sprintf('%s/%d/CEDS_%s_0.1x0.1_%d.nc',CEDSDir,SimYr,Spec{sp},SimYr);
        for sc = 1:length(Sectors)
            if sc == 1
                emi_so2 = ncread(Infname,sprintf('%s_%s',Spec{sp},Sectors{sc})); % lat lon time
                emi_so2(isnan(emi_so2)) = 0;
            else
                tspec0 = ncread(Infname,sprintf('%s_%s',Spec{sp},Sectors{sc})); % lat lon time
                tspec0(isnan(tspec0)) = 0;
                emi_so2 = emi_so2 + tspec0; clear tspec0
            end
        end

        emi_so2 = permute(emi_so2,[2 1 3]);
        lat = ncread(Infname,'lat');
        lon = ncread(Infname,'lon');

    end
    sfname = sprintf('%s/%s_seasonal_allsectors.mat',SaveDir,EmissionName);
    save(sfname,'emi_so2','tunit','lat','lon')
    fprintf('%s saved\n',sfname)
    clear emi_so2
end

%% CEDS-2021 
if S_CEDS_2021 == 1

    EmissionName = 'CEDS-2021';
    CEDSDir = './CEDS-2021';
    Sectors = {'agr','ene','ind','rco',   'shp','slv','tra','wst'};
    
    Spec = {'SO2','NO','N2O'}; % NOx + SO2
    
    for sp = 1%:length(Spec)

        Infname = sprintf('%s/%d/%s-em-anthro_CMIP_CEDS_%d.nc',CEDSDir,SimYr,Spec{sp},SimYr);
        for sc = 1:length(Sectors)
            if sc == 1
                emi_so2 = ncread(Infname,sprintf('%s_%s',Spec{sp},Sectors{sc})); % lat lon time
                emi_so2(isnan(emi_so2)) = 0;
            else
                tspec0 = ncread(Infname,sprintf('%s_%s',Spec{sp},Sectors{sc})); % lat lon time
                tspec0(isnan(tspec0)) = 0;
                emi_so2 = emi_so2 + tspec0; clear tspec0
            end
        end

        emi_so2 = permute(emi_so2,[2 1 3]);

    end
    % regrid to 0.1
    Lat = ncread(Infname,'lat');
    Lon = ncread(Infname,'lon');
    fname = sprintf('%s/%s_seasonal_allsectors.mat',SaveDir,'CEDS-2023');
    load(fname,'lat','lon');
    emi_so2 = interp3(Lon,Lat',1:12, emi_so2, lon, lat',1:12,'nearest');

    sfname = sprintf('%s/%s_seasonal_allsectors.mat',SaveDir,EmissionName);
    save(sfname,'emi_so2','tunit','lat','lon')
    fprintf('%s saved\n',sfname)
    clear emi_so2

end
%% EDGARv6.1
if S_EDGARv61 == 1
    
    EmissionName = 'EDGARv6.1';
    files = dir('./EDGARv6.1/*_nc');

    for ff = 1:length(files)
        filename = files(ff).name;
        if filename(1) ~= '.'

            if contains(filename,num2str(SimYr)) % sectors with monthly data

                for Mn = 1:12
                    file2 = dir(sprintf('./EDGARv6.1/%s/*_%d_*.nc',filename,Mn));
                    for mm = 1:length(file2)
                        filename2 = file2(mm).name;
                        if filename2(1) ~= '.'

                            file2read = sprintf('./EDGARv6.1/%s/%s',filename,filename2);
                            fprintf('Reading %s\n',file2read)

                            % initialize
                            if ~exist('emi_so2','var')
                                lat = ncread(file2read,'lat');
                                lon = ncread(file2read,'lon');
                                emi_so2 = zeros(length(lat),length(lon),12);
                            end

                            temi = ncread(file2read,'emi_so2')';
                            temi(isnan(temi)) = 0;
                            emi_so2(:,:,Mn) = emi_so2(:,:,Mn) + temi;

                        end
                    end
                end

            else % sectors with only annual mean
                file2 = dir(sprintf('./EDGARv6.1/%s/*%s*.nc',filename,num2str(SimYr)));
                for mm = 1:length(file2)
                    filename2 = file2(mm).name;
                    if filename2(1) ~= '.'
                        file2read = sprintf('./EDGARv6.1/%s/%s',filename,filename2);
                        fprintf('Reading %s\n',file2read)

                        % initialize
                        if ~exist('emi_so2','var')
                            lat = ncread(file2read,'lat');
                            lon = ncread(file2read,'lon');
                            emi_so2 = zeros(length(lat),length(lon),12);
                        end
                        temi = ncread(file2read,'emi_so2')';
                        temi(isnan(temi)) = 0;
                        for Mn = 1:12
                            emi_so2(:,:,Mn) = emi_so2(:,:,Mn) + temi; % [May 13: remove '/12' since this is annul mean emission flux]
                        end
                    end
                end
            end
        end
    end

    lon(lon>180) = lon(lon>180)-360;
    [lon,I ] = sort(lon);
    emi_so2 = emi_so2(:,I,:);
    disp(size(emi_so2))

    sfname = sprintf('%s/%s_seasonal_allsectors.mat',SaveDir,EmissionName);
    save(sfname,'emi_so2','tunit','lat','lon')
    fprintf('%s saved\n',sfname)
    clear emi_so2
  

end


%% HTAPv3
if S_HTAPv3 == 1
    EmissionName = 'HTAPv3';
    Sectors = {'AGR','ENE','IND','RCO',   'SHP','SLV','TRA','WST'};

    file2read = sprintf('./HTAPv3/2018/HTAPv3_SO2_0.1x0.1_%s.nc',num2str(SimYr));
    
    fprintf('Reading %s\n',file2read)
    for sc = 1:length(Sectors)
        if sc == 1
            emi_so2 = ncread(file2read,sprintf('SO2_%s',Sectors{sc})); % lat lon time
            emi_so2(isnan(emi_so2)) = 0;
        else
            tspec0 = ncread(file2read,sprintf('SO2_%s',Sectors{sc})); % lat lon time
            tspec0(isnan(tspec0)) = 0;
            emi_so2 = emi_so2 + tspec0; clear tspec0
        end
    end

    emi_so2 = permute(emi_so2,[2 1 3]);
    lat = ncread(file2read,'lat');
    lon = ncread(file2read,'lon');

    sfname = sprintf('%s/%s_seasonal_allsectors.mat',SaveDir,EmissionName);
    save(sfname,'emi_so2','tunit','lat','lon')
    fprintf('%s saved\n',sfname)
end

%% Load data for maps and plots
if (S_maps + S_seasonal_plots) > 0

    fname = sprintf('%s/CEDS-2023_seasonal_allsectors.mat',SaveDir);
    ceds = load(fname);
    fname = sprintf('%s/EDGARv6.1_seasonal_allsectors.mat',SaveDir);
    edgar = load(fname);
    fname = sprintf('%s/HTAPv3_seasonal_allsectors.mat',SaveDir);
    htap = load(fname);

end

%% Making maps
if S_maps == 1
    % annual mean
    cedsann = mean(ceds.emi_so2,3,'omitnan');
    edgarann = mean(edgar.emi_so2,3,'omitnan');
    htapann = mean(htap.emi_so2,3,'omitnan');

    % make maps
    compare_emission_map(htap.lat,htap.lon, cedsann,'CEDS SO_2',edgarann,'EDGAR SO_2',htapann,'HTAP SO_2',tunit)

    % save maps
    sfname = sprintf('%s/Map_compare_ceds_edgar_htap_ann.png',SaveDir);
    saveas(gcf,sfname)
    fprintf('saving %s\n',sfname)

    % jan to jul ratio
    if S_global_jan2jul_ratio == 1

        mapdata1 = ceds.emi_so2(:,:,1)./ceds.emi_so2(:,:,7);
        mapdata2 = edgar.emi_so2(:,:,1)./edgar.emi_so2(:,:,7);
        mapdata3 = htap.emi_so2(:,:,1)./htap.emi_so2(:,:,7);
        % make figure for the ratio
        compare_emission_map(htap.lat,htap.lon, mapdata1,'CEDS ratio',mapdata2,'EDGAR ratio',mapdata3,'HTAP ratio',tunit)
        sfname = sprintf('%s/Map_compare_ceds_edgar_htap_jan2july.png',SaveDir);
        saveas(gcf,sfname)
        fprintf('saving %s\n',sfname)

        % screen out non-1 area and see what remains
        mapdata1(mapdata1==1) = nan;
        mapdata2(mapdata2==1) = nan;
        mapdata3(mapdata3==1) = nan;
        compare_emission_map(htap.lat,htap.lon, mapdata1,'CEDS ratio',mapdata2,'EDGAR ratio',mapdata3,'HTAP ratio',tunit)
        sfname = sprintf('%s/Map_compare_ceds_edgar_htap_jan2july_no_ones.png',SaveDir);
        saveas(gcf,sfname)
        fprintf('saving %s\n',sfname)

    end

end

%% Making seasonal plots
if S_seasonal_plots == 1

    [lat, lon, regional_mapdata] = making4seasonmaps(lat,lon,mapdata,Spec,Unit,rng);

    % making scatter plot for 12-month variation for the Middle East
    % and for each SPARTAN site in the area
    for Mn = 1:12
        site_em(:,Mn) = interp2(rlon,rlat',regional_emission(:,:,Mn),longitudes,latitudes);
        [region_pwm(Mn),region_std(Mn),~] = Pop_weighted_STD(regional_emission(:,:,Mn),rlat,rlon);
    end
    Site_cities2 = Site_cities;
    Site_cities2(isnan(site_em(:,1))) = [];
    site_em(isnan(site_em(:,1)),:) = [];
    disp(size(site_em))

    % make plot
    make_seasonality_plot(site_em,Site_cities2,region_pwm,region_std,tunit)
    % save the 4 season figure
    sfname = sprintf('%s/Plot_12-Month_%s_%s_%s_%d.png',SaveDir,EmissionName,Spec{sp},Sectors{sc},SimYr);
    saveas(gcf,sfname)
    fprintf('saving %s\n',sfname)
    clear site_em
end





%% Functions
function    make_seasonality_plot(site_em,Site_cities,region_pwm,region_std,tunit)
fz = 16;
figure('Position',[10 10 1200 600])

% first plot
subplot(1,2,1)
for ii = 1:length(Site_cities)
    plot(1:12, site_em(ii,:),'DisplayName',Site_cities{ii})
    hold on
end
plot(1:12, region_pwm,'DisplayName','Middle East Mean')
hold on
% adding error bar
er = errorbar(1:12,region_pwm,region_std,region_std);
er.LineStyle = 'none';
hold on

legend('Location','best')
ylabel(tunit)
xticks(1:12)
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'fontsize',fz)

% second plot - normalized
site_em = site_em./mean(site_em,2);
region_pwm = region_pwm./mean(region_pwm);
region_std = region_std./mean(region_pwm);

subplot(1,2,2)
for ii = 1:length(Site_cities)
    plot(1:12, site_em(ii,:),'DisplayName',Site_cities{ii})
    hold on
end
plot(1:12, region_pwm,'DisplayName','Middle East Mean')
hold on
% adding error bar
er = errorbar(1:12,region_pwm,region_std,region_std);
er.LineStyle = 'none';
hold on

legend('Location','best')
ylabel('Normalized emission')
ylim([0.5 1.5])
xticks(1:12)
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'fontsize',fz)

end

function  season_ratio_maps(lat,lon, mapdata,Spec)
[latmin,latmax,lonmin,lonmax,~,~,~,~,PaperSize] = Map4_FigPara('global');
fz = 20;
rng = [0.9 1.1];

% load mask file
load('NewLandMask-0.05.mat') % 'LAT','LON'
wMASK = double(MASKp); clear MASKp
wMASK_land = wMASK(:,:,2); clear wMASK
wMASK_land(wMASK_land==0) = NaN;
wMASK_land(abs(wMASK_land)>0) = 1;
Lat = double(LAT); Lon = double(LON);

wMASK_land = interp2(Lon,Lat',wMASK_land,lon ,lat','nearest');

fig = figure('WindowState','maximized','Position',PaperSize);

% making figure
worldmap([latmin latmax],[lonmin lonmax]);
surfm(lat,lon, mapdata.*wMASK_land)
setm(gca,'Grid','off','MapProjection','miller','parallellabel','off','meridianlabel','off')
tightmap

bordersm
load coastlines
plotm(coastlat,coastlon)
hold on

if rng(1) ==0
    cm = cbrewer('div','RdYlBu',100,'spline');
    cm(cm>1) = 1;
    cm(cm<0) = 0;
    colormap(gca,flipud(cm));
elseif rng(1) < 0
    cm = cbrewer('div','RdBu',100,'spline');
    cm(cm>1) = 1;
    cm(cm<0) = 0;
    colormap(gca,flipud(cm));
end
set(gca,'clim',rng); %  colorscale = log or linear


cb2=colorbar('vertical','fontsize',fz,'fontweight', 'bold');
cblb = sprintf('Jan to Jul Ratio');
set(get(cb2,'YLabel'),'string',cblb,'fontsize',fz,'fontweight','bold','FontName','Helvetica');
hold on


end

function  compare_emission_map(lat,lon, mapdata1,note1, mapdata2, note2,mapdata3,note3,tunit)

[latmin,latmax,lonmin,lonmax,Positions,cbvserpos,notepos,PaperSize] = Map6_FigPara('global');

fz = 20;
rng = [0.1 10].*mean(mean(mapdata1,'omitnan'),'omitnan');
disp(rng)

% load mask file
load('NewLandMask-0.05.mat') % 'LAT','LON'
wMASK = double(MASKp); clear MASKp
wMASK_land = wMASK(:,:,2); clear wMASK
wMASK_land(wMASK_land==0) = NaN;
wMASK_land(abs(wMASK_land)>0) = 1;
Lat = double(LAT); Lon = double(LON);

wMASK_land = interp2(Lon,Lat',wMASK_land,lon ,lat','nearest');

fig = figure('WindowState','maximized','Position',PaperSize);

% making figure
ps = 2;
submap_emi(ps, mapdata1,wMASK_land, lat,lon,latmin,latmax,lonmin,lonmax,Positions,cbvserpos,notepos,rng,tunit,note1)
ps = 4;
submap_emi(ps, mapdata2,wMASK_land, lat,lon,latmin,latmax,lonmin,lonmax,Positions,cbvserpos,notepos,rng,tunit,note2)
ps = 6;
submap_emi(ps, mapdata3,wMASK_land, lat,lon,latmin,latmax,lonmin,lonmax,Positions,cbvserpos,notepos,rng,tunit,note3)

end


function submap_emi(ps, mapdata,wMASK_land, lat,lon,latmin,latmax,lonmin,lonmax,Positions,cbvserpos,notepos,rng,tunit,note)
subplot('Position',Positions(ps,:));
worldmap([latmin latmax],[lonmin lonmax]);
surfm(lat,lon, mapdata.*wMASK_land)
% surfm(lat,lon, mapdata)
setm(gca,'Grid','off','MapProjection','miller','parallellabel','off','meridianlabel','off')
tightmap

fz = 14; 

bordersm
load coastlines
plotm(coastlat,coastlon)
hold on

if rng(1) ==0
    cm = cbrewer('div','RdYlBu',100,'spline');
    cm(cm>1) = 1;
    cm(cm<0) = 0;
    colormap(gca,flipud(cm));
elseif rng(1) < 0
    cm = cbrewer('div','RdBu',100,'spline');
    cm(cm>1) = 1;
    cm(cm<0) = 0;
    colormap(gca,flipud(cm));
end
set(gca,'clim',rng); %  colorscale = log or linear

cb2=colorbar('vertical','fontsize',fz,'fontweight', 'bold','Position',cbvserpos(ps,:));
cblb = sprintf('so2 emission %s',tunit);
set(get(cb2,'YLabel'),'string',cblb,'fontsize',fz,'fontweight','bold','FontName','Helvetica');
hold on

% calc population weighted mean
[M,STDW,~] = Pop_weighted_STD(mapdata,lat,lon);
specmeantexts = sprintf('%s\nPWM = %6.2e %s \n%s = %6.2e',note,M, tunit, '\mu',STDW./M);
% print note on figure
text(max(xlim)-notepos(1)*(max(xlim)-min(xlim)), min(ylim)+notepos(2)*(max(ylim)-min(ylim)),...
    specmeantexts ,'fontsize',fz+2,'fontweight','bold','Horiz','right', 'Vert','bottom');

end

function [lat, lon, mapdata] = making4seasonmaps(lat,lon,mapdata,Spec,Unit,rng)
[latmin,latmax,lonmin,lonmax,Positions,~,cbvserpos,notepos,PaperSize] = Map4_FigPara('middle east');
fz = 16;

% output lat lon
latind = find(lat>=latmin & lat <= latmax);
lonind = find(lon>=lonmin & lon <= lonmax);
mapdata = mapdata(latind,lonind,:);
lat = lat(latind);
lon = lon(lonind);

% load mask file
load('NewLandMask-0.05.mat') % 'LAT','LON'
wMASK = double(MASKp); clear MASKp
wMASK_land = wMASK(:,:,2); clear wMASK
wMASK_land(wMASK_land==0) = NaN;
wMASK_land(abs(wMASK_land)>0) = 1;
Lat = double(LAT); Lon = double(LON);

wMASK_land = interp2(Lon,Lat',wMASK_land,lon ,lat','nearest');

Mns = [1 4 7 10];
figlabel = {'Jan','Apr','Jul','Oct'};

for ps = 1:4
    if ps == 1
        fig = figure('WindowState','maximized','Position',PaperSize);
    end
    % making figure
    subplot('Position',Positions(ps,:));
    worldmap([latmin latmax],[lonmin lonmax]);
    surfm(lat,lon, mapdata(:,:,Mns(ps)).*wMASK_land)
    setm(gca,'Grid','off','MapProjection','miller','parallellabel','off','meridianlabel','off')
    tightmap
    
    bordersm
    load coastlines
    plotm(coastlat,coastlon)
    hold on
    
    if rng(1) ==0
        cm = cbrewer('div','RdYlBu',100,'spline');
        cm(cm>1) = 1;
        cm(cm<0) = 0;
        colormap(gca,flipud(cm));
    elseif rng(1) < 0
        cm = cbrewer('div','RdBu',100,'spline');
        cm(cm>1) = 1;
        cm(cm<0) = 0;
        colormap(gca,flipud(cm));
    end
    set(gca,'clim',rng); %  colorscale = log or linear


     if mod(ps,2) == 0
        cb2=colorbar('vertical','Position',cbvserpos(ps,:),'fontsize',fz,'fontweight', 'bold');
        cblb = sprintf('%s (%s)',Spec,Unit);
        set(get(cb2,'YLabel'),'string',cblb,'fontsize',fz,'fontweight','bold','FontName','Helvetica');
        hold on      
     end

     % calc population weighted mean
     [M,STDW,~] = Pop_weighted_STD(mapdata(:,:,Mns(ps)),lat,lon);
     specmeantexts = sprintf('%s\nPWM = %6.2e %s \n%s = %6.2e',figlabel{ps},M, Unit, '\mu',STDW./M);

     % print note on figure
     text(max(xlim)-notepos(1)*(max(xlim)-min(xlim)), min(ylim)+notepos(2)*(max(ylim)-min(ylim)),...
         specmeantexts ,'fontsize',fz+2,'fontweight','bold','Horiz','right', 'Vert','bottom');


end

end