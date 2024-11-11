% July-17-2024
% Haihui Zhu

% Read SO2, SO4, BC, PM25 from GCHP ACAG and ConcAboveSfc collections
% SO2 is sampled at TROPOMIT overpass time.
clear 
close all

% ---- Input & Switches --------------
YEARS = 2018; % year range of interest, currently means year of simulation
SimName  = 'edgar_2018'; % change here
InDir  = '/nobackup/hzhu3/3.edgar_2018/OutputDir/'; % change here
OutDir = '/nobackup/hzhu3/3.edgar_2018/processed/'; % change here

Mns = 1:12;
% ---- End of input ----------------------
daysinmonth = [31 28 31 30 31 30 31 31 30 31 30 31]; 

%% prep
fname = sprintf('%s/GEOSChem.ACAG.%d%.2d%.2d_%.2d30z.nc4',InDir,YEARS(1),01,01,0);
lat=double(ncread(fname,'lats'));
lon=double(ncread(fname,'lons'));
lon(lon>180) = lon(lon>180) - 360;

tLAT = -89.5 :1: 89.5;
tLON = -179.5:1:179.5;
[xLON, xLAT] = meshgrid(tLON,tLAT);

R = 8.314; % m3⋅Pa⋅K−1⋅mol−1
Na = 6.02e23;
SO4MW = 96;
BCMW = 12;

RH_growth_sia = 1.10; % The SIA recommendation at 35% RH is less certain 
% since it depends on the efflorescence RH of the SIA in the aerosol mixture 
% under the variable conditions of the instruments, collection media, and 
% laboratories involved. Given knowledge gaps about the aerosol phase at 
% low RH, the proposed growth factor of 1.1 assumes that half of the 
% particles are aqueous (growth factor of 1.19 for Kappa-Kohler) and the 
% other half are crystalline (growth factor of unity).

% all RH growth factor:
%       35%       50% 
% SIA 1.1000    1.3498
% BC  1.0000    1.0000 % BC does not grow until RH 80
% OC  1.0490    1.0706
% BrC 1.1630    1.2393
% SSa 1.8559    2.4236
% SSc 1.8689    2.4472

%% reading sim output
for Yr = YEARS
    so2y = zeros(length(tLAT),length(tLON));
    pm25y = zeros(length(tLAT),length(tLON));
    so4y = zeros(length(tLAT),length(tLON));
    bcy = zeros(length(tLAT),length(tLON));
    
    for Mn = Mns
        so2m = zeros(length(tLAT),length(tLON));
        pm25m = zeros(length(tLAT),length(tLON));
        so4m = zeros(length(tLAT),length(tLON));
        bcm = zeros(length(tLAT),length(tLON));

        for Dy = 1:daysinmonth(Mn)
            tic
            sfname = sprintf('%s/GCHP_SO2_SO4_BC_PM25_%s_%.2d%.2d.nc',OutDir,SimName,Mn,Dy);
            if exist(sfname,'file') % for re-run without repeating work
                so2 = ncread(sfname,'so2');
                so4 = ncread(sfname,'so4');
                pm25 = ncread(sfname,'pm25');
                bc = ncread(sfname,'bc');
            
            else
                so2d = zeros(length(tLAT),length(tLON));
                pm25d = zeros(length(tLAT),length(tLON));
                so4d = zeros(length(tLAT),length(tLON));
                bcd = zeros(length(tLAT),length(tLON));

                for Hr = 0:23
                    fname = sprintf('%s/GEOSChem.ACAG.%d%.2d%.2d_%.2d30z.nc4',InDir,Yr,Mn,Dy,Hr);
                    fname2 = sprintf('%s/GEOSChem.ConcAboveSfc.%d%.2d%.2d_%.2d00z.nc4',InDir,Yr,Mn,Dy,Hr);
                    so2h = double(ncread(fname,'SpeciesConcVV_SO2')); % mol/mol dry air
                    pm25h = double(ncread(fname,'PM25')); % ug/m3
                    so4h = double(ncread(fname2,'SpeciesConcALT1_SO4')); % secondary species conc at 2m measurement height
                    bch = double(ncread(fname,'SpeciesConcVV_BCPI')) + double(ncread(fname,'SpeciesConcVV_BCPO')); % mol/mol dry air
                    
                    % SO2:
                    % unit conversion [mol/mol dry air] to [molecule per cm2]:
                    % mol/mol dry air => mol/m3 dry air => mol/m2 dry air =>
                    % mol/m2 surface area => moelc/m2 => molec/cm2
                    hBxH = double(ncread(fname,'Met_BXHEIGHT')); % m
                    tT = double(ncread(fname,'Met_T')); % K
                    tP = 100*double(ncread(fname,'Met_PMID')); % Pa
                    vol = 1.*R.*tT./tP; % m3

                    so2h = sum(so2h./vol.*hBxH, 4); % mol/m2
                    so2h = so2h.* Na .* 1e-4; % molec/cm2
                    clear vol hBxH

                    % SO4: 
                    % mixing ratio to ug/m3
                    unitconv = 1e6 / 8.314 .* (tP./tT); clear tP tT
                    so4h = so4h .* unitconv(:,:,:,1) .* SO4MW .* RH_growth_sia;

                    % BC:
                    bch = bch .* unitconv(:,:,:,1) .* BCMW ; % no growth

                    % Regrid to lat-lon
                    so2h = cube2latlon2d(so2h,lat,lon,xLAT,xLON);
                    pm25h = cube2latlon2d(pm25h(:,:,:,1),lat,lon,xLAT,xLON);
                    so4h = cube2latlon2d(so4h,lat,lon,xLAT,xLON);
                    bch = cube2latlon2d(bch(:,:,:,1),lat,lon,xLAT,xLON);

                    % ==== sampling at TROPOMI overpass time [Algorithm 2] ======
                    % this has the same problem as algorithm 1 but is easier to
                    % read, more concise.
                    % see Sat_overpass_time_cheatsheet.xlsx (on onedrive) for reasonings
                    loc = (13-Hr)*15;
                    if loc > 180
                        loc = -180+(loc-180);
                    end
                    LonInd = find(tLON >= loc-7.5 & tLON < loc+7.5);
                    so2d(:,LonInd,:) = so2d(:,LonInd,:) + so2h(:,LonInd,:);
                    pm25d = pm25h + pm25d;
                    so4d = so4h + so4d;
                    bcd = bch + bcd;
                end

                % save daily data
                so2 = so2d; clear so2d
                pm25 = pm25d./24; clear pm25d
                so4 = so4d./24; clear so4d
                bc = bcd./24; clear bcd
                savenc(sfname,so2,so4,pm25,bc,tLAT, tLON)
                fprintf('%s saved\n',sfname)
            end
            % Calc monthly mean
            so2m = so2m + so2;
            so4m = so4m + so4;
            pm25m = pm25m + pm25;
            bcm = bcm + bc;

            toc
        end
        so2 = so2m ./ Dy; clear so2m
        so4 = so4m ./ Dy; clear so4m
        pm25 = pm25m ./ Dy; clear pm25m
        bc = bcm ./ Dy; clear bcm
        % save monthly data
        sfname = sprintf('%s/GCHP_SO2_SO4_BC_PM25_%s_%.2d.nc',OutDir,SimName,Mn);
        savenc(sfname,so2,so4,pm25,bc,tLAT, tLON)
        fprintf('%s saved\n',sfname)

        so2y = so2y + so2;
        so4y = so4y + so4;
        pm25y = pm25y + pm25;
        bcy = bcy + bc;

    end
    % Calc annual mean
    so2 = so2y./Mn; clear so2y
    so4 = so4y./Mn; clear so4y
    pm25 = pm25y./Mn; clear pm25y
    bc = bcy./Mn; clear bcy
    % save monthly data
    sfname = sprintf('%s/GCHP_SO2_SO4_BC_PM25_%s_annual.nc',OutDir,SimName);
    savenc(sfname,so2,so4,pm25,bc,tLAT, tLON)
    fprintf('%s saved\n',sfname)

end

disp('Done')


%% FUNCTION
function var = cube2latlon2d(var1,lat,lon,xLAT,xLON)

        var2 = scatteredInterpolant(lon(:),lat(:),var1(:),'linear');
        var  = var2(xLON,xLAT);

        var(var<0) = 0;

end


function savenc (sfname,so2,so4,pm25, bc, tLAT, tLON)

if exist(sfname,'file')
    delete(sfname)
end
% initiate
nccreate(sfname, 'latitude', 'Dimensions', {'lat', length(tLAT)});
nccreate(sfname, 'longitude', 'Dimensions', {'lon', length(tLON)});
nccreate(sfname, 'so2', 'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});
nccreate(sfname, 'so4', 'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});
nccreate(sfname, 'pm25', 'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});
nccreate(sfname, 'bc', 'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});
% write data 
ncwrite(sfname, 'latitude', tLAT);
ncwrite(sfname, 'longitude', tLON);
ncwrite(sfname, 'so2', so2);
ncwrite(sfname, 'so4', so4);
ncwrite(sfname, 'pm25', pm25);
ncwrite(sfname, 'bc', bc);

% adding attribute
ncwriteatt(sfname, 'latitude', 'units', 'degrees_north');
ncwriteatt(sfname, 'longitude', 'units', 'degrees_east');

ncwriteatt(sfname, 'so2', 'units', 'molec/cm2');  
ncwriteatt(sfname, 'so2', 'long_name', 'Sulfur Dioxide Vertical Column Density at S5P overpass time');

ncwriteatt(sfname, 'so4', 'units', 'ug/m3');  
ncwriteatt(sfname, 'so4', 'long_name', 'Surface sulfate concentration at RH 35%');

ncwriteatt(sfname, 'pm25', 'units', 'ug/m3');  
ncwriteatt(sfname, 'pm25', 'long_name', 'Surface PM25 concentration at RH 35%');

ncwriteatt(sfname, 'bc', 'units', 'ug/m3');  
ncwriteatt(sfname, 'bc', 'long_name', 'Surface black carbon concentration at RH 35%');

end



function fig = mapx3(LAT_GC,LON_GC,SO2_GC,LAT_TROPOMI,LON_TROPOMI,SO2_TROPOMI)

[latmin,latmax,lonmin,lonmax,Positions,cbvserpos,notepos,PaperSize] = Map6_FigPara('global');
fz = 14;

mgc = mean(mean(SO2_GC,'omitnan'),'omitnan');
rng = [0.1 5].*mgc;

% % load mask file
% load('NewLandMask-0.05.mat') % 'LAT','LON'
% wMASK = double(MASKp); clear MASKp
% wMASK_land = wMASK(:,:,2); clear wMASK
% wMASK_land(wMASK_land==0) = NaN;
% wMASK_land(abs(wMASK_land)>0) = 1;
% Lat = double(LAT); Lon = double(LON);
% 
% wMASK_land = interp2(Lon,Lat',wMASK_land,lon ,lat','nearest');

fig = figure('WindowState','maximized','Position',PaperSize);

% making subplot 1
ps = 2;
subplot('Position',Positions(ps,:))
worldmap([latmin latmax],[lonmin lonmax]);
surfm(LAT_GC,LON_GC, SO2_GC)
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

cb2=colorbar('vertical','fontsize',fz,'fontweight', 'bold','Position',cbvserpos(ps,:));
cblb = sprintf('GC SO_2 [molec/cm^2]');
set(get(cb2,'YLabel'),'string',cblb,'fontsize',fz,'fontweight','bold','FontName','Helvetica');
hold on



% making subplot 2
ps = 4;
subplot('Position',Positions(ps,:))
worldmap([latmin latmax],[lonmin lonmax]);
surfm(LAT_TROPOMI,LON_TROPOMI, SO2_TROPOMI)
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

cb2=colorbar('vertical','fontsize',fz,'fontweight', 'bold','Position',cbvserpos(ps,:));
cblb = sprintf('TROPOMI SO_2 [molec/cm^2]');
set(get(cb2,'YLabel'),'string',cblb,'fontsize',fz,'fontweight','bold','FontName','Helvetica');
hold on




% making subplot 3 - raletive difference
rng = [0.5 1.5];
mapdata = interp2(LON_TROPOMI,LAT_TROPOMI',SO2_TROPOMI,LON_GC,LAT_GC');
mapdata = SO2_GC./mapdata;

ps = 6;
subplot('Position',Positions(ps,:))
worldmap([latmin latmax],[lonmin lonmax]);
surfm(LAT_GC,LON_GC, mapdata)
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

cb2=colorbar('vertical','fontsize',fz,'fontweight', 'bold','Position',cbvserpos(ps,:));
cblb = sprintf('%s SO_2 [unitless]','\Delta');
set(get(cb2,'YLabel'),'string',cblb,'fontsize',fz,'fontweight','bold','FontName','Helvetica');
hold on

end

 = 1;
    cm(cm<0) = 0;
    colormap(gca,flipud(cm));
end
set(gca,'clim',rng); %  colorscale = log or linear

cb2=colorbar('vertical','fontsize',fz,'fontweight', 'bold','Position',cbvserpos(ps,:));
cblb = sprintf('%s SO_2 [unitless]','\Delta');
set(get(cb2,'YLabel'),'string',cblb,'fontsize',fz,'fontweight','bold','FontName','Helvetica');
hold on

end



