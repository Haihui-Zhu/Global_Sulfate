% Aug-4-2024
% Haihui Zhu

% Read SO2 below PBL, sampled at TROPOMIT overpass time.
clear 
close all
addpath('/storage1/fs1/rvmartin/Active/haihuizhu/1.code') % for mapping bordersm

% ---- Input & Switches --------------
YEARS = 2019; % year range of interest, currently means year of simulation
SimName  = 'ceds_2019'; % change here
% InDir  = '/nobackup/hzhu3/1.ceds_2018/OutputDir/'; % change here
% OutDir = '/nobackup/hzhu3/1.ceds_2018/processed/'; % change here
InDir  = '~/my-projects2/5.GEOS-Chem/7.GCHP-13.4.0/4.C90_GFEDd_2m_CEDS_vert_diel/OutputDir_ACAG/'; % change here
OutDir = '~/my-projects2/5.GEOS-Chem/7.GCHP-13.4.0/4.C90_GFEDd_2m_CEDS_vert_diel/Processed/'; % change here

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
    
    for Mn = Mns
        so2m = zeros(length(tLAT),length(tLON));

        for Dy = 1:daysinmonth(Mn)
            tic
            sfname = sprintf('%s/GCHP_SO2_PBL_%s_%.2d%.2d.nc',OutDir,SimName,Mn,Dy);
            if exist(sfname,'file') % for re-run without repeating work
                so2 = ncread(sfname,'so2');
            
            else
                so2d = zeros(length(tLAT),length(tLON));

                for Hr = 0:23
                    fname = sprintf('%s/GEOSChem.ACAG.%d%.2d%.2d_%.2d30z.nc4',InDir,Yr,Mn,Dy,Hr);
                    
                    % SO2:
%                     so2h = double(ncread(fname,'SpeciesConcVV_SO2')); % mol/mol dry air
                    so2h = double(ncread(fname,'SpeciesConc_SO2')); % mol/mol dry air
                    
                    % zero out values above PBLH:
                    pblh = double(ncread(fname,'Met_PBLH')); % m
                    hBxH = double(ncread(fname,'Met_BXHEIGHT')); % m
                    BxH_cs = cumsum(hBxH,4); % cummulated sum of box height

                    for ii = 1:size(BxH_cs,1)
                        for jj =  1:size(BxH_cs,2)
                            for ff = 1:size(BxH_cs,3)
                                ind = find( BxH_cs(ii,jj,ff,:)>pblh(ii,jj,ff) );
                                so2h(ii,jj,ff,ind) = 0;
                            end
                        end
                    end
%                     % print out average # of layers within PBL
%                     layers = zeros(size(so2h));
%                     layers(so2h>0)=1;
%                     layers = sum(layers,4);
%                     quickmap(pblh,BxH_cs,layers,lat,lon,sprintf('layers-pbl-%.2d-%.2d.png',Dy,Hr));


                    % unit conversion [mol/mol dry air] to [molecule per cm2]:
                    % mol/mol dry air => mol/m3 dry air => mol/m2 dry air =>
                    % mol/m2 surface area => moelc/m2 => molec/cm2
                    tT = double(ncread(fname,'Met_T')); % K
                    tP = 100*double(ncread(fname,'Met_PMID')); % Pa
                    vol = 1.*R.*tT./tP; % m3

                    so2h = sum(so2h./vol.*hBxH, 4); % mol/m2
                    so2h = so2h.* Na .* 1e-4; % molec/cm2
                    clear vol hBxH

                    % Regrid to lat-lon
                    so2h = cube2latlon2d(so2h,lat,lon,xLAT,xLON);

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
                end

                % save daily data
                so2 = so2d; clear so2d
                savenc(sfname,so2,tLAT, tLON)
                fprintf('%s saved\n',sfname)
            end
            % Calc monthly mean
            so2m = so2m + so2;

            toc
        end
        so2 = so2m ./ Dy; clear so2m
        % save monthly data
        sfname = sprintf('%s/GCHP_SO2_PBL_%s_%.2d.nc',OutDir,SimName,Mn);
        savenc(sfname,so2,tLAT, tLON)
        fprintf('%s saved\n',sfname)

        so2y = so2y + so2;

    end
    % Calc annual mean
    so2 = so2y./Mn; clear so2y
    % save monthly data
    sfname = sprintf('%s/GCHP_SO2_PBL_%s_annual.nc',OutDir,SimName);
    savenc(sfname,so2,tLAT, tLON)
    fprintf('%s saved\n',sfname)

end

disp('Done')


%% FUNCTION
function var = cube2latlon2d(var1,lat,lon,xLAT,xLON)

        var2 = scatteredInterpolant(lon(:),lat(:),var1(:),'linear');
        var  = var2(xLON,xLAT);

        var(var<0) = 0;

end

function savenc (sfname,so2, tLAT, tLON)

if exist(sfname,'file')
    delete(sfname)
end
% initiate
nccreate(sfname, 'latitude', 'Dimensions', {'lat', length(tLAT)});
nccreate(sfname, 'longitude', 'Dimensions', {'lon', length(tLON)});
nccreate(sfname, 'so2', 'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});
% write data 
ncwrite(sfname, 'latitude', tLAT);
ncwrite(sfname, 'longitude', tLON);
ncwrite(sfname, 'so2', so2);

% adding attribute
ncwriteatt(sfname, 'latitude', 'units', 'degrees_north');
ncwriteatt(sfname, 'longitude', 'units', 'degrees_east');

ncwriteatt(sfname, 'so2', 'units', 'molec/cm2');  
ncwriteatt(sfname, 'so2', 'long_name', 'GCHP SO2 density below PBL at S5P overpass time');

end


function quickmap(pblh, csbxh, layer,lat,lon,fname)
csbxh = csbxh(:,:,:,4); % see if the 4th layer is hgiher than pblh (should be around 6th)
pblh  = permute(pblh,[2,1,3]);
csbxh = permute(csbxh,[2,1,3]);
layer = permute(layer,[2,1,3]);

figure('WindowState','maximized')

subplot(3,1,1)
worldmap('World')
for ff = 1:6
    surfm(lat(:,:,ff),lon(:,:,ff),pblh(:,:,ff))
end
colorbar
bordersm
set(gca,'clim',[-500 1000])
title('pblh')

subplot(3,1,2)
worldmap('World')
for ff = 1:6
    surfm(lat(:,:,ff),lon(:,:,ff),csbxh(:,:,ff))
end
colorbar
bordersm
set(gca,'clim',[-500 1000])
title('cs bxh')

subplot(3,1,3)
worldmap('World')
for ff = 1:6
    surfm(lat(:,:,ff),lon(:,:,ff),layer(:,:,ff))
end
colorbar
bordersm
set(gca,'clim',[-2 4])
title('# of layers')


saveas(gcf,fname)
fprintf('%s saved\n',fname)
close 

end
