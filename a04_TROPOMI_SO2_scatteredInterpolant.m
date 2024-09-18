% This script reads the raw TROPOMI SO2 VCD from NASA L2 offline product

% clear % uncomment for interactive jobs
close all
addpath('/storage1/fs1/rvmartin/Active/haihuizhu/1.code')
addpath('/storage1/fs1/rvmartin2/Active/haihuizhu/10.SpatialETA')
RootDir = '/storage1/fs1/rvmartin/Active';

TROPOMI_SO2_indir = '/storage1/fs1/rvmartin/Active/haihuizhu/4.SPARTAN_SO4/02.TROPOMI_SO2_Ref/NASA';

SaveDir = sfmkdir(sprintf('./03.TROPOMI_SO2_Result/04.NASA_L2_RPRO_VCD_CF03'));

% regrid settings
olat_res  = 0.1;
olon_res  = 0.1;
finegrid_str   = '0.1x0.1deg';
reftime = datetime([2010, 01,01,0,0,0]); % reference time in the raw data
tess_nvar_finegrid = 1; % number of variables (SO2)

daily = 1;
% Mns = 1; % uncomment for interactive jobs

% other constants
daysinmonth = [31 28 31 30 31 30 31 31 30 31 30 31]; 

Av = 6.02e23; % avogadro's number

% now start processing
tlat = -90+olat_res/2 :olat_res: 90-olat_res/2; % grid center
tlon = -180+olon_res/2:olon_res:180-olon_res/2; % grid center
[mlon,mlat] = meshgrid(tlon,tlat);
clat = -90  :olat_res: 90; % grid corner
clon = -180 :olon_res:180; % grid corner

for yr = 2019
    if daily == 1
    for mn = Mns
        for dy = 1:daysinmonth(mn)
            tic

            files = dir(sprintf('%s/S5P_RPRO_L2__SO2__*%d%.2d%.2dT*_*.nc',TROPOMI_SO2_indir,yr,mn,dy)); 
            
            temp = [];
            templat = [];
            templon = [];
            % to debug
%             tempsza = [];
%             tempqa = [];

            for ff = 1:length(files)
                tfile = files(ff).name;
                if tfile(1) ~= '.'

                    tfullpath = sprintf('%s/%s',TROPOMI_SO2_indir,tfile);
                    
                    % Samping method:
                    % Time in the title is the utc time and matches the time var in the file.
                    % So collect all pixels with UTC time on the desired date.
                    % This would lead to minor discrepancy vs GCHP when comparing at
                    % a daily scale. So this would only recommend for monthly
                    % analysis or at coarser temporal resolutions.

                    % Compare time in the title and time in the file:
                    fprintf('reading %s\n',tfile);
                    ttime = ncread(tfullpath,'PRODUCT/time')./3600/24;
                    delta_time = ncread(tfullpath,'PRODUCT/delta_time')./3600/24/1000;
                    ttime = reftime+delta_time+ttime;
                    
                    % NOTE: This is how GCHP SO2 was sampled.
                    % loc = (13-Hr)*15;
                    % if loc > 180
                    %     loc = -180+(loc-180);
                    % end
                    % LonInd = find(tLON >= loc-7.5 & tLON < loc+7.5);
                    % so2d(:,LonInd,:) = so2d(:,LonInd,:) + so2h(:,LonInd,:);

                    % Now reading to get the VCD
                    vcd = ncread(tfullpath,'PRODUCT/sulfurdioxide_total_vertical_column'); % in mol m-2

                    % quality assurance value filter set to 0.75 to avoid
                    % large coverage of cloud, snow/ice 
                    qa_value = ncread(tfullpath,'PRODUCT/qa_value');

                    % cloud fraction
                    cf = ncread(tfullpath,'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_fraction_intensity_weighted');


                    % Solar Zenith Angle below 60 degree
                    sza = ncread(tfullpath,'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle');

                    % read lat lon and corner lat lon
                    lat = ncread(tfullpath,'PRODUCT/latitude');
                    lon = ncread(tfullpath,'PRODUCT/longitude');

                    % find the spots that satisfy the quality check and
                    % fall in the same day


                    good_idx1 = find(...
                        qa_value >0.75 & ...
                        cf  < 0.3 & ...
                        sza < 60 & ...
                        ttime>=datetime([yr,mn,dy,0,0,0]));
                    lon(~good_idx1) = nan;
                    lat(~good_idx1) = nan;
                    vcd(~good_idx1) = nan;

                    % remove nan or inf in lat and lon
                    ind = find(isnan(lat) | isnan(lon) | isinf(lat) | isinf(lon));
                    lat(ind) = [];
                    lon(ind) = [];
                    vcd(ind) = [];

                    temp = [temp vcd];
                    templat = [templat lat];
                    templon = [templon lon];

                    % to debug:
                    %{
                    sza(~good_idx1) = nan;
                    qa_value(~good_idx1) = nan;
                    sza(ind) = [];
                    qa_value(ind) = [];
                    tempsza = [tempsza sza];
                    tempqa = [tempqa qa_value];
                    %}
                end

            end
            % to debug:
            %{
            disp(size(temp))
            disp(size(tempsza))
            disp(size(tempqa))
            
            sza = scatter_regrid(mlat, mlon,tempsza,templat, templon);
            datestr = sprintf('%d%.2d%.2d_sza',yr,mn,dy);
            rng = [0,100];
            submap(sza,tlat,tlon,SaveDir,datestr,rng)

            qa = scatter_regrid(mlat, mlon,tempqa,templat, templon);
            datestr = sprintf('%d%.2d%.2d_qa',yr,mn,dy);
            rng = [0,1];
            submap(qa,tlat,tlon,SaveDir,datestr,rng)
            %}
            % regrid
            so2 = scatter_regrid(mlat, mlon,temp,templat, templon);
            so2 = so2*Av*1e-4; % convert mol m-2 to molec cm-2
            datestr = sprintf('%d%.2d%.2d',yr,mn,dy);
            sfname = [SaveDir, '/', 'SO2_',datestr,'_',finegrid_str,'.mat'];
            save(sfname,'so2','tlat','tlon');
            fprintf('%s saved\n',sfname)
            toc

%             % making a quick map
%             rng = [0,1];% convert to DU
%             so2 = so2/2.69e16;
%             submap(so2,tlat,tlon,SaveDir,datestr,rng)

        end
    end
    end

    for mn = Mns
        tic
        for dy = 1:daysinmonth(mn)
            datestr = sprintf('%d%.2d%.2d',yr,mn,dy);
            fname = [SaveDir, '/', 'SO2_',datestr,'_',finegrid_str,'.mat'];
            load(fname,'so2','tlat','tlon');
            if dy ==1
                mso2 = nan(length(tlat), length(tlon),daysinmonth(mn));
            end
            mso2(:,:,dy) = so2;
            
            clear so2
        end
        so2 = mean(mso2,3,'omitnan');
        datestr = sprintf('%d%.2d',yr,mn);
        sfname = [SaveDir, '/', 'SO2_',datestr,'_',finegrid_str,'.mat'];
        save(sfname,'so2','tlat','tlon');
        fprintf('%s saved\n',sfname)

        % making a quick map
        rng = [0,1];% convert to DU
        so2 = so2/2.69e16;
        submap(so2,tlat,tlon,SaveDir,datestr,rng)
        toc
    end


    for mn = 1:12
        datestr = sprintf('%d%.2d',yr,mn);
        fname = [SaveDir, '/', 'SO2_',datestr,'_',finegrid_str,'.mat'];
        load(fname,'so2','tlat','tlon');
        if mn ==1
            mso2 = nan(length(tlat), length(tlon),12);
        end
        mso2(:,:,mn) = so2;

        clear so2
    end
    so2 = mean(mso2,3,'omitnan');
    datestr = sprintf('%d',yr);
    sfname = [SaveDir, '/', 'SO2_',datestr,'_',finegrid_str,'.mat'];
    save(sfname,'so2','tlat','tlon');
    fprintf('%s saved\n',sfname)

    % making a quick map
    rng = [0,1];% convert to DU
    so2 = so2/2.69e16;
    submap(so2,tlat,tlon,SaveDir,datestr,rng)


end
%% FUNCTIONS
function so2s = scatter_regrid(mlat, mlon,temp,templat, templon)
    var2 = scatteredInterpolant(templon(:),templat(:),temp(:),'linear');
    so2s  = var2(mlon,mlat);
end

function submap(SO2,lat,lon,OutDir,datestr,rng)


figure('windowstate','maximized')
subplot(2,1,1)
worldmap('world')
surfm(lat,lon,SO2)
    
cm = cbrewer('seq','Reds',30,'spline');
cm(cm>1) = 1;
cm(cm<0) = 0;
colormap(gca,cm);

set(gca,'clim',rng)
bordersm
cb = colorbar;
cblb = sprintf('SO2 regrid [DU]');
set(get(cb,'YLabel'),'string',cblb,'fontweight','bold','FontName','Helvetica');

saveas(gcf,[OutDir,'/SO2_regrid_',datestr,'.png'])
end

