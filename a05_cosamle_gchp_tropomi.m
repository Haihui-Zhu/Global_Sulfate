% July-22-2024
% Haihui Zhu

% Cosample gchp so2 with tropomi so2

clear 
close all

% ---- Input & Switches --------------
year = 2021; % year of tropomi data
SimName  = 'htap_2018'; % change here
gchp_dir  = './05.GCHP_outputs/2.htap_2018'; % change here
qcstr = 'CF03-SZA50-QA75';
trpomi_dir = ['/storage1/fs1/rvmartin2/Active/haihuizhu/02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_' qcstr];


Mns = 1:12;
% ---- End of input ----------------------
daysinmonth = [31 28 31 30 31 30 31 31 30 31 30 31];

% Now start looping:
for Mn = Mns
    tic
    for Dy = 1:daysinmonth(Mn)

        tro_fname = sprintf('%s/Tropomi_Regrid_%d%.2d%.2d_%s.nc',trpomi_dir,year,Mn,Dy,qcstr);
        % search if tropomi so2 exists
        if exist(tro_fname, 'file')

            % load tropomi so2 lat lon
            tro_so2 = ncread(tro_fname, 'so2');

            if ~exist('tlat', 'var')
                tlat = ncread(tro_fname, 'lat');
                tlon = ncread(tro_fname, 'lon');
            end

            % load gchp so2 lat lon
            sim_fname = sprintf('%s/GCHP_SO2_SO4_BC_PM25_%s_%.2d%.2d.nc', gchp_dir, SimName, Mn, Dy);
            % sim_fname = sprintf('%s/GCHP_SO2_SO4_PM25_%s_%.2d%.2d.nc', gchp_dir, SimName, Mn, Dy);
            sim_so2 = ncread(sim_fname, 'so2');

            if ~exist('sim_lat', 'var')
                sim_lat = ncread(sim_fname, 'latitude');
                sim_lon = ncread(sim_fname, 'longitude');
            end

            % create mask and a mat of N tracking number of days sampled.
            mask = nan(size(tro_so2));
            mask(abs(tro_so2) > 0) = 1;
            % check for any 0:
            ind = find(tro_so2 == 0);

            if ~isempty(ind)
                warning('%s boxes = 0\n', length(ind))
            end

            % regridding gchp so2 to tropomi
            sim_so2 = interp2(sim_lon, sim_lat', sim_so2, tlon, tlat');

            % appy mask
            sim_so2 = sim_so2 .* mask;

            % add it to existing gchp_regrid_tro_so2
            if ~exist('so2_sim_out', 'var')
                so2_sim_out = sim_so2;
                so2_tro_out = tro_so2;
                N = mask;
            else
                so2_sim_out(:, :, 2) = sim_so2;
                so2_tro_out(:, :, 2) = tro_so2;
                N(:, :, 2) = mask;

                so2_sim_out = sum(so2_sim_out, 3, 'omitnan');
                so2_tro_out = sum(so2_tro_out, 3, 'omitnan');
                N = sum(N, 3, 'omitnan'); % sum of nan = 0;
            end

            clear sim_so2 tro_so2
        else
            fprintf('%s not found\n', tro_fname)
        end
    end
    toc
    % gchp_regrid_so2 ./ N to get to cosampled mean
    N(N==0) = nan;
    so2_sim = so2_sim_out ./ N;
    so2_tro = so2_tro_out./N;

    % save monthly data
    sfname = sprintf('%s/gchp_so2_cosampled_tropomi_%s_%.2d.nc',trpomi_dir,SimName,Mn);
    savenc(sfname, so2_sim, so2_tro, tlat, tlon)
    fprintf('%s saved\n',sfname)
    
    clear so2_*
end
%}

% A separate process to calc annual mean (do not want to load too many
% matrices at the same time)
n=0;
for Mn = Mns
    % load monthly data
    fname = sprintf('%s/gchp_so2_cosampled_tropomi_%s_%.2d.nc',trpomi_dir,SimName,Mn);
    if exist(fname, 'file')
        fprintf('%s added\n',fname)
        n = n+1;
        mso2_sim = ncread(fname, 'so2_gchp');
        mso2_tro = ncread(fname, 'so2_tro');
        % initiate or add to existing vars
        if ~exist('so2_gchp', 'var') == 1
            so2_gchp = mso2_sim;
            so2_tro = mso2_tro;
            tlat = ncread(fname,'lat');
            tlon = ncread(fname,'lon');
        else
            so2_gchp(:, :, n) = mso2_sim;
            so2_tro(:, :, n) = mso2_tro;
        end
    end

end
% save
so2_sim = mean(so2_gchp, 3, 'omitnan');
so2_tro = mean(so2_tro, 3, 'omitnan');
sfname = sprintf('%s/gchp_so2_cosampled_tropomi_%s_annual.nc', trpomi_dir, SimName); 
savenc(sfname, so2_sim, so2_tro, tlat, tlon)
fprintf('%s saved\n',sfname)

%% FUNCTION
function savenc (sfname, so2_sim, so2_tro, tLAT, tLON)

if exist(sfname,'file')
    delete(sfname)
end
% initiate
nccreate(sfname, 'lat', 'Dimensions', {'lat', length(tLAT)});
nccreate(sfname, 'lon', 'Dimensions', {'lon', length(tLON)});
nccreate(sfname, 'so2_gchp', 'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});
nccreate(sfname, 'so2_tro',  'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});

% write data 
ncwrite(sfname, 'lat', tLAT);
ncwrite(sfname, 'lon', tLON);
ncwrite(sfname, 'so2_gchp', so2_sim);
ncwrite(sfname, 'so2_tro', so2_tro);

% adding attribute
ncwriteatt(sfname, 'lat', 'units', 'degrees_north');
ncwriteatt(sfname, 'lat', 'long_name', 'latitude')

ncwriteatt(sfname, 'lon', 'units', 'degrees_east');
ncwriteatt(sfname, 'lon', 'long_name', 'longitude')

ncwriteatt(sfname, 'so2_gchp', 'units', 'molec/cm2');  
ncwriteatt(sfname, 'so2_gchp', 'long_name', 'GCHP simulated sulfur dioxide vertical column density at S5P overpass time; Cosampled with TROPOMI'); 

ncwriteatt(sfname, 'so2_tro', 'units', 'molec/cm2');
ncwriteatt(sfname, 'so2_tro', 'long_name', 'TROPOMI sulfur dioxide vertical column density at S5P overpass time; Reduced noise using GCHP background SO2 VCD');
end

