% July-22-2024
% Haihui Zhu

% Cosample gchp so2 with tropomi so2

clear 
close all

% ---- Input & Switches --------------
year = 2021; % year of tropomi data
SimName  = 'ceds_2021'; % change here
gchp_dir  = './05.GCHP_outputs/4.ceds_2021'; % change here
qcstr = 'CF03-SZA50-QA75';
trpomi_dir = ['./02.TROPOMI_SO2_Ref/NASA_SO2_Tesellation_' qcstr];

unitconverter = 2.69e16;

Mns = 1:12;
% ---- End of input ----------------------
daysinmonth = [31 28 31 30 31 30 31 31 30 31 30 31];

% read the Pacific mask
fname = './utils/pacific_mask.nc';
pmask= ncread(fname, 'mask')'; % transpose to make it lat x lon
plat = ncread(fname, 'latitude');
plon = ncread(fname, 'longitude');

% Now start looping:

for Mn = Mns
    tic
    for Dy = 1:daysinmonth(Mn)

        tro_fname = sprintf('%s/Tropomi_Regrid_%d%.2d%.2d_%s.nc', trpomi_dir, year, Mn, Dy, qcstr);
        tro_o3_fname = sprintf('%s/Regrid_tropomi_O3_%d%.2d%.2d_%s.nc', trpomi_dir, year, Mn, Dy, qcstr);
        % search if tropomi so2 exists
        if exist(tro_fname, 'file')

            % load tropomi so2 lat lon
            tro_so2 = ncread(tro_fname, 'so2');
            % load tropomi o3
            if exist(tro_o3_fname, 'file')
                o3 = ncread(tro_o3_fname,'O3')./unitconverter;
            else
                o3 = ones(size(tro_so2));
            end

            if ~exist('tlat', 'var')
                tlat = ncread(tro_fname, 'lat');
                tlon = ncread(tro_fname, 'lon');
                pmask = pmask((plat <= max(tlat) & plat >= min(tlat) - .01), (plon <= max(tlon)+0.01 & plon >= min(tlon)));
                % disp('pmask size:')
                % disp(size(pmask))
            end

            % load gchp so2 lat lon
            switch SimName
                case 'ceds_2021'
                    sim_fname = sprintf('%s/GCHP_SO2_SO4_PM25_%s_%.2d%.2d.nc', gchp_dir, SimName, Mn, Dy);
                otherwise
                    sim_fname = sprintf('%s/GCHP_SO2_SO4_BC_PM25_%s_%.2d%.2d.nc', gchp_dir, SimName, Mn, Dy);
            end
            sim_so2 = ncread(sim_fname, 'so2');

            if ~exist('sim_lat', 'var')
                sim_lat = ncread(sim_fname, 'latitude');
                sim_lon = ncread(sim_fname, 'longitude');
            end

            % create mask and a mat of N tracking number of days sampled.
            mask = nan(size(tro_so2));
            mask(abs(tro_so2) > 0) = 1;
            mask(o3>400) == nan;
            % check for any 0:
            ind = find(tro_so2 == 0);

            if ~isempty(ind)
                warning('%s boxes = 0\n', length(ind))
            end

            % regridding gchp so2 to tropomi
            sim_so2 = interp2(sim_lon, sim_lat', sim_so2, tlon, tlat');

            % appy mask
            sim_so2 = sim_so2 .* mask;

            % find background signal over the Pacific
            bg_gchp = mean(sim_so2.*pmask, 2, 'omitnan'); % take the average along longitude
            bg_tro  = mean(tro_so2.*pmask, 2, 'omitnan');
            noise = bg_tro - bg_gchp;
            tro_so2_nr = tro_so2 - noise; % lead to negative values

            % add it to existing gchp_regrid_tro_so2
            if ~exist('so2_sim_out', 'var')
                so2_sim_out = sim_so2;
                so2_tro_out = tro_so2;
                so2_tro_nr_out = tro_so2_nr;
                bg_gchp_out = bg_gchp;
                bg_tro_out = bg_tro;
                % disp('size of bg_gchp_out:')
                % disp(size(bg_gchp_out))
                N = mask;
            else
                so2_sim_out(:, :, 2) = sim_so2;
                so2_tro_out(:, :, 2) = tro_so2;
                so2_tro_nr_out(:, :, 2) = tro_so2_nr;
                bg_gchp_out(:, Dy) = bg_gchp;
                bg_tro_out(:, Dy) = bg_tro;
                N(:, :, 2) = mask;

                so2_sim_out = sum(so2_sim_out, 3, 'omitnan');
                so2_tro_out = sum(so2_tro_out, 3, 'omitnan');
                so2_tro_nr_out = sum(so2_tro_nr_out, 3, 'omitnan');
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
    so2_tro = so2_tro_out ./ N;
    so2_tro_nr = so2_tro_nr_out ./ N;
    bg_gchp = mean(bg_gchp_out, 2, 'omitnan');
    bg_tro  = mean(bg_tro_out, 2, 'omitnan');

    % save monthly data
    sfname = sprintf('%s/gchp_so2_cosampled_tropomi_%s_noisereduced_wO3_%.2d.nc',trpomi_dir,SimName,Mn);
    savenc(sfname, so2_sim, so2_tro, so2_tro_nr, bg_gchp, bg_tro, tlat, tlon)
    fprintf('%s saved\n',sfname)
    
    clear so2_*
end
%}

% A separate process to calc annual mean (do not want to load too many
% matrices at the same time)
n=0;
for Mn = 1:12
    % load monthly data
    fname = sprintf('%s/gchp_so2_cosampled_tropomi_%s_noisereduced_wO3_%.2d.nc',trpomi_dir,SimName,Mn);
    if exist(fname, 'file')
        fprintf('%s added\n',fname)
        n = n+1;
        mso2_sim = ncread(fname, 'so2_gchp');
        mso2_tro = ncread(fname, 'so2_tro');
        mso2_tro_nr = ncread(fname, 'so2_tro_nr');
        mbg_sim = ncread(fname, 'so2_bg_gchp');
        mbg_tro = ncread(fname, 'so2_bg_tro');
        % initiate or add to existing vars
        if ~exist('so2_gchp', 'var') == 1
            so2_gchp = mso2_sim;
            so2_tro = mso2_tro;
            so2_tro_nr = mso2_tro_nr;
            bg_gchp = mbg_sim;
            bg_tro  = mbg_tro;
            tlat = ncread(fname,'lat');
            tlon = ncread(fname,'lon');
        else
            so2_gchp(:, :, n) = mso2_sim;
            so2_tro(:, :, n) = mso2_tro;
            so2_tro_nr(:, :, n) = mso2_tro_nr;
            bg_gchp(:, n) = mbg_sim;
            bg_tro(:, n) = mbg_tro;
        end
    end

end
% save
so2_sim = mean(so2_gchp, 3, 'omitnan');
so2_tro = mean(so2_tro, 3, 'omitnan');
so2_tro_nr = mean(so2_tro_nr, 3, 'omitnan');
bg_gchp = mean(bg_gchp, 2, 'omitnan');
bg_tro = mean(bg_tro, 2, 'omitnan');
sfname = sprintf('%s/gchp_so2_cosampled_tropomi_%s_noisereduced_wO3_annual.nc', trpomi_dir, SimName); 
savenc(sfname, so2_sim, so2_tro, so2_tro_nr, bg_gchp, bg_tro, tlat, tlon)
fprintf('%s saved\n',sfname)

%% FUNCTION
function savenc (sfname, so2_sim, so2_tro, so2_tro_nr, bg_gchp, bg_tro, tLAT, tLON)

if exist(sfname,'file')
    delete(sfname)
end
% initiate
nccreate(sfname, 'lat', 'Dimensions', {'lat', length(tLAT)});
nccreate(sfname, 'lon', 'Dimensions', {'lon', length(tLON)});
nccreate(sfname, 'so2_gchp', 'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});
nccreate(sfname, 'so2_tro', 'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});
nccreate(sfname, 'so2_tro_nr', 'Dimensions', {'lat', length(tLAT), 'lon', length(tLON)});
nccreate(sfname, 'so2_bg_gchp', 'Dimensions', {'lat', length(tLAT)});
nccreate(sfname, 'so2_bg_tro',  'Dimensions', {'lat', length(tLAT)});
% write data 
ncwrite(sfname, 'lat', tLAT);
ncwrite(sfname, 'lon', tLON);
ncwrite(sfname, 'so2_gchp', so2_sim);
ncwrite(sfname, 'so2_tro', so2_tro);
ncwrite(sfname, 'so2_tro_nr', so2_tro_nr);
ncwrite(sfname, 'so2_bg_gchp', bg_gchp);
ncwrite(sfname, 'so2_bg_tro', bg_tro);

% adding attribute
ncwriteatt(sfname, 'lat', 'units', 'degrees_north');
ncwriteatt(sfname, 'lat', 'long_name', 'latitude')

ncwriteatt(sfname, 'lon', 'units', 'degrees_east');
ncwriteatt(sfname, 'lon', 'long_name', 'longitude')

ncwriteatt(sfname, 'so2_gchp', 'units', 'molec/cm2');  
ncwriteatt(sfname, 'so2_gchp', 'long_name', 'GCHP simulated sulfur dioxide vertical column density at S5P overpass time; Cosampled with TROPOMI'); 

ncwriteatt(sfname, 'so2_tro', 'units', 'molec/cm2');
ncwriteatt(sfname, 'so2_tro', 'long_name', 'TROPOMI sulfur dioxide vertical column density at S5P overpass time');

ncwriteatt(sfname, 'so2_tro_nr', 'units', 'molec/cm2');
ncwriteatt(sfname, 'so2_tro_nr', 'long_name', 'Similar to so2_tro, with reduced noise using GCHP background SO2 VCD');

ncwriteatt(sfname, 'so2_bg_gchp', 'units', 'molec/cm2');
ncwriteatt(sfname, 'so2_bg_gchp', 'long_name', 'GCHP simulated BACKGROUND sulfur dioxide vertical column density at S5P overpass time; Cosampled with TROPOMI; Background definition see pacific_mask.nc');

ncwriteatt(sfname, 'so2_bg_tro', 'units', 'molec/cm2');
ncwriteatt(sfname, 'so2_bg_tro', 'long_name', 'TROPOMI BACKGROUND sulfur dioxide vertical column density at S5P overpass time; Background definition see pacific_mask.nc');
end

