%% FUNCTION
function mapout = MapRegion(mapdata,lon,lat,fz,RegionName,maparea)
% addpath('/storage1/fs1/rvmartin/Active/haihuizhu/1.code')
% ===== Data Processing =====
% making fiure
figure('Position',[10 10 800 600])
subplot('Position',[0.1 0.2 0.6 0.6])
if ismember(maparea,'world')
    mapout = worldmap([-55 70],[-170 170]);
else
    mapout = worldmap(maparea);
end

% apply country border & coastlines
bordersm
load coastlines
plotm(coastlat,coastlon)
setm(gca,'Grid','off','MapProjection','miller','parallellabel','off','meridianlabel','off')

if ismember(maparea,'china')
    cn = load('china.province.mat');
    plotm(cn.lat,cn.long)
end
if ismember(maparea,'us')
    bordersm('continental us')
end

surfm(lat,lon,mapdata)

cm = cbrewer('qual','Paired', numel(RegionName) , 'spline');
% cm = cbrewer('div','Spectral', numel(RegionName) , 'spline');
cm(cm>1) = 1;
cm(cm<0) = 0;
colormap(gca, cm)
set(gca, 'clim',[1-0.5 numel(RegionName)+0.5]) % exclude other country

% color bar definition
cb=colorbar('vertical','fontsize',fz,'Direction','reverse');
set(get(cb,'YLabel'),'string','Regions','fontsize',fz,'FontName','Helvetica','fontweight','bold');
cb.Ticks = (1:numel(RegionName));
cb.TickLabels = RegionName;

fprintf('Done mask map - %s.\n',maparea)

end
