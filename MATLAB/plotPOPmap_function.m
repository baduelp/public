function plotPOPmap_function(PopArray, LatArray, LongArray, PopIND, PopVals, PopVals_name,...
        shapeMarker, shapeLegend, shapeSize, shapeColor, contourbool, printdir, titlename, figname, printbool)

% requires landareas.shp, worldrivers.shp, World_Countries.shp from WORLDDATAMAP package and white_red_10_cmap.mat, white_green_10_cmap.mat files to be downloaded in Matlab directory

% PopArray : vector of marker names
% LatArray : vector of marker latitudes
% LonArray : vector of marker longitudes 
% PopIND : vector of marker indexes for marker properties such as defined in shapeMarker e.g. [1 2 1 1 1] 
% PopVals : vector of marker values (if positive and max >1 from white to red, if negative and positive from green to red, if not indexed by shapeColor)
% PopVals_name : name of the marker values e.g. 'shape'
% shapeMarker : cell array of shapes to be used based on PopIND e.g. {'o' 's'} 
% shapeLegend : cell array of legends to be used based on PopIND e.g. {'circle' 'square'}
% shapeSize : vector of marker sizes based on PopIND e.g. [5 6]
% shapeColor : cell array of marker colors e.g. {[0 0 0] [1 0 0]}
% contourbool : if set to 1 country borders will be drawn
% printdir : directory for printing figure
% titlename : title of the map
% figname : name of the figure file
% printbool : if set to 1 a figure will be printed

beep off
close gcf
figure1 = figure('Color',[1 1 1]);
set(figure1,'Resize', 'on','PaperOrientation', 'landscape','PaperUnits',...
    'points','PaperSize', [900 400], 'PaperPosition', [0 0 900 400],...
    'InvertHardcopy','off');
hold on
set(gca, 'XColor', 'k', 'YColor', 'k');
lgdname=''; contourname='';

PositionOffset=zeros(length(PopArray),1); 
LatOffset=0.2*ones(length(PopArray),1); 

for i=1:length(PopArray)
    if LongArray(i)>15
        PositionAlignment{i}='center';
        VertAlignment{i}='bottom';
    else
        PositionAlignment{i}='center';
        VertAlignment{i}='bottom';
    end
end

latlim = [floor(min(LatArray)/5)*5 ceil(max(LatArray)/5)*5]; 
lonlim = [floor(min(LongArray)/5)*5 ceil(max(LongArray)/5)*5]; 

axesm('mercator','Frame','on', 'FEdgeColor','k')
setm(gca,'MapProjection','mercator', 'FLatLimit',[0 8],'MapLatLimit', latlim, 'MapLonLimit', lonlim) %
gridm('Gcolor', [0 0 0])
tightmap
lat_int=min(ceil((latlim(2)-latlim(1))/10/5)*5,10);
lon_int=min(ceil((lonlim(2)-lonlim(1))/10/5)*5,20);
setm(gca,'MeridianLabel','on', 'ParallelLabel','on', 'MLabelParallel',floor(min(LatArray)/10)*10,'MLabelLocation',lon_int,'PLabelMeridian',floor(min(LongArray)/10)*10,'PLabelLocation',lat_int)
setm(gca,'Grid','on', 'Gcolor','k', 'mlinelocation',lon_int,'plinelocation',lat_int)
   

landareas = shaperead('landareas.shp','UseGeoCoords', true,'BoundingBox',[lonlim',latlim']);
% 

totLat=[];
totLong=[];
for i=1:length(landareas) %14
    totLat=vertcat(totLat,landareas(i).Lat');
    totLong=vertcat(totLong,landareas(i).Lon');
    geoshow (landareas(i).Lat,landareas(i).Lon, 'DisplayType', 'polygon','FaceColor',[0.9 0.9 0.9]);
end

geoshow (totLat,totLong,  'Color', 'k');

geoshow(flipud(totLat),flipud(totLong), 'DisplayType', 'polygon','FaceColor',[195/256 216/256 1])

geoshow('worldrivers.shp', 'Color', 'k')

if contourbool
    geoshow('World_Countries.shp', 'faceColor', 'none', 'edgeColor', 'k')
end
    
PosVec3=[0 0 0 0];

h=[];
noemptyshapeind=find(~cellfun(@isempty,shapeMarker)); 
for s=noemptyshapeind
    if ~isempty(find(strcmp(shapeMarker{s}, 's')))
    h(end+1)=plot(NaN, NaN, 'Marker',shapeMarker{s},'MarkerSize', shapeSize(s), 'MarkerEdgeColor','k','MarkerFaceColor',shapeColor{s}, 'LineStyle','none');
    else
    h(end+1)=plot(NaN, NaN, 'Marker',shapeMarker{s},'MarkerSize', shapeSize(s), 'MarkerEdgeColor',shapeColor{s},'MarkerFaceColor',shapeColor{s}, 'LineStyle','none');
    end
end
if ~isempty(shapeLegend)
    legd=legend(h,shapeLegend(noemptyshapeind),'location', 'NorthWest', 'interpreter', 'none');
    set(legd,'EdgeColor',[0 0 0],'Color',[1 1 1]);
    PosVec3=get(legd, 'Position');
end

if max(PopVals)>1 & min(PopVals)>=0
    
    load(['white_red_10_cmap.mat']);
    colormap(cmap);
    
    maxVals=max(abs(PopVals));
%     maxVals=quantile(PopVals,0.95);
    PopVals_norm=PopVals./maxVals;
    PopVals_norm(find(PopVals_norm>1))=1;
    ticknames=cellstr(num2str(floor([min(PopVals):(maxVals-min(PopVals))/4:maxVals])'));
%     ticknames{end}=['>' ticknames{end}];
    
    for i=1:length(PopArray)
        if ~isnan(PopIND(i))
            if ~isempty(find(strcmp(shapeMarker{PopIND(i)}, 's')))
            geoshow(LatArray(i), LongArray(i), 'Marker',shapeMarker{PopIND(i)},'MarkerSize', shapeSize(PopIND(i)), 'MarkerEdgeColor','k',...
                'MarkerFaceColor',cmap(floor(PopVals_norm(i)*(length(cmap)-1)+1),:));%cmap(PopVals(i),:));%cmap(colorIND(reIND(i)),:));
            else
            geoshow(LatArray(i), LongArray(i), 'Marker',shapeMarker{PopIND(i)},'MarkerSize', shapeSize(PopIND(i)), 'MarkerEdgeColor',shapeColor{PopIND(i)},...
                'MarkerFaceColor',cmap(floor(PopVals_norm(i)*(length(cmap)-1)+1),:));%cmap(PopVals(i),:));%cmap(colorIND(reIND(i)),:));
            end
        end
    end
    
    cb=colorbar('location', 'eastoutside', 'color', 'k', 'Ticks',[0:0.25:1], 'TickLength', 0, 'TickLabels', ticknames); %[2600:600:5000]
    cb.Label.String=PopVals_name;
    
elseif min(PopVals)<0
    
    load(['white_red_10_cmap.mat']);
    cmap1=cmap;
    load(['white_green_10_cmap.mat']);
    cmap2=cmap(end:-1:1,:);
    cmap=vertcat(cmap2, cmap1(2:end,:));
    colormap(cmap)
    
    maxVal=max(abs(PopVals));
%     maxVal=5;
    PopVals_norm=PopVals./(maxVal);
    PopVals_norm(find(PopVals_norm>1))=1;
    PopVals_norm(find(PopVals_norm<-1))=-1;
    PopVals_norm=PopVals_norm*9+10;
    for i=1:length(PopArray)
        if ~isnan(PopIND(i))
        geoshow(LatArray(i), LongArray(i), 'Marker',shapeMarker{PopIND(i)},'MarkerSize', shapeSize(PopIND(i)), 'MarkerEdgeColor',shapeColor{PopIND(i)},...
            'MarkerFaceColor',cmap(round(PopVals_norm(i)),:));
        end
    end
    cb_ticks=0:-maxVal*2/4:-maxVal;
    cb_ticks=[cb_ticks(end:-1:2) 0:maxVal*2/4:maxVal]
    if maxVal>0.1
        cb=colorbar('location', 'eastoutside', 'color', 'k', 'Ticks',[0:0.25:1], 'TickLength', 0, 'TickLabels', num2str(cb_ticks', '%1.1f')); %[2600:600:5000]
    else
        cb=colorbar('location', 'eastoutside', 'color', 'k', 'Ticks',[0:0.25:1], 'TickLength', 0, 'TickLabels', num2str(cb_ticks','%1.0e')); %[2600:600:5000]
    end
    cb.Label.String=strrep(PopVals_name, '_','-');
elseif max(PopVals)>0
    background_inds=find(PopIND==1);
    for i=background_inds' 
        geoshow(LatArray(i), LongArray(i), 'Marker',shapeMarker{PopIND(i)},'MarkerSize', shapeSize(PopIND(i)), 'MarkerEdgeColor', shapeColor{PopIND(i)},... 
            'MarkerFaceColor',shapeColor{PopIND(i)});
    end
    
    indOI=setdiff(1:length(PopArray),background_inds);
    for i=indOI 
        if ~isempty(find(strcmp(shapeMarker{PopIND(i)}, 's')))
        geoshow(LatArray(i), LongArray(i), 'Marker',shapeMarker{PopIND(i)},'MarkerSize', shapeSize(PopIND(i)), 'MarkerEdgeColor','k',... 
            'MarkerFaceColor',shapeColor{PopIND(i)});
        else
        geoshow(LatArray(i), LongArray(i), 'Marker',shapeMarker{PopIND(i)},'MarkerSize', shapeSize(PopIND(i)), 'MarkerEdgeColor',shapeColor{PopIND(i)},... 
            'MarkerFaceColor',shapeColor{PopIND(i)});
        end
    end
end


if ~isempty(shapeLegend)
    legd=legend(h,shapeLegend(noemptyshapeind),'location', 'EastOutside', 'interpreter', 'none');
    set(legd,'EdgeColor',[0 0 0],'Color',[1 1 1]);
end


if ~isempty(titlename)
    title(titlename, 'interpreter', 'none')
    figname=[figname ' - ' titlename];
end

if ~isempty(PopVals_name)
    PopVals_name=strrep(PopVals_name,'\', '');
    figname=[figname ' - ' PopVals_name];
end

if ~isempty(contourname)
    figname=[figname contourname];
end

figtitle=[figname ' LON' num2str(lonlim(1)) '-' num2str(lonlim(2)) ' LAT'  num2str(latlim(1)) '-' num2str(latlim(2)) ' ' date]

if printbool
    currdir=pwd;
    cd(printdir)
    print(gcf,  '-painters', '-dpdf',[figtitle ' ' date '.pdf']);
    print(gcf,  '-painters', '-dpng','-r1200', [figtitle ' ' date '.png']);
    cd(currdir);
end
