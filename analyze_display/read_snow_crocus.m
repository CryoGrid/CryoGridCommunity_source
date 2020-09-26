%% Read out and display snow variables of interest
% R. Zweigel Sept. 2019

% Load your file

clear result

depth = 1; % depth of soil to be assessed
height = 3.5; % height above ground to be assessed
spacing = 0.01; % Vertical Grid resolution

out = self;

alt = out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.upperPos;

grid = [alt-depth:spacing:alt+height]';

result.T = [];
result.waterIce = [];
result.water = [];
result.d = [];
result.gs = [];
result.s = [];

variableList = fieldnames(result);
numberOfVariables = size(variableList,1);


%%
for i = 1:size(out.STRATIGRAPHY,2)
    altitudeLowestCell = out.STRATIGRAPHY{1,i}{end,1}.STATVAR.lowerPos;
    layerThick = [];
    area=[];
    
    for j = 1:size(out.STRATIGRAPHY{1,i},1)
        layerThick = [layerThick; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick];
        area=[area; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.area];

    end
    depths = [0; cumsum(layerThick)];
    depths = -(depths-depths(end,1));
    depths = (depths(1:end-1,1)+depths(2:end,1))./2 + altitudeLowestCell;
    
    temp = NaN(size(layerThick,1), numberOfVariables);
    pos = 1;
    for j = 1:size(out.STRATIGRAPHY{1,i},1)
        fieldLength = size(out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick,1);
        for k = 1:numberOfVariables
            if any(strcmp(fieldnames(out.STRATIGRAPHY{1,i}{j,1}.STATVAR), variableList{k,1}))
                temp(pos:pos+fieldLength-1,k) = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.(variableList{k,1});
            end
        end
        pos = pos+fieldLength;
    end
    
    
    for k=1:numberOfVariables
        if strcmp(variableList{k,1}, 'water') || strcmp(variableList{k,1}, 'waterIce')
            temp(:,k) = temp(:,k)./layerThick ./ area;
        end
        result.(variableList{k,1}) = [result.(variableList{k,1}) interp1(depths, temp(:,k), grid)];
    end
    
end

result.waterIce = result.waterIce.*1000;
result.gs = result.gs.*1000;

%% PLOTTING

for i=1:numberOfVariables
    figure
    imagesc(out.TIMESTAMP, grid, result.(variableList{i,1}))
    axis xy
    datetick
    if strcmp(variableList{i,1}, 'T')
        hold on
        contour(out.TIMESTAMP, grid, result.(variableList{i,1}),[0 0],'k')
    end
    if strcmp(variableList{i,1}, 'waterIce')
        title('Density')
    else
        title(variableList{i,1})
    end
    colorbar
    ylabel('m.a.s.l.')
end