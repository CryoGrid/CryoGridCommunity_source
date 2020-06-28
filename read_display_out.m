%function result = read_display_out()

%interpolates the given target variables to the variable "new_grid" and
%displays them as color plots

%define the grid to which variables are interpolated
target_cell_size = 0.02;
new_grid = [395:target_cell_size:402]';
threshold = 0.5;  %cell displayed as bulk value, no interpolation

%define target variables
result.T = [];
result.waterIce =[];
result.water = [];
result.ice = [];
result.Xice = [];
result.Xwater = [];
result.XwaterIce = [];
result.waterPotential = [];
%result.saltConc = [];
%result.salt_c_brine=[];
result.class_number = [];
%must all have the same dimensions in all fields

variableList = fieldnames(result);
numberOfVariables = size(variableList,1);

%out=self;

for i=1:size(out.STRATIGRAPHY,2)
    
    altitudeLowestCell = out.STRATIGRAPHY{1,i}{end,1}.STATVAR.lowerPos;
    
    
    %read out and accumulate over all classes
    
    layerThick=[];
    area=[];
    
    for j=1:size(out.STRATIGRAPHY{1,i},1)
        layerThick=[layerThick; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick];
        area=[area; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.area];
    end
    
    %     depths = [0; cumsum(layerThick)];
    %     depths = -(depths-depths(end,1));
    %     depths = (depths(1:end-1,1)+depths(2:end,1))./2 + altitudeLowestCell;
    
    
    temp=repmat(NaN, size(layerThick,1), numberOfVariables);
    pos=1;
    for j = 1:size(out.STRATIGRAPHY{1,i},1)
        fieldLength = size(out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick,1);
        for k=1:numberOfVariables-1
            if any(strcmp(fieldnames(out.STRATIGRAPHY{1,i}{j,1}.STATVAR), variableList{k,1}))
                temp(pos:pos+fieldLength-1,k) = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.(variableList{k,1});
            end
        end
        temp(pos:pos+fieldLength-1,numberOfVariables) = zeros(fieldLength,1) + size(out.STRATIGRAPHY{1,i},1)+1-j; %assigna class number starting with 1 from the bottom
        pos = pos+fieldLength;
    end
    
    %compute targate variables
    for k=1:numberOfVariables
        if strcmp(variableList{k,1}, 'saltConc')
            pos_waterIce = find(strcmp(variableList, 'waterIce'));
            temp(:,k) = temp(:,k)./ (temp(:,pos_waterIce) ./layerThick./area);  %divide by total water content
        end
    end
    
    for k=1:numberOfVariables
        if strcmp(variableList{k,1}, 'water') || strcmp(variableList{k,1}, 'ice') || strcmp(variableList{k,1}, 'waterIce') || strcmp(variableList{k,1}, 'XwaterIce') || strcmp(variableList{k,1}, 'Xwater') || strcmp(variableList{k,1}, 'Xice') || strcmp(variableList{k,1}, 'saltConc')
            temp(:,k) = temp(:,k)./layerThick./area;
        end
        %result.(variableList{k,1}) = [result.(variableList{k,1}) interp1(depths, temp(:,k), new_grid)];
    end
    
    %duplicate large "bulk cells"
    layerThick_dummy = layerThick;
    for k=size(layerThick_dummy,1):-1:1
        if layerThick_dummy(k,1) > threshold
            layerThick =[layerThick(1:k-1,1); target_cell_size/2; layerThick(k,1)-target_cell_size; target_cell_size/2; layerThick(k+1:end,1)];
            area = [area(1:k-1,:); area(k,:); area(k,:); area(k,:); area(k+1:end,:)];
            temp = [temp(1:k-1,:); temp(k,:); temp(k,:); temp(k,:); temp(k+1:end,:)];
        end
    end
    
    %add top and bottom half cell
    %layerThick = [layerThick(1,1)/2; layerThick];
    area = [area(1,1); area];
    temp = [temp(1,:); temp];
    
%     layerThick = [ layerThick; layerThick(end,1)/2];
%     area = [area; area(end,1)];
%     temp = [temp; temp(end,:)];
    
    
    
    %interpolate to new grid
    depths = [0; cumsum(layerThick)];
    depths = -(depths-depths(end,1));
    %depths = (depths(1:end-1,1)+depths(2:end,1))./2 + altitudeLowestCell;
    depths = depths + altitudeLowestCell;
    
    for k=1:numberOfVariables
        result.(variableList{k,1}) = [result.(variableList{k,1}) interp1(depths, temp(:,k), new_grid)];
    end
    
end

%plot
for i=1:numberOfVariables
    figure
    imagesc(out.TIMESTAMP, new_grid, result.(variableList{i,1}))
    axis xy
    datetick
    title(variableList{i,1})
    colorbar
end
