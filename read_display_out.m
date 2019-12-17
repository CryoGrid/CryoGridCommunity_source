%function result = read_display_out()

%interpolates the given target variables to the variable "new_grid" and
%displays them as color plots

%define the grid to which variables are interpolated
new_grid = [395:0.02:400]';

%define target variables
result.T = [];
result.waterIce =[];
result.water = [];
result.ice = [];
result.saltConc = [];
result.salt_c_brine=[];
%must all have the same dimensions in all fields

variableList = fieldnames(result);
numberOfVariables = size(variableList,1); 


for i=1:size(out.STRATIGRAPHY,2)

    altitudeLowestCell = out.STRATIGRAPHY{1,i}{end,1}.STATVAR.lowerPos;
    
    
    %accumulate over all classes

    layerThick=[];
    
    for j=1:size(out.STRATIGRAPHY{1,i},1)
        layerThick=[layerThick; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick];
    end
    
    depths = [0; cumsum(layerThick)];
    depths = -(depths-depths(end,1));
    depths = (depths(1:end-1,1)+depths(2:end,1))./2 + altitudeLowestCell;

    
    temp=repmat(NaN, size(layerThick,1), numberOfVariables); 
    pos=1;
    for j = 1:size(out.STRATIGRAPHY{1,i},1)
         fieldLength = size(out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick,1);
         for k=1:numberOfVariables
             if any(strcmp(fieldnames(out.STRATIGRAPHY{1,i}{j,1}.STATVAR), variableList{k,1}))
                temp(pos:pos+fieldLength-1,k) = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.(variableList{k,1});
             end
         end
         pos = pos+fieldLength;
     end
     
     for k=1:numberOfVariables
         if strcmp(variableList{k,1}, 'saltConc')
             pos_waterIce = find(strcmp(variableList, 'waterIce'));
             temp(:,k) = temp(:,k)./ temp(:,pos_waterIce);  %divide by total water content
         end
     end
     for k=1:numberOfVariables
         if strcmp(variableList{k,1}, 'water') || strcmp(variableList{k,1}, 'ice') || strcmp(variableList{k,1}, 'waterIce')
             temp(:,k) = temp(:,k)./layerThick;
         end
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
