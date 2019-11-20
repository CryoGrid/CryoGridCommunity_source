%% Construct control sheet
clear all

addpath(genpath('modules'))

filename = 'test_crocus';
cols = 0;
rows = 0;

out_class = {'OUT_parallel'}; % alternatively 'OUT_all'
lateral_class = {'LATERAL_snow'};
structural_classes = {'FORCING_seb'; 'GRID_user_defined'; 'STRAT_layers'; 'STRAT_classes'; 'STRAT_linear'};
ground_classes = {'GROUND_freeW_bucketW_seb_snow';
                'GROUND_freeW_seb';
                'SNOW_crocus_no_inheritance'};
classlist = [out_class; lateral_class; structural_classes; ground_classes];

for i = 1:size(classlist,1)
    temp_class = eval(classlist{i,1});
    temp{1,i} = write_excel(temp_class); 
    cols = [cols, size(temp{1,i},2)];
    rows = [rows, size(temp{1,i},1)];
end

empty = cell(1,max(cols(:)));
sheet = cell(1,max(cols(:)));
 
for i = 1:size(classlist,1)
    sheet = [ sheet;  empty; [temp{1,i}, cell(rows(i+1),max(cols(:))-cols(i+1))] ];
end

for i = 1:size(classlist,1)
    if strcmp(classlist(i,1),'STRAT_classes') == 1
        pos = sum(rows(1:i)) +i +2;
    end
end

sheet{pos+5,2} = ground_classes{1,1};
sheet{pos+6,2} = ground_classes{2,1};
sheet{pos+9,2} = ground_classes{end,1};

xlswrite([ filename '.xlsx'],[sheet])
