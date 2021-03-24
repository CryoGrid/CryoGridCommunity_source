number_of_cores = 4;
my_core = 4;
%file_folder = [forcing.PARA.forcing_ground_folder tile.RUN_INFO.PARA.run_name '/'];
file_folder = 'H:\projects\GlobPermafrost\CCI_Permafrost\test_data_newCryoGrid\test\test_Broegger/';
run_name = 'test_Broegger';

variables = {'ERA_melt_bare'; 'ERA_melt_forest'; 'ERA_snowfall_downscaled'; 'ERA_T_downscaled'; 'final_av_T'; 'final_MODIS_weight'};
for i=1:size(variables,1)  %make this only for a certain part
    %file_info = dir([forcing.PARA.forcing_ground_folder tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variables{i,1} '*.nc']);
    file_info = dir([file_folder run_name '_' variables{i,1}  '*.nc'])
    file_info = struct2cell(file_info);
    file_info = file_info(1,:)';
    file_info = cell2mat(file_info);
    start_time = file_info(:,end-10:end-3); %file_info(:,end-19:end-12);
    start_time = datenum(start_time, 'yyyymmdd');
    [start_time, pos]= sort(start_time, 'ascend');
    
    forcing.DATA.(variables{i,1}) = [];
    forcing.DATA.timestamp = [];
    
    segment_size = ceil(size(file_info,1)./number_of_cores);
    
    %for j=1:size(file_info,1) %make this dependent
    for j=(my_core-1).*segment_size+1:min(my_core.*segment_size, size(file_info,1))
        filename = file_info(pos(j),:);
        
        data = ncread([file_folder filename], variables{i,1}, [1 1], [Inf Inf], [1 1]);
        forcing.DATA.(variables{i,1}) = [ forcing.DATA.(variables{i,1}) data];
        %make the 8 dependent on the forcing
        forcing.DATA.timestamp = [forcing.DATA.timestamp start_time(pos(j),1) + 3.5 + [0:8:8.*(size(data,2)-1)]];
    end
    
    
    
    %split the files in
    
    for j=1:number_of_cores
        
        segment_size = ceil(size(forcing.DATA.(variables{i,1}),1)./number_of_cores);
        
        start_pos = (j-1).*segment_size+1
        end_pos = min(j.*segment_size, size(forcing.DATA.(variables{i,1}),1));
        
        nccreate([file_folder 'forcing_intermediate' num2str(my_core) '_' num2str(j) '.nc'], variables{i,1}, 'datatype', 'single', 'Dimensions', ...
            {'x',end_pos-start_pos+1,'y', size(forcing.DATA.(variables{i,1}),2)}, 'FillValue', -9999);
        
        ncwrite([file_folder 'forcing_intermediate' num2str(my_core) '_' num2str(j) '.nc'], ...
            variables{i,1}, forcing.DATA.(variables{i,1})(start_pos:end_pos,:), [1 1], [1 1]);
        
        if i==1

             nccreate([file_folder 'forcing_intermediate' num2str(my_core) '_' num2str(j) '.nc'], 'timestamp', 'datatype', 'single', 'Dimensions', ...
                 {'y', size(forcing.DATA.timestamp,2)}, 'FillValue', -9999);
            
            ncwrite([file_folder 'forcing_intermediate' num2str(my_core) '_' num2str(j) '.nc'], ...
                'timestamp', forcing.DATA.timestamp, [1], [1]);
        end
    end
    
    forcing.DATA.(variables{i,1}) = [];

end

labBarrier()
% 
variables = [variables; 'timestamp'];

% 
for i=1:size(variables,1)-1
    forcing.DATA.(variables{i,1}) = [];
    for j=1:number_of_cores
          data = ncread([file_folder 'forcing_intermediate' num2str(j) '_' num2str(my_core) '.nc'], variables{i,1}, [1 1], [Inf Inf], [1 1]);
          forcing.DATA.(variables{i,1}) = [forcing.DATA.(variables{i,1}) data];
    end
    
    
    %change filename
    nccreate([file_folder run_name 'forcing_' num2str(my_core) '.nc'], variables{i,1}, 'datatype', 'single', 'Dimensions', ...
        {'x',size(forcing.DATA.(variables{i,1}),1),'y',size(forcing.DATA.(variables{i,1}),2)}, 'FillValue', -9999);
    
    ncwrite([file_folder run_name 'forcing_' num2str(my_core) '.nc'], ...
        variables{i,1}, forcing.DATA.(variables{i,1}), [1 1], [1 1]);

    
end

%timestamp
i=size(variables,1);
forcing.DATA.(variables{i,1}) = [];
for j=1:number_of_cores
    data = ncread([file_folder run_name 'forcing_intermediate' num2str(j) '_' num2str(my_core) '.nc'], variables{i,1}, [1], [Inf], [ 1]);
    forcing.DATA.(variables{i,1}) = [forcing.DATA.(variables{i,1}) data'];
end


%change filename
nccreate([file_folder run_name 'forcing_' num2str(my_core) '.nc'], variables{i,1}, 'datatype', 'single', 'Dimensions', ...
    {'y',size(forcing.DATA.(variables{i,1}),2)}, 'FillValue', -9999);

ncwrite([file_folder run_name 'forcing_' num2str(my_core) '.nc'], ...
    variables{i,1}, forcing.DATA.(variables{i,1}), [1], [1]);


