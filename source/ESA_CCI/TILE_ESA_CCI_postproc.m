classdef TILE_ESA_CCI_postproc < matlab.mixin.Copyable
    %only works for yearly output

    properties
        RUN_INFO
        PARA
        STATVAR
        CONST
        TEMP

    end
    
    methods
        
        
        function tile = provide_PARA(tile)
            tile.PARA.input_file_name = [];
            tile.PARA.input_folder = [];
            tile.PARA.check_folder = [];
            tile.PARA.out_folder = [];
            

        end
        
        function tile = provide_CONST(tile)
        end
        
        function tile = provide_STATVAR(tile)

        end
        
        
        function tile = finalize_init(tile)
                      
            if ~(exist([tile.PARA.out_folder tile.RUN_INFO.PARA.run_name])==7)
                mkdir([tile.PARA.out_folder tile.RUN_INFO.PARA.run_name]);
            end
            tile.PARA.out_folder = [tile.PARA.out_folder tile.RUN_INFO.PARA.run_name '/'];
            
            %get the number of variables when called for the first time
            start_range = (tile.PARA.run_index(1,1)-1).*tile.RUN_INFO.PARA.number_of_cells_per_tile + 1;
            end_range = min(tile.PARA.run_index(1,1).*tile.RUN_INFO.PARA.number_of_cells_per_tile, tile.RUN_INFO.PARA.total_number_of_cells);
            file_name = [tile.PARA.input_folder tile.RUN_INFO.PARA.run_name '/' tile.PARA.input_file_name '_' num2str(start_range) '_' num2str(end_range) '.nc'];
            nc_info = ncinfo(file_name);
            tile.PARA.variable_names ={};
            for i=1:size(nc_info.Variables,2)
                tile.PARA.variable_names = [tile.PARA.variable_names; nc_info.Variables(i).Name];
            end
            
            tile.PARA.check_folder = [tile.PARA.check_folder tile.RUN_INFO.PARA.run_name '/'];
            if ~(exist(tile.PARA.check_folder)==7)
                mkdir(tile.PARA.check_folder)
            end
        end
        
        function tile = run_model(tile) 


            core_list = ceil([0:size(tile.PARA.variable_names,1)./tile.RUN_INFO.PARA.num_ranks:size(tile.PARA.variable_names,1)])';

            core_list = core_list(2:end,1);
            variable_start_index = core_list(tile.RUN_INFO.PARA.my_rank,1);
            permission_list = zeros(size(tile.PARA.variable_names,1),1);
            dummy = [0; core_list];
            if dummy(tile.RUN_INFO.PARA.my_rank+1,1) ~= dummy(tile.RUN_INFO.PARA.my_rank,1)
                permission_list(variable_start_index,1)=1;
%                 if tile.RUN_INFO.PARA.my_rank == tile.RUN_INFO.PARA.num_ranks
%                     permission_list(variable_start_index:size(tile.PARA.variable_names,1),1) = 1;
%                     permission_list(1:core_list(1,1)-1,1) = 1;
%                 else
%                     permission_list(variable_start_index:core_list(tile.RUN_INFO.PARA.my_rank+1,1)-1,1) = 1;
%                 end
                if tile.RUN_INFO.PARA.my_rank == 1
                    permission_list(1:variable_start_index,1) = 1;
                    permission_list(core_list(end,1)+1:end,1) = 1;
                else
                    permission_list(core_list(tile.RUN_INFO.PARA.my_rank-1,1)+1:variable_start_index,1) = 1;
                end
            end


             for i = [[variable_start_index:-1:1] [size(tile.PARA.variable_names, 1):-1:variable_start_index+1]]

                start_permission = permission_list(i,1);
                [i start_permission]

                datapackage = [];
                for j=1:size(tile.PARA.run_index,1)
                    
                    start_range = (tile.PARA.run_index(j,1)-1).*tile.RUN_INFO.PARA.number_of_cells_per_tile + 1;
                    if j==1
                        tile.PARA.start_range = start_range;
                    end
                    end_range = min(tile.PARA.run_index(j,1).*tile.RUN_INFO.PARA.number_of_cells_per_tile, tile.RUN_INFO.PARA.total_number_of_cells);

                    file_name = [tile.PARA.input_folder tile.RUN_INFO.PARA.run_name '/' tile.PARA.input_file_name '_' num2str(start_range) '_' num2str(end_range) '.nc'];
                    
                    datapackage = cat(1, datapackage, ncread(file_name, tile.PARA.variable_names{i,1}));
                end
                
                out_file_name = [tile.PARA.out_folder tile.PARA.variable_names{i,1} '.nc'];

                
                tag = ['read_' tile.PARA.variable_names{i,1} '_'];

                prev = tile.RUN_INFO.PARA.my_rank - 1;
                if prev == 0
                    prev = tile.RUN_INFO.PARA.num_ranks;
                end
               
                done = 0;
                while done == 0
                    %[tile.PARA.check_folder tag num2str(next) '.check']
                    if start_permission == 1 || exist([tile.PARA.check_folder tag num2str(prev) '.check'])==2
                        
                        if ~(exist(out_file_name)==2)
                            nccreate(out_file_name, tile.PARA.variable_names{i,1}, 'datatype', 'single', 'Dimensions', ...
                                {'x', tile.RUN_INFO.PARA.total_number_of_cells, 'y', size(datapackage,2), 'z', size(datapackage,3)}, 'FillValue', -9999);
                        end
                        
                        ncwrite(out_file_name, tile.PARA.variable_names{i,1}, datapackage, [tile.PARA.start_range 1 1]);
                        
                        done = 1;
                        start_permission = 0;
                        dlmwrite([tile.PARA.check_folder tag num2str(tile.RUN_INFO.PARA.my_rank) '.check'], []);
                    else
                        pause(0.5)
%                         if (exist([tile.PARA.check_folder tag num2str(tile.RUN_INFO.PARA.my_rank) '.check']))==2
%                             done=1;
%                         end
                    end
                    
                end
            end
            
            
            
        end

    end
end

