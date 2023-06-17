%========================================================================
%works together with the RUN_INFO..._NMPI class
%
% S. Westermann, Dec 2022
%========================================================================

classdef ESA_CCI_global_run_MODIS_deg < matlab.mixin.Copyable

    properties
        PARA
        CONST
    end
    
    methods
        function update = provide_PARA(update)

            update.PARA.deg_tile_list_file = [];
            update.PARA.deg_tile_list_folder = [];
            update.PARA.MODIS_tile_number = []; %if empty select automatically depending on previously simulated tiles
        end
        
        function update = provide_STATVAR(update)

        end
        
        function update = provide_CONST(update)
            
        end
        
        function update = finalize_init(update)

        end
        
        function update = assign_run_info_properties(update, run_info)
            load([update.PARA.deg_tile_list_folder update.PARA.deg_tile_list_file]);
            if isempty(update.PARA.MODIS_tile_number) || isnan(update.PARA.MODIS_tile_number)
                still2do = find(progress_list(:,1) == 0);
                rng('shuffle') %different random number sequence each time
                pos = max(1, min(length(still2do), ceil(rand(1) .* length(still2do))));
                update.PARA.MODIS_tile_number = progress_list(still2do(pos));
            end
            tile_str = ['_' num2str(MODIS_deg_list(update.PARA.MODIS_tile_number,1)) '_' num2str(MODIS_deg_list(update.PARA.MODIS_tile_number,2)) ...
                '_' num2str(MODIS_deg_list(update.PARA.MODIS_tile_number,3)) '_' num2str(MODIS_deg_list(update.PARA.MODIS_tile_number,4))];
            run_info.PARA.run_name = [run_info.PARA.run_name tile_str];
            run_info.PPROVIDER.CLASSES.COORDINATES_FROM_FILE_CCI{1,1}.PARA.proj_file_name = ['MODIS' tile_str '_2019.nc'];
            run_info.PPROVIDER.CLASSES.merge_MODIS_ERA{1,1}.PARA.MODIS_file = ['MODIS' tile_str];
            if run_info.PARA.worker_number ==  1
                progress_list(update.PARA.MODIS_tile_number,1) = 1;
                save([update.PARA.deg_tile_list_folder update.PARA.deg_tile_list_file], 'MODIS_deg_list', 'progress_list');
            end
        end
        
        function update = finalize_progress_list(update, run_info)
            load([update.PARA.deg_tile_list_folder update.PARA.deg_tile_list_file]);
            if run_info.PARA.worker_number ==  1
                progress_list(update.PARA.MODIS_tile_number,2) = progress_list(update.PARA.MODIS_tile_number,2) + 1;
                if progress_list(update.PARA.MODIS_tile_number,2) == run_info.PARA.number_of_cores
                    progress_list(update.PARA.MODIS_tile_number,2) = -1;
                end
                    
                save([update.PARA.deg_tile_list_folder update.PARA.deg_tile_list_file], 'MODIS_deg_list', 'progress_list');
            end
        end
        
        
        
%         %-------------param file generation-----
%         function update = param_file_info(update)
%             update = provide_PARA(update);
% 
%             update.PARA.STATVAR = [];
%             update.PARA.class_category = 'ASSIGN_TILE_PROPERTIES';
%             update.PARA.default_value = [];
%             
%             update.PARA.comment.class_name =  {'name of class in which the variable is changed'};
%             update.PARA.options.class_name.name =  'H_LIST';
%             update.PARA.options.class_name.entries_x = {'TILE_1D_standard' 'TILE_1D_standard' 'TILE_1D_standard'};
%             
%             update.PARA.comment.class_index =  {'index of class in which the variable is changed'};
%             update.PARA.options.class_index.name =  'H_LIST';
%             update.PARA.options.class_index.entries_x = {'1' '1' '1'};
% 
%             update.PARA.comment.variable =  {'name of variable which is changed (must be assigned by DATA PROVIDER class)'};
%             update.PARA.options.variable.name =  'H_LIST';
%             update.PARA.options.variable.entries_x = {'latitude' 'longitude' 'altitude'};
% 
%         end
        
    end
end

