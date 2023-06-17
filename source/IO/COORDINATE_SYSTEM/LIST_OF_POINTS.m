%========================================================================
% CryoGrid SPATIAL_REFERENCE class LIST_OF_POINTS
% S. Westermann, Dec 2022
%========================================================================


classdef LIST_OF_POINTS < matlab.mixin.Copyable

    properties
        RUN_INFO
        PARA
        CONST
        STATVAR
        TEMP
        ACTION
    end
    
    methods
        function proj = provide_PARA(proj)
            proj.PARA.info_file_folder = []; %first row contains variable names
            proj.PARA.info_file_name = []; 
            
            proj.PARA.assign_tile_properties_class = [];
            proj.PARA.assign_tile_properties_class_index = [];
        end
        
        function proj = provide_STATVAR(proj)

        end
        
        function proj = provide_CONST(proj)
            
        end
        
        function proj = finalize_init(proj)
            
            [~,~,raw] = xlsread([proj.PARA.info_file_folder proj.PARA.info_file_name]);
            
            for i=1:size(raw,2)
                if strcmp(raw{2,i}, 'num')
                    proj.STATVAR.(raw{1,i}) = cell2mat(raw(3:end,i));
                else
                    proj.STATVAR.(raw{1,i}) = raw(3:end,i);
                end
            end
            
            variables = fieldnames(proj.STATVAR);
            lat_yes=0;
            lon_yes=0;
            alt_yes=0;
            for i=1:size(variables,1)
                if strcmp(variables{i,1}, 'latitude')
                    lat_yes = 1;
                end
                if strcmp(variables{i,1}, 'longitude')
                    lat_yes = 1;
                end
                if strcmp(variables{i,1}, 'altitude')
                    lat_yes = 1;
                end
            end
            if ~lat_yes
                proj.STATVAR.latitude = repmat(70, size(proj.STATVAR.(variables{1,1}),1),1);
            end
            if ~lon_yes
                proj.STATVAR.longitude = repmat(70, size(proj.STATVAR.(variables{1,1}),1),1);
            end
            if ~alt_yes
                proj.STATVAR.altitde = repmat(70, size(proj.STATVAR.(variables{1,1}),1),1);
            end
            proj.STATVAR.key = [1:size(proj.STATVAR.latitude,1)]';
            
%             %apply masks before data sets
%             proj.STATVAR.mask = logical(proj.STATVAR.longitude.*1);
%             for i=1:size(proj.PARA.mask_class_index,1)
%                 mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.mask_class{i,1}){proj.PARA.mask_class_index(i,1),1});
%                 mask_class.PARENT = proj;
%                 mask_class = finalize_init(mask_class);
%                 mask_class = apply_mask(mask_class); %can be additive or subtractive
%             end
%             
%             %reduce the list to the ones inside the masks
%             mask = proj.STATVAR.mask;
%             fn = fieldnames(proj.STATVAR);
%             for i=1:size(fn,1)  
%                     proj.STATVAR.(fn{i,1})(~mask) = [];
%             end
%             
% 
%             %load data sets
%             for i=1:size(proj.PARA.data_class,1)
%                 data_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_class{i,1}){proj.PARA.data_class_index(i,1),1});
%                 data_class.PARENT = proj;
%                 data_class = finalize_init(data_class);
%                 data_class = load_data(data_class); %can be additive or subtractive
%             end
%             
%             for i=1:size(proj.PARA.data_mask_class_index,1)
%                 mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_mask_class{i,1}){proj.PARA.data_mask_class_index(i,1),1});
%                 mask_class.PARENT = proj;
%                 mask_class = finalize_init(mask_class);
%                 mask_class = apply_mask(mask_class); %can be additive or subtractive
%             end
%             
%             mask = proj.STATVAR.mask;
%             fn = fieldnames(proj.STATVAR);
%             for i=1:size(fn,1)  
%                     proj.STATVAR.(fn{i,1})(~mask) = [];
%             end
            
            for i=1:size(proj.PARA.assign_tile_properties_class,1)
                proj.ACTION{i,1} = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.assign_tile_properties_class{i,1}){proj.PARA.assign_tile_properties_class_index(i,1),1});
                proj.ACTION{i,1} = finalize_init(proj.ACTION{i,1});
                proj.ACTION{i,1}.PROJ = proj;
            end

        end
        
 
%         %-------------param file generation-----
%         function proj = param_file_info(proj)
%             proj = provide_PARA(proj);
%             
%             proj.PARA.STATVAR = [];
%             proj.PARA.class_category = 'SPATIAL_REFERENCE';
%             proj.PARA.default_value = [];
%             
%             proj.PARA.comment.max_lat = {'maximum latitude of model domain'}; 
%             
%             proj.PARA.comment.min_lat = {'minimum latitude of model domain'};
%             
%             proj.PARA.comment.lat_grid_cell_size = {'latitudinal spacing of grid'};
%             
%             proj.PARA.comment.max_lon = {'maximum longitude of model domain'};
%             
%             proj.PARA.comment.min_lon = {'minimum longitude of model domain'};
%             
%             proj.PARA.comment.lon_grid_cell_size = {'longitudinal spacing of grid'};
%             
%             proj.PARA.comment.mask_class = {'list of mask classes, constains the region of interest based on the coordinates'};
%             proj.PARA.options.mask_class.name = 'H_LIST';
%             proj.PARA.options.mask_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.data_class = {'list of data provider classes, provide data for each target location'};
%             proj.PARA.options.data_class.name = 'H_LIST';
%             proj.PARA.options.data_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.data_mask_class = {'list of data mask classes, constrains the region of interest based on the data provided by data provider classes'};
%             proj.PARA.options.data_mask_class.name = 'H_LIST';
%             proj.PARA.options.data_mask_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.assign_tile_properties_class = {'translates the data to changes to the classes used in the simulations for each point, basically cutomizing the "parameter file" for each taret location'};
%             proj.PARA.options.assign_tile_properties_class.name = 'H_LIST';
%             proj.PARA.options.assign_tile_properties_class.entries_x = {'update_one2one' 'tag_out_w_run_number'};            
%         end
        
    end
end

