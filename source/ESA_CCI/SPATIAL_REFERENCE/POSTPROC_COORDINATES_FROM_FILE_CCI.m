%========================================================================
% CryoGrid SPATIAL_REFERENCE class COORDINATES_FROM_FILE
% Defines the locations of the target points as the grid coordinates of an 
% existing raster file. The region of interest can be selected by optional
% mask classes. For each point, data (e.g. calculated from DEM, or
% other sources) are read by data provder classes. Based on these data, 
% data mask classes can be used to further constrain the modeled domain. 
% NOTE: currently only fucntional for UTM coordinates  
%
% S. Westermann, Dec 2022
%========================================================================

%defines a regular grid in geographical coordinates, with fixed resolution

classdef POSTPROC_COORDINATES_FROM_FILE_CCI < COORDINATES_FROM_FILE_CCI 

    properties

        RESULT
    end
    
    methods
        
        function proj = assign_offset_scale_factor(proj)
            proj.PARA.scale_factor = [0.002;0.002;0.002; 0.002;0.002;0.002; 0.01; 1; 1; 0.01; 1; 0.002;0.002;0.002;0.002;0.002;0.002];
            proj.PARA.scale_offset = [-70;-70;-70;-70;-70;-70; 0; -10; 0; 0; 0; -70;-70;-70;-70;-70;-70];
        end
        
        
        function proj = load_spatial_reference(proj, spatial_reference_file)
            load(spatial_reference_file);
            proj.PARA = spatial_reference.PARA;
            proj.STATVAR = spatial_reference.STATVAR;
            proj.TEMP = spatial_reference.TEMP;
            proj.PARA.spatial_reference_file = spatial_reference_file;
        end
        
        
        function proj = finalize_init(proj)
            spatial_reference_file = proj.PARA.proj_file_name;
            dividers = find(proj.PARA.proj_file_name(1,:)=='_');
            proj.PARA.max_lat = str2num(spatial_reference_file(1, dividers(2)+1:dividers(3)-1));
            proj.PARA.min_lat = str2num(spatial_reference_file(1, dividers(1)+1:dividers(2)-1));
            proj.PARA.max_lon = str2num(spatial_reference_file(1, dividers(4)+1:dividers(5)-1));
            proj.PARA.min_lon = str2num(spatial_reference_file(1, dividers(3)+1:dividers(4)-1));
            lat = [proj.PARA.max_lat:-proj.PARA.delta_lat:proj.PARA.min_lat];
            lon = [proj.PARA.min_lon:proj.PARA.delta_lon:proj.PARA.max_lon];
            lon = (lon(1:end-1)+lon(2:end))/2;
            lat = (lat(1:end-1)+lat(2:end))/2;
            [lon, lat] = meshgrid(lon,lat);
            proj.RESULT.latitude = lat;
            proj.RESULT.longitude = lon;
        end
        
        function proj = assign2grid(proj, variable)
            proj.RESULT.(variable) = proj.RESULT.latitude.*NaN;
            proj.RESULT.(variable)(proj.STATVAR.key) = proj.STATVAR.(variable);
        end
        
        function proj = assign2grid_multiIndex(proj, variable)
            proj.RESULT.(variable) = [];
            for i=1:size(proj.STATVAR.(variable),2)
                data = proj.RESULT.latitude.*NaN;
                data(proj.STATVAR.key) = proj.STATVAR.(variable)(:,i);
                proj.RESULT.(variable) = cat(3, proj.RESULT.(variable), data);
            end
        end
        
        function proj = accumulate_single_variable(proj, variable)
            dummy = ncread([proj.PARA.nc_file_folder proj.PARA.nc_file_name '_1_' num2str(min(size(proj.STATVAR.latitude,1), proj.PARA.number_of_slices)) '.nc'], variable); 
            data = zeros(size(proj.STATVAR.latitude,1), size(dummy,2),size(dummy,3)) .* NaN;
            
            for i=1:proj.PARA.number_of_slices:size(data,1)
                end_index = min(size(data,1), i+proj.PARA.number_of_slices-1);
                data(i:end_index, :,:) = ncread([proj.PARA.nc_file_folder proj.PARA.nc_file_name '_' num2str(i) '_' num2str(end_index) '.nc'], variable); 
            end
            proj.STATVAR.(variable) = data;
        end
        
        function proj = accumulate_and_save_all(proj)
            info = ncinfo([proj.PARA.nc_file_folder proj.PARA.nc_file_name '_1_' num2str(min(size(proj.STATVAR.latitude,1), proj.PARA.number_of_slices)) '.nc']);
            for i=1:size(info.Variables,2)
                variable = info.Variables(i).Name;
                disp(variable)
                proj = accumulate_single_variable(proj, variable);
                data = double(proj.STATVAR.(variable)) .* proj.PARA.scale_factor(i) + proj.PARA.scale_offset(i);
                save([proj.PARA.target_folder variable '_' num2str(proj.PARA.min_lat) '_' num2str(proj.PARA.max_lat) '_' num2str(proj.PARA.min_lon) '_' num2str(proj.PARA.max_lon) '.mat'], 'data')
                proj.STATVAR.(variable) = [];
            end
        end
        
%         function proj = provide_PARA(proj)
%             proj.PARA.nc_file_folder = '../CryoGridCommunity_results/test_CCI_saga/';
%             proj.PARA.nc_file_name = 'out';
%             
%             proj.PARA.number_of_slices = 250;
%             
%         end

        
        
%         max_lat = 80;
% min_lat = 75;
% min_lon = 0;
% max_lon = 180;
% lat = [max_lat:-0.01:min_lat];
% lon = [min_lon:0.04:max_lon];
% [lon, lat] = meshgrid(lon,lat);
% lon = (lon(:,1:end-1)+lon(:, 2:end))/2;
% lat = [max_lat:-0.01:min_lat];
% lon = [min_lon:0.04:max_lon];
% lon = (lon(1:end-1)+lon(2:end))/2;
% lat = (lat(1:end-1)+lat(2:end))/2;
% [lon, lat] = meshgrid(lon,lat);
% data = lon.*NaN;
% data(spatial_reference.STATVAR.key) = spatial_reference.STATVAR.altitude;
% imagesc(data)
% data = lon.*NaN;
% data(spatial_reference.STATVAR.key) = spatial_reference.STATVAR.geothermal;
% imagesc(data)
        
%         function proj = provide_PARA(proj)
%             proj.PARA.proj_file_folder = [];
%             proj.PARA.proj_file_name = [];
%             
%             proj.PARA.mask_class = []; %acts on the entire 2d matirx
%             proj.PARA.mask_class_index = [];
%             
%             proj.PARA.data_class = [];
%             proj.PARA.data_class_index = [];
%             
%             proj.PARA.data_mask_class = []; %
%             proj.PARA.data_mask_class_index = [];
%             
%             proj.PARA.assign_tile_properties_class = [];
%             proj.PARA.assign_tile_properties_class_index = [];
%         end
        
        
        
%         function proj = finalize_init(proj)
%             
%             proj.STATVAR.latitude = ncread([proj.PARA.proj_file_folder proj.PARA.proj_file_name], 'latitude');
%             proj.STATVAR.longitude = ncread([proj.PARA.proj_file_folder proj.PARA.proj_file_name], 'longitude');
%             proj.STATVAR.key = ncread([proj.PARA.proj_file_folder proj.PARA.proj_file_name], 'key');
%             proj.PARA.delta_lat = ncread([proj.PARA.proj_file_folder proj.PARA.proj_file_name], 'delta_lat');
%             proj.PARA.delta_lon = ncread([proj.PARA.proj_file_folder proj.PARA.proj_file_name], 'delta_lon');
%             proj.STATVAR.key_internal = [1:size(proj.STATVAR.latitude,1)]';
% 
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
%             %load data sets
%             for i=1:size(proj.PARA.data_class,1)
%                 data_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_class{i,1}){proj.PARA.data_class_index(i,1),1});
%                 data_class.PARENT = proj;
%                 data_class = finalize_init(data_class);
%                 data_class = load_data(data_class); %can be additive or subtractive
%             end
%             
%             %remove all cells with NaN values
%             fn = fieldnames(proj.STATVAR);
%             for i=1:size(fn,1)
%                     proj.STATVAR.mask(logical(sum(isnan(proj.STATVAR.(fn{i,1})),2))) = 0;
%             end
%             mask = proj.STATVAR.mask;
%             fn = fieldnames(proj.STATVAR);
%             for i=1:size(fn,1)  
%                     proj.STATVAR.(fn{i,1})(~mask(:,1), :) = [];
%             end
%             
%             %apply data masks
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
% 
%             for i=1:size(proj.PARA.assign_tile_properties_class,1)
%                 proj.ACTION{i,1} = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.assign_tile_properties_class{i,1}){proj.PARA.assign_tile_properties_class_index(i,1),1});
%                 proj.ACTION{i,1} = finalize_init(proj.ACTION{i,1});
%                 proj.ACTION{i,1}.PROJ = proj;
%             end

%         end
        
        
        
%         %-------------param file generation-----
%         function proj = param_file_info(proj)
%             proj = provide_PARA(proj);
%             
%             proj.PARA.STATVAR = [];
%             proj.PARA.class_category = 'SPATIAL_REFERENCE';
%             proj.PARA.default_value = [];
%             
%             proj.PARA.comment.proj_file_folder = {'folder where data set providing the spatial reference is located'};
% 
%             proj.PARA.comment.proj_file_name = {'file name of data set providing the spatial reference'};
%             
%             proj.PARA.comment.mask_class = {'list of mask classes, constains the region of interest based on the coordinates'};
%             proj.PARA.options.mask_class.name = 'H_LIST';
% %             proj.PARA.options.mask_class.entries_x = {};
%             proj.PARA.options.mask_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.data_class = {'list of data provider classes, provide data for each target location'};
%             proj.PARA.options.data_class.name = 'H_LIST';
% %             proj.PARA.options.data_class.entries_x = {''};
%             proj.PARA.options.data_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.data_mask_class = {'list of data mask classes, constrains the region of interest based on the data provided by data provider classes'};
%             proj.PARA.options.data_mask_class.name = 'H_LIST';
% %             proj.PARA.options.data_mask_class.entries_x = {''};
%             proj.PARA.options.data_mask_class_index.name = 'H_LIST';
% 
%             proj.PARA.comment.assign_tile_properties_class = {'translates the data to changes to the classes used in the simulations for each point, basically cutomizing the "parameter file" for each taret location'};
%             proj.PARA.options.assign_tile_properties_class.name = 'H_LIST';
%             proj.PARA.options.assign_tile_properties_class.entries_x = {'update_one2one' 'tag_out_w_run_number'};   
%             
%             proj.PARA.options.assign_tile_properties_class_index.name = 'H_LIST';
%             proj.PARA.options.assign_tile_properties_class_index.entries_x = {'1' '1'};   
%         end
  
    end
end

