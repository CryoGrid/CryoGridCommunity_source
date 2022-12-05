%defines a regular grid in geographical coordinates, with fixed resolution

classdef LAT_LON < matlab.mixin.Copyable

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
            proj.PARA.max_lat = [];
            proj.PARA.min_lat = [];
            proj.PARA.lat_grid_cell_size = [];
            proj.PARA.max_lon = [];
            proj.PARA.min_lon = [];
            proj.PARA.lon_grid_cell_size = [];
            
            proj.PARA.mask_class = []; %acts on the entire 2d matirx
            proj.PARA.mask_class_index = [];
            
            proj.PARA.data_class = [];
            proj.PARA.data_class_index = [];
            
            proj.PARA.data_mask_class = []; %
            proj.PARA.data_mask_class_index = [];
            
            proj.PARA.assign_tile_properties_class = [];
            proj.PARA.assign_tile_properties_class_index = [];
        end
        
        function proj = provide_STATVAR(proj)

        end
        
        function proj = provide_CONST(proj)
            
        end
        
        function proj = finalize_init(proj)
            %span up grid and initialize matrices for lat and lon
            proj.STATVAR.latitude = [proj.PARA.min_lat:proj.PARA.lat_grid_cell_size:proj.PARA.max_lat]';
            proj.STATVAR.longitude = [proj.PARA.min_lon:proj.PARA.lon_grid_cell_size:proj.PARA.max_lon]';
            [proj.STATVAR.longitude, proj.STATVAR.latitude] = meshgrid(proj.STATVAR.longitude, proj.STATVAR.latitude);
            
            %reduce to list 
            proj.STATVAR.latitude = proj.STATVAR.latitude(:);
            proj.STATVAR.longitude = proj.STATVAR.longitude(:);
            proj.STATVAR.key = [1:size(proj.STATVAR.latitude,1)]';
            
            %apply masks before data sets
            proj.STATVAR.mask = logical(proj.STATVAR.longitude.*1);
            for i=1:size(proj.PARA.mask_class_index,1)
                mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.mask_class{i,1}){proj.PARA.mask_class_index(i,1),1});
                mask_class.PARENT = proj;
                mask_class = finalize_init(mask_class);
                mask_class = apply_mask(mask_class); %can be additive or subtractive
            end
            
            %reduce the list to the ones inside the masks
            mask = proj.STATVAR.mask;
            fn = fieldnames(proj.STATVAR);
            for i=1:size(fn,1)  
                    proj.STATVAR.(fn{i,1})(~mask) = [];
            end
            

            %load data sets
            for i=1:size(proj.PARA.data_class,1)
                data_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_class{i,1}){proj.PARA.data_class_index(i,1),1});
                data_class.PARENT = proj;
                data_class = finalize_init(data_class);
                data_class = load_data(data_class); %can be additive or subtractive
            end
            
            for i=1:size(proj.PARA.data_mask_class_index,1)
                mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_mask_class{i,1}){proj.PARA.data_mask_class_index(i,1),1});
                mask_class.PARENT = proj;
                mask_class = finalize_init(mask_class);
                mask_class = apply_mask(mask_class); %can be additive or subtractive
            end
            
            mask = proj.STATVAR.mask;
            fn = fieldnames(proj.STATVAR);
            for i=1:size(fn,1)  
                    proj.STATVAR.(fn{i,1})(~mask) = [];
            end
            
            for i=1:size(proj.PARA.assign_tile_properties_class,1)
                proj.ACTION{i,1} = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.assign_tile_properties_class{i,1}){proj.PARA.assign_tile_properties_class_index(i,1),1});
                proj.ACTION{i,1} = finalize_init(proj.ACTION{i,1});
                proj.ACTION{i,1}.PROJ = proj;
            end

        end
        
 
        
    end
end

