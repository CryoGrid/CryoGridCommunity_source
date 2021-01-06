
classdef RUN_ESA_CCI < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        TILE
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)

            run_info.PARA.run_name = [];
            
            run_info.PARA.tile_preproc_class = [];
            run_info.PARA.tile_preproc_class_index = [];

            run_info.PARA.number_of_cells_per_tile = [];
            
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            run_info.PARA.number_of_runs = [];
            
            run_info.PARA.base_projection_class = [];
            run_info.PARA.base_projection_class_index = [];
            
            run_info.PARA.mask_class = []; %generally several classes - array
            run_info.PARA.mask_class_index = [];
            
            run_info.PARA.DEM_class = [];
            run_info.PARA.DEM_class_index = [];
            
            run_info.PARA.landcover_class = [];
            run_info.PARA.landcover_class_index = [];
            
            run_info.PARA.geothermal_class = [];
            run_info.PARA.geothermal_class_index = [];
            
            
            

            
        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
        
        function run_info = finalize_init(run_info)
%             domain_size %size of computation domain (see DOMAIN class)
%             domain_info %contain start index and size of a domain, see function decompose_in_domains
%             ensemble_size %size of ensemble (not yet active)
%             
%             %needs to go to MODIS class
%             threshold_T_GlobPF_map % only pixels with temperature lower than a threshold in a GlobPF product are included
%             mask %specific to MODIS tiles; 1200x1200 array with logicals, 1 include cell, 0 not include cell

            %provided by coordinate system MODIS LST classes
            run_info.STATVAR.spatial = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.base_projection_class){run_info.PARA.base_projection_class_index,1});
            run_info.STATVAR.spatial.RUN_INFO = run_info;
            run_info.STATVAR.spatial = finalize_init(run_info.STATVAR.spatial);
            
            run_info.STATVAR.latitude = run_info.STATVAR.spatial.STATVAR.latitude;
            run_info.STATVAR.longitude = run_info.STATVAR.spatial.STATVAR.longitude;
            run_info.STATVAR.key = run_info.STATVAR.spatial.STATVAR.key;
            run_info.STATVAR.list_of_MODIS_tiles = run_info.STATVAR.spatial.STATVAR.list_of_MODIS_tiles; %use to loop over several tiles
            run_info.PARA.total_number_of_cells = size(run_info.STATVAR.key,1);
            
            run_info.PARA.dem = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.DEM_class){run_info.PARA.DEM_class_index,1});
            run_info.PARA.dem = finalize_init(run_info.PARA.dem);
            run_info.STATVAR.altitude = get_altitudes(run_info.PARA.dem, run_info);

            run_info.PARA.geothermal = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.geothermal_class){run_info.PARA.geothermal_class_index,1});
            run_info.PARA.geothermal = finalize_init(run_info.PARA.geothermal);
            run_info.STATVAR.geothermal = get_geothermal_heatflux(run_info.PARA.geothermal, run_info);
            
            run_info.PARA.landcover = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.landcover_class){run_info.PARA.landcover_class_index,1});
            run_info.PARA.landcover = finalize_init(run_info.PARA.landcover);
            run_info.STATVAR.landcover = get_landcover(run_info.PARA.landcover, run_info);

%             %provided by distributed STRATIGRAPHY class
%             run_info.PARA.stratigraphy = [repmat(1, 80, 1); repmat(2, 20, 1)]; %vector of stratigraphies
%             run_info.PARA.landcover = run_info.PARA.stratigraphy;           %vector of landcover information
%             run_info.PARA.urban = repmat(0, 100, 1);

% 
     %        run_info.PARA.k_mineral = repmat(3, 100, 1); %vector with thermal conductivity mineral fraction

% for parallel runs
%         num_ranks
%         my_rank


        end
        
        function [run_info, tile] = run_preproc(run_info)
            for i=1:size(run_info.PARA.tile_preproc_class,1)
                tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{i,1}){run_info.PARA.tile_preproc_class_index(i,1),1});
                tile.RUN_INFO = run_info;
                tile = finalize_init(tile);
                tile = run_model(tile);
            end
        end
        
        
        function [run_info, tile] = run_model(run_info)
            
            number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
            for run_index = 1:number_of_tiles
                
                disp(['running range index ' num2str(run_index)])
                
                range = [(run_index-1).*run_info.PARA.number_of_cells_per_tile+1:min(run_index.*run_info.PARA.number_of_cells_per_tile, run_info.PARA.total_number_of_cells)]';
                
                for i=1:1%size(run_info.PARA.tile_class,1)
                    disp(['running tile number ' num2str(i)])
                    for j=1:1%run_info.PARA.number_of_runs(i,1)
                        disp(['running round ' num2str(j)])
                        
                        %load the next tile from the PROVIDER
                        tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                        tile.PARA.number_of_realizations = size(range,1);
                        tile.PARA.range = range;
                        
                        tile.PARA.geothermal = run_info.STATVAR.geothermal(range,1);
                        %REMOVE
                        %tile.PARA.k_mineral = run_info.PARA.k_mineral(range,1);
                        [~, pos] = max(run_info.STATVAR.landcover(range,1),[], 2);
                        tile.PARA.stratigraphy = pos(1,1) .*0 +1;
                        %REMOVE
                        
                        %tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});
                        tile.RUN_INFO = run_info;
                        
                        tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
                       
                        %tile = run_model(tile);
                    end
                end
            end
            
            
            %do the run(s) 
            %
            %  %time integration
            
        end
 
        
        function run_info = customize(run_info)

            
        end
        
        
        
        
        
    end
end



