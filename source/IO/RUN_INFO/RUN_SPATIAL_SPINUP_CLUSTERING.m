
classdef RUN_SPATIAL_SPINUP_CLUSTERING < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        SPATIAL
        CLUSTER
        TILE
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)
            
            run_info.PARA.number_of_cores = [];
            
            %run_info.PARA.number_of_cells_per_tile = [];
            
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            run_info.PARA.number_of_runs_per_tile = []; %vector
            
            run_info.PARA.projection_class = [];
            run_info.PARA.projection_class_index = [];
            
            run_info.PARA.clustering_class = [];
            run_info.PARA.clustering_class_index = [];
        end
        
        function run_info = provide_CONST(run_info)
            
        end
        
        function run_info = provide_STATVAR(run_info)
            
        end
        
        
        function run_info = finalize_init(run_info)
            
            disp('get spatial data')
            %provided by coordinate system MODIS LST classes
            run_info.SPATIAL = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.projection_class){run_info.PARA.projection_class_index,1});
            run_info.SPATIAL.RUN_INFO = run_info;
            run_info.SPATIAL = finalize_init(run_info.SPATIAL);
            
            run_info.CLUSTER = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.clustering_class){run_info.PARA.clustering_class_index,1});
            run_info.CLUSTER.RUN_INFO = run_info;
            run_info.CLUSTER = finalize_init(run_info.CLUSTER);
            run_info.CLUSTER = compute_clusters(run_info.CLUSTER);
        end
        
        
        function [run_info, tile] = run_model(run_info)
            
            
            tile = 0;
            if run_info.PARA.number_of_cores > 1 %parallelized
                parpool(run_info.PARA.number_of_cores)
                spmd
                    worker_number = labindex;
                    for sample_number = worker_number:run_info.PARA.number_of_cores:size(run_info.CLUSTER.STATVAR.sample_centroid_index,1)
                        
                        run_number = run_info.CLUSTER.STATVAR.sample_centroid_index(sample_number,1);
                        for ai=1:size(run_info.SPATIAL.ACTION,1)
                            run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class
                        end
                        
                        disp(['running grid cell ' num2str(run_number)])
                        %as normal 1D run
                        for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
                            disp(['running tile number ' num2str(i)])
                            for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                                disp(['running round ' num2str(j)])
                                
                                new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                                new_tile.RUN_INFO = run_info;
                                new_tile = finalize_init(new_tile);
                                tile = new_tile;
                                run_info.TILE = tile;
                                
                                tile = run_model(tile);  %time integration
                            end
                        end
                    end
                end
                delete(gcp('nocreate'));
            else
                for sample_number = 1:size(run_info.CLUSTER.STATVAR.sample_centroid_index,1)
                    
                    run_number = run_info.CLUSTER.STATVAR.sample_centroid_index(sample_number,1);
                    for i=1:size(run_info.SPATIAL.ACTION,1)
                        run_info.SPATIAL.ACTION{i,1} = assign_tile_properties(run_info.SPATIAL.ACTION{i,1}, run_number); %writes the provider class
                    end
                    
                    disp(['running grid cell ' num2str(run_number)])
                    %as normal 1D run
                    for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
                        disp(['running tile number ' num2str(i)])
                        for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                            disp(['running round ' num2str(j)])
                            
                            new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                            new_tile.RUN_INFO = run_info;
                            new_tile = finalize_init(new_tile);
                            tile = new_tile;
                            run_info.TILE = tile;
                            
                            tile = run_model(tile);  %time integration
                        end
                    end
                end
            end
            
        end
        
        
        
    end
end



