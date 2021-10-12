% base class build a model tile

classdef RUN_ENSEMBLE < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        TILE
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)

            
            run_info.PARA.number_of_tiles = []; 

            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];

        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        

        
        function run_info = finalize_init(run_info)
 
        end
        
        
        
        function [run_info, tile] = run_model(run_info)

            parpool(run_info.PARA.number_of_tiles)
            spmd
                run_info.PARA.worker_number = labindex;
                                
                %read worker-specific parameter file
                run_info.PPROVIDER = read_parameters(run_info.PPROVIDER);
                                
                %tile = run_info.PPROVIDER.FUNCTIONAL_CLASSES.TILE{run_info.PARA.tile_number,1};
                tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});

                tile.RUN_INFO = run_info;
                run_info.TILE = tile;

                tile = finalize_init(tile);
                
                tile = run_model(tile);  %time integration
            end
        end
 
        
        
    end
end



