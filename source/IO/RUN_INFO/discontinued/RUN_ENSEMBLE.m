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
                
                tile.PARA.run_name = [tile.PARA.run_name '_' num2str(run_info.PARA.worker_number)];
                
                
                tile = run_model(tile);  %time integration
            end
        end
 
        
        %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'RUN_INFO';
            
            out.PARA.default_value.number_of_tiles = {30};
            out.PARA.comment.number_of_tiles = {'number of ensemble members/cores'};
            
            out.PARA.default_value.tile_class = {'TILE_1D_standard'};
            out.PARA.comment.tile_class = {'TILE class'};
            
            out.PARA.default_value.tile_class_index = {1};
            out.PARA.comment.tile_class_index = {'TILE class index'};
            
        end
        
    end
end



