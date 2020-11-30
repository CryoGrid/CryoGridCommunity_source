% base class build a model tile

classdef RUN_1D_STANDARD < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        TILE
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)

            run_info.PARA.coordinates = [];
            run_info.PARA.crs = [];
            run_info.PARA.height_system = [];

            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
        function run_info = initialize_excel(run_info)
            
        end
        
        function run_info = finalize_init(run_info)
 
        end
        
        
        
        function [run_info, tile] = run(run_info)
            %this could first open spmd and assign run_number depending on
            %worker, then do another round of pprovider
            %it could also do a loop over different tile representing
            %different sections of the run, e.g. initial inial init, spin-up, actual run 
            
            tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});
            tile.RUN_INFO = run_info;
            run_info.TILE = tile;

            
            %do the run(s) - this is blanked out, and the code is handled in the main run file - for more complicated or opertaional runs, it should be handled here. 
            %assemble the stratigraphy, etc.
            tile = finalize_init(tile);
            tile = run(tile);  %time integration

        end
 
        
        function run_info = customize(run_info)
            %FUNCTION TO BE EDITED BY USER - when inheriting from this class
            %here customizations can be done by directly writing pprovider
            %one can for example change parameters in the different
            %subsurface classes
            %run_info.PPROVIDER.FUNCTIONAL_CLASSES.XX = YY;
            %run_info.PPROVIDER.FUNCTIONAL_CLASSES.TILE{run_info.PARA.tile_number,1}.PARA.lat = f(id)
            
        end
        
        
        
        
        
    end
end



