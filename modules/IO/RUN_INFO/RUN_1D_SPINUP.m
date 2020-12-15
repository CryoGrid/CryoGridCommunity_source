% base class build a model tile

classdef RUN_1D_SPINUP < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        TILE
        PARA
        CONST
        STATVAR
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)

            run_info.PARA.coordinates = [];
            run_info.PARA.crs = [];
            run_info.PARA.height_system = [];
            
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            run_info.PARA.number_of_runs_per_tile = []; %vector
        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
        function run_info = initialize_excel(run_info)
            
        end
        
        function run_info = finalize_init(run_info)
 
        end
        
        
        
        function [run_info, tile] = run_model(run_info)
            %this could first open spmd and assign run_number depending on
            %worker, then do another round of pprovider
            %it could also do a loop over different tile representing
            %different sections of the run, e.g. initial inial init, spin-up, actual run 
            
            
            run_info = customize(run_info);
            
            tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{1,1}){run_info.PARA.tile_class_index(1,1),1});

            tile.RUN_INFO = run_info;
            run_info.TILE = tile;
            
            %do the run(s) - this is blanked out, and the code is handled in the main run file - for more complicated or opertaional runs, it should be handled here. 
            %assemble the stratigraphy, etc.
            tile = finalize_init(tile);
            
            
            tile = run_model(tile);  %time integration
            
            for i=2:size(run_info.PARA.tile_class,1)
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
 
        
        function run_info = customize(run_info)
            %FUNCTION TO BE EDITED BY USER - when inheriting from this class
            %here customizations can be done by directly writing pprovider
            %one can for example change parameters in the different
            %subsurface classes
            %run_info.PPROVIDER.FUNCTIONAL_CLASSES.XX = YY;
            
            
        end
        
        
        
        
        
    end
end



