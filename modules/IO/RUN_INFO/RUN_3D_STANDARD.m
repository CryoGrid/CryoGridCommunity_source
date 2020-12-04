% base class build a model tile

classdef RUN_3D_STANDARD < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        TILE
        %pprovider
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)

            run_info.PARA.coordinates = [];
            run_info.PARA.crs = [];
            run_info.PARA.height_system = [];
            
            run_info.PARA.number_of_tiles = []; %3;
            run_info.PARA.param_file_number = []; %[1;2;3];

            run_info.PARA.connected = [];
            run_info.PARA.contact_length = [];
            run_info.PARA.distance = [];

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
            parpool(run_info.PARA.number_of_tiles)
            spmd
                run_info.PARA.worker_number = labindex;
                %update the worker-specific name of the parameter file and the run
                run_info.PPROVIDER = update_parameter_file(run_info.PPROVIDER, run_info.PARA.param_file_number(run_info.PARA.worker_number,1));
                run_info.PPROVIDER = update_run_name(run_info.PPROVIDER, run_info.PARA.worker_number);
                %read worker-specific parameter file
                run_info.PPROVIDER = read_parameters(run_info.PPROVIDER);
                
                run_info = customize(run_info);
                
                %tile = run_info.PPROVIDER.FUNCTIONAL_CLASSES.TILE{run_info.PARA.tile_number,1};
                tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});

                tile.RUN_INFO = run_info;
                run_info.TILE = tile;

                tile = finalize_init(tile);
                
                tile = run(tile);  %time integration
            end
        end
 
        
        function run_info = customize(run_info)
            %FUNCTION TO BE EDITED BY USER - when inheriting from this class
            %here customizations can be done by directly writing pprovider
            %one can for example change parameters in the different
            %subsurface classes
            %run_info.PPROVIDER.FUNCTIONAL_CLASSES.XX = YY;
            %
            %use run_info.PARA.worker_number for customization

        end
        
        
        
        
        
    end
end



