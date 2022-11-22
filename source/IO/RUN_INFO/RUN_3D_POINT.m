%========================================================================
% CryoGrid RUN_INFO class RUN_3D_STANDARD
% RUN_INFO class designed to run several TILE class in parallel, including
% the possibility for laterally coupled TILE classes with the LATERAL class
% LATERAL_3D. The topological relationships between the TILE classes are
% provided in this class
% RUN_3D_STANDARD can also be used for parallel runs of several independent
% TILE classes, e.g. for a sensitivity analysis. Note that RUN_3D_STANDARD
% does not offer the possibility to combine parallel and sequential
% simulations of TILE classes, which would be required to run a large number 
% of TILE classes
% Matlab parallel computing toolbox is required!

% S. Westermann, Jan 2021
%========================================================================

classdef RUN_3D_POINT < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        SPATIAL
        TILE

    end
    
    methods
        
        function run_info = provide_PARA(run_info)

            run_info.PARA.3D_point_class = [];
            run_info.PARA.3D_point_class_index = [];
            
%             run_info.PARA.number_of_tiles = []; %3;
%             run_info.PARA.param_file_number = []; %[1;2;3];            
%             run_info.PARA.connected = [];
%             run_info.PARA.contact_length = [];
%             run_info.PARA.distance = [];

            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];

        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
        
        function run_info = finalize_init(run_info)
 
            run_info.SPATIAL = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.point_class){run_info.PARA.point_class_index,1});
            run_info.SPATIAL.RUN_INFO = run_info;
            run_info.SPATIAL = finalize_init(run_info.SPATIAL);
            
        end
        
        
        
        function [run_info, tile] = run_model(run_info)
            %this could first open spmd and assign run_number depending on
            %worker, then do another round of pprovider
            %it could also do a loop over different tile representing
            %different sections of the run, e.g. initial inial init, spin-up, actual run 
            
            parpool(run_info.SPATIAL.PARA.number_of_tiles)
            spmd
                run_info.PARA.worker_number = labindex; %->NEEDS TWO DIFFERENT VARIABLES, ONE THE ACTUAL WORKER NUMBER AND THE OTHER ONE THE NUMBER IN THE 3D_POINT CLASS
                [run_info, tile] = setup_run(run_info);
                
                tile = run_model(tile);  %time integration
            end
        end
 

        function [run_info, tile] = setup_run(run_info)
            %update the worker-specific name of the parameter file and the run
            run_info.PPROVIDER = update_parameter_file(run_info.PPROVIDER, run_info.SPATIAL.PARA.param_file_number(run_info.PARA.worker_number,1));
            run_info.PPROVIDER = update_run_name(run_info.PPROVIDER, run_info.PARA.worker_number);
                            
            %read worker-specific parameter file
            run_info.PPROVIDER = read_parameters(run_info.PPROVIDER);
                                        
            tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});

            tile.RUN_INFO = run_info;
            run_info.TILE = tile;
            
            fn = fieldnames(run_info.SPATIAL.STATVAR);
            for i=1:size(fn,1) 
                if ~isempty(run_info.SPATIAL.STATVAR.(fn{i,1}))
                    tile.PARA.(fn{i,1}) = run_info.SPATIAL.STATVAR.(fn{i,1});
                end
            end
            
            tile = finalize_init(tile);
        end


        
        
        %-------------param file generation-----
        function run_info = param_file_info(run_info)
            run_info = provide_PARA(run_info);

            run_info.PARA.STATVAR = [];
            run_info.PARA.class_category = 'RUN_INFO';
            
            run_info.PARA.default_value.number_of_tiles = {3};
            run_info.PARA.comment.number_of_tiles = {'number of tiles/cores'};
            
            run_info.PARA.options.param_file_number.name =  'H_LIST';
            run_info.PARA.options.param_file_number.entries_x = {1 2 3};
            

            run_info.PARA.options.connected.name = 'MATRIX';
            run_info.PARA.options.connected.entries_matrix = {'0' '1' '0'; '1' '0' '1'; '0' '1' '0'};
            run_info.PARA.comment.connected = {'1 if connected'};
            
            run_info.PARA.options.contact_length.name = 'MATRIX';
            run_info.PARA.options.contact_length.entries_matrix = {'0' '1' '0'; '1' '0' '1'; '0' '1' '0'};
            run_info.PARA.comment.contact_length = {'lateral contact length between tiles'};

            run_info.PARA.options.distance.name = 'MATRIX';
            run_info.PARA.options.distance.entries_matrix = {'0' '1' '0'; '1' '0' '1'; '0' '1' '0'};
            run_info.PARA.comment.distance = {'distance between tiles'};
            
            run_info.PARA.default_value.tile_class = {'TILE_1D_standard'};
            run_info.PARA.comment.tile_class = {'TILE class'};
            
            run_info.PARA.default_value.tile_class_index = {1};
            run_info.PARA.comment.tile_class_index = {'TILE class index'};
        end
        
        

        
    end
end



