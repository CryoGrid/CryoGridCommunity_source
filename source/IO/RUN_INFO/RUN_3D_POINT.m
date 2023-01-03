%========================================================================
% CryoGrid RUN_INFO class RUN_3D_POINT
% RUN_INFO class designed to run several TILE class in parallel, including
% the possibility for laterally coupled TILE classes with the LATERAL class
% LATERAL_3D. 
% Matlab parallel computing toolbox is required!
%
% S. Westermann, Jan 2021
% S. Westermann, Dec 2022
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

            run_info.PARA.point_3D_class = [];
            run_info.PARA.point_3D_class_index = [];

            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];

        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
        
        function run_info = finalize_init(run_info)
 
            run_info.SPATIAL = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.point_3D_class){run_info.PARA.point_3D_class_index,1});
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
            
%             %this is double and in principle unnecessary, potentially
%             %remove this, then all references to these fields in LATERAL need to be changed!! 
%             run_info.PARA.connected = run_info.SPATIAL.PARA.connected;
%             run_info.PARA.contact_length = run_info.SPATIAL.PARA.contact_length;
%             run_info.PARA.distance = run_info.SPATIAL.PARA.distance;
%             %end remove!
            
            fn = fieldnames(run_info.SPATIAL.STATVAR);
            for i=1:size(fn,1) 
                if ~isempty(run_info.SPATIAL.STATVAR.(fn{i,1}))
                    tile.PARA.(fn{i,1}) = run_info.SPATIAL.STATVAR.(fn{i,1})(run_info.PARA.worker_number,1);
                end
            end
            
            tile = finalize_init(tile);
        end


        
        
        %-------------param file generation-----
        function run_info = param_file_info(run_info)
            run_info = provide_PARA(run_info);

            run_info.PARA.STATVAR = [];
            run_info.PARA.class_category = 'RUN_INFO';
            
            run_info.PARA.default_value.tile_class = {'TILE_1D_standard'};
            run_info.PARA.comment.tile_class = {'TILE class'};
            
            run_info.PARA.default_value.tile_class_index = {1};
            run_info.PARA.comment.tile_class_index = {'TILE class index'};
            
            run_info.PARA.default_value.point_3D_class = {'POINT_3D_SIMPLE'};
            run_info.PARA.comment.point_3D_class = {'3D point class containing information on locations and topological relationships between tiles'};
            
            run_info.PARA.default_value.point_3D_class_index = {1};
        end
        
    end
end



