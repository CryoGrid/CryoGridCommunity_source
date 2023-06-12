% base class build a model tile

classdef RUN_ENSEMBLE < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        TILE
        SPATIAL
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)

            
            run_info.PARA.ensemble_size = []; 

            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];

            run_info.PARA.point_class = [];
            run_info.PARA.point_class_index = [];
        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        

        
        function run_info = finalize_init(run_info)
            run_info.SPATIAL.STATVAR.latitude = 60;
            run_info.SPATIAL.STATVAR.longitude = 10;
            run_info.SPATIAL.STATVAR.altitude = 0;
            run_info.SPATIAL.STATVAR.area = 1;
            run_info.SPATIAL.STATVAR.slope_angle = 0;
            run_info.SPATIAL.STATVAR.aspect = 0;
            run_info.SPATIAL.STATVAR.skyview_factor = 0;
            run_info.SPATIAL.STATVAR.horizon_bins = 0;
            run_info.SPATIAL.STATVAR.horizon_angles = 0;
            
            if ~isempty(run_info.PARA.point_class) && sum(isnan(run_info.PARA.point_class))==0
                run_info.SPATIAL = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.point_class){run_info.PARA.point_class_index,1});
                run_info.SPATIAL.RUN_INFO = run_info;
                run_info.SPATIAL = finalize_init(run_info.SPATIAL);
            end
        end
        
        
        
        function [run_info, tile] = run_model(run_info)

            parpool(run_info.PARA.ensemble_size)
            spmd
                run_info.PARA.worker_number = labindex;
                                
                %read worker-specific parameter file
                %run_info.PPROVIDER = read_parameters(run_info.PPROVIDER);
                                
                tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});
                tile.RUN_INFO = run_info;
                run_info.TILE = tile;
                
                fn = fieldnames(run_info.SPATIAL.STATVAR);
                for i=1:size(fn,1)  %be careful, does not work if empty array (and not NaN) is willingly assigned to a parameter
                    if ~isempty(run_info.SPATIAL.STATVAR.(fn{i,1}))
                        tile.PARA.(fn{i,1}) = run_info.SPATIAL.STATVAR.(fn{i,1});
                    end
                end
                
                tile.PARA.worker_number = run_info.PARA.worker_number;
                tile.PARA.ensemble_size = run_info.PARA.ensemble_size;
                tile = finalize_init(tile);
                
                tile.PARA.run_name = [tile.PARA.run_name '_' num2str(run_info.PARA.worker_number)];
                
                tile = run_model(tile);  %time integration
            end
        end
 
        
        %-------------param file generation-----
%         function out = param_file_info(out)
%             out = provide_PARA(out);
% 
%             out.PARA.STATVAR = [];
%             out.PARA.options = [];
%             out.PARA.class_category = 'RUN_INFO';
%             
%             out.PARA.default_value.ensemble_size = {30};
%             out.PARA.comment.ensemble_size = {'number of ensemble members/cores'};
%             
%             out.PARA.default_value.tile_class = {'TILE_1D_standard'};
%             out.PARA.comment.tile_class = {'TILE class'};
%             
%             out.PARA.default_value.tile_class_index = {1};
%             out.PARA.comment.tile_class_index = {'TILE class index'};
%             
%         end
        
    end
end



