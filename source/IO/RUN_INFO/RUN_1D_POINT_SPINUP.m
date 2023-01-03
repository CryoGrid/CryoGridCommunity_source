%========================================================================
% CryoGrid RUN_INFO class RUN_1D_POINT_SPINUP
% RUN_INFO class which runs several TILE class sequentially
% can be used for model spin-up when using TILE_BUILDER classes which
% initialize subsequent TILE classesbased on the results of the previous
% TILE
% can also be used for seuqntial runs of several independent TILE classes,
% e.g. for a sensitivity analysis

% S. Westermann, Jan 2021
% S. westermann, Dec 2022
%========================================================================

classdef RUN_1D_POINT_SPINUP < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        TILE
        PARA
        CONST
        SPATIAL
        STATVAR
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)
            
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            run_info.PARA.number_of_runs_per_tile = []; %vector
            
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

            
            for i=1:size(run_info.PARA.tile_class,1)
                disp(['running tile number ' num2str(i)])
                for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                    disp(['running round ' num2str(j)])
                                        
                    new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                    fn = fieldnames(run_info.SPATIAL.STATVAR);
                    for k=1:size(fn,1)  %be careful, does not work if empty array (and not NaN) is willingly assigned to a parameter
                        if ~isempty(run_info.SPATIAL.STATVAR.(fn{k,1}))
                            new_tile.PARA.(fn{k,1}) = run_info.SPATIAL.STATVAR.(fn{k,1});
                        end
                    end
                    new_tile.RUN_INFO = run_info;
                    new_tile = finalize_init(new_tile);
                    
                    tile = new_tile;
                    run_info.TILE = tile;
                    
                    tile = run_model(tile);  %time integration
                end
            end
            
        end
 
        
        
        %-------------param file generation-----
        function run_info = param_file_info(run_info)
            run_info = provide_PARA(run_info);

            run_info.PARA.STATVAR = [];
            run_info.PARA.class_category = 'RUN_INFO';
            run_info.PARA.default_value = [];
            run_info.PARA.comment = [];
            
            run_info.PARA.options.tile_class.name =  'H_LIST';
            run_info.PARA.options.tile_class.entries_x = {'TILE_1D_standard' 'TILE_1D_standard'};
            
            run_info.PARA.options.tile_class_index.name =  'H_LIST'; 
            run_info.PARA.options.tile_class_index.entries_x = {1 2};
            
            run_info.PARA.options.number_of_runs_per_tile.name =  'H_LIST'; % 
            run_info.PARA.options.number_of_runs_per_tile.entries_x = {1 1};
            
            run_info.PARA.comment.point_class = {'point class providing information on the location; if empty, no location, altitude = 0m and area = 1m2 is assumed'};
            
        end
        
        
    end
end



