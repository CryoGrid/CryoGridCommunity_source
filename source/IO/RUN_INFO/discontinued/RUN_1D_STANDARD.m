%========================================================================
% CryoGrid RUN_INFO class RUN_1D_STANDARD
% RUN_INFO class which runs a single TILE class only
% used for simplest possible simulation setup with only a single model run 

% S. Westermann, Jan 2021
%========================================================================
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
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
%         function run_info = initialize_excel(run_info)
%             
%         end
        
        function run_info = finalize_init(run_info)
 
        end
        

        function [run_info, tile] = setup_run(run_info)
            %run_info = customize(run_info)
            
            tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});
            tile.RUN_INFO = run_info;
            run_info.TILE = tile;
            
            
            %do the run(s) 
            tile = finalize_init(tile);
        end
 
        
        
        function [run_info, tile] = run_model(run_info)
            %run_info = customize(run_info)
            [run_info, tile] = setup_run(run_info);
            tile = run_model(tile);  %time integration
        end
 

        function run_info = customize(run_info)
            %FUNCTION TO BE EDITED BY USER - when inheriting from this class
            %here customizations can be done by directly writing pprovider
            %one can for example change parameters in the different
            %subsurface classes
%             run_info.PPROVIDER.CLASSES.GROUND_freeW_seb.PARA.albedo = i/100;
%             run_info.PPROVIDER.CLASSES.FROCING_seb.filename = ['myfocing_ ' num2str(i) '.mat'];
%             run_info.PPROVIDER.FUNCTIONAL_CLASSES.XX = YY;
%             run_info.PPROVIDER.FUNCTIONAL_CLASSES.TILE{run_info.PARA.tile_number,1}.PARA.lat = f(id)
            
        end
        
        
        %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'RUN_INFO';
            
            out.PARA.default_value.tile_class = {'TILE_1D_standard'};
            out.PARA.comment.tile_class = {'TILE class'};
            
            out.PARA.default_value.tile_class_index = {1};
            out.PARA.comment.tile_class_index = {'TILE class index'};
            
        end
        
        
    end
end



