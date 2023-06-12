%========================================================================
% CryoGrid STRATIGRAPHY_STATVAR class STRAT_linear 
% STRAT_linear defines the initial stratigraphy of model
% state variables by linearly interpolating between values at depths
% provided. Depths must be given as depth below the ground surface, and the
% final depth value must extend below the depth of the model domain.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef STRAT_linear2 < matlab.mixin.Copyable
    
    properties
		strat_linear_index
        STATVAR
		PARA
        CONST
        depth
        variable_names
        variable_values
        variable_gridded
    end
    
    methods


		function stratigraphy = provide_PARA(stratigraphy)
			stratigraphy.PARA.points = [];
        end
        
        function stratigraphy = provide_CONST(stratigraphy)

        end
        
        function stratigraphy = provide_STATVAR(stratigraphy)

        end 

        
        function stratigraphy = finalize_init(stratigraphy, tile)

            variables = fieldnames(stratigraphy.PARA.points);
            depth = stratigraphy.PARA.points.depth;
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'depth')
                    stratigraphy.STATVAR.(variables{i,1}) = interp1(depth, stratigraphy.PARA.points.(variables{i,1}), tile.GRID.STATVAR.MIDPOINTS, 'linear');
                end
            end           
		end


        
         %-------------param file generation-----
         function stratigraphy = param_file_info(stratigraphy)
             stratigraphy = provide_PARA(stratigraphy);
             %default
             stratigraphy.PARA.default_value = [];
             stratigraphy.PARA.STATVAR = [];
             stratigraphy.PARA.comment = [];
             stratigraphy.PARA.class_category = 'STRATIGRAPHY_STATVAR';
             stratigraphy.PARA.options.points.name = 'STRAT_MATRIX';
             stratigraphy.PARA.options.points.entries_y = {0; 10; 5000};
             stratigraphy.PARA.options.points.is_statvar_matrix = 1; %fill with STATVAR that are identified for initialization in GROUND classes
         end
        
    end
    
end