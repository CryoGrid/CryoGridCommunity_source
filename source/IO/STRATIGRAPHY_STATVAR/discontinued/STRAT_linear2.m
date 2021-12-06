%========================================================================
% CryoGrid STRATIGRAPHY class STRAT_linear defines the initial stratigraphy
% state variables by linearly interpolating between values at depths
% provided. Depths must be given as depth below the ground surface, and the
% final depth value must extend below the depth of the model domain.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef STRAT_linear2 < matlab.mixin.Copyable
    
    properties
		strat_linear_index
		PARA
        CONST
        depth
        variable_names
        variable_values
        variable_gridded
    end
    
    methods
        

		


		function self = provide_PARA(self)
			self.PARA.points = [];

        end
        
        function self = provide_CONST(self)

        end
        
        function self = provide_STATVAR(self)

        end 
        

        function strat = finalize_init(strat, tile)

            
            variables = fieldnames(strat.PARA.points);
            depth = strat.PARA.points.depth;
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'depth')
                    if ~isfield(tile.GRID.STATVAR, variables{i,1})
                        tile.GRID.STATVAR.(variables{i,1}) = tile.GRID.STATVAR.layerThick .* 0;
                    end
                    store_var = interp1(depth, strat.PARA.points.(variables{i,1}), tile.GRID.STATVAR.MIDPOINTS, 'linear');
                    ind = find(tile.PARA.stratigraphy >= strat.PARA.class_number);
                    store_var = repmat(store_var, 1, tile.PARA.number_of_realizations);
                    tile.GRID.STATVAR.(variables{i,1})(:,ind) = store_var(:,ind);
                end
            end
		end
        
       


        function xls_out = write_excel(self)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'STRATIGRAPHY','index';'STRAT_linear',1;NaN,NaN;'depth','T';'[m]','[degree C]';'TOP',NaN;0,1;1,0;10,-5;100,0;5000,20;'BOTTOM',NaN;'STRATIGRAPHY_END',NaN};
        end
        

        
    end
    
end