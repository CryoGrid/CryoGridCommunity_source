%========================================================================
% CryoGrid STRATIGRAPHY_STATVAR class STRAT_linear 
% STRAT_linear defines the initial stratigraphy of model
% state variables by linearly interpolating between values at depths
% provided. Depths must be given as depth below the ground surface, and the
% final depth value must extend below the depth of the model domain.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef STRAT_linear < matlab.mixin.Copyable
    
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
                    tile.GRID.STATVAR.(variables{i,1}) = interp1(depth, stratigraphy.PARA.points.(variables{i,1}), tile.GRID.STATVAR.MIDPOINTS, 'linear');
                end
            end           
		end

        function stratigraphy = finalize_init_GROUND_multi_tile(stratigraphy, GRID)
            variables = fieldnames(stratigraphy.PARA.points);
            depth = stratigraphy.PARA.points.depth;
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'depth')
                    GRID.STATVAR.(variables{i,1}) = interp1(depth, stratigraphy.PARA.points.(variables{i,1}), GRID.STATVAR.MIDPOINTS, 'linear');
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

% 		
% 		function self = populate_variables(self, pprovider)
% 			% POPULATE_VARIABLES  Updates the PARA structure with values from pprovider. Assigns values from the PARA structure to the corresponding class properties.
%             %
%             %   ARGUMENTS:
%             %   pprovider:  instance of PARAMETER_PROVIDER class
%             
% 			self.PARA = pprovider.populate_struct(self.PARA, 'STRAT_linear', mfilename('class'), self.strat_linear_index);
% 			
% 			fn_substruct = fieldnames(self.PARA.initial_cond);
% 			p = properties(self);
% 			for i = 1:size(fn_substruct, 1)
% 				if any(strcmp(p, fn_substruct(i)))
% 					index = find(strcmp(p, fn_substruct{i}));
% 					self.(p{index}) = self.PARA.initial_cond.(fn_substruct{i});
% 				else
% 					self.variable_names = [self.variable_names fn_substruct(i)];
% 					self.variable_values = [self.variable_values self.PARA.initial_cond.(fn_substruct{i})];
% 				end
% 			end
% 			self.depth = cell2mat(self.depth);
% 			self.variable_values = cell2mat(self.variable_values);					
% 		end		
% 		
% 
% 
%         function xls_out = write_excel(self)
% 			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
% 			
%             xls_out = {'STRATIGRAPHY','index';'STRAT_linear',1;NaN,NaN;'depth','T';'[m]','[degree C]';'TOP',NaN;0,1;1,0;10,-5;100,0;5000,20;'BOTTOM',NaN;'STRATIGRAPHY_END',NaN};
%         end
%         

        
    end
    
end