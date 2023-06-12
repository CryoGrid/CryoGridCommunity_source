
%========================================================================
% CryoGrid STRATIGRAPHYSTATVAR class STRAT_layers 
% STRAT_layers defines the initial stratigraphy of model state variables 
% as seperate layers with constant values. The depth below the surface 
% for the top position of each layer must be provided, and the last 
% layer is assumed to reach until the bottom of the model domain.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef STRAT_layers2 < matlab.mixin.Copyable
    
    properties
		strat_layers_index
        STATVAR
		PARA
        CONST
        depth
        variable_names
        variable_values
        variable_gridded
    end
    
    methods
        

		function stratigraphy = provide_PARA(stratigraphy)  %this is problematic, these fields must be completely arbitray, it changes accoriding to the requirements of the classes
			
			stratigraphy.PARA.layers = [];
        end
        
        function stratigraphy = provide_CONST(stratigraphy)

        end
        
        function stratigraphy = provide_STATVAR(stratigraphy)

        end 
        
        
        function stratigraphy = finalize_init(stratigraphy, tile)
			
            variables = fieldnames(stratigraphy.PARA.layers);
            depth = stratigraphy.PARA.layers.depth;
            depth = [depth; Inf];
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'depth')
                    stratigraphy.STATVAR.(variables{i,1}) = tile.GRID.STATVAR.MIDPOINTS .* 0;
                    for j=1:size(depth,1)-1
                        range = tile.GRID.STATVAR.MIDPOINTS > depth(j,1) & tile.GRID.STATVAR.MIDPOINTS <= depth(j+1,1);
                        stratigraphy.STATVAR.(variables{i,1})(range) = stratigraphy.PARA.layers.(variables{i,1})(j,1);
                    end
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
             stratigraphy.PARA.options.layers.name = 'STRAT_MATRIX';
             stratigraphy.PARA.options.layers.entries_y = {0; 10};
             stratigraphy.PARA.options.layers.is_statvar_matrix = 1; %fill with STATVAR that are identified for initialization in GROUND classes
         end
            
%             variable_gridded = repmat(tile.GRID.MIDPOINTS .* 0, 1, size(self.variable_values,2));
%             for j=1:size(self.variable_values,1)-1
%                 range = grid.MIDPOINTS > self.depth(j,1) & grid.MIDPOINTS <= self.depth(j+1,1);
%                 for i=1:size(self.variable_values,2)
%                     self.variable_gridded(range,i) = self.variable_gridded(range,i) + self.variable_values(j,i);
%                 end
%             end
%             range = grid.MIDPOINTS > self.depth(end,1);
%             for i=1:size(self.variable_values,2)
%                     self.variable_gridded(range,i) = self.variable_gridded(range,i) + self.variable_values(end,i);
%             end
		
        
        
% 		
% 		function self = populate_variables(self, pprovider)
% 			% POPULATE_VARIABLES  Updates the PARA structure with values from pprovider. Assigns values from the PARA structure to the corresponding class properties.
%             %
%             %   ARGUMENTS:
%             %   pprovider:  instance of PARAMETER_PROVIDER class
% 			self.PARA = pprovider.populate_struct(self.PARA, 'STRAT_layers', mfilename('class'), self.strat_layers_index);
% 			
% 			fn_substruct = fieldnames(self.PARA.layers);
% 			p = properties(self);
% 			for i = 1:size(fn_substruct, 1)
% 				if any(strcmp(p, fn_substruct(i)))
% 					index = find(strcmp(p, fn_substruct{i}));
% 					self.(p{index}) = self.PARA.layers.(fn_substruct{i});
% 				else
% 					self.variable_names = [self.variable_names fn_substruct(i)];
% 					self.variable_values = [self.variable_values self.PARA.layers.(fn_substruct{i})];
% 				end
% 			end
% 			self.depth = cell2mat(self.depth);
% 			self.variable_values = cell2mat(self.variable_values);		
%             
%         end
%         
%         
% 		
% 
% 
%         function xls_out = write_excel(self)
% 			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
%             
% 			xls_out = {'STRATIGRAPHY','index';'STRAT_linear',1;NaN,NaN;'depth','T';'[m]','[degree C]';'TOP',NaN;0,1;1,0;10,-5;100,0;5000,20;'BOTTOM',NaN;'STRATIGRAPHY_END',NaN};
%         end
        
%         % ==========================================
%         % DEPRECATED METHODS
%         % to be deleted when new implementation
%         % is validated for backwards compatibility
%         % ==========================================
%         
%         function self = initalize_from_file(self, section)
%             st = dbstack;
%             warning(['DEPRECATION WARNING: Method ' st.name '() is deprecated and will be removed.' newline,...
%                 'Use PARAMETER_PROVIDER class to obtain parameter values.']);
%             
%             pos_list = get_range_TOP_BOTTOM(section);
%             self.depth = cell2mat(section(pos_list(1,1):pos_list(1,2), 1));
%             self.variable_names={};
%             self.variable_values=[];
%             i=2;
%             field=cell2mat(section(pos_list(1,1)-3, i));
%             
%             while i<=size(section,2) && ~isnan(field(1))
%                 field=cell2mat(section(pos_list(1,1)-3, i));
%                 self.variable_names=[self.variable_names section{pos_list(1,1)-3, i}];
%                 self.variable_values = [self.variable_values cell2mat(section(pos_list(1,1):pos_list(1,2), i))];
%                 i=i+1;
%             end
%         end
%         
%         function self = initialize_from_table(self, table)
%             % INITIALIZE_FROM_TABLE  Initializes the variables from the output table of the csv parser.
%             
%             %	ARGUMENTS:
%             %	table:	Matlab output table from the csv parser
%             
%             st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                 'Code should be moved to new PARAMETER_PROVIDER class ',...
%                 'to streamline file access and the population of parameters.']);
%             
%             self.depth = table.depth;
%             self.variable_names = {};
%             self.variable_values = [];
%             i = 2;
%             field = table.Properties.VariableNames{i};
%             % Are the commented lines really necessary ? Seems to function
%             % without it ...
%             while i<=length(table.Properties.VariableNames) %&& ~isnan(field(1))
%                 self.variable_names = [self.variable_names table.Properties.VariableNames{i}];
%                 self.variable_values = [self.variable_values table.(self.variable_names{i-1})];
%                 i=i+1;
%                 %field = table.Properties.VariableNames{i};
%             end
%         end
%         
%         function self = interpolate_to_grid(self, grid)
%             st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                 'Parameter initialization should be finalized in the ' mfilename('class') '.finalize_setup() ']);
%             
%             self.variable_gridded = repmat(grid.MIDPOINTS .* 0, 1, size(self.variable_values,2));
%             for j=1:size(self.variable_values,1)-1
%                 range = grid.MIDPOINTS > self.depth(j,1) & grid.MIDPOINTS <= self.depth(j+1,1);
%                 for i=1:size(self.variable_values,2)
%                     self.variable_gridded(range,i) = self.variable_gridded(range,i) + self.variable_values(j,i);
%                 end
%             end
%             range = grid.MIDPOINTS > self.depth(end,1);
%             for i=1:size(self.variable_values,2)
%                 self.variable_gridded(range,i) = self.variable_gridded(range,i) + self.variable_values(end,i);
%             end
%         end
        
    end
    
end