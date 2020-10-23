%========================================================================
% CryoGrid STRATIGRAPHY class STRAT_linear defines the initial stratigraphy
% state variables by linearly interpolating between values at depths
% provided. Depths must be given as depth below the ground surface, and the
% final depth value must extend below the depth of the model domain.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef STRAT_linear
    
    properties
		strat_linear_index
		PARA
        depth
        variable_names
        variable_values
        variable_gridded
    end
    
    methods
        
        %constructor
		function self = STRAT_linear(varargin)               % Temporary definition, to allow old code to run
        %function self = STRAT_linear(index, pprovider, grid)      % Definition to be used when old code is no longer supported
            % CONSTRUCTOR for STRAT_linear
            %   Reads in variable profile from the specified file.
            %
            %   ARGUMENTS:
            %   index:      user defined class index
            %   pprovider:  instance of PARAMETER_PROVIDER class
			%	grid:		instance of GRID class
            
            % The following is only needed to allow legacy code to run
            % May be removed when deprecated functions are removed
            if nargin==3
                index = varargin{1};
                pprovider = varargin{2};
				grid = varargin{3};
            else
                st = dbstack;
                warning(['DEPRECATION WARNING: Instantiating ' st.name '() with no arguments is deprecated.' newline,...
                         'You should update your code to take advantage of new IO interface.']);
                return
            end
            % End allow legacy code
            
            self.strat_linear_index = index;
            self = self.initialize();
			self = self.populate_variables(pprovider);
            self = self.finalize_setup(grid);
        end
		
		function self = initialize(self)
            % INITIALIZE  Initializes all properties needed by the class.

            self.depth = [];
			self.variable_names = {};
			self.variable_values = [];
			self.variable_gridded = [];
			self = self.initialize_PARA ();
        end	

		function self = initialize_PARA(self)
			% INITIALIZE_PARA  Initializes initial conditions (initial conditions) in the PARA structure.
			
			self.PARA.initial_cond.depth = [];
			self.PARA.initial_cond.T = [];
		end
		
		function self = populate_variables(self, pprovider)
			% POPULATE_VARIABLES  Updates the PARA structure with values from pprovider. Assigns values from the PARA structure to the corresponding class properties.
            %
            %   ARGUMENTS:
            %   pprovider:  instance of PARAMETER_PROVIDER class
            
			self.PARA = pprovider.populate_struct(self.PARA, 'STRAT_linear', mfilename('class'), self.strat_linear_index);
			
			fn_substruct = fieldnames(self.PARA.initial_cond);
			p = properties(self);
			for i = 1:size(fn_substruct, 1)
				if any(strcmp(p, fn_substruct(i)))
					index = find(strcmp(p, fn_substruct{i}));
					self.(p{index}) = self.PARA.initial_cond.(fn_substruct{i});
				else
					self.variable_names = [self.variable_names fn_substruct(i)];
					self.variable_values = [self.variable_values self.PARA.initial_cond.(fn_substruct{i})];
				end
			end
			self.depth = cell2mat(self.depth);
			self.variable_values = cell2mat(self.variable_values);					
		end		
		
		function self = finalize_setup(self, grid)
			% FINALIZE_SETUP  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
			
			%	ARGUMENTS:
			%	grid:	instance of GRID class
			
            for i=1:size(self.variable_values,2)
                self.variable_gridded = [self.variable_gridded; interp1(self.depth, self.variable_values(:,i), grid.MIDPOINTS, 'linear')];
            end
		end

        function xls_out = write_excel(self)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'STRATIGRAPHY','index';'STRAT_linear',1;NaN,NaN;'depth','T';'[m]','[degree C]';'TOP',NaN;0,1;1,0;10,-5;100,0;5000,20;'BOTTOM',NaN;'STRATIGRAPHY_END',NaN};
        end
        
% 		% ==========================================
%         % DEPRECATED METHODS
%         % to be deleted when new implementation
%         % is validated for backwards compatibility
%         % ==========================================
% 		
%         function self = initalize_from_file(self, section)
% 			st = dbstack;
% 			warning(['DEPRECATION WARNING: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Use PARAMETER_PROVIDER class to obtain parameter values.']);
% 					 
%             pos_list = get_range_TOP_BOTTOM(section);
%             self.depth = cell2mat(section(pos_list(1,1):pos_list(1,2), 1));
%             self.variable_names={};
%             self.variable_values=[];
%             i=2;
%             field=cell2mat(section(pos_list(1,1)-3, i));
% 			
%             while i<=size(section,2) && ~isnan(field(1))
%                 self.variable_names=[self.variable_names section{pos_list(1,1)-3, i}];
%                 self.variable_values = [self.variable_values cell2mat(section(pos_list(1,1):pos_list(1,2), i))];
%                 i=i+1;
%                 field=cell2mat(section(pos_list(1,1)-3, i));
%             end
%         end
% 		
% 		function self = initialize_from_table(self, table)
% 			% INITIALIZE_FROM_TABLE  Initializes the variables from the output table of the csv parser.
% 			
% 			%	ARGUMENTS:
% 			%	table:	Matlab output table from the csv parser
% 			
% 			st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Code should be moved to new PARAMETER_PROVIDER class ',...
%                      'to streamline file access and the population of parameters.']);   
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
% 			st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Parameter initialization should be finalized in the ' mfilename('class') '.finalize_setup() ']);
%             
% 			self.variable_gridded = [];
%             size(self.variable_values,2);
%             for i=1:size(self.variable_values,2)
%                 self.variable_gridded = [self.variable_gridded; interp1(self.depth, self.variable_values(:,i), grid.MIDPOINTS, 'linear')];
%             end
%         end
        
    end
    
end

