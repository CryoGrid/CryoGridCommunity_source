%========================================================================
% CryoGrid STRATIGRAPHY class STRAT_linear defines the initial stratigraphy
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
        
        %constructor
% 		function self = STRAT_linear(varargin)               % Temporary definition, to allow old code to run
%         %function self = STRAT_linear(index, pprovider, grid)      % Definition to be used when old code is no longer supported
%             % CONSTRUCTOR for STRAT_linear
%             %   Reads in variable profile from the specified file.
%             %
%             %   ARGUMENTS:
%             %   index:      user defined class index
%             %   pprovider:  instance of PARAMETER_PROVIDER class
% 			%	grid:		instance of GRID class
%             
%             % The following is only needed to allow legacy code to run
%             % May be removed when deprecated functions are removed
%             if nargin==3
%                 index = varargin{1};
%                 pprovider = varargin{2};
% 				grid = varargin{3};
%             else
%                 st = dbstack;
%                 warning(['DEPRECATION WARNING: Instantiating ' st.name '() with no arguments is deprecated.' newline,...
%                          'You should update your code to take advantage of new IO interface.']);
%                 return
%             end
%             % End allow legacy code
%             
%             self.strat_linear_index = index;
%             self = self.initialize();
% 			self = self.populate_variables(pprovider);
%             self = self.finalize_setup(grid);
%         end
		
		function self = initialize(self)
            % INITIALIZE  Initializes all properties needed by the class.

            self.depth = [];
			self.variable_names = {};
			self.variable_values = [];
			self.variable_gridded = [];
			self = self.initialize_PARA ();
        end	

		function self = provide_PARA(self)
			% INITIALIZE_PARA  Initializes initial conditions (initial conditions) in the PARA structure.
			self.PARA.points = [];
% 			self.PARA.initial_cond.depth = [];
% 			self.PARA.initial_cond.T = [];
        end
        
        function self = provide_CONST(self)

        end
        
        function self = provide_STATVAR(self)

        end 
        
        function self = initialize_excel(self)
            
        end
        
        function strat = finalize_init(strat, tile)
			% FINALIZE_SETUP  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
            variables = fieldnames(strat.PARA.points);
            depth = strat.PARA.points.depth;
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'depth')
                    tile.GRID.STATVAR.(variables{i,1}) = interp1(depth, strat.PARA.points.(variables{i,1}), tile.GRID.STATVAR.MIDPOINTS, 'linear');
                end
            end
			
            %conversion of variables, make this a dedicated class??
            for i=1:size(variables,1)
                if strcmp(variables{i,1}, 'waterIce') || strcmp(variables{i,1}, 'mineral') || strcmp(variables{i,1}, 'organic')
                    tile.GRID.STATVAR.(variables{i,1}) = tile.GRID.STATVAR.(variables{i,1}) .* tile.GRID.STATVAR.layerThick .* tile.PARA.area;
                 %ADD CONVERSION OF OTHER VARIABLES HERE
                elseif strcmp(variables{i,1}, 'Xice')
                    
                end
            end
            
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
		


        function xls_out = write_excel(self)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'STRATIGRAPHY','index';'STRAT_linear',1;NaN,NaN;'depth','T';'[m]','[degree C]';'TOP',NaN;0,1;1,0;10,-5;100,0;5000,20;'BOTTOM',NaN;'STRATIGRAPHY_END',NaN};
        end
        

        
    end
    
end