%========================================================================
% CryoGrid GRID class  for defining the compute grid
% GRID_user_defined defines the compute grid as layers with constant grid spacing
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef GRID_user_defined

    properties
        grid_index
		PARA
        GRID
        MIDPOINTS
        LAYERTHICK
        variable_names
        variable_gridded
    end
    
    methods
        function self = GRID_user_defined (varargin)               % Temporary definition, to allow old code to run
        %function self = GRID_user_defined(index, pprovider, forcing)      % Definition to be used when old code is no longer supported
            % CONSTRUCTOR for GRID_user_defined
            %   Reads in grid parameters from the specified file.
            %
            %   ARGUMENTS:
            %   index:      user defined class index
            %   pprovider:  instance of PARAMETER_PROVIDER class
			%	forcing:	instance of FORCING class
            
            % The following is only needed to allow legacy code to run
            % May be removed when deprecated functions are removed
            if nargin==3
                index = varargin{1};
                pprovider = varargin{2};
                forcing = varargin{3};
            else
                st = dbstack;
                warning(['DEPRECATION WARNING: Instantiating ' st.name '() with no arguments is deprecated.' newline,...
                         'You should update your code to take advantage of new IO interface.']);
                return
            end
            % End allow legacy code
            
            self.grid_index = index;
            self = self.initialize();
            self = self.populate_GRID(pprovider);
            self = self.finalize_setup(forcing);
        end
        
        function self = initialize(self)
			% INITIALIZE  Initializes all properties needed by the class.

			self.MIDPOINTS = [];
			self.LAYERTHICK = [];
			self.variable_names = [];
			self.variable_gridded = [];
			self.GRID = [];
			self = self.initialize_PARA();
        end
		
		function self = initialize_PARA(self)
			% INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  
			
			self.PARA.grid.upper = [];
			self.PARA.grid.spacing = [];
			self.PARA.grid.lower = [];
		end
		
		function self = populate_GRID(self, pprovider)
			% POPULATE_GRID  Updates the PARA structure with values from pprovider. Assigns values from the PARA structure to the corresponding class properties.
            %
            %   ARGUMENTS:
            %   pprovider:  instance of PARAMETER_PROVIDER class
			
            self.PARA = pprovider.populate_struct(self.PARA, 'GRID', mfilename('class'), self.grid_index);
			
			fn_substruct = fieldnames(self.PARA.grid);
			grid_breaks = [];
			for i = 1:size(fn_substruct, 1)
				grid_breaks = [grid_breaks self.PARA.grid.(fn_substruct{i})];
			end
			
			grid_breaks = cell2mat(grid_breaks);
            grid_breaks(2:end,1) = grid_breaks(2:end,1) + grid_breaks(2:end,2);
            for i=1:size(grid_breaks)
				self.GRID = [self.GRID; [grid_breaks(i,1):grid_breaks(i,2):grid_breaks(i,3)]'];
            end
		end
		
		function self = finalize_setup(self, forcing)
			% FINALIZE_SETUP  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
			
			%	ARGUMENTS:
			%	forcing:	instance of FORCING class
			
			self.GRID(self.GRID > forcing.PARA.domain_depth)=[]; 
            self.MIDPOINTS = (self.GRID(2:end,1) + self.GRID(1:end-1,1))./2;
            self.LAYERTHICK = (self.GRID(2:end,1) - self.GRID(1:end-1,1));
		end

        function xls_out = write_excel(self)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'GRID','index',NaN;'GRID_user_defined',1,NaN;NaN,NaN,NaN;'upper','spacing','lower';'[m]','[m]','[m]';'TOP',NaN,NaN;0,0.0500000000000000,2;2,0.100000000000000,10;10,0.500000000000000,30;30,1,50;50,5,100;100,10,150;150,50,500;500,100,1000;'BOTTOM',NaN,NaN;'GRID_END',NaN,NaN};
        end
        
		
		% ==========================================
        % DEPRECATED METHODS
        % to be deleted when new implementation
        % is validated for backwards compatibility
        % ==========================================
		
		
%         function self = initalize_from_file(self, section)
% 			st = dbstack;
%             warning(['DEPRECATION WARNING: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Use PARAMETER_PROVIDER class to obtain parameter values.']);
%             
% 			pos_list = get_range_TOP_BOTTOM(section);
%             grid_breaks = cell2mat(section(pos_list(1,1):pos_list(1,2), 1:3));
%             grid_breaks(2:end,1) = grid_breaks(2:end,1) + grid_breaks(2:end,2);
%             self.GRID =[];
%             for i=1:size(grid_breaks)
%                 self.GRID = [self.GRID; [grid_breaks(i,1):grid_breaks(i,2):grid_breaks(i,3)]'];
%             end
%         end
%         
%         function self = initialize_from_table(self, table)
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
%             grid_breaks = table2array(table);
%             grid_breaks(2:end,1) = grid_breaks(2:end,1) + grid_breaks(2:end,2);
%             self.GRID =[];
%             for i=1:size(grid_breaks)
%                 self.GRID = [self.GRID; [grid_breaks(i,1):grid_breaks(i,2):grid_breaks(i,3)]'];
%             end
%         end
% 		
%         function self = reduce_grid(self, forcing)
% 			% REDUCE_GRID  Adapts grid to domain depth specified during the forcing initialization.
% 			
% 			%	ARGUMENTS:
% 			%	forcing:	instance of FORCING class
% 			
% 			st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Parameter initialization should be finalized in the ' mfilename('class') '.finalize_setup() ']);
% 					 
%             self.GRID(self.GRID > forcing.PARA.domain_depth)=[]; 
%             self.MIDPOINTS = (self.GRID(2:end,1) + self.GRID(1:end-1,1))./2;
%             self.LAYERTHICK = (self.GRID(2:end,1) - self.GRID(1:end-1,1));
%         end
    end
    
end

