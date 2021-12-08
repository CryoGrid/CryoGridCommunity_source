%========================================================================
% CryoGrid GRID class  for defining the compute grid
% GRID_user_defined defines the compute grid as layers with constant grid spacing
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef GRID_user_defined < matlab.mixin.Copyable

    properties

		PARA
        CONST
        STATVAR

    end
    
    methods
%         function self = GRID_user_defined (varargin)               % Temporary definition, to allow old code to run
%         %function self = GRID_user_defined(index, pprovider, forcing)      % Definition to be used when old code is no longer supported
%             % CONSTRUCTOR for GRID_user_defined
%             %   Reads in grid parameters from the specified file.
%             %
%             %   ARGUMENTS:
%             %   index:      user defined class index
%             %   pprovider:  instance of PARAMETER_PROVIDER class
% 			%	forcing:	instance of FORCING class
%             
%             % The following is only needed to allow legacy code to run
%             % May be removed when deprecated functions are removed
%             if nargin==3
%                 index = varargin{1};
%                 pprovider = varargin{2};
%                 forcing = varargin{3};
%             else
%                 st = dbstack;
%                 warning(['DEPRECATION WARNING: Instantiating ' st.name '() with no arguments is deprecated.' newline,...
%                          'You should update your code to take advantage of new IO interface.']);
%                 return
%             end
%             % End allow legacy code
%             
%             self.grid_index = index;
%             self = self.initialize();
%             self = self.populate_GRID(pprovider);
%             self = self.finalize_setup(forcing);
%         end
        

		
        
		function self = provide_PARA(self)
			% INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  
			self.PARA.grid =[];
        end
        
        function self = provide_CONST(self)

        end
        
        function self = provide_STATVAR(self)

        end
        
        %called before TILE exists        
%         function self = initialize_excel(self)
%             self.STATVAR.GRID = [];
%             for i=1:size(self.PARA.grid.upper,1)-1
%                 self.STATVAR.GRID = [self.STATVAR.GRID; [self.PARA.grid.upper(i,1):self.PARA.grid.spacing(i,1):self.PARA.grid.upper(i+1,1)-self.PARA.grid.spacing(i,1)]'];
%             end
%             self.STATVAR.GRID = [self.STATVAR.GRID; [self.PARA.grid.upper(end,1):self.PARA.grid.spacing(end,1):self.PARA.grid.lower(end,1)]'];
%             
%         end
        
       
        %called when TILE exists
        function self = finalize_init(self, tile)
            
            self.STATVAR.GRID = [];
            for i=1:size(self.PARA.grid.upper,1)-1
                self.STATVAR.GRID = [self.STATVAR.GRID; [self.PARA.grid.upper(i,1):self.PARA.grid.spacing(i,1):self.PARA.grid.upper(i+1,1)-self.PARA.grid.spacing(i,1)]'];
            end
            self.STATVAR.GRID = [self.STATVAR.GRID; [self.PARA.grid.upper(end,1):self.PARA.grid.spacing(end,1):self.PARA.grid.lower(end,1)]'];
            
            %---
            
            %delete grid points below domain_depth
            self.STATVAR.GRID(self.STATVAR.GRID > tile.PARA.domain_depth)=[];
            
            self.STATVAR.MIDPOINTS = (self.STATVAR.GRID(2:end,1) + self.STATVAR.GRID(1:end-1,1))./2;
            self.STATVAR.layerThick = (self.STATVAR.GRID(2:end,1) - self.STATVAR.GRID(1:end-1,1));
            self.STATVAR.layerDistance = (self.STATVAR.MIDPOINTS(2:end,1) - self.STATVAR.MIDPOINTS(1:end-1,1))./2;
        end
        
        
        
                
% 		function self = populate_GRID(self, pprovider)
% 			% POPULATE_GRID  Updates the PARA structure with values from pprovider. Assigns values from the PARA structure to the corresponding class properties.
%             %
%             %   ARGUMENTS:
%             %   pprovider:  instance of PARAMETER_PROVIDER class
% 			
%             self.PARA = pprovider.populate_struct(self.PARA, 'GRID', mfilename('class'), self.grid_index);
% 			
% 			fn_substruct = fieldnames(self.PARA.grid);
% 			grid_breaks = [];
% 			for i = 1:size(fn_substruct, 1)
% 				grid_breaks = [grid_breaks self.PARA.grid.(fn_substruct{i})];
% 			end
% 			
% 			grid_breaks = cell2mat(grid_breaks);
%             grid_breaks(2:end,1) = grid_breaks(2:end,1) + grid_breaks(2:end,2);
%             for i=1:size(grid_breaks)
% 				self.GRID = [self.GRID; [grid_breaks(i,1):grid_breaks(i,2):grid_breaks(i,3)]'];
%             end
%         end
		
        
% 		function self = finalize_init(self, forcing)
% 			% FINALIZE_SETUP  Performs all additional property
%             %   initializations and modifications. Checks for some (but not
%             %   all) data validity.
% 			
% 			%	ARGUMENTS:
% 			%	forcing:	instance of FORCING class
% 			
% 			self.GRID(self.GRID > forcing.PARA.domain_depth)=[]; 
%             self.MIDPOINTS = (self.GRID(2:end,1) + self.GRID(1:end-1,1))./2;
%             self.LAYERTHICK = (self.GRID(2:end,1) - self.GRID(1:end-1,1));
%         end

        
        function xls_out = write_excel(self)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'GRID','index',NaN;'GRID_user_defined',1,NaN;NaN,NaN,NaN;'upper','spacing','lower';'[m]','[m]','[m]';'TOP',NaN,NaN;0,0.0500000000000000,2;2,0.100000000000000,10;10,0.500000000000000,30;30,1,50;50,5,100;100,10,150;150,50,500;500,100,1000;'BOTTOM',NaN,NaN;'GRID_END',NaN,NaN};
        end
        

    end
    
end

