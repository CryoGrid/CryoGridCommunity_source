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

		function self = provide_PARA(self)
			self.PARA.grid =[];
        end
        
        function self = provide_CONST(self)

        end
        
        function self = provide_STATVAR(self)

        end
        
       
        %called when TILE exists
        function self = finalize_init(self, tile)
            
            self.STATVAR.GRID = [];
            for i=1:size(self.PARA.grid.upper,1)-1
                self.STATVAR.GRID = [self.STATVAR.GRID; [self.PARA.grid.upper(i,1):self.PARA.grid.spacing(i,1):self.PARA.grid.upper(i+1,1)-self.PARA.grid.spacing(i,1)]'];
            end
            self.STATVAR.GRID = [self.STATVAR.GRID; [self.PARA.grid.upper(end,1):self.PARA.grid.spacing(end,1):self.PARA.grid.lower(end,1)]'];

            %delete grid points below domain_depth
            self.STATVAR.GRID(self.STATVAR.GRID > tile.PARA.domain_depth)=[];
            
            self.STATVAR.MIDPOINTS = (self.STATVAR.GRID(2:end,1) + self.STATVAR.GRID(1:end-1,1))./2;
            self.STATVAR.layerThick = (self.STATVAR.GRID(2:end,1) - self.STATVAR.GRID(1:end-1,1));
            self.STATVAR.layerDistance = (self.STATVAR.MIDPOINTS(2:end,1) - self.STATVAR.MIDPOINTS(1:end-1,1))./2;
        end
        

        function self = finalize_init_GROUND_multi_tile(self, tile)
            
            self.STATVAR.GRID = [];
            for i=1:size(self.PARA.grid.upper,1)-1
                self.STATVAR.GRID = [self.STATVAR.GRID; [self.PARA.grid.upper(i,1):self.PARA.grid.spacing(i,1):self.PARA.grid.upper(i+1,1)-self.PARA.grid.spacing(i,1)]'];
            end
            self.STATVAR.GRID = [self.STATVAR.GRID; [self.PARA.grid.upper(end,1):self.PARA.grid.spacing(end,1):self.PARA.grid.lower(end,1)]'];

            %delete grid points below domain_depth            
            self.STATVAR.MIDPOINTS = (self.STATVAR.GRID(2:end,1) + self.STATVAR.GRID(1:end-1,1))./2;
            self.STATVAR.layerThick = (self.STATVAR.GRID(2:end,1) - self.STATVAR.GRID(1:end-1,1));
            self.STATVAR.layerDistance = (self.STATVAR.MIDPOINTS(2:end,1) - self.STATVAR.MIDPOINTS(1:end-1,1))./2;
        end

        
        %-------------param file generation-----
         function stratigraphy = param_file_info(stratigraphy)
             stratigraphy = provide_PARA(stratigraphy);
             %default
             stratigraphy.PARA.default_value = [];
             stratigraphy.PARA.STATVAR = [];
             stratigraphy.PARA.comment = [];
             stratigraphy.PARA.class_category = 'GRID';
             stratigraphy.PARA.options.grid.name = 'V_MATRIX';
             stratigraphy.PARA.options.grid.entries_x = {'upper' 'spacing' 'lower'}; %fill with STATVAR that are identified for initialization in GROUND classes
             stratigraphy.PARA.options.grid.entries_matrix = {'0' '5e-2' '1'; '1' '1e-1' '2'; '2' '2e-1' '5'; '5' '5e-1' '30'; '30' '5' '100'};
 
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

        
%         function xls_out = write_excel(self)
% 			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
% 			
%             xls_out = {'GRID','index',NaN;'GRID_user_defined',1,NaN;NaN,NaN,NaN;'upper','spacing','lower';'[m]','[m]','[m]';'TOP',NaN,NaN;0,0.0500000000000000,2;2,0.100000000000000,10;10,0.500000000000000,30;30,1,50;50,5,100;100,10,150;150,50,500;500,100,1000;'BOTTOM',NaN,NaN;'GRID_END',NaN,NaN};
%         end
        

    end
    
end

