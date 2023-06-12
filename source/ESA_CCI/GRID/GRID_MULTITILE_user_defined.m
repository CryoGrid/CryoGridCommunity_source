%========================================================================
% CryoGrid GRID class  for defining the compute grid
% GRID_user_defined defines the compute grid as layers with constant grid spacing
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef GRID_MULTITILE_user_defined < matlab.mixin.Copyable

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
                        
            %repeat for all cells/ensemble members
            self.STATVAR.depth = repmat(self.STATVAR.GRID, 1, tile.PARA.tile_size);
            self.STATVAR.midPoints = repmat(self.STATVAR.MIDPOINTS, 1, tile.PARA.tile_size);
            self.STATVAR.layerThick = repmat(self.STATVAR.layerThick, 1, tile.PARA.tile_size);
            self.STATVAR.layerDistance = repmat(self.STATVAR.layerDistance, 1, tile.PARA.tile_size);
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
                

    end
    
end

