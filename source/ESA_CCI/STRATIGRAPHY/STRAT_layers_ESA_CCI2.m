

classdef STRAT_layers_ESA_CCI2 < matlab.mixin.Copyable
    
    properties
        
		PARA
        CONST
        STATVAR
    end
    
    methods
        

		function strat = provide_PARA(strat) 
            strat.PARA.stratigraphy_index = [];
            strat.PARA.stratigraphy_fraction = [];
                        
            strat.PARA.strat_statvar_class = []; %no H_LIST, just a single class
            strat.PARA.strat_statvar_class_index_start = [];
            strat.PARA.strat_statvar_class_index_end = [];
        end
        
        function strat = provide_CONST(strat)

        end
        
        function strat = provide_STATVAR(strat)

        end 
        
        
        function strat = finalize_init(strat, tile)
			
            %make list of ID's so that strat.PARA.stratigraphy_index can be
            %matched
            class_index_list = [];
            for j = strat.PARA.strat_statvar_class_index_start:strat.PARA.strat_statvar_class_index_end
                class_index_list = [class_index_list; tile.RUN_INFO.PPROVIDER.CLASSES.(strat.PARA.strat_statvar_class){j,1}.PARA.ID];
            end
            
            %first grid cell/ensemble member
            tile.GRID.STATVAR.MIDPOINTS = tile.GRID.STATVAR.midPoints(:,1);
            tile.GRID.STATVAR.GRID = tile.GRID.STATVAR.depth(:,1);
            class_index = find(strat.PARA.stratigraphy_index(1,1) == class_index_list(:,1));
            strat_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(strat.PARA.strat_statvar_class){class_index,1});
            strat_class.PARA.stratigraphy_fraction = strat.PARA.stratigraphy_fraction(1,1); %overwrite local variable
            strat_class = finalize_init(strat_class, tile);
            variables = fieldnames(strat_class.STATVAR);
            for j = 1:size(variables,1)
                tile.GRID.STATVAR.(variables{j,1}) = strat_class.STATVAR.(variables{j,1});
            end

            for i=2:size(strat.PARA.stratigraphy_index,2)
                tile.GRID.STATVAR.MIDPOINTS = tile.GRID.STATVAR.midPoints(:,i);
                tile.GRID.STATVAR.GRID = tile.GRID.STATVAR.depth(:,i);
                class_index = find(strat.PARA.stratigraphy_index(1,i) == class_index_list(:,1));
                strat_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(strat.PARA.strat_statvar_class){class_index,1});
                strat_class.PARA.stratigraphy_fraction = strat.PARA.stratigraphy_fraction(1,i); %overwrite local variable
                strat_class = finalize_init(strat_class, tile);
                for j = 1:size(variables,1)
                    tile.GRID.STATVAR.(variables{j,1}) = [tile.GRID.STATVAR.(variables{j,1}) strat_class.STATVAR.(variables{j,1})];
                end
            end

        end
        
    end
    
end