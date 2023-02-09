

classdef STRAT_layers_ESA_CCI < matlab.mixin.Copyable
    
    properties
        
		PARA
        CONST
        STATVAR
    end
    
    methods
        

		function strat = provide_PARA(strat) 

			strat.PARA.ID = [];
            strat.PARA.name = [];
            
			strat.PARA.layers_organic_stratigraphy = [];
            strat.PARA.layers_mineral_stratigraphy = [];
        end
        
        function strat = provide_CONST(strat)

        end
        
        function strat = provide_STATVAR(strat)

        end 
        
        
        function strat = finalize_init(strat, tile)
			
            variables = fieldnames(strat.PARA.layers_organic_stratigraphy);
            depth1 = strat.PARA.layers_organic_stratigraphy.depth;
            depth1 = [depth1; Inf];
            depth2 = strat.PARA.layers_mineral_stratigraphy.depth;
            depth2 = [depth2; Inf];
            
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'depth')
%                     if ~isfield(tile.GRID.STATVAR, variables{i,1})
%                         tile.GRID.STATVAR.(variables{i,1}) = tile.GRID.STATVAR.layerThick .* 0;
%                     end
                    store_var1 = tile.GRID.STATVAR.GRID(1:end-1,1) .*0;
                    store_var2 = store_var1;
                    for j=1:size(depth1,1)-1
                        range = tile.GRID.STATVAR.MIDPOINTS > depth1(j,1) & tile.GRID.STATVAR.MIDPOINTS <= depth1(j+1,1);
                        store_var1(range,1) = strat.PARA.layers_organic_stratigraphy.(variables{i,1})(j,1);
                        range = tile.GRID.STATVAR.MIDPOINTS > depth2(j,1) & tile.GRID.STATVAR.MIDPOINTS <= depth2(j+1,1);
                        store_var2(range,1) = strat.PARA.layers_mineral_stratigraphy.(variables{i,1})(j,1);
                    end
                    strat.STATVAR.(variables{i,1}) = [store_var1 store_var2];
                   % strat.STATVAR.depth = [depth1(1:end-1,1) depth2(1:end-1,1)];
                    
%                     ind = find(tile.PARA.stratigraphy >= strat.PARA.class_number);
%                     store_var = repmat(store_var, 1, tile.PARA.number_of_realizations);
%                     tile.GRID.STATVAR.(variables{i,1})(:,ind) = store_var(:,ind);
                end
            end
            
        end
        
            

        
    end
    
end