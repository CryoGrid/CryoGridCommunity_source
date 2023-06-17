

classdef STRAT_ensemble_two_confining_stratigraphies < matlab.mixin.Copyable
    
    properties
        
		PARA
        CONST
        STATVAR
    end
    
    methods
        

		function strat = provide_PARA(strat) 
            strat.PARA.stratigraphy_fraction = [];
                        
          %  strat.PARA.ensemble_size = 1;
            strat.PARA.strat_statvar_class1 = []; 
            strat.PARA.strat_statvar_class1_index = [];
            strat.PARA.strat_statvar_class2 = [];
            strat.PARA.strat_statvar_class2_index = [];
        end
        
        function strat = provide_CONST(strat)

        end
        
        function strat = provide_STATVAR(strat)

        end 
        
        
        function strat = finalize_init(strat, tile)
			
%             if size(strat.PARA.stratigraphy_fraction,2) == 1
%                 strat.PARA.stratigraphy_fraction = repmat(strat.PARA.stratigraphy_fraction, 1, strat.PARA.ensemble_size);
%             end
            
            %first grid cell/ensemble member
%             tile.GRID.STATVAR.MIDPOINTS = tile.GRID.STATVAR.midPoints(:,1);
%             tile.GRID.STATVAR.GRID = tile.GRID.STATVAR.depth(:,1);
            strat_class1 = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(strat.PARA.strat_statvar_class1){strat.PARA.strat_statvar_class1_index,1});
            strat_class2 = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(strat.PARA.strat_statvar_class2){strat.PARA.strat_statvar_class2_index,1});
            strat_class1 = finalize_init(strat_class1, tile);
            strat_class2 = finalize_init(strat_class2, tile);            
            variables = fieldnames(strat_class1.STATVAR);
            for j = 1:size(variables,1)
                tile.GRID.STATVAR.(variables{j,1}) = strat.PARA.stratigraphy_fraction(1,1) .* strat_class1.STATVAR.(variables{j,1}) + ...
                    (1-strat.PARA.stratigraphy_fraction(1,1)) .* strat_class2.STATVAR.(variables{j,1});
                if strcmp(variables{j,1}, 'soil_type')
                    tile.GRID.STATVAR.(variables{j,1}) = round(tile.GRID.STATVAR.(variables{j,1}));
                end
            end

%             for i=2:strat.PARA.ensemble_size
%                 tile.GRID.STATVAR.MIDPOINTS = tile.GRID.STATVAR.midPoints(:,i);
%                 tile.GRID.STATVAR.GRID = tile.GRID.STATVAR.depth(:,i);
%                 strat_class1 = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(strat.PARA.strat_statvar_class1){strat.PARA.strat_statvar_class1_index,1});
%                 strat_class2 = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(strat.PARA.strat_statvar_class2){strat.PARA.strat_statvar_class2_index,1});
%                 strat_class1 = finalize_init(strat_class1, tile);
%                 strat_class2 = finalize_init(strat_class2, tile);
%                 for j = 1:size(variables,1)
%                     tile.GRID.STATVAR.(variables{j,1}) = [tile.GRID.STATVAR.(variables{j,1}) strat.PARA.stratigraphy_fraction(1,i) .* strat_class1.STATVAR.(variables{j,1}) + ...
%                     (1-strat.PARA.stratigraphy_fraction(1,i)) .* strat_class2.STATVAR.(variables{j,1})];
%                     if strcmp(variables{j,1}, 'soil_type')
%                         tile.GRID.STATVAR.(variables{j,1}) = round(tile.GRID.STATVAR.(variables{j,1}));
%                     end
%                 end
%             end

        end
        
    end
    
end