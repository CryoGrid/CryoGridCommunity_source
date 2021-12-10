%========================================================================
% CryoGrid INTERACTION (IA) class for multiTile class above a normal class
% S. Westermann, October 2020
%========================================================================

classdef IA_MULTI_TILE10 < IA_BASE
    
    properties
       IA_MULTI_TILE 
    end
    
    
    methods
        
        function finalize_init(ia_multiTile, tile)
            for i=1:size(ia_multiTile.PREVIOUS.STATVAR.SUB_TILES_BOTTOM,1)
                ia_multiTile.IA_MULTI_TILE{i,1} = get_IA_class(class(ia_multiTile.PREVIOUS.STATVAR.SUB_TILES_BOTTOM{i,1}.PREVIOUS), class(ia_multiTile.NEXT));
                ia_multiTile.PREVIOUS.STATVAR.SUB_TILES_BOTTOM{i,1}.PREVIOUS.IA_NEXT = ia_multiTile.IA_MULTI_TILE{i,1};
                ia_multiTile.IA_MULTI_TILE{i,1}.PREVIOUS = ia_multiTile.PREVIOUS.STATVAR.SUB_TILES_BOTTOM{i,1}.PREVIOUS;
                ia_multiTile.IA_MULTI_TILE{i,1}.NEXT = ia_multiTile.NEXT;
            end
        end
        
        function get_boundary_condition_m(ia_multiTile, tile)
           for i = 1:size(ia_multiTile.IA_MULTI_TILE,1)
               ia_multiTile.NEXT.IA_PREVIOUS = ia_multiTile.IA_MULTI_TILE{i,1};
               get_boundary_condition_m(ia_multiTile.IA_MULTI_TILE{i,1}, tile);
           end
           ia_multiTile.NEXT.IA_PREVIOUS = ia_multiTile;
        end
        
    end
end