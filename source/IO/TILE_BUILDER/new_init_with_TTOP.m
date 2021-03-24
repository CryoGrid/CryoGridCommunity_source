classdef new_init_with_TTOP
    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_new_init_TTOP_GlobPermafrost(builder.TILE);            
        end
        
    end
end

