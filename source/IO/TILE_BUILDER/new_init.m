classdef new_init
    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_new_init(builder.TILE);            
        end
    end
end

