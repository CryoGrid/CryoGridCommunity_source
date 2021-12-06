classdef new_init_steady_state
    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_new_init_steady_state(builder.TILE);            
        end
        
    end
end

