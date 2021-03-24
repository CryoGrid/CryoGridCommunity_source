classdef update_forcing_out

    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_update_forcing_out(builder.TILE);            
        end
    end
end

