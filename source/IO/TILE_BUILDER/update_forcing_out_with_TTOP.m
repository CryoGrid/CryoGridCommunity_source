classdef update_forcing_out_with_TTOP

    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_update_forcing_out_TTOP(builder.TILE);            
        end
    end
end

