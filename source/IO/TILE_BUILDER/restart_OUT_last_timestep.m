classdef restart_OUT_last_timestep
    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_restart_OUT_last_timestep(builder.TILE);            
        end
    end
end

