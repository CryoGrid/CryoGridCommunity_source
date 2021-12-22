classdef restart_from_OUT_update_classes
    
    properties
        TILE
    end
    
    methods
        
        function build_tile(builder)
            builder.TILE = build_tile_restart_OUT_update_classes(builder.TILE);            
        end

    end
end

