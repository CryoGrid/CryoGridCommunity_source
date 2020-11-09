classdef TILE

    
    properties
        FORCING
        CONST
        t
        timestep
        next_break_time
        LATERAL
        TOP
        BOTTOM
        RUN_NUMBER
        RESULT_PATH
        OUT
    end
    
    methods

        function tile = interpolate_forcing_tile(tile)
             tile.FORCING = interpolate_forcing(tile.t, tile.FORCING);
        end

        function tile = interact_lateral(tile)
            tile.LATERAL = interact(tile.LATERAL, tile);
        end
        
        function tile = store_OUT_tile(tile)
            %tile.OUT = store_OUT(tile.OUT, tile.t, tile.TOP, tile.BOTTOM, tile.FORCING, tile.RUN_NUMBER, tile.timestep, tile.RESULT_PATH,);
            tile.OUT = store_OUT(tile.OUT, tile);
        end
        
    end
end

