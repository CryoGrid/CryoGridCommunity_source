%========================================================================
% CryoGrid TILE_BUILDER class new_init_steady_state
% CryoGrid TLE_BUILDER class used to initialize the TILE class with a
% steady-state temperature profile computed from the geothermal heat flux 
% and a give temperature at a defined depth, as well as profiles of soil 
% properties provieded by the user. The temperature/depth can either be 
% provided by the user, or calculated by an INIT_STEADY_STATE class, e.g. 
% using the model forcing or the output of a previous simulation. 
% Note that the code (i.e. what the TILE_BUILDER class actually does in in
% the respective TILE classes which the TILE_BUILDER class is compatible
% with
% S. Westermann, Aug 2021
%========================================================================

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

