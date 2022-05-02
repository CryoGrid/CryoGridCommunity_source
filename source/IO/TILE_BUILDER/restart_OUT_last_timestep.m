%========================================================================
% CryoGrid TILE_BUILDER class restart_OUT_last_timestep
% CryoGrid TLE_BUILDER class used to initialize the TILE class with a
% CryoGrid stratigraphy written by the OUT class OUT_last_timestep. In 
% order to change the model FORCING or OUT classes, use in sequential 
% simulation of TILE classes and use TILE_BUILDER update_forcing_out for 
% the following TILE class.
% Note that the code (i.e. what the TILE_BUILDER class actually does in in
% the respective TILE classes which the TILE_BUILDER class is compatible
% with
% S. Westermann, Jan 2021
%========================================================================


classdef restart_OUT_last_timestep
    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_restart_OUT_last_timestep(builder.TILE);            
        end
%         function builder = provide_PARA(builder)
%             builder.PARA.restart_file_path = [];
%             builder.PARA.restart_file_name = [];
%         end
    end
end

