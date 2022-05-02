%========================================================================
% CryoGrid TLE_BUILDER class update_forcing_out
% CryoGrid TLE_BUILDER class to be used in the context of sequential 
% simulations of several TILE classes after each other. Reads the final 
% state of the previous TILE class for initialization of the CryoGrid
% stratigraphy, but overwrites FORCING and OUT classes with user-rpovided
% ones.
% Note that the code (i.e. what the TILE_BUILDER class actually does in in
% the respective TILE classes which the TILE_BUILDER class is compatible
% with
% S. Westermann, Jan 2021
%========================================================================

classdef update_forcing_out

    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_update_forcing_out(builder.TILE);            
        end
% 
%         function builder = provide_PARA(builder)
%             
%             builder.PARA.forcing_class = [];
%             builder.PARA.forcing_class_index = [];
%             
%             builder.PARA.out_class = [];
%             builder.PARA.out_class_index = [];
%             
%         end
    end
end

