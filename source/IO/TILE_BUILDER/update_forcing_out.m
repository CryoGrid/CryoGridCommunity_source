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

