%========================================================================
% CryoGrid TLE_BUILDER class new_init
% CryoGrid TLE_BUILDER class which accomplishes the initialization of the
% TILE class based only on user-defined profiles of mode state avraibles. 
% Note that the code (i.e. what the TILE_BUILDER class actually does in in
% the respective TILE classes which the TILE_BUILDER class is compatible
% with
% S. Westermann, Jan 2021
%========================================================================


classdef new_init
    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_new_init(builder.TILE);            
        end
        
%         function builder = provide_PARA(builder)
%             
%             builder.PARA.latitude = [];
%             builder.PARA.longitude = [];
%             builder.PARA.altitude = [];
%             builder.PARA.domain_depth = [];
%             builder.PARA.area = [];
%             
%             builder.PARA.forcing_class = [];
%             builder.PARA.forcing_class_index = [];
%             builder.PARA.grid_class = [];
%             builder.PARA.grid_class_index = [];
%             builder.PARA.out_class = [];
%             builder.PARA.out_class_index = [];
%             builder.PARA.strat_classes_class = [];
%             builder.PARA.strat_classes_class_index = [];
%             builder.PARA.strat_statvar_class = [];
%             builder.PARA.strat_statvar_class_index = [];
%             builder.PARA.lateral_class = [];
%             builder.PARA.lateral_class_index = [];
%             builder.PARA.lateral_IA_classes = [];
%             builder.PARA.lateral_IA_classes_index = [];
%         end
    end
end

