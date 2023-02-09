%========================================================================
% CryoGrid ASSIGN_TILE_PROPERTIES class update_one2one
% updates variables that have the same name in the CryoGrid class as
% provided by the DATA_PROVODER class
%
% S. Westermann, Dec 2022
%========================================================================

classdef assign_range < matlab.mixin.Copyable

    properties
        PARA
        CONST
        PROJ
    end
    
    methods
        function update = provide_PARA(update)

            update.PARA.class_name = [];
            update.PARA.class_index = [];
            
        end
        
        function update = provide_STATVAR(update)

        end
        
        function update = provide_CONST(update)
            
        end
        
        function update = finalize_init(update)
 
        end
        
        function update = assign_tile_properties(update, range)
            for i=1:size(update.PARA.class_name,1)
                update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.class_name{i,1}){update.PARA.class_index(i,1),1}.PARA.range = range;
                 update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.class_name{i,1}){update.PARA.class_index(i,1),1}.PARA.number_of_realizations = size(range,1);
            end
        end
        
        %-------------param file generation-----
        function update = param_file_info(update)
            update = provide_PARA(update);

            update.PARA.STATVAR = [];
            update.PARA.class_category = 'ASSIGN_TILE_PROPERTIES';
            update.PARA.default_value = [];
            
            update.PARA.comment.class_name =  {'name of class in which the variable is changed'};
            update.PARA.options.class_name.name =  'H_LIST';
            update.PARA.options.class_name.entries_x = {'TILE_1D_standard' 'TILE_1D_standard' 'TILE_1D_standard'};
            
            update.PARA.comment.class_index =  {'index of class in which the variable is changed'};
            update.PARA.options.class_index.name =  'H_LIST';
            update.PARA.options.class_index.entries_x = {'1' '1' '1'};

            update.PARA.comment.variable =  {'name of variable which is changed (must be assigned by DATA PROVIDER class)'};
            update.PARA.options.variable.name =  'H_LIST';
            update.PARA.options.variable.entries_x = {'latitude' 'longitude' 'altitude'};

        end
        
    end
end

