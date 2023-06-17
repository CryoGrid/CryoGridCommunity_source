%========================================================================
% CryoGrid ASSIGN_TILE_PROPERTIES class tag_out_w_run_number
% add the run_number as tag to the OUT class
%
% S. Westermann, Dec 2022
%========================================================================

classdef tag_out_w_run_number < matlab.mixin.Copyable

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
        
        function update = assign_tile_properties(update, run_number)
            for i=1:size(update.PARA.class_name,1)
                update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.class_name{i,1}){update.PARA.class_index(i,1),1}.PARA.tag = num2str(run_number);
            end
       
        end
        
        
        %-------------param file generation-----
        function update = param_file_info(update)
            update = provide_PARA(update);

            update.PARA.STATVAR = [];
            update.PARA.class_category = 'ASSIGN_TILE_PROPERTIES';
            update.PARA.default_value = [];
            
            update.PARA.comment.class_name =  {'name of OUT class in which tag is to be added'};
            update.PARA.options.class_name.name =  'H_LIST';
            
            update.PARA.comment.class_index =  {'index of OUT class in which tag is to be added'};
            update.PARA.options.class_index.name =  'H_LIST';

        end
    end
end

