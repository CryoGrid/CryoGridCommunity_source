%defines a regular grid in geographical coordinates, with fixed resolution

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
    end
end

