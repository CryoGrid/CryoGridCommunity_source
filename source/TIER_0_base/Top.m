%========================================================================
% CryoGrid TOP class 
% top confining class of stratigraphy, storage for sleeping GROUND classes, pointer to lateral class 
% S. Westermann, Oct 2020
%========================================================================

classdef Top  < matlab.mixin.Copyable 

    properties
        NEXT
        IA_NEXT
        STORE
        LATERAL
%         FORCING
%         TIME
%         RUN_NUMBER
%         RESULT_PATH
    end
    
    methods
        function obj = init_top(obj, top_class)
            obj.NEXT = top_class;
        end
        
    end
end

