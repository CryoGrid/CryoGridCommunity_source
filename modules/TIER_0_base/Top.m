classdef Top  < matlab.mixin.Copyable 

    properties
        NEXT
        STORE
        LATERAL
    end
    
    methods
        function obj = init_top(obj, top_class)
            obj.NEXT = top_class;
        end
        
    end
end

