classdef Top  < matlab.mixin.Copyable 
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NEXT
    end
    
    methods
        function obj = init_top(obj, top_class)
            obj.NEXT = top_class;
        end
        
    end
end

