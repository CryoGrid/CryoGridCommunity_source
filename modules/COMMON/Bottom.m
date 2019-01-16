classdef Bottom  < matlab.mixin.Copyable 
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        PREVIOUS
    end
    
    methods
        function obj = init_bottom(obj, bottom_class)
            obj.PREVIOUS = bottom_class;
        end
        
    end
end

