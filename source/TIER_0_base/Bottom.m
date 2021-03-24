%========================================================================
% CryoGrid Bottom class 
% bottom confining class of stratigraphy
% S. Westermann, Oct 2020
%========================================================================

classdef Bottom  < matlab.mixin.Copyable 
    
    properties
        PREVIOUS
    end
    
    methods
        function obj = init_bottom(obj, bottom_class)
            obj.PREVIOUS = bottom_class;
        end
        
    end
end

