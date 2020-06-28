classdef BASE_LATERAL < matlab.mixin.Copyable
    
    properties
        CONST %constants
        PARA %external service parameters, all other
        STATVAR  %energy, water content, etc.
        TEMP  %derivatives in prognostic timestep and optimal timestep
        PARENT
    end
    
end

