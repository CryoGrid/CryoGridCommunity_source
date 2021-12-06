%========================================================================
% CryoGrid BASE class for all LATERAL classes 
% contains variables and initialization routines
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, Oct 2020
%========================================================================


classdef BASE_LATERAL < matlab.mixin.Copyable
    
    properties
        class_index
        CONST    %constants
        PARA     %external service parameters, all other
        STATVAR  %energy, water content, etc.
        TEMP     %derivatives in prognostic timestep and optimal timestep
        PARENT
    end
    
    methods
        

    end
end

