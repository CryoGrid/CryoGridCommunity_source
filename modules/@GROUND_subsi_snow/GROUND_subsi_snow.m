classdef GROUND_subsi_snow < GROUND_base_class
    %GROUND_SUBSIDENCE Summary of this class goes here
    %   Detailed explanation goes here
    methods
        
        %mandatory functions for each class

        function snow = provide_variables(snow)  %initializes the subvariables as empty arrays
            
            snow = provide_variables@GROUND_base_class(snow); %call function of the base class
            snow = provide_PARA(snow); %add additional variables
%             ground = provide_CONST(ground);
%             ground = provide_STATVAR(ground);
        end
    end
end

