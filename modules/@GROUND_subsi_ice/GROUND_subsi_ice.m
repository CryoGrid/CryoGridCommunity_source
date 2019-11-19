classdef GROUND_subsi_ice < GROUND_base_class
    %GROUND_SUBSIDENCE Summary of this class goes here
    %   Detailed explanation goes here
    methods
        
        %mandatory functions for each class

        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            
            ground = provide_variables@GROUND_base_class(ground); %call function of the base class
            ground = provide_PARA(ground); %add additional variables
%             ground = provide_CONST(ground);
%             ground = provide_STATVAR(ground);
        end
    end
end

