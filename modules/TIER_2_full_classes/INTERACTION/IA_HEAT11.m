%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction between two GROUND
% classes without water cycle
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11 <  IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat, tile)
            get_boundary_condition_HEAT_m(ia_heat);
        end
        
        %SNOW
        function get_IA_CHILD_boundary_condition_u(ia_heat, tile)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat);
        end
           
    end
end