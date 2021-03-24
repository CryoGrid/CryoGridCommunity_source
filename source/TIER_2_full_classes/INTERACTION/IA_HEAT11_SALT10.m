%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction between two GROUND
% classes and zero salt flux boundary conditions for GROUND class above 
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_SALT10 <  IA_HEAT & IA_SALT
    
    methods
        
        function get_boundary_condition_m(ia_heat, tile)
            get_boundary_condition_HEAT_m(ia_heat);
            get_boundary_condition_ZEROFLUX_SALT_PREVIOUS_m(ia_heat);
        end
           
    end
end