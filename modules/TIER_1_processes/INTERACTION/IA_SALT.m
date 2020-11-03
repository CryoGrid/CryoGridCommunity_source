%========================================================================
% CryoGrid TIER1 INTERACTION (IA) class for functions related to salt
% diffusion
% NOTE: at this pooint, only zero flux boundary conditions are implemented
% S. Westermann, October 2020
%========================================================================

classdef IA_SALT < IA_BASE
    
    methods
        
        function get_boundary_condition_ZEROFLUX_SALT_NEXT_m(ia_heat_water) %coupling between classes without (PREVIOUS) and with (NEXT) salt balance
            ia_heat_water.NEXT.TEMP.F_ub_salt = 0;
            ia_heat_water.NEXT.TEMP.d_salt(1) = ia_heat_water.NEXT.TEMP.d_salt(1) + 0;
        end
        
        function get_boundary_condition_ZEROFLUX_SALT_PREVIOUS_m(ia_heat_water) %coupling between classes with (PREVIOUS) and without (NEXT) salt balance
            ia_heat_water.PREVIOUS.TEMP.F_lb_salt = 0;
            ia_heat_water.PREVIOUS.TEMP.d_salt(end) = ia_heat_water.PREVIOUS.TEMP.d_salt(end) + 0;
        end
  
    end
end


           
        
      