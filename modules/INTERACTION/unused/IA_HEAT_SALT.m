classdef IA_HEAT_SALT < IA_HEAT 
    %transition from salty to non-salty class, zero flux boundary
    %condition imposed at both boundaries, this should be changd depnding on context of use

    

    methods
        
        function get_boundary_condition_m(ia_heat_salt)
            get_boundary_condition_m@IA_HEAT(ia_heat_salt)
            
            stratigraphy1 = ia_heat_salt.PREVIOUS;
            stratigraphy2 = ia_heat_salt.NEXT;
            
            stratigraphy1.TEMP.F_lb_salt = 0;
            stratigraphy2.TEMP.F_ub_salt = 0; 
        end
    end
end
   
   
