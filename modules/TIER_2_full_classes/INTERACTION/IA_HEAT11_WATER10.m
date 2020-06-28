classdef IA_HEAT11_WATER10 < IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_ZEROFLUX_PREVIOUS_m(ia_heat_water);
        end
        
        function get_IA_CHILD_boundary_condition_u(ia_heat_water)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat_water);
            get_boundary_condition_ZEROFLUX_PREVIOUS_m(ia_heat_water);
        end
        
        
    end
end