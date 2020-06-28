classdef IA_HEAT11_WATER10_LAKE <  IA_HEAT & IA_WATER
    
    methods
        
        function get_boundary_condition_m(ia_heat_water)
            get_boundary_condition_HEAT_LAKE_m(ia_heat_water);
            get_boundary_condition_ZEROFLUX_PREVIOUS_m(ia_heat_water);
        end
        
        %----------
        
        function  get_boundary_condition_HEAT_LAKE_m(ia_heat_water)
            get_boundary_condition_HEAT_LAKE_m@IA_HEAT(ia_heat_water);
        end
        
        function get_boundary_condition_ZEROFLUX_PREVIOUS_m(ia_heat_water)
             get_boundary_condition_ZEROFLUX_PREVIOUS_m@IA_WATER(ia_heat_water);
        end
           
    end
end