classdef IA_HEAT11_WATER11_RichardsEq_SNOW_XICE < IA_HEAT11_WATER11_SNOW_XICE 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_RichardsEq_Xice_SNOW_m(ia_heat_water);
        end

        
        function get_IA_CHILD_boundary_condition_u(ia_heat_water)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat_water);
            get_boundary_condition_RichardsEq_Xice_SNOW_m(ia_heat_water);
        end
        
%          remove_excessWater_CHILD(ia_heat_water) handled in IA_HEAT11_WATER11_SNOW_XICE
        
        
    end
end