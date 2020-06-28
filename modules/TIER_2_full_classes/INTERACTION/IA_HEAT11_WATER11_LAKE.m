classdef IA_HEAT11_WATER11_LAKE < IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water)
            get_boundary_condition_HEAT_LAKE_m(ia_heat_water);
            get_boundary_condition_BUCKET_LAKE_m(ia_heat_water);
        end
           
    end
end