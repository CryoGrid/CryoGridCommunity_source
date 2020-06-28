classdef IA_HEAT11 <  IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat)
            get_boundary_condition_HEAT_m(ia_heat);
        end
        
        %SNOW
        function get_IA_CHILD_boundary_condition_u(ia_heat)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat);
        end
           
    end
end