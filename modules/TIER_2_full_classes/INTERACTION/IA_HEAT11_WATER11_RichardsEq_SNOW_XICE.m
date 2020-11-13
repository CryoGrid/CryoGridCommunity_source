%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction and water fluxes between a SNOW class
% and a GROUND class with Richards equation water and excess ice
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_WATER11_RichardsEq_SNOW_XICE < IA_HEAT11_WATER11_SNOW_XICE 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water, tile)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_RichardsEq_Xice_SNOW_m(ia_heat_water);
        end

        
        function get_IA_CHILD_boundary_condition_u(ia_heat_water, tile)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat_water);
            get_boundary_condition_RichardsEq_Xice_SNOW_m(ia_heat_water);
        end
        
%          remove_excessWater_CHILD(ia_heat_water) handled in IA_HEAT11_WATER11_SNOW_XICE
        
        
    end
end