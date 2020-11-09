%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction and water fluxes between a LAKE class
% and a GROUND class with Richards equation water and excess ice
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_WATER11_RichardsEq_LAKE_XICE < IA_HEAT11_WATER11_LAKE_XICE 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water, tile)
            get_boundary_condition_HEAT_LAKE_m(ia_heat_water);
            get_boundary_condition_RichardsEq_Xice_LAKE_m(ia_heat_water);
        end
        
        %trigger functions handled in IA_HEAT11_WATER11_LAKE_XICE
                
    end
end