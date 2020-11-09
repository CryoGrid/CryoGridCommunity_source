%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction and water fluxes between a SNOW class
% and a GROUND class with bucket water
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_WATER11_SNOW < IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water, tile)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_m(ia_heat_water);
        end
        
        %SNOW
        function get_IA_CHILD_boundary_condition_u(ia_heat_water, tile)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_m(ia_heat_water);
        end

    end
end