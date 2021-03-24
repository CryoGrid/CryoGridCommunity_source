%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction and water fluxes between a SNOW class
% and a LAKE class 
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_WATER11_SNOW_LAKE < IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water, tile)
            get_boundary_condition_HEAT_m(ia_heat_water);
            %MODIFY
            %ia_heat_water.PREVIOUS.TEMP.F_lb_water = 0;
            %ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = 0;
            %ia_heat_water.NEXT.TEMP.F_ub_water = 0;
            %ia_heat_water.NEXT.TEMP.F_ub_water_energy = 0;
            %get_boundary_condition_BUCKET_SNOW_m(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_LAKE_m(ia_heat_water);
        end
        
        %SNOW
        function get_IA_CHILD_boundary_condition_u(ia_heat_water, tile)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat_water);
            %ia_heat_water.PREVIOUS.TEMP.F_lb_water = 0;
            %ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = 0;
            %ia_heat_water.NEXT.TEMP.F_ub_water = 0;
            %ia_heat_water.NEXT.TEMP.F_ub_water_energy = 0;
            %get_boundary_condition_BUCKET_SNOW_m(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_LAKE_m(ia_heat_water);
        end
        
        function remove_excessWater_CHILD(ia_heat_water) %move excessWater from SNOW to water of PARENT, 0 energy transfer since meltwater must be zero
            ia_heat_water.NEXT.STATVAR.waterIce(1) = ia_heat_water.NEXT.STATVAR.waterIce(1) + ia_heat_water.PREVIOUS.STATVAR.excessWater;
            ia_heat_water.NEXT.STATVAR.water(1) = ia_heat_water.NEXT.STATVAR.water(1) + ia_heat_water.PREVIOUS.STATVAR.excessWater;
            ia_heat_water.NEXT.STATVAR.layerThick(1) = ia_heat_water.NEXT.STATVAR.layerThick(1) + ia_heat_water.PREVIOUS.STATVAR.excessWater ./ ia_heat_water.NEXT.STATVAR.area(1);
            ia_heat_water.PREVIOUS.STATVAR.excessWater = 0;
        end

    end
end