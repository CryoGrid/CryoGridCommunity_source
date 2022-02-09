%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction and water fluxes between a SNOW class
% and a GROUND class with bucket water and excess ice
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_WATER11_SNOW_XICE < IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water, tile)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_XICE_m(ia_heat_water);
        end
        
        %SNOW
        function get_IA_CHILD_boundary_condition_u(ia_heat_water, tile)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_XICE_m(ia_heat_water);
        end
        
        function remove_excessWater_CHILD(ia_heat_water) %move excessWater from SNOW to water and Xwater from PARENT, 0 energy transfer since meltwater must be zero
            
            space_left = max(0,ia_heat_water.NEXT.STATVAR.layerThick(1) .* ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.mineral(1) ...
                - ia_heat_water.NEXT.STATVAR.organic(1) - ia_heat_water.NEXT.STATVAR.waterIce(1) - ia_heat_water.NEXT.STATVAR.XwaterIce(1)); 
            water_in = min(space_left, ia_heat_water.PREVIOUS.STATVAR.excessWater);
            Xwater_in = max(0, ia_heat_water.PREVIOUS.STATVAR.excessWater - water_in);
            
            ia_heat_water.NEXT.STATVAR.waterIce(1) = ia_heat_water.NEXT.STATVAR.waterIce(1) + water_in;
            ia_heat_water.NEXT.STATVAR.XwaterIce(1) = ia_heat_water.NEXT.STATVAR.XwaterIce(1) + Xwater_in;
            ia_heat_water.NEXT.STATVAR.Xwater(1) = ia_heat_water.NEXT.STATVAR.Xwater(1) + Xwater_in;
            ia_heat_water.NEXT.STATVAR.layerThick(1) = ia_heat_water.NEXT.STATVAR.layerThick(1) + Xwater_in ./ ia_heat_water.NEXT.STATVAR.area(1);
            ia_heat_water.NEXT.STATVAR.energy(1,1) = ia_heat_water.NEXT.STATVAR.energy(1,1) + ia_heat_water.PREVIOUS.STATVAR.excessWater_energy;
            ia_heat_water.PREVIOUS.STATVAR.excessWater = 0;
            ia_heat_water.PREVIOUS.STATVAR.excessWater_energy = 0;
        end
        
    end
end