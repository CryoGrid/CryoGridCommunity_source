classdef IA_HEAT11_WATER11_SNOW_XICE < IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_XICE_m(ia_heat_water);
        end
        
        %SNOW
        function get_IA_CHILD_boundary_condition_u(ia_heat_water)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat_water);
            get_boundary_condition_BUCKET_SNOW_XICE_m(ia_heat_water);
        end
        
        function remove_excessWater_CHILD(ia_heat_water) %move excessWater from SNOW to Xwater from PARENT, 0 energy transfer since meltwater must be zero
            ia_heat_water.NEXT.STATVAR.XwaterIce(1) = ia_heat_water.NEXT.STATVAR.XwaterIce(1) + ia_heat_water.PREVIOUS.STATVAR.excessWater;
            ia_heat_water.NEXT.STATVAR.Xwater(1) = ia_heat_water.NEXT.STATVAR.Xwater(1) + ia_heat_water.PREVIOUS.STATVAR.excessWater;
            ia_heat_water.NEXT.STATVAR.layerThick(1) = ia_heat_water.NEXT.STATVAR.layerThick(1) + ia_heat_water.PREVIOUS.STATVAR.excessWater ./ ia_heat_water.NEXT.STATVAR.area(1);
            ia_heat_water.PREVIOUS.STATVAR.excessWater = 0;
        end
        
    end
end