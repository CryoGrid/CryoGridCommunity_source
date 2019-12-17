classdef IA_HEAT_WATER_SNOW_GROUND < matlab.mixin.Copyable %zero water flux between classes, lower class does not have water balance 

    
     properties
        PREVIOUS
        NEXT
    end
    
    methods

        function  get_boundary_condition_m(ia_heat_water)
            stratigraphy1 = ia_heat_water.PREVIOUS;
            stratigraphy2 = ia_heat_water.NEXT;
            flux = (stratigraphy1.STATVAR.T(end) - stratigraphy2.STATVAR.T(1)) .* stratigraphy1.STATVAR.thermCond(end) .* stratigraphy2.STATVAR.thermCond(1) ./...
                (stratigraphy1.STATVAR.thermCond(end).* stratigraphy2.STATVAR.layerThick(1)./2 + stratigraphy2.STATVAR.thermCond(1).* stratigraphy1.STATVAR.layerThick(end)./2 );
            
            stratigraphy1.TEMP.F_lb = -flux;
            stratigraphy2.TEMP.F_ub = flux;
        end
        
        function get_boundary_condition_water_m(ia_heat_water)
            saturation = ia_heat_water.PREVIOUS.STATVAR.water(end,1) ./ max(1e-12, ia_heat_water.PREVIOUS.STATVAR.layerThick(end,1) - ia_heat_water.PREVIOUS.STATVAR.ice(end,1));
            waterMobile = double(saturation > ia_heat_water.PREVIOUS.PARA.field_capacity);
            ia_heat_water.PREVIOUS.TEMP.d_water_out(end,1) = waterMobile .* ia_heat_water.PREVIOUS.PARA.hydraulicConductivity .* ia_heat_water.PREVIOUS.STATVAR.water(end,1) ./ ia_heat_water.PREVIOUS.STATVAR.layerThick(end,1);
            ia_heat_water.NEXT.TEMP.F_ub_water = ia_heat_water.PREVIOUS.TEMP.d_water_out(end,1);
        end
        
        function finalize_boundary_condition_water_m(ia_heat_water, timestep)
             %limit outflow to field capacity
            ia_heat_water.PREVIOUS.TEMP.d_water_out(end,1)  = min(ia_heat_water.PREVIOUS.TEMP.d_water_out(end,1), max(0, ia_heat_water.PREVIOUS.STATVAR.water(end,1) - ia_heat_water.PREVIOUS.PARA.field_capacity .* (ia_heat_water.PREVIOUS.STATVAR.layerThick(end,1) - ia_heat_water.PREVIOUS.STATVAR.ice(end,1))));
            ia_heat_water.NEXT.TEMP.d_water_in(1,1) = ia_heat_water.PREVIOUS.TEMP.d_water_out(end,1);
            %limit inflow so that unity is not exceeded
            ia_heat_water.NEXT.TEMP.d_water_in(1,1) = min(ia_heat_water.NEXT.TEMP.d_water_in(1,1), ground.STATVAR.layerThick(1,1) - ground.STATVAR.mineral(1,1) - ground.STATVAR.organic(1,1) - ground.STATVAR.waterIce(1,1));
            ia_heat_water.PREVIOUS.TEMP.d_water_out(end,1) = ia_heat_water.NEXT.TEMP.d_water_in(1,1);
            
            ia_heat_water.NEXT.TEMP.d_water_in(1,1) = ia_heat_water.NEXT.TEMP.d_water_in(1,1) ./timestep;

        end
    end
end
