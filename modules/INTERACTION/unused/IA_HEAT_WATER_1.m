classdef IA_HEAT_WATER_1 < matlab.mixin.Copyable %zero water flux between classes, lower class does not have water balance 

    
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
            ia_heat_water.PREVIOUS.TEMP.d_water_out(end,1) = 0; %zero flux out for upper class, lower class does not have water balance
        end
        
        function finalize_boundary_condition_water_m(ia_heat_water, timestep)
            %do nothing, water is not routed further down
        end
    end
end
