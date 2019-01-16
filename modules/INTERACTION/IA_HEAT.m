classdef IA_HEAT < matlab.mixin.Copyable 

     properties
        PREVIOUS
        NEXT
    end
    
    methods

        function  get_boundary_condition_m(ia_heat)
            stratigraphy1 = ia_heat.PREVIOUS;
            stratigraphy2 = ia_heat.NEXT;
            flux = (stratigraphy1.STATVAR.T(end) - stratigraphy2.STATVAR.T(1)) .* stratigraphy1.STATVAR.thermCond(end) .* stratigraphy2.STATVAR.thermCond(1) ./...
                (stratigraphy1.STATVAR.thermCond(end).* stratigraphy2.STATVAR.layerThick(1)./2 + stratigraphy2.STATVAR.thermCond(1).* stratigraphy1.STATVAR.layerThick(end)./2 );
            
            stratigraphy1.TEMP.F_lb = -flux;
            stratigraphy2.TEMP.F_ub = flux;
        end
    end
end
