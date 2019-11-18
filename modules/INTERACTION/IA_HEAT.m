%> IA_HEAT heat-only interaction between two classes
%> heat is transfered between two classes via a heatflux proportional to
%> the temperature difference at the border
%> no other state variables are touched
classdef IA_HEAT < matlab.mixin.Copyable 

    properties
        PREVIOUS
        NEXT
    end
    
    methods

        function  get_boundary_condition_m(ia_heat)
            % get the two classes involved
            stratigraphy1 = ia_heat.PREVIOUS;
            stratigraphy2 = ia_heat.NEXT;
            
            % calculate the heat flux based on the temperature difference
            % at the border
            flux = (stratigraphy1.STATVAR.T(end) - stratigraphy2.STATVAR.T(1)) .* stratigraphy1.STATVAR.thermCond(end) .* stratigraphy2.STATVAR.thermCond(1) ./...
                (stratigraphy1.STATVAR.thermCond(end).* stratigraphy2.STATVAR.layerThick(1)./2 + stratigraphy2.STATVAR.thermCond(1).* stratigraphy1.STATVAR.layerThick(end)./2 );
            
            % assign the flux to the lower/upper boundary of the
            % upper/lower class, respectively
            stratigraphy1.TEMP.heatFlux_lb = -flux;
            stratigraphy2.TEMP.heatFlux_ub = -flux;
            
            % add temperature at upper boundary of the lower class 
            % for conductivity calculations
            stratigraphy2.TEMP.T_ub = stratigraphy1.STATVAR.T(end);
        end
    end
end
