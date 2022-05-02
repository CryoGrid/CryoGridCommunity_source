%========================================================================
% CryoGrid TIER1 INTERACTION (IA) class for functions related to heat conduction
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT < IA_BASE

    methods
        %heat conduction between two normal GROUND classes
        function ia_heat = get_boundary_condition_HEAT_m(ia_heat)
            stratigraphy1 = ia_heat.PREVIOUS;
            stratigraphy2 = ia_heat.NEXT;
            flux = (stratigraphy1.STATVAR.T(end) - stratigraphy2.STATVAR.T(1)) .* stratigraphy1.STATVAR.thermCond(end) .* stratigraphy2.STATVAR.thermCond(1) ./...
                (stratigraphy1.STATVAR.thermCond(end).* stratigraphy2.STATVAR.layerThick(1)./2 + stratigraphy2.STATVAR.thermCond(1).* stratigraphy1.STATVAR.layerThick(end)./2 );
            
            %flux = flux .* (stratigraphy1.STATVAR.area(end) + stratigraphy2.STATVAR.area(1)) ./ 2;
            flux = flux .* min(stratigraphy1.STATVAR.area(end), stratigraphy2.STATVAR.area(1));
            
            stratigraphy1.TEMP.F_lb = -flux;
            stratigraphy2.TEMP.F_ub = flux;
            stratigraphy1.TEMP.d_energy(end) = stratigraphy1.TEMP.d_energy(end) - flux;
            stratigraphy2.TEMP.d_energy(1) = stratigraphy2.TEMP.d_energy(1) + flux;
        end
        
        %heat conduction between LAKE and GROUND class
        function  ia_heat = get_boundary_condition_HEAT_LAKE_m(ia_heat) %coupling between LAKE (PREVIOUS) and GROUND (NEXT) when unfrozen and the lake is a big cell
            stratigraphy1 = ia_heat.PREVIOUS;
            stratigraphy2 = ia_heat.NEXT;
            flux = (stratigraphy1.STATVAR.T(end) - stratigraphy2.STATVAR.T(1))  .* stratigraphy2.STATVAR.thermCond(1) ./ (stratigraphy2.STATVAR.layerThick(1) ./ 2);
            
            %flux = flux .* (stratigraphy1.STATVAR.area(end) + stratigraphy2.STATVAR.area(1)) ./ 2;
            flux = flux .* min(stratigraphy1.STATVAR.area(end), stratigraphy2.STATVAR.area(1));
            
            stratigraphy1.TEMP.F_lb = -flux;
            stratigraphy2.TEMP.F_ub = flux;
            stratigraphy1.TEMP.d_energy(end) = stratigraphy1.TEMP.d_energy(end) - flux;
            stratigraphy2.TEMP.d_energy(1) = stratigraphy2.TEMP.d_energy(1) + flux;
        end
        
        %heat conduction between SNOW in CHILD phase and GROUND class
        function get_boundary_condition_HEAT_IA_CHILD(ia_heat)
            stratigraphy1 = ia_heat.PREVIOUS; %snow
            stratigraphy2 = ia_heat.NEXT; %ground
            flux = (stratigraphy1.STATVAR.T - stratigraphy2.STATVAR.T(1)) .* stratigraphy1.STATVAR.thermCond .* stratigraphy2.STATVAR.thermCond(1) ./...
                (stratigraphy1.STATVAR.thermCond.* stratigraphy2.STATVAR.layerThick(1)./2 + stratigraphy2.STATVAR.thermCond(1).* stratigraphy1.STATVAR.layerThick./2 );
            
            flux = flux .* stratigraphy1.STATVAR.area(end); %SNOW area
            
            stratigraphy1.TEMP.F_lb = -flux;
            stratigraphy2.TEMP.F_ub = stratigraphy2.TEMP.F_ub + flux;  %ground already has F_ub from surface nergy balance
            stratigraphy1.TEMP.d_energy(end) = stratigraphy1.TEMP.d_energy(end) - flux;
            stratigraphy2.TEMP.d_energy(1) = stratigraphy2.TEMP.d_energy(1) + flux;
        end
        
        %heat conduction with GROUND_TTOP_SIMPLE2
        function ia_heat = get_boundary_condition_HEAT_TTOP_m(ia_heat) 
            stratigraphy1 = ia_heat.PREVIOUS;
            stratigraphy2 = ia_heat.NEXT;
            flux = (stratigraphy1.STATVAR.MAGT - stratigraphy2.STATVAR.T(1)) .* stratigraphy2.STATVAR.thermCond(1) ./stratigraphy2.STATVAR.layerThick(1) ;
            
            flux = flux .* stratigraphy2.STATVAR.area(1);
            
%             stratigraphy1.TEMP.F_lb = -flux;
%             stratigraphy2.TEMP.F_ub = flux;
%             stratigraphy1.TEMP.d_energy(end) = stratigraphy1.TEMP.d_energy(end) - flux;
            stratigraphy2.TEMP.d_energy(1) = stratigraphy2.TEMP.d_energy(1) + flux;
        end
        
    end
end
