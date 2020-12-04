classdef IA_HEAT_GROUND_SNOW_VEGETATION < matlab.mixin.Copyable
    
    %> Simone: Interaction class between Vegetation and Ground
    %> Soil heat flux (W/m^2) up and down
    
    properties
        PREVIOUS
        NEXT
        IA_PARENT_GROUND
        IA_CHILD_SNOW
        FRACTIONAL_SNOW_COVER
        STATUS
    end
    
    methods
        
        % 0 no child
        % 1 child initialized / growing
        % 2 Child existing fractionally or melting
        % -1 child existing
        
        
        function ia_heat_ground_snow_vegetation = get_boundary_condition_m(ia_heat_ground_snow_vegetation)
            ground = ia_heat_ground_snow_vegetation.NEXT;
            snow = ia_heat_ground_snow_vegetation.PREVIOUS;

                flux = (snow.STATVAR.T(end) - ground.STATVAR.T(1)) .* snow.STATVAR.thermCond(end) .* ground.STATVAR.thermCond(1) ./...
                    (snow.STATVAR.thermCond(end).* ground.STATVAR.layerThick(1)./2 + ground.STATVAR.thermCond(1).* snow.STATVAR.layerThick(end)./2 );
                
                snow.TEMP.F_lb = -flux;
                ground.TEMP.F_ub = flux;

        end
    end
end