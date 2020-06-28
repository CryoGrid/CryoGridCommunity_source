classdef HEAT_FLUXES_LATERAL < BASE


    methods

        
        %---------3D-coupled-fluxes-------

        function ground = lateral3D_pull_heat_simple(ground, lateral)
    
            
            if isempty(lateral.PARENT.STATVAR.depths_heat)
                depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            else
                depths = ground.STATVAR.upperPos - cumsum(ground.STATVAR.layerThick);
            end
            lateral.PARENT.STATVAR.depths_heat = [lateral.PARENT.STATVAR.depths_heat; depths];
            lateral.PARENT.STATVAR.thermCond = [lateral.PARENT.STATVAR.thermCond; ground.STATVAR.thermCond];
            lateral.PARENT.STATVAR.T_heat = [lateral.PARENT.STATVAR.T_heat; ground.STATVAR.T];
        end

        
        function ground = lateral3D_push_heat_simple(ground, lateral)

            ground.STATVAR.energy = ground.STATVAR.energy + lateral.PARENT.STATVAR.heat_flux(1:size(ground.STATVAR.energy,1),1);
            lateral.PARENT.STATVAR.heat_flux(1:size(ground.STATVAR.energy,1),:) = [];
            
        end

    end
end

