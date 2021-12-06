%========================================================================
% CryoGrid TIER1 library class for functions for lateral fluxes of heat
% contains push and pull functions used in lateral heat flux classes, e.g. LAT3D_HEAT 
% S. Westermann, October 2020
%========================================================================


classdef HEAT_FLUXES_LATERAL < BASE


    methods

        function ground = lateral_push_heat_simple(ground, lateral)
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            depths = (depths(1:end-1,1) + depths(2:end,1))./2;
            
            cross_section = lateral.PARA.heatReservoir_contact_length .* ground.STATVAR.layerThick .* double(depths >= lateral.PARA.lowerElevation & depths <= lateral.PARA.upperElevation);
            
            flux = ground.STATVAR.thermCond .* (lateral.PARA.reservoir_T - ground.STATVAR.T)./ lateral.PARA.distance_heatReservoir .* cross_section;
            flux = flux .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec;
            
            ground.STATVAR.energy = ground.STATVAR.energy + flux;
        end
        
        %read information from GROUND class and send it to the LATERAL class   
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

         % add lateral heat flux to STATVAR energy 
        function ground = lateral3D_push_heat_simple(ground, lateral)
            ground.STATVAR.energy = ground.STATVAR.energy + lateral.PARENT.STATVAR.heat_flux(1:size(ground.STATVAR.energy,1),1);
            lateral.PARENT.STATVAR.heat_flux(1:size(ground.STATVAR.energy,1),:) = [];
        end

    end
end

