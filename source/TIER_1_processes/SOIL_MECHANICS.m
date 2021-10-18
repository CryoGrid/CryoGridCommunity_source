%========================================================================
% CryoGrid TIER1 library class, functions related to soil mechanics
% J. Schmidt, May 2021
%========================================================================

classdef SOIL_MECHANICS < BASE

    methods
        
        %Calculate settlement out of water flow, as well as overburden pressure and bearing capacity for the next timestep
        function ground = soil_mechanics(ground)
            
            %Porosity, void ratio and saturation is updated as this changed for saturated cells in advanced prognostics
            ground.STATVAR.porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ground.STATVAR.layerThick./ground.STATVAR.area;
            ground.STATVAR.void_ratio = ground.STATVAR.porosity ./ (1 - ground.STATVAR.porosity);
            ground.STATVAR.saturation = (ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.porosity;
            
            %Find all unsaturated cells --> will be updated in this function
            unsaturated = ground.STATVAR.saturation <= 1-1e-6;

            %Calculate the overburden pressure for each cell
            overburden_pressure_per_cell_total = (ground.STATVAR.water .* ground.CONST.density_water + ground.STATVAR.ice .* ground.CONST.density_ice + ground.STATVAR.mineral .* ground.CONST.density_mineral + ...
                ground.STATVAR.organic .* ground.CONST.density_organic) ./ ground.STATVAR.area .*ground.CONST.g; %[Pa]
            obp_old = ground.STATVAR.overburden_pressure;
            ground.STATVAR.overburden_pressure = cumsum(overburden_pressure_per_cell_total)-overburden_pressure_per_cell_total./2; %[Pa]
            
            %Add external pressure
            ground.STATVAR.overburden_pressure = ground.STATVAR.overburden_pressure + ground.STATVAR.external_pressure;
            
            %Make mean of current overburden pressure and 9x old overburden pressure --> needed for stable calculations
            ground.STATVAR.overburden_pressure = (ground.STATVAR.overburden_pressure + 9.*obp_old)./10; 
            
            %Update porosity and layerThick for unsaturated cells
            ground.STATVAR.porosity(unsaturated) = (ground.STATVAR.initial_voidRatio(unsaturated) - ground.STATVAR.compression_index(unsaturated) .* log10(ground.STATVAR.overburden_pressure(unsaturated)./ground.STATVAR.reference_pressure(unsaturated))) ./ ...
                (1 + ground.STATVAR.initial_voidRatio(unsaturated) - ground.STATVAR.compression_index(unsaturated) .* log10(ground.STATVAR.overburden_pressure(unsaturated)./ground.STATVAR.reference_pressure(unsaturated)));
            %do not allow higher porosities than initial void ration
            ground.STATVAR.porosity(unsaturated) = min(ground.STATVAR.porosity(unsaturated), ground.STATVAR.initial_voidRatio(unsaturated) ./(1 + ground.STATVAR.initial_voidRatio(unsaturated)));
            ground.STATVAR.layerThick(unsaturated) = max(((ground.STATVAR.mineral(unsaturated) + ground.STATVAR.organic(unsaturated)) ./ (1 - ground.STATVAR.porosity(unsaturated))) ./ ground.STATVAR.area(unsaturated), ...
                (ground.STATVAR.mineral(unsaturated) + ground.STATVAR.organic(unsaturated)+ ground.STATVAR.waterIce(unsaturated))./ ground.STATVAR.area(unsaturated));
            
            %Recalculate porosity, void ratio and saturation for all cells            
            ground.STATVAR.porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ground.STATVAR.layerThick./ground.STATVAR.area;
            ground.STATVAR.saturation = (ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.porosity;
            ground.STATVAR.void_ratio = ground.STATVAR.porosity ./ (1 - ground.STATVAR.porosity);

            %Calculate bearing capacity (= effective stress on soil)
            bc_old = ground.STATVAR.bearing_capacity;
            ground.STATVAR.bearing_capacity = (10.^((ground.STATVAR.initial_voidRatio - ground.STATVAR.porosity - ground.STATVAR.initial_voidRatio .* ground.STATVAR.porosity) ./ (ground.STATVAR.compression_index - ground.STATVAR.porosity .* ground.STATVAR.compression_index))) .* ground.STATVAR.reference_pressure;
            ground.STATVAR.bearing_capacity = max(ground.STATVAR.reference_pressure, ground.STATVAR.bearing_capacity); %generate "Xice" (void ratio becomes smaller than initial void ratio)
            ground.STATVAR.bearing_capacity = (ground.STATVAR.bearing_capacity + 9.*bc_old)./10;  
            
            %Split waterIce up in XwaterIce and waterIce
            porosity_equilibrium = ground.STATVAR.initial_voidRatio ./ (1 + ground.STATVAR.initial_voidRatio);
            waterIce_equilibrium = (ground.STATVAR.mineral + ground.STATVAR.organic) ./ (1 - porosity_equilibrium) - ground.STATVAR.mineral - ground.STATVAR.organic;
            
            ground.STATVAR.XwaterIce = ground.STATVAR.waterIce.*0;
            excessIce = ground.STATVAR.saturation >=1-1e-6 & waterIce_equilibrium < ground.STATVAR.waterIce;
            ground.STATVAR.XwaterIce(excessIce) = ground.STATVAR.waterIce(excessIce) - waterIce_equilibrium(excessIce);
            ground.STATVAR.waterIce(excessIce) = waterIce_equilibrium(excessIce);
        end
    end
end

