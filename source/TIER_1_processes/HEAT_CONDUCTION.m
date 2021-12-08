%========================================================================
% CryoGrid TIER1 library class for functions related heat conduction
% NOTE: this class also contains code related to the free water freeze curve, 
% as well as functions computing thermal conductivity  
% S. Westermann, October 2020
%========================================================================

classdef HEAT_CONDUCTION < BASE
    
    methods
        
        %-----derivatives----------
        %conductive heat flux between grid cells
        function ground = get_derivative_energy(ground)
            fluxes = (ground.STATVAR.T(1:end-1) - ground.STATVAR.T(2:end)) .* ground.STATVAR.thermCond(1:end-1) .* ground.STATVAR.thermCond(2:end) ./...
                (ground.STATVAR.thermCond(1:end-1).* ground.STATVAR.layerThick(2:end)./2 +  ground.STATVAR.thermCond(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            
            fluxes = fluxes .*(ground.STATVAR.area(1:end-1) + ground.STATVAR.area(2:end))./2; %[J/sec]
            
            d_energy=ground.STATVAR.energy.*0;
            
            d_energy(1) =  - fluxes(1);
            d_energy(2:end-1) = fluxes(1:end-1) - fluxes(2:end);
            d_energy(end) =  fluxes(end);
            
            ground.TEMP.d_energy = ground.TEMP.d_energy + d_energy;
        end
        
        %-----------timesteps----------
        %limit maximum energy change between timesteps
        function timestep = get_timestep_heat_coduction(ground)  
            timestep = ground.PARA.dE_max ./ (max(abs(ground.TEMP.d_energy) ./ ground.STATVAR.layerThick./ ground.STATVAR.area));
        end
        
        %----diagnostic functions---------
        %free water freeze curve
        function ground = get_T_water_freeW(ground)
            
            Lf = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            
            E_frozen = -Lf.*ground.STATVAR.waterIce;
            
            ground.STATVAR.T = double(ground.STATVAR.energy < E_frozen) .* (ground.STATVAR.energy - E_frozen) ./ (c_i.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic) + ...
                double(ground.STATVAR.energy >0) .* ground.STATVAR.energy ./ (c_w.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic);
            ground.STATVAR.ice = double(ground.STATVAR.energy <= E_frozen) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > E_frozen & ground.STATVAR.energy < 0) .* ground.STATVAR.energy ./ (-Lf);
            ground.STATVAR.water = double(ground.STATVAR.energy >= 0) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > - Lf.*ground.STATVAR.waterIce & ground.STATVAR.energy < 0) .* (ground.STATVAR.energy + Lf.*ground.STATVAR.waterIce) ./ Lf;
        end
        
        %calculate energy from temperature and water contents, free water
        %freeze curve
        function ground = get_E_freeW(ground) %required for initialization

            T = ground.STATVAR.T;
            mineral= ground.STATVAR.mineral;
            organic = ground.STATVAR.organic;
            waterIce = ground.STATVAR.waterIce;
            layerThick = ground.STATVAR.layerThick;
            area = ground.STATVAR.area;

            energy = T.*(mineral .* ground.CONST.c_m + organic .* ground.CONST.c_o + double(T>=0).*(waterIce .* ground.CONST.c_w) + ...
                double(T<0).*(waterIce .* ground.CONST.c_i )) - double(T<0) .* (waterIce) .* ground.CONST.L_f;
            
            ground.STATVAR.waterIce = waterIce .* layerThick .* area; % [m3]
            ground.STATVAR.mineral = mineral .* layerThick .* area; % [m3]
            ground.STATVAR.organic = organic .* layerThick .* area; % [m3]
            ground.STATVAR.energy = energy .* layerThick .* area;  % [J]
            
            ground.STATVAR.water = double(T>=0) .* waterIce .* layerThick .* area;  % [m3]
            ground.STATVAR.ice = double(T<0) .* waterIce .* layerThick .* area; %[m3]
            ground.STATVAR.air = (1-mineral-organic-waterIce) .* layerThick .* area;  % [m3]
            ground = conductivity(ground);
            
        end
        
        %---thermal conductivities--------------
        function ground = conductivity_mixing_squares(ground)
            
            water = ground.STATVAR.water./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ice = ground.STATVAR.ice ./ ground.STATVAR.layerThick./ ground.STATVAR.area;
            mineral = ground.STATVAR.mineral./ground.STATVAR.layerThick./ ground.STATVAR.area;
            organic = ground.STATVAR.organic./ground.STATVAR.layerThick./ ground.STATVAR.area;
            air = 1 - mineral - organic - water - ice;
            
            ground.STATVAR.thermCond = (water.* ground.CONST.k_w.^0.5 + ice.* ground.CONST.k_i.^0.5 ...
                + mineral.* ground.CONST.k_m.^0.5 + organic.* ground.CONST.k_o.^0.5 + air.* ground.CONST.k_a.^0.5).^2;
        end
        
        function ground = conductivity_mixing_squares_Xice(ground)
            water = (ground.STATVAR.water + ground.STATVAR.Xwater)./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ice = (ground.STATVAR.ice + ground.STATVAR.Xice)./ ground.STATVAR.layerThick./ ground.STATVAR.area;
            mineral = ground.STATVAR.mineral./ground.STATVAR.layerThick./ ground.STATVAR.area;
            organic = ground.STATVAR.organic./ground.STATVAR.layerThick./ ground.STATVAR.area;
            air = 1 - mineral - organic - water - ice;
            
            ground.STATVAR.thermCond = (water.* ground.CONST.k_w.^0.5 + ice.* ground.CONST.k_i.^0.5 ...
                + mineral.* ground.CONST.k_m.^0.5 + organic.* ground.CONST.k_o.^0.5 + air.* ground.CONST.k_a.^0.5).^2;
            
        end
        
        function ground = thermalConductivity_CLM4_5(ground)
            
            k_dry_organic = 0.05; %slightly nonsense...
                        
            waterIce = ground.STATVAR.waterIce./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            water = ground.STATVAR.water./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ice = ground.STATVAR.ice ./ ground.STATVAR.layerThick./ ground.STATVAR.area;
            mineral = ground.STATVAR.mineral./ground.STATVAR.layerThick./ ground.STATVAR.area;
            organic = ground.STATVAR.organic./ground.STATVAR.layerThick./ ground.STATVAR.area;
            porosity = 1 - mineral - organic;
            organic_fraction = organic ./ (mineral + organic);
            saturation = waterIce./porosity;
            
            k_solids = organic_fraction .* ground.CONST.k_o + (1- organic_fraction) .* ground.CONST.k_m;
            k_sat = k_solids.^(1-porosity) .* ground.CONST.k_w .^(water./waterIce.* porosity) .* ground.CONST.k_i .^(ice./waterIce.* porosity);
            
            bulk_density = 2700 .* (1-porosity);
            k_dry_mineral = (0.135 .* bulk_density + 64.7) ./ (2700 - 0.947 .* bulk_density);
            k_dry = organic_fraction .* k_dry_organic + (1- organic_fraction) .* k_dry_mineral;

            Kersten_number = double(ground.STATVAR.T>=0) .* max(0, log(saturation) ./ log(10) + 1) +  double(ground.STATVAR.T<0) .* saturation;
            
            ground.STATVAR.thermCond = Kersten_number .* k_sat + (1- Kersten_number) .* k_dry;
            
        end
        
        function snow = conductivity_snow_Yen(snow)
            
            ki = 2.2196 - 0.0062489 .* snow.STATVAR.T + 0.00010154.*snow.STATVAR.T.^2;
            snow.STATVAR.thermCond = ki.*(snow.STATVAR.waterIce./snow.STATVAR.layerThick./ snow.STATVAR.area).^1.88;
        end
        
    end
end

