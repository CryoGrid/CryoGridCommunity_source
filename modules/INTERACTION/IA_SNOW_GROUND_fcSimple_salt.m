classdef IA_SNOW_GROUND_fcSimple_salt
    
    properties
        IA_PARENT_GROUND
        IA_CHILD_SNOW
        STATUS
    end
    
    methods
        function ia_snow_ground = get_boundary_condition_u(ia_snow_ground, forcing)
            if ia_snow_ground.STATUS == 1 || forcing.TEMP.snowfall > 0 %non-zero SWE, but snow is still a child, or snowfall occurring
                ia_snow_ground.IA_CHILD_SNOW.STATVAR.Lstar = ia_snow_ground.IA_PARENT_GROUND.STATVAR.Lstar;  %assign L_star from parent
                ia_snow_ground.IA_CHILD_SNOW = get_boundary_condition_u(ia_snow_ground.IA_CHILD_SNOW, forcing); %call the native function for the snow class
                ia_snow_ground.STATUS = 1;
            end
        end
        
        function ia_snow_ground = get_derivative_energy(ia_snow_ground)
            if ia_snow_ground.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_snow_ground.IA_CHILD_SNOW;
                ground = ia_snow_ground.IA_PARENT_GROUND;
                fraction_snow = snow.STATVAR.waterIce./(snow.PARA.swe_per_cell./2);
                fraction_snow = min(1, fraction_snow);
                
                %change parent variable
                ground.TEMP.F_ub = (1-fraction_snow) .* ground.TEMP.F_ub + fraction_snow.* snow.TEMP.F_ub; %mix the surface heat fluxes from snow and ground
                
                ground.STATVAR.Lout = (1-fraction_snow) .* ia_snow_ground.IA_PARENT_GROUND.STATVAR.Lout + fraction_snow.* snow.STATVAR.Lout; %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Sout = (1-fraction_snow) .* ia_snow_ground.IA_PARENT_GROUND.STATVAR.Sout + fraction_snow.* snow.STATVAR.Sout; %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Qh = (1-fraction_snow) .* ground.STATVAR.Qh + fraction_snow.* snow.STATVAR.Qh; %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Qe = (1-fraction_snow) .* ground.STATVAR.Qe + fraction_snow.* snow.STATVAR.Qe; %mix the surface heat fluxes from snow and ground
                %mixing of Qh and Qe will impact atmospheric stability
                
                layerThick = ground.STATVAR.layerThick;
                layerThick(1,1) = layerThick(1,1) + snow.STATVAR.layerThick; %add thickness of snow cell for calculation of effective conductivities
                fluxes = (ground.STATVAR.T(1:end-1) - ground.STATVAR.T(2:end)) .* ground.STATVAR.thermCond(1:end-1) .* ground.STATVAR.thermCond(2:end) ./...
                    (ground.STATVAR.thermCond(1:end-1).* layerThick(2:end)./2 +  ground.STATVAR.thermCond(2:end).* layerThick(1:end-1)./2 );
                
                d_energy=ground.STATVAR.energy.*0;
                d_energy(1) = ground.TEMP.F_ub - fluxes(1);
                d_energy(2:end-1) = fluxes(1:end-1) - fluxes(2:end);
                d_energy(end) = ground.TEMP.F_lb + fluxes(end);
                
                ground.TEMP.d_energy = d_energy;

            else
                ia_snow_ground.IA_PARENT_GROUND = get_derivative_energy(ia_snow_ground.IA_PARENT_GROUND); %call normal function
            end
        end
        
        function ia_snow_ground = advance_prognostic(ia_snow_ground, timestep)
            if ia_snow_ground.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_snow_ground.IA_CHILD_SNOW;
                
                snow.STATVAR.energy = snow.STATVAR.energy + timestep .* snow.TEMP.snow_energy + timestep .* snow.TEMP.rain_energy; %only the energy from new snow
                snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* (snow.TEMP.snowfall + snow.TEMP.rainfall);
                
                snow.STATVAR.target_density = (snow.STATVAR.ice + timestep .* snow.TEMP.snowfall)./ (snow.STATVAR.layerThick + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./snow.CONST.rho_w));
                
                %snow.STATVAR.target_density = snow.STATVAR.ice_fraction; %remove later
                
                %snow.STATVAR.layerThick = snow.STATVAR.layerThick + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./snow.CONST.rho_w);
            end
        end
        
        function ia_snow_ground = compute_diagnostic(ia_snow_ground)
            if ia_snow_ground.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_snow_ground.IA_CHILD_SNOW;
                ground = ia_snow_ground.IA_PARENT_GROUND;
                
               
                [ground, snow] = get_T_water_salt_FreezeDepress_Xice(ia_snow_ground, ground, snow); %function below
                
                %snow.STATVAR.layerThick = snow.STATVAR.ice ./snow.STATVAR.ice_fraction; %conserving old density
                snow.STATVAR.layerThick = snow.STATVAR.ice ./snow.STATVAR.target_density; %conserving old density
                %runoff = max(0, snow.STATVAR.water - snow.PARA.field_capacity .* (1-snow.STATVAR.ice_fraction) .* snow.STATVAR.layerThick);  %could be used to reroute in the ground
                runoff = max(0, snow.STATVAR.water - snow.PARA.field_capacity .* (1-snow.STATVAR.target_density) .* snow.STATVAR.layerThick);  %could be used to reroute in the ground
                %snow.STATVAR.water = min(snow.STATVAR.water, snow.PARA.field_capacity .* (1-snow.STATVAR.ice_fraction) .* snow.STATVAR.layerThick);
                snow.STATVAR.water = min(snow.STATVAR.water, snow.PARA.field_capacity .* (1-snow.STATVAR.target_density) .* snow.STATVAR.layerThick);
                snow.STATVAR.waterIce = snow.STATVAR.water + snow.STATVAR.ice;
                
                
                if snow.STATVAR.waterIce == 0;
                    ia_snow_ground.STATUS = 0;
                    snow = initialize_zero_snow(snow, ground); %set all variables to zero
                end
            end
        end
        
        function ia_snow_ground = mix_conductivity(ia_snow_ground)
            if ia_snow_ground.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_snow_ground.IA_CHILD_SNOW;
                ground = ia_snow_ground.IA_PARENT_GROUND;
                
                snow = conductivity(snow);
                mixed_cond = (snow.STATVAR.layerThick + ground.STATVAR.layerThick(1,1)) .* snow.STATVAR.thermCond .* ground.STATVAR.thermCond(1,1) ./ ...
                    (snow.STATVAR.layerThick .* ground.STATVAR.thermCond(1,1) + ground.STATVAR.layerThick(1,1) .* snow.STATVAR.thermCond);
                ground.STATVAR.thermCond(1,1) = mixed_cond;
            end
        end
        
        function ia_snow_ground = check_trigger(ia_snow_ground)
            if ia_snow_ground.STATUS == 1 && ia_snow_ground.IA_CHILD_SNOW.STATVAR.waterIce > ia_snow_ground.IA_CHILD_SNOW.PARA.swe_per_cell./2
                
                ia_snow_ground.IA_CHILD_SNOW.PREVIOUS = ia_snow_ground.IA_PARENT_GROUND.PREVIOUS;
                ia_snow_ground.IA_CHILD_SNOW.PREVIOUS.NEXT = ia_snow_ground.IA_CHILD_SNOW;
                
                ia_snow_ground.IA_CHILD_SNOW.NEXT = ia_snow_ground.IA_PARENT_GROUND;
                ia_snow_ground.IA_PARENT_GROUND.PREVIOUS = ia_snow_ground.IA_CHILD_SNOW;
                
                
                temp_store = ia_snow_ground.IA_CHILD_SNOW; %snow class
                
                ia_class = get_IA_class(class(temp_store), class(temp_store.NEXT)); %put snow on top
                temp_store.IA_NEXT = ia_class;
                temp_store.IA_NEXT.PREVIOUS = temp_store;
                temp_store.IA_NEXT.NEXT = temp_store.NEXT;
                temp_store.NEXT.IA_PREVIOUS = temp_store.IA_NEXT;
                
                ia_snow_ground = [];  % destroy the class
            end
        end
        
        function [ground, snow] = get_T_water_salt_FreezeDepress_Xice(ia_snow_ground, ground, snow)

            c_w = ground.CONST.c_w; c_i = ground.CONST.c_i; c_o = ground.CONST.c_o; c_m = ground.CONST.c_m;
            R = ground.CONST.R; Tmfw = ground.CONST.Tmfw; L_f = ground.CONST.L_f;
            
            deltaT = ground.STATVAR.deltaT(1,1);

            energy = ground.STATVAR.energy(1) + snow.STATVAR.energy; %total energy of combined cell
            
            freeWaterIce = snow.STATVAR.waterIce; 
            
            waterIce = ground.STATVAR.waterIce(1,1);
            mineral = ground.STATVAR.mineral(1,1);
            organic = ground.STATVAR.organic(1,1);

            N = ground.STATVAR.saltConc(1,1) ./ ground.STATVAR.waterIce(1);
            
            A = 1 + c_i.*deltaT.*freeWaterIce./(waterIce.*L_f);
            A1 = c_w.*waterIce + c_m.*mineral + c_o.*organic;
            A2 = c_i.*waterIce + c_m.*mineral + c_o.*organic;
            A3 = (c_w-c_i).*waterIce;
            A4 = waterIce .* L_f;
            B = -L_f./(R.* Tmfw.^2);
            
            %quadratic equation in Tm, the onset of freezing T, a*Tm^2 + b*Tm + c = 0
            %zero-th order terms
            c = - N .* A2 .* deltaT - N .* A4;
            
            %first order terms
            b = - N .* A3 + B .*(energy ./ A + L_f .*freeWaterIce ./ A + A2 .* deltaT + A4);
            
            %second-order terms
            a = B.* (-c_i.*freeWaterIce./A - (c_i .* deltaT .* freeWaterIce .*A1) ./ (A .*A4) - A2);
            
            %Tm_2 = (-b + sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
            Tm_1 = (-b - sqrt(b.^2 - 4.*a.*c)) ./ (2.*a); %this is the right branch!!
            
            thresh1 = 0;
            thresh2 = - L_f.*freeWaterIce;
            thresh3 = - L_f .* freeWaterIce + (c_w.*waterIce + c_m.*mineral + c_o.*organic + c_i .*freeWaterIce) .* (-R.* Tmfw.^2 ./L_f).* N;
            
            T = double(energy >= thresh1) .* energy ./(c_w .*waterIce + c_w.*freeWaterIce + c_m.*mineral + c_o.*organic) + ...
                double(energy < thresh2 & energy > thresh3) .* (energy + L_f.*freeWaterIce) ./  (c_w .*waterIce + c_i.*freeWaterIce + c_m.*mineral + c_o.*organic) + ...
                double(energy <= thresh3) .* (Tm_1 - deltaT + deltaT.* (energy./A - c_i.*freeWaterIce./A .*Tm_1 - c_i.*deltaT.*freeWaterIce./(A.*A4) .* A1.*Tm_1 + L_f.*freeWaterIce ./ A ...
                - (A2 .*(Tm_1-deltaT) - A4)) ./ (A3.*Tm_1 + A2.*deltaT + A4));
            
            
            ice = double(energy <= thresh3) .* waterIce .*(1-(energy./A - c_i.*freeWaterIce./A .*Tm_1 - c_i.*deltaT.*freeWaterIce./(A.*A4) .* A1.*Tm_1 + L_f.*freeWaterIce ./ A ...
                - (A2 .*(Tm_1-deltaT) - A4)) ./ (A3.*Tm_1 + A2.*deltaT + A4));
            water = waterIce - ice;
            salt_c_brine = ground.STATVAR.saltConc(1,1) ./ water;
            salt_c_brine(isnan(salt_c_brine))=0;

            freeIce = double(energy >= thresh2 & energy < thresh1) .* energy ./thresh2 .* freeWaterIce + double(energy < thresh2) .* freeWaterIce;
            freeWater = freeWaterIce - freeIce;
            
            
            ground.STATVAR.T(1,1) = T;
            ground.STATVAR.salt_c_brine(1,1) = salt_c_brine;
            
            ground.STATVAR.ice(1,1) = ice;
            ground.STATVAR.water(1,1) = water;

            snow.STATVAR.T = T;
            snow.STATVAR.ice = freeIce;
            snow.STATVAR.water = freeWater;

            
            snow.STATVAR.energy = snow.STATVAR.ice .* (c_i.*snow.STATVAR.T - L_f); %redistribute energy between snow and ground cell
            ground.STATVAR.energy(1) = energy - snow.STATVAR.energy;
        end
        
    end
    
end

