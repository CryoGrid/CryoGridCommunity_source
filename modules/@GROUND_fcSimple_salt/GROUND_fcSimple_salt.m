%base class for a GROUND object with a simple frteeze curve (constant freezing range between 0 and -deltaT degree C), including diffusion of salt - upper
%boundary condition is not specified; %superclass that cannot be run alone

classdef GROUND_fcSimple_salt < GROUND_base_class

        
    methods
        %mandatory functions for each class
        
        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_variables@GROUND_base_class(ground);
            ground = provide_PARA(ground);
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        
        function variable = initialize_from_file(ground, variable, section)
            variable = initialize_from_file@GROUND_base_class(ground, variable, section);
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_base_class(ground, forcing);
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_base_class(ground, grid, forcing, depths);
           
            %overwrite energy values assigned in base class
            ground.STATVAR.saltConc = ground.STATVAR.saltConc .* ground.STATVAR.waterIce; %total number of moles
            ground = get_E_water_salt_FreezeDepress_Xice(ground); %energy, water, ice, salt_c_brine
            ground = conductivity(ground);
            ground = diffusivity_salt(ground); % [m2/sec]
%             ground = initialize_STATVAR_from_file@GROUND_subsi(ground, forcing, grid);
        end
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            %assign upper boundary for heat in subclass 
            %assign upper boundary for salt in subclass
        end
        
        function ground = get_boundary_condition_l(ground) 
            ground = get_boundary_condition_l@GROUND_base_class(ground);
            %assign lower boundary for salt in subclass
        end
        
        function ground = get_derivatives_prognostic(ground)
            ground = get_derivative_energy(ground);
            ground = get_derivative_salt(ground); 
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables, in this case only check energy derivative
            timestep = get_timestep@GROUND_base_class(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground = advance_prognostic@GROUND_base_class(ground, timestep);
            ground.STATVAR.saltConc = max(0,ground.STATVAR.saltConc + timestep .*ground.TEMP.d_salt);
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
             %assigned in subclass
        end
        
        function ground = compute_diagnostic(ground, forcing) %function from base class fully overwritten
            ground = get_T_water_salt_FreezeDepress_Xice(ground);
            ground = conductivity(ground);
            ground = diffusivity_salt(ground); 
        end
        

        
        %non-mandatory functions -> required here so that they are usable
        %in subclasses
        
        function ground = get_derivative_salt(ground)
            fluxes = (ground.STATVAR.salt_c_brine(1:end-1) - ground.STATVAR.salt_c_brine(2:end)) .* ground.STATVAR.diffusivitySalt(1:end-1) .* ground.STATVAR.diffusivitySalt(2:end) ./...
                (ground.STATVAR.diffusivitySalt(1:end-1).* ground.STATVAR.layerThick(2:end)./2 +  ground.STATVAR.diffusivitySalt(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            %unit mole/m2, so flux through 1m2 unit cross section between cells!
            d_salt=ground.STATVAR.energy.*0;
            d_salt(1) = ground.TEMP.F_ub_salt - fluxes(1);
            d_salt(2:end-1) = fluxes(1:end-1) - fluxes(2:end);
            d_salt(end) = ground.TEMP.F_lb_salt + fluxes(end);
            
            ground.TEMP.d_salt = d_salt;
        end
        
        
        function ground = get_T_water_salt_FreezeDepress_Xice(ground)
            
            L_f = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            R = ground.CONST.R;
            Tmfw = ground.CONST.Tmfw;
            
            deltaT = ground.STATVAR.deltaT;
            %
            freeWaterIce = 0; %[m]  %change later for snow
            
            waterIce = ground.STATVAR.waterIce;
            mineral = ground.STATVAR.mineral;
            organic = ground.STATVAR.organic;
            %N = ground.STATVAR.saltConc;
            
            energy = ground.STATVAR.energy;
            
            
%             waterIce = ground.STATVAR.waterIce./ground.STATVAR.layerThick;
%             mineral = ground.STATVAR.mineral./ground.STATVAR.layerThick;
%             organic = ground.STATVAR.organic./ground.STATVAR.layerThick;
            N = ground.STATVAR.saltConc./ground.STATVAR.waterIce;
            
%             energy = ground.STATVAR.energy./ground.STATVAR.layerThick;
            
            A = 1 + c_i.*deltaT.*freeWaterIce./(waterIce.*L_f);
            
            A1 = c_w.*waterIce + c_m.*mineral + c_o.*organic;
            
            A2 = c_i.*waterIce + c_m.*mineral + c_o.*organic;
            
            A3 = (c_w-c_i).*waterIce;
            
            A4 = waterIce .* L_f;
            
            B = -L_f./(R.* Tmfw.^2);
            
            %quadratic equation in Tm, the onset of freezing T
            %a*Tm^2 + b*Tm + c = 0
            
            %zero-th order terms
            c = - N .* A2 .* deltaT - N .* A4;
            
            %first order terms
            b = - N .* A3 + B .*(energy ./ A + L_f .*freeWaterIce ./ A + A2 .* deltaT + A4);
            
            %second-order terms
            a = B.* (-c_i.*freeWaterIce./A - (c_i .* deltaT .* freeWaterIce .*A1) ./ (A .*A4) - A2);

            %Tm_2 = (-b + sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
            Tm_1 = (-b - sqrt(b.^2 - 4.*a.*c)) ./ (2.*a); %thius is the right branch!!

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
            
            freeIce = double(energy >= thresh2 & energy < thresh1) .* energy ./thresh2 .* freeWaterIce + double(energy < thresh2) .* freeWaterIce;
            freeWater = freeWaterIce - freeIce;

            ground.STATVAR.T=T;
            ground.STATVAR.ice=ice;
            ground.STATVAR.water=water;
            
%             ground.STATVAR.ice=ice .* ground.STATVAR.layerThick;
%             ground.STATVAR.water=water.* ground.STATVAR.layerThick;
            
            salt_c_brine = ground.STATVAR.saltConc ./ water;
            salt_c_brine(isnan(salt_c_brine))=0;
            ground.STATVAR.salt_c_brine = salt_c_brine;
        end
        
        
        function ground = get_E_water_salt_FreezeDepress_Xice(ground)
            
            L_f = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            R = ground.CONST.R;
            Tmfw = ground.CONST.Tmfw;
            
            deltaT = ground.STATVAR.deltaT;
            %
            freeWaterIce = 0; %[m]  %no Xice here
            
            waterIce = ground.STATVAR.waterIce;
            mineral = ground.STATVAR.mineral;
            organic = ground.STATVAR.organic;
%             N = ground.STATVAR.saltConc;
            
%             waterIce = ground.STATVAR.waterIce./ground.STATVAR.layerThick;
%             mineral = ground.STATVAR.mineral./ground.STATVAR.layerThick;
%             organic = ground.STATVAR.organic./ground.STATVAR.layerThick;
            N = ground.STATVAR.saltConc./ground.STATVAR.waterIce;
            
            T = ground.STATVAR.T;
            
            A = 1 + c_i.*deltaT.*freeWaterIce./(waterIce.*L_f);
            A1 = c_w.*waterIce + c_m.*mineral + c_o.*organic;
            A2 = c_i.*waterIce + c_m.*mineral + c_o.*organic;
            A3 = (c_w-c_i).*waterIce;
            A4 = waterIce .* L_f;
            B = -L_f./(R.* Tmfw.^2);
            
            %quadratic equation in Tm, the onset of freezing T
            %a*Tm^2 + b*Tm + c = 0
            
            %zero-th order terms
            c = N.*R.*Tmfw.^2./L_f.*deltaT;
            
            %first order terms
            b=T+deltaT;
            
            %second-order terms
            a=-1;
            
            Tm_1 = (-b + sqrt(b.^2 - 4.*a.*c)) ./ (2.*a); %this is the right branch!!
            %Tm_2 = (-b - sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
            
            thresh1 = 0;
            thresh2 = (-R.* Tmfw.^2 ./L_f).* N;
            
            E1 = A2 .* (Tm_1-deltaT) - A4;
            E2 = A1 .*Tm_1;
            
            thirdTerm= (N./(-L_f./(R.*Tmfw.^2).*Tm_1) .*(E2-E1) + E1);
            thirdTerm(isnan(thirdTerm) | thirdTerm==Inf)=0;
            
            energy = double(T>= thresh1) .* (c_w .*waterIce + c_w.*freeWaterIce + c_m.*mineral + c_o.*organic) .*T + ...
                double(T < thresh1 & T >= thresh2) .* ((c_w .*waterIce + c_i.*freeWaterIce + c_m.*mineral + c_o.*organic) .*T - freeWaterIce .*L_f) + ...
                double(T < thresh2) .* thirdTerm;
            
            ice = double(T <= thresh2) .* waterIce .* (1-(energy - E1)./(E2-E1));
            water = waterIce - ice;
            
            ground.STATVAR.energy=energy;
            ground.STATVAR.water = water;
            ground.STATVAR.ice = ice;
            
%             ground.STATVAR.energy=energy .* ground.STATVAR.layerThick;
%             ground.STATVAR.water = water .* ground.STATVAR.layerThick;
%             ground.STATVAR.ice = ice .* ground.STATVAR.layerThick;
            
            salt_c_brine = ground.STATVAR.saltConc ./ water;
            salt_c_brine(isnan(salt_c_brine))=0;
            ground.STATVAR.salt_c_brine = salt_c_brine;
        end
        
        
        function ground = get_derivative_energy(ground)
            ground = get_derivative_energy@GROUND_base_class(ground);
        end
        
        
        function ground = conductivity(ground)
            ground = conductivity@GROUND_base_class(ground);
        end
        
        
        
    end
end