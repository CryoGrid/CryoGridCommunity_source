classdef SALT < BASE
    
    methods
        
        function ground = get_boundary_condition_u_ZERO_SALT(ground)
            ground.TEMP.F_ub_salt=0; %zero flux bc, could be changed to Dirichlet BC, but this makes only sense when temperature bc is also Dirichlet
        end
        
        function ground = get_boundary_condition_l_ZERO_SALT(ground)
            ground.TEMP.F_lb_salt=0; %zero flux bc
        end
        

        %-----derivatives----------
        
        function ground = get_derivative_salt(ground)
            fluxes = (ground.STATVAR.salt_c_brine(1:end-1) - ground.STATVAR.salt_c_brine(2:end)) .* ground.STATVAR.diffusivitySalt(1:end-1) .* ground.STATVAR.diffusivitySalt(2:end) ./...
                (ground.STATVAR.diffusivitySalt(1:end-1).* ground.STATVAR.layerThick(2:end)./2 +  ground.STATVAR.diffusivitySalt(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            %unit mol/m2, so flux through 1m2 unit cross section between cells!
            d_salt=ground.STATVAR.energy.*0;
            d_salt(1) =  - fluxes(1);
            d_salt(2:end-1) = fluxes(1:end-1) - fluxes(2:end);  
            d_salt(end) =  + fluxes(end);
            
            d_salt = d_salt.*ground.STATVAR.area;  %multiply by area, in [mol/sec] 
            
            ground.TEMP.d_salt = ground.TEMP.d_salt + d_salt;
        end
        
        %-----------timesteps----------
        function timestep = get_timestep_salt(ground)
            timestep = ground.PARA.dt_max; %change later
        end
        
        %----diagnostic functions---------
        
        function ground = get_T_water_salt_fcSimple_Xice(ground)
            
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
            
            waterIce = ground.STATVAR.waterIce ./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
            organic = ground.STATVAR.organic ./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
            energy = ground.STATVAR.energy ./ ground.STATVAR.area ./ ground.STATVAR.layerThick ;
            N = ground.STATVAR.saltConc./ ground.STATVAR.area ./ground.STATVAR.layerThick; %concentration in grid cell mol/m3, not concentration in waterIce phase
            
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
            
            ground.STATVAR.T = T;
            ground.STATVAR.ice = ice .* ground.STATVAR.layerThick .*  ground.STATVAR.area;
            ground.STATVAR.water = water .* ground.STATVAR.area .* ground.STATVAR.layerThick;
            
            salt_c_brine = N ./ water;
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
            freeWaterIce = 0; %[m]  %no Xice here            

            deltaT = ground.STATVAR.deltaT;
            waterIce = ground.STATVAR.waterIce; %still in [-]
            mineral = ground.STATVAR.mineral;
            organic = ground.STATVAR.organic;
            area = ground.STATVAR.area;
            T = ground.STATVAR.T;
            
            N = ground.STATVAR.saltConc.*waterIce;  %this should be saltConc.*waterIce?  in the initialization, it should be the concentration in the brine water?
            %N = ground.STATVAR.saltConc./waterIce;  the original
            %formulation, this does not make sense at all for me
            
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
            
            ground.STATVAR.energy = energy .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.water = water .*  ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.ice = ice .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            
            ground.STATVAR.waterIce = ground.STATVAR.waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.organic = ground.STATVAR.organic .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.air = ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.waterIce - ground.STATVAR.mineral - ground.STATVAR.organic;
            
            ground.STATVAR.saltConc = N .* ground.STATVAR.layerThick .* ground.STATVAR.area;  %unit mol
            
            salt_c_brine = N ./ water;  %unit mol/m3
            salt_c_brine(isnan(salt_c_brine))=0;
            ground.STATVAR.salt_c_brine = salt_c_brine;
        end
        
        
        
        %---diffusivity salt--------------
        function ground = diffusivity_salt(ground)
            
            water = ground.STATVAR.water./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            
            D0 = ((6.06 + 9.60)/2  + max(ground.STATVAR.T, 0) .* (0.297  + 0.438)/2) .* 1e-10; %from Boudreau, B., 1997, Diagenetic Models and thier implementation, Springer, Berlin.
            %average between values for Na+ and Cl-
            
            ground.STATVAR.diffusivitySalt = D0 .* water ./ground.PARA.tortuosity.^2;
        end
        
        
    end
end

