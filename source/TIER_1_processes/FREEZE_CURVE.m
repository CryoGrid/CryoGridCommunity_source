%========================================================================
% CryoGrid TIER1 library class for functions related to the freeze curve
% caclulated by the freezing = drying assumption, Dall'Amico et al., 2011
% S. Westermann, October 2020
%========================================================================


classdef FREEZE_CURVE < BASE
    
    
    properties
        LUT   %additional variable LUT for LookUpTable 
    end
        
    methods
        
        
        %----diagnostic functions---------
        function ground = get_T_water_freezeC(ground)
            
            L_sl = ground.CONST.L_f ./ ground.CONST.rho_w; %---------------
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            Tmfw = ground.CONST.Tmfw;
            
            energy = ground.STATVAR.energy ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            waterIce = ground.STATVAR.waterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            soil_type = ground.STATVAR.soil_type;
            
            n = ground.STATVAR.n;
            alpha = ground.STATVAR.alpha;
            thetaRes = ground.STATVAR.thetaRes;
        
            m=1-1./n;
            porosity = 1-mineral-organic;
            
            waterPotZero = -1./alpha .*(((waterIce - thetaRes)./(porosity - thetaRes)).^(-1./m) -1).^(1./n);
            waterPotZero = real(waterPotZero);

            Tstar =  g .* Tmfw ./ L_sl .* waterPotZero;
            
            unfrozen_positive = energy >= 0;
            unfrozen_negative = energy < 0 & energy >= Tstar .* (waterIce .* c_i + mineral .* c_m + organic .* c_o);
            linear_range = energy >= Tstar .* (waterIce .* c_i + mineral .* c_m + organic .* c_o);
            
            ground.STATVAR.waterPotential(linear_range)  = waterPotZero(linear_range);
            ground.STATVAR.water(linear_range) = waterIce(linear_range);
            ground.STATVAR.T(unfrozen_positive) = energy(unfrozen_positive) ./ (waterIce(unfrozen_positive) .* c_w + mineral(unfrozen_positive) .* c_m + organic(unfrozen_positive) .* c_o);
            ground.STATVAR.T(unfrozen_negative) = energy(unfrozen_negative) ./ (waterIce(unfrozen_negative) .* c_i + mineral(unfrozen_negative) .* c_m + organic(unfrozen_negative) .* c_o);
            
            %interpolate in look-up tables
            Tstar = Tstar(~linear_range) + Tmfw;
            a = g .* Tstar ./ L_sl;
            A = Tstar - Tmfw - a .*waterPotZero(~linear_range);
            b = waterIce(~linear_range) .* c_i + mineral(~linear_range) .* c_m + organic(~linear_range) .* c_o;
            B = rho_w .* L_sl .* (waterIce(~linear_range) - thetaRes(~linear_range));
            C = rho_w .* L_sl .* (porosity(~linear_range) - thetaRes(~linear_range));
            c = -1 ./ alpha(~linear_range);
            soil_type = ground.STATVAR.soil_type(~linear_range)-1;
            size_LUT = ground.PARA.LUT_size_waterIce .* ground.PARA.LUT_size_T;
            
            energy_shifted_and_scaled = real((energy(~linear_range) - b.*A + B) ./ (a.*b.*c));
            scale_factor = real(C ./ (a.*b.*c));
            
            sf_matrix = ground.LUT.sf_matrix;
            LUT = ground.LUT.lut_energy_T_water;
            
            pos_sf = (scale_factor - ground.LUT.min_sf(soil_type + 1,1)) ./ (ground.LUT.max_sf(soil_type + 1,1) - ground.LUT.min_sf(soil_type + 1,1)) .* ground.PARA.LUT_size_waterIce;
            pos_sf(pos_sf<1) = 1;
            pos_sf(pos_sf > ground.PARA.LUT_size_waterIce -1) = ground.PARA.LUT_size_waterIce - 1;
            fraction_sf = pos_sf - floor(pos_sf);
            pos_sf = floor(pos_sf);
            
            min_sf_interp = sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,2) + ...
                fraction_sf .* (sf_matrix(pos_sf + 1 + soil_type.*ground.PARA.LUT_size_waterIce,2) - sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,2));
            max_sf_interp = sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,3) + ...
                fraction_sf .* (sf_matrix(pos_sf + 1 + soil_type.*ground.PARA.LUT_size_waterIce,3) - sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,3));

            pos_energy = (energy_shifted_and_scaled - min_sf_interp) ./ (max_sf_interp - min_sf_interp) .* ground.PARA.LUT_size_T;
            pos_energy(pos_energy<1) = 1;
            pos_energy(pos_energy > ground.PARA.LUT_size_T -1) = ground.PARA.LUT_size_T-1;
            fraction_energy = pos_energy - floor(pos_energy);
            pos_energy = floor(pos_energy);

            left = LUT(pos_sf + (pos_energy-1).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) + ...
                fraction_sf .* (LUT(pos_sf+1 + (pos_energy-1).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) - ...
                LUT(pos_sf + (pos_energy-1).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT));
            right = LUT(pos_sf + (pos_energy).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) + ...
                fraction_sf .* (LUT(pos_sf+1 + (pos_energy).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) - ...
                LUT(pos_sf + (pos_energy).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT));

            X_interp = left + fraction_energy .*(right - left);
            
            waterPot_interp = X_interp./-alpha(~linear_range);
            T_interp = (waterPot_interp-waterPotZero(~linear_range)).* Tstar.*g ./L_sl + Tstar - Tmfw;
            
            water_interp  = (thetaRes(~linear_range) + (porosity(~linear_range) - thetaRes(~linear_range)) .* (1+(X_interp).^n(~linear_range)).^-m(~linear_range));
            water_interp = real(water_interp);
            
            ground.STATVAR.waterPotential(~linear_range)  = waterPot_interp;
            ground.STATVAR.water(~linear_range) = water_interp;
            ground.STATVAR.T(~linear_range) = T_interp;
           
            ground.STATVAR.water(ground.STATVAR.water > waterIce) = waterIce(ground.STATVAR.water > waterIce); %elminates small rounding errors and avoids negative ice contents 

            ground.STATVAR.water = ground.STATVAR.water .* ground.STATVAR.layerThick .* ground.STATVAR.area;            
            ground.STATVAR.ice = ground.STATVAR.waterIce - ground.STATVAR.water;
            
%             %TEST SEBAS
%             ground.STATVAR.waterPotentialExp = waterPotZero;
%             ground.STATVAR.waterPotentialExp(~linear_range) = ground.STATVAR.waterPotentialExp(~linear_range) + waterPot_interp;
%             
%             ground.STATVAR.waterPotential = ground.STATVAR.waterPotentialExp;
% %             saturation_water = ground.STATVAR.water ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
% %             saturation_ice = ground.STATVAR.ice ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
% %             %satuaration_water(~linear_range) ./ (1- satuaration_ice(~linear_range));
% %             ground.STATVAR.waterPotentialExp(~linear_range) = real(-1./alpha(~linear_range) .*((saturation_water(~linear_range) ./ (1- saturation_ice(~linear_range))).^(-1./m(~linear_range)) -1).^(1./n(~linear_range)));
% %             
% %             ground.STATVAR.waterPotentialExp2 = double(ground.STATVAR.T<0) .* (ground.STATVAR.T)./273.15 .* ground.CONST.L_f ./ ground.CONST.g ./ ground.CONST.rho_w;
% %             ground.STATVAR.waterPotentialExp3 = ground.STATVAR.waterPotential;
% %             ground.STATVAR.waterPotentialExp3(~linear_range) = ground.STATVAR.waterPotentialExp3(~linear_range) - ground.STATVAR.waterPotentialExp(~linear_range);
%             
        end
        
        
        
        function ground = get_T_water_freezeC_Xice(ground)
            
            L_sl = ground.CONST.L_f ./ ground.CONST.rho_w;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            Tmfw = ground.CONST.Tmfw;
            
            %waterPotZero and T_star must be known here
            energy = ground.STATVAR.energy;
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;  %volumetric fraction of the total volume incl Xice 
            organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            waterIce = ground.STATVAR.waterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            XwaterIce = ground.STATVAR.XwaterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            
            n = ground.STATVAR.n;
            alpha = ground.STATVAR.alpha;
            thetaRes = ground.STATVAR.thetaRes;
            
            m=1-1./n;
            porosity = 1 - mineral - organic - XwaterIce;
            waterPotZero = real(-1./alpha .*(((waterIce - thetaRes)./(porosity - thetaRes)).^(-1./m)-1).^(1./n));
            Tstar =  g .* Tmfw ./ L_sl .* waterPotZero;            
            
            %distinguish four different regimes
            %1. linear positive unfrozen
            unfrozen_positive = (energy >= 0);
            ground.STATVAR.waterPotential(unfrozen_positive)  = waterPotZero(unfrozen_positive);
            ground.STATVAR.water(unfrozen_positive) = ground.STATVAR.waterIce(unfrozen_positive);
            ground.STATVAR.Xwater(unfrozen_positive) = ground.STATVAR.XwaterIce(unfrozen_positive);
            ground.STATVAR.T(unfrozen_positive) = energy(unfrozen_positive) ./ ((ground.STATVAR.waterIce(unfrozen_positive) + ground.STATVAR.XwaterIce(unfrozen_positive)) .* c_w + ...
                ground.STATVAR.mineral(unfrozen_positive) .* c_m + ground.STATVAR.organic(unfrozen_positive) .* c_o);

            %2. Xice melt regime, T=0, calculate Xice and Xwater, "normal soil" unfrozen
            Xice_melt = (energy <0 & energy >= - ground.CONST.L_f .* ground.STATVAR.XwaterIce);
            ground.STATVAR.waterPotential(Xice_melt)  = waterPotZero(Xice_melt);
            ground.STATVAR.water(Xice_melt) = ground.STATVAR.waterIce(Xice_melt);
            ground.STATVAR.Xwater(Xice_melt) = ground.STATVAR.XwaterIce(Xice_melt) .*(1 + energy(Xice_melt)./ (ground.CONST.L_f .* ground.STATVAR.XwaterIce(Xice_melt)));
            ground.STATVAR.T(Xice_melt) = 0;
           
            energy = energy + ground.CONST.L_f .* ground.STATVAR.XwaterIce; %subtract the latent part of the Xice energy (-Lf * XwaterIce) from energy
            %and treat the Xice part in exactly the same fashion as mineral and organic
            
            %3. linear negative,  XwaterIce fully frozen
            unfrozen_negative = (energy < 0 & energy >= Tstar .* ((ground.STATVAR.waterIce + ground.STATVAR.XwaterIce) .* c_i + ...
                ground.STATVAR.mineral .* c_m + ground.STATVAR.organic .* c_o));
            ground.STATVAR.waterPotential(unfrozen_negative)  = waterPotZero(unfrozen_negative);
            ground.STATVAR.water(unfrozen_negative) = ground.STATVAR.waterIce(unfrozen_negative);
            ground.STATVAR.Xwater(unfrozen_negative) = 0;
            ground.STATVAR.T(unfrozen_negative) = energy(unfrozen_negative) ./ ((ground.STATVAR.waterIce(unfrozen_negative) + ground.STATVAR.XwaterIce(unfrozen_negative)) .* c_i + ...
                ground.STATVAR.mineral(unfrozen_negative) .* c_m + ground.STATVAR.organic(unfrozen_negative) .* c_o);
            
            %4. freeze curve region, XwaterIce fully frozen
            freeze_curve = energy < Tstar .* ((ground.STATVAR.waterIce + ground.STATVAR.XwaterIce) .* c_i + ground.STATVAR.mineral .* c_m + ground.STATVAR.organic .* c_o);
            
            %transform to [J/m3], at this point XwaterIce is treated like mineral/organic, latent heat of XwaterIce is already subtracted
            energy = energy(freeze_curve) ./ ground.STATVAR.layerThick(freeze_curve) ./ ground.STATVAR.area(freeze_curve);

            %interpolate in look-up tables
            Tstar = Tstar(freeze_curve) + Tmfw;
            a = g .* Tstar ./ L_sl;
            A = Tstar - Tmfw - a .*waterPotZero(freeze_curve);
            b = (waterIce(freeze_curve) + XwaterIce(freeze_curve)) .* c_i + mineral(freeze_curve) .* c_m + organic(freeze_curve) .* c_o;
            B = rho_w .* L_sl .* (waterIce(freeze_curve) - thetaRes(freeze_curve));
            C = rho_w .* L_sl .* (porosity(freeze_curve) - thetaRes(freeze_curve));
            c = -1 ./ alpha(freeze_curve);
            soil_type = ground.STATVAR.soil_type(freeze_curve) - 1;
            size_LUT = ground.PARA.LUT_size_waterIce .* ground.PARA.LUT_size_T;
            
            energy_shifted_and_scaled = real((energy - b.*A + B) ./ (a.*b.*c));
            scale_factor = real(C ./ (a.*b.*c));

            sf_matrix = ground.LUT.sf_matrix;
            LUT = ground.LUT.lut_energy_T_water;
            
            pos_sf = (scale_factor - ground.LUT.min_sf(soil_type + 1,1)) ./ (ground.LUT.max_sf(soil_type + 1,1) - ground.LUT.min_sf(soil_type + 1,1)) .* ground.PARA.LUT_size_waterIce;
            pos_sf(pos_sf<1) = 1;
            pos_sf(pos_sf > ground.PARA.LUT_size_waterIce -1) = ground.PARA.LUT_size_waterIce - 1;
            fraction_sf = pos_sf - floor(pos_sf);
            pos_sf = floor(pos_sf);
            
            min_sf_interp = sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,2) + ...
                fraction_sf .* (sf_matrix(pos_sf + 1 + soil_type.*ground.PARA.LUT_size_waterIce,2) - sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,2));
            max_sf_interp = sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,3) + ...
                fraction_sf .* (sf_matrix(pos_sf + 1 + soil_type.*ground.PARA.LUT_size_waterIce,3) - sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,3));

            pos_energy = (energy_shifted_and_scaled - min_sf_interp) ./ (max_sf_interp - min_sf_interp) .* ground.PARA.LUT_size_T;
            pos_energy(pos_energy<1) = 1;
            pos_energy(pos_energy > ground.PARA.LUT_size_T -1) = ground.PARA.LUT_size_T-1;
            fraction_energy = pos_energy - floor(pos_energy);
            pos_energy = floor(pos_energy);

            left = LUT(pos_sf + (pos_energy-1).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) + ...
                fraction_sf .* (LUT(pos_sf+1 + (pos_energy-1).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) - ...
                LUT(pos_sf + (pos_energy-1).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT));
            right = LUT(pos_sf + (pos_energy).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) + ...
                fraction_sf .* (LUT(pos_sf+1 + (pos_energy).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) - ...
                LUT(pos_sf + (pos_energy).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT));

            X_interp = left + fraction_energy .*(right - left);
            
            waterPot_interp = X_interp./-alpha(freeze_curve);
            T_interp = (waterPot_interp-waterPotZero(freeze_curve)).* Tstar.*g ./L_sl + Tstar - Tmfw;
            water_interp  = real(thetaRes(freeze_curve) + (porosity(freeze_curve) - thetaRes(freeze_curve)) .* (1+(X_interp).^n(freeze_curve)).^-m(freeze_curve));
            
            ground.STATVAR.waterPotential(freeze_curve)  = waterPot_interp;
            ground.STATVAR.water(freeze_curve) = water_interp .* ground.STATVAR.layerThick(freeze_curve) .* ground.STATVAR.area(freeze_curve);
            ground.STATVAR.water(ground.STATVAR.water > ground.STATVAR.waterIce) = ground.STATVAR.waterIce(ground.STATVAR.water > ground.STATVAR.waterIce); %elminates small rounding errors and avoids negative ice contents 
         
            ground.STATVAR.Xwater(freeze_curve) = 0;
            ground.STATVAR.T(freeze_curve) = T_interp;
                    
            ground.STATVAR.ice = ground.STATVAR.waterIce - ground.STATVAR.water;
            ground.STATVAR.Xice = ground.STATVAR.XwaterIce - ground.STATVAR.Xwater;
        end
        
        
        
        
        % get energy from temeprature and water contents, normally part of initializiation
        
        function ground = get_E_freezeC(ground) %required for initialization
        
            L_sl = ground.CONST.L_f ./ ground.CONST.rho_w ;%--------------;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            Tmfw = ground.CONST.Tmfw;
            
            T = ground.STATVAR.T;
            mineral= ground.STATVAR.mineral;
            organic = ground.STATVAR.organic;
            waterIce = ground.STATVAR.waterIce;
            layerThick = ground.STATVAR.layerThick;
            area = ground.STATVAR.area;
            soil_type = ground.STATVAR.soil_type;
            
            n = double(soil_type == 1) .* ground.CONST.vanGen_n(1) + double(soil_type == 2) .* ground.CONST.vanGen_n(2) + double(soil_type == 3) .* ground.CONST.vanGen_n(3) + double(soil_type == 4) .* ground.CONST.vanGen_n(4) + double(soil_type == 5) .* ground.CONST.vanGen_n(5);
            ground.STATVAR.n = n;
            alpha = double(soil_type == 1) .* ground.CONST.vanGen_alpha(1) + double(soil_type == 2) .* ground.CONST.vanGen_alpha(2) + double(soil_type == 3) .* ground.CONST.vanGen_alpha(3) + double(soil_type == 4) .* ground.CONST.vanGen_alpha(4) + + double(soil_type == 5) .* ground.CONST.vanGen_alpha(5);
            ground.STATVAR.alpha = alpha;
            thetaRes = double(soil_type == 1) .* ground.CONST.vanGen_residual_wc(1) + double(soil_type == 2) .* ground.CONST.vanGen_residual_wc(2) + double(soil_type == 3) .* ground.CONST.vanGen_residual_wc(3) + double(soil_type == 4) .* ground.CONST.vanGen_residual_wc(4) +  double(soil_type == 5) .* ground.CONST.vanGen_residual_wc(5);
            ground.STATVAR.thetaRes = thetaRes;
                        
            porosity = 1-mineral-organic;
            m=1-1./n;
            
            waterPotZero = real( -1./alpha .*(((waterIce - thetaRes)./(porosity - thetaRes)).^(-1./m)-1).^(1./n));
            Tstar = Tmfw + g .* Tmfw ./ L_sl .* waterPotZero;
            
            waterPot = waterPotZero + (L_sl./g./Tstar .* (T - Tstar + Tmfw)).*double(T < Tstar - Tmfw);
            waterC  = double(T < Tstar - Tmfw) .* real(thetaRes + (porosity - thetaRes) .* (1+(-alpha.*waterPot).^n).^-m) + double(T >= Tstar-Tmfw) .* waterIce;
            
            energy = double (T>=0) .* (waterIce .* c_w + mineral .* c_m + organic .* c_o) .*T + double(T<0) .* ((waterIce .* c_i + mineral .* c_m + organic .* c_o) .*T - rho_w .*L_sl .* (waterIce - waterC));
                  
            ground.STATVAR.water = waterC .* layerThick .* area;  % [m3]
            ground.STATVAR.ice = (waterIce - waterC) .*  layerThick .* area; %[m3]          
            
            ground.STATVAR.waterIce = waterIce .* layerThick .* area; % [m3]
            ground.STATVAR.mineral = mineral .* layerThick .* area; % [m3]
            ground.STATVAR.organic = organic .* layerThick .* area; % [m3]
            ground.STATVAR.energy = energy .* layerThick .* area;  % [J]

            ground.STATVAR.air = (1-mineral-organic-waterIce) .* layerThick .* area;  % [m3]
            
            ground.STATVAR.waterPotential = waterPot;
        end
        
        
        
        function ground = get_E_freezeC_Xice(ground) %required for initialization
        
            L_sl = ground.CONST.L_f ./ ground.CONST.rho_w;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            Tmfw = ground.CONST.Tmfw;
            
            T = ground.STATVAR.T;
            mineral= ground.STATVAR.mineral; %properties of the matrix when Xice is removed!!!
            organic = ground.STATVAR.organic;
            waterIce = ground.STATVAR.waterIce;
            layerThick = ground.STATVAR.layerThick;
            area = ground.STATVAR.area;
            soil_type = ground.STATVAR.soil_type;
            
            Xice = ground.STATVAR.Xice .* double(T <= 0); % Xice initially only possible when frozen, provided in multiples of the "matrix", 1 means 50 vol% Xice, 50 % normal soil

            n = double(soil_type == 1) .* ground.CONST.vanGen_n(1) + double(soil_type == 2) .* ground.CONST.vanGen_n(2) + double(soil_type == 3) .* ground.CONST.vanGen_n(3) + double(soil_type == 4) .* ground.CONST.vanGen_n(4) + double(soil_type == 5) .* ground.CONST.vanGen_n(5);
            ground.STATVAR.n = n;
            alpha = double(soil_type == 1) .* ground.CONST.vanGen_alpha(1) + double(soil_type == 2) .* ground.CONST.vanGen_alpha(2) + double(soil_type == 3) .* ground.CONST.vanGen_alpha(3) + double(soil_type == 4) .* ground.CONST.vanGen_alpha(4) + + double(soil_type == 5) .* ground.CONST.vanGen_alpha(5);
            ground.STATVAR.alpha = alpha;
            thetaRes = double(soil_type == 1) .* ground.CONST.vanGen_residual_wc(1) + double(soil_type == 2) .* ground.CONST.vanGen_residual_wc(2) + double(soil_type == 3) .* ground.CONST.vanGen_residual_wc(3) + double(soil_type == 4) .* ground.CONST.vanGen_residual_wc(4) +  double(soil_type == 5) .* ground.CONST.vanGen_residual_wc(5);
            ground.STATVAR.thetaRes = thetaRes;
            
            porosity = 1-mineral-organic;
            m=1-1./n;
            
            waterPotZero = real(-1./alpha .*(((waterIce - thetaRes)./(porosity - thetaRes)).^(-1./m)-1).^(1./n));
            Tstar = Tmfw + g .* Tmfw ./ L_sl .* waterPotZero;
            
            waterPot = waterPotZero + (L_sl./g./Tstar .* (T - Tstar + Tmfw)).*double(T < Tstar - Tmfw);
            waterC  = double(T < Tstar - Tmfw) .* real(thetaRes + (porosity - thetaRes) .* (1+(-alpha.*waterPot).^n).^-m) + double(T >= Tstar-Tmfw) .* waterIce;
            
            energy = double (T>=0) .* (waterIce .* c_w + mineral .* c_m + organic .* c_o) .*T + double(T<0) .* ((waterIce .* c_i + mineral .* c_m + organic .* c_o) .*T - rho_w .*L_sl .* (waterIce - waterC));
                  
            ground.STATVAR.water = waterC .* layerThick ./ (1 + Xice) .* area;  % [m3]
            ground.STATVAR.ice = (waterIce - waterC) .*  layerThick ./ (1 + Xice) .* area; %[m3]          
            
            ground.STATVAR.waterIce = waterIce .* layerThick ./ (1 + Xice) .* area; % [m3]
            ground.STATVAR.mineral = mineral .* layerThick ./ (1 + Xice) .* area; % [m3]
            ground.STATVAR.organic = organic .* layerThick ./ (1 + Xice) .* area; % [m3]
            ground.STATVAR.energy = (energy + Xice .* (ground.STATVAR.T .* ground.CONST.c_i - ground.CONST.L_f)) ./ (1 + Xice) .* layerThick .* area;  % [J]
            ground.STATVAR.XwaterIce = Xice ./ (1 + Xice) .* layerThick .* area;
            ground.STATVAR.Xice = ground.STATVAR.XwaterIce .*double(T<=0);
            ground.STATVAR.Xwater = ground.STATVAR.XwaterIce .* double(T>0);

            ground.STATVAR.air = (1-mineral-organic-waterIce) .* layerThick ./ (1 + Xice) .* area;  % [m3]
            
            ground.STATVAR.waterPotential = waterPot;
        end
        
        
        
        %---look-up tables initialization-----------
        function ground = create_LUT_freezeC(ground) % creates lookup table LUT
            
            disp('creating look-up tables')
            
            L_sl = ground.CONST.L_f ./ ground.CONST.rho_w ; %---------------
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            Tmfw = ground.CONST.Tmfw;
            
            ground.LUT.min_sf = [];
            ground.LUT.max_sf = [];
            ground.LUT.sf_matrix = [];
            ground.LUT.lut_energy_T_water = [];
            
            for soil_type = 1:size(ground.CONST.vanGen_alpha,2)
                
                alpha = ground.CONST.vanGen_alpha(1, soil_type);
                n = ground.CONST.vanGen_n(1, soil_type);
                thetaRes = ground.CONST.vanGen_residual_wc(1, soil_type);
                
                store =[];
                min_sf = 1e20;
                max_sf = -1e20;
                for j=0:1  %max and min matrix fill
                    for waterIce = ground.PARA.min_waterIce:0.005:ground.PARA.max_waterIce

                        mineral = double(j).*(1 - waterIce) + double(~j) .* ground.PARA.min_mineral_organic;
                        organic = 0;
                        porosity = 1-mineral-organic;
                        
                        m=1-1./n;
                        
                        T=[ground.PARA.min_T 0];
                        
                        
                        waterPotZero = -1./alpha .*(((waterIce - thetaRes)./(porosity - thetaRes)).^(-1./m)-1).^(1./n);
                        if j==1 || waterIce + mineral == 1 %correct small rounding errors for waterPotZero = 0
                            waterPotZero=0;
                        end
                        Tstar = Tmfw + g .* Tmfw ./ L_sl .* waterPotZero;
                        
                        
                        waterPot = waterPotZero + (L_sl./g./Tstar .* (T - Tstar + Tmfw)).*double(T < Tstar - Tmfw);
                        waterC  = double(T < Tstar - Tmfw) .* (thetaRes + (porosity - thetaRes) .* (1+(-alpha.*waterPot).^n).^-m) + double(T >= Tstar-Tmfw) .* waterIce;
                        
                        energy = double (T>=0) .* (waterIce .* c_w + mineral .* c_m + organic .* c_o) .*T + double(T<0) .* ((waterIce .* c_i + mineral .* c_m + organic .* c_o) .*T - rho_w .*L_sl .* (waterIce - waterC));
                        
                        X = real(-alpha .*waterPot);
                        
                        a = g .* Tstar ./ L_sl;
                        A = Tstar - 273.15 - a .*waterPotZero;
                        b = waterIce .* c_i + mineral .* c_m + organic .* c_o;
                        B = rho_w .* L_sl .* (waterIce - thetaRes);
                        C = rho_w .* L_sl .* (porosity - thetaRes);
                        c = -1 ./ alpha;
                        
                        energy_shifted_and_scaled = real((energy - b.*A + B) ./ (a.*b.*c));
                        scale_factor = real(C ./ (a.*b.*c));  %abbreviated sf, this is the variable that contains the water content
                        
                        store = [store ; [waterIce T scale_factor energy_shifted_and_scaled X]];
                        
                        min_sf = min(min_sf, scale_factor);
                        max_sf = max(max_sf, scale_factor);
                        
                    end
                end
                
                %delete duplicates
                i=1;
                while i<=size(store,1)
                    if sum(store(i,4)-store(:,4)==0)>1
                        store(i,:)=[];
                        i=1;
                    end
                    i=i+1;
                end

                sf_matrix = [];
                
                for sf = linspace(min_sf, max_sf, ground.PARA.LUT_size_waterIce)
                    %sf_matrix = [sf_matrix; [sf interp1(store(:,4), store(:,6), sf) interp1(store(:,4), store(:,5), sf) interp1(store(:,4), store(:,8), sf) interp1(store(:,4), store(:,7), sf)]];
                    sf_matrix = [sf_matrix; [sf sf interp1(store(:,4), store(:,5), sf) 0 interp1(store(:,4), store(:,7), sf)]];

                end
                
                %first cell scale factor
                %second cell: minimum of energy_shifted_and_scaled
                %third cell: maximum of energy_shifted_and_scaled
                %forth cell: minimum of target variable X, related to T and matrix potential
                %fifth cell: maximum of target variable X, related to T and matrix potential
                
                LUT=[];
                
                for i=1:size(sf_matrix,1)
                    
                    X_values = linspace(sf_matrix(i,4), sf_matrix(i,5)+100, 10000);
                    energy_shifted_and_scaled_values = real(X_values + sf_matrix(i,1) .*(1 + X_values.^n).^-m);
                    
                    LUT=[LUT; interp1(energy_shifted_and_scaled_values, X_values, linspace(sf_matrix(i,2), sf_matrix(i,3), ground.PARA.LUT_size_T))];
                    
                end
                
                ground.LUT.min_sf = cat(1, ground.LUT.min_sf, min_sf);
                ground.LUT.max_sf = cat(1, ground.LUT.max_sf, max_sf);
                ground.LUT.sf_matrix = cat(1, ground.LUT.sf_matrix, sf_matrix);
                ground.LUT.lut_energy_T_water = cat(1, ground.LUT.lut_energy_T_water, LUT(:));
                
            end
            
        end

    end
end

