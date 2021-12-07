%========================================================================
% CryoGrid TIER1 library class for functions related to the freeze curve
% of Karra and Painter 2014
% S. Westermann, December 2020
%========================================================================


classdef FREEZE_CURVE_KarraPainter < BASE
    
    
    properties
        LUT   %additional variable LUT for LookUpTable
    end
    
    methods
        
        
        %----diagnostic functions---------

        
        function ground = get_T_water_freezeC(ground)
            
            sat_waterIce_min = 0.005;
            LUT_size_gamma = 2.^9;
            LUT_size_T = 2.^10;
            porosity_max = 0.95;
            porosity_min = 0.05;
            
            beta_interface = 2.2;  %make this a constant!
            
            L_sl = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            T0 = ground.CONST.Tmfw;
            
            energy = ground.STATVAR.energy ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            waterIce = ground.STATVAR.waterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            soil_type = ground.STATVAR.soil_type;
            
            n = ground.STATVAR.n;
            alpha = ground.STATVAR.alpha ./ rho_w ./ g;
            
            m=1-1./n;
            porosity = 1-mineral-organic;
            sat_waterIce = waterIce ./ porosity;
            
            %matric potential (Pa)
            mwp0 = real(1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n));
            
            linear_range = energy >= 0;
            
            %linear range 
            ground.STATVAR.waterPotential(linear_range)  = -mwp0(linear_range) ./ rho_w ./ g;
            ground.STATVAR.water(linear_range) = waterIce(linear_range);
            ground.STATVAR.T(linear_range) = energy(linear_range) ./ (waterIce(linear_range) .* c_w + mineral(linear_range) .* c_m + organic(linear_range) .* c_o);
            %done!
            
            if sum(double(~linear_range)) > 0
                %non-linear range
                alpha = alpha(~linear_range);
                n= n(~linear_range);
                m=m(~linear_range);
                mwp0 = mwp0(~linear_range);
                porosity = porosity(~linear_range);
                sat_waterIce = sat_waterIce(~linear_range);
                
                C0 = mineral(~linear_range) .* c_m + organic(~linear_range) .* c_o + porosity .* sat_waterIce .* c_i; %vector
                X = -L_sl./ T0 .* beta_interface;
                L0 = porosity.* L_sl; %vector
                gamma = C0 ./ alpha ./ X ./ L0 ; %vector
                
                soil_type = ground.STATVAR.soil_type(~linear_range); %-1;
                
                E_prime = energy(~linear_range) ./ L0 + sat_waterIce + C0 .* mwp0 ./ X./ L0; %takes soil type into account, alpha is assigned correctly
                
                T_nonlinear = E_prime.*0;
                satWater_nonlinear = E_prime.*0;
                matric_pot_nonlinear = E_prime.*0;
                
                %find correct gamma
                depth_gamma_LUT = log(LUT_size_gamma)./log(2);
                gamma_left_index = (soil_type-1) .* (LUT_size_gamma+1) + 1 ; %start with first index that belongs to correct soil code
                for i = depth_gamma_LUT-1:-1:0
                    gamma_left_index = gamma_left_index + 2.^i .* double(gamma >= ground.LUT.gamma(gamma_left_index + 2.^i,1));
                end
                scale_factor_gamma = (gamma - ground.LUT.gamma(gamma_left_index, 1)) ./ (ground.LUT.gamma(gamma_left_index+1, 1) - ground.LUT.gamma(gamma_left_index, 1));
                
                %             ground.STATVAR.gamma = gamma;
                %             ground.STATVAR.gamma_left_index = gamma_left_index;
                %             ground.STATVAR.scale_factor_gamma = scale_factor_gamma;
                
                %left side fit E_prime
                depth_T_LUT = log(LUT_size_T)./log(2);
                %T_left_index = (soil_type-1) .* (LUT_size_gamma+1) .* (LUT_size_T+1) + (gamma_left_index-1) .* (LUT_size_T+1) + 1; %soil_type offset + gamma offset
                T_left_index = (gamma_left_index-1) .* (LUT_size_T+1) + 1; %soil_type offset already contained in gamma offset
                for i = depth_T_LUT-1:-1:0
                    T_left_index = T_left_index + 2.^i .* double(E_prime >= ground.LUT.lut_E_prime(T_left_index + 2.^i,1));
                end
                scale_factor_T = (E_prime - ground.LUT.lut_E_prime(T_left_index, 1)) ./ (ground.LUT.lut_E_prime(T_left_index+1, 1) - ground.LUT.lut_E_prime(T_left_index, 1));
                T_prime_left = ground.LUT.lut_T_prime(T_left_index, 1) + scale_factor_T .* (ground.LUT.lut_T_prime(T_left_index+1, 1) - ground.LUT.lut_T_prime(T_left_index, 1));
                
                %             ground.STATVAR.T_left_index = T_left_index;
                %             ground.STATVAR.T_prime_left = T_prime_left;
                
                %right side fit E_prime
                %T_right_index = (soil_type-1) .* (LUT_size_gamma+1) .* (LUT_size_T+1) + (gamma_left_index+1-1) .* (LUT_size_T+1) + 1; %soil_type offset + gamma offset for gamma right side
                T_right_index = (gamma_left_index + 1 - 1) .* (LUT_size_T+1) + 1; %soil_type offset already contained in gamma offset
                for i = depth_T_LUT-1:-1:0
                    T_right_index = T_right_index + 2.^i .* double(E_prime >= ground.LUT.lut_E_prime(T_right_index + 2.^i,1));
                end
                scale_factor_T = (E_prime - ground.LUT.lut_E_prime(T_right_index, 1)) ./ (ground.LUT.lut_E_prime(T_right_index+1, 1) - ground.LUT.lut_E_prime(T_right_index, 1));
                T_prime_right = ground.LUT.lut_T_prime(T_right_index, 1) + scale_factor_T .* (ground.LUT.lut_T_prime(T_right_index+1, 1) - ground.LUT.lut_T_prime(T_right_index, 1));
                
                %             ground.STATVAR.T_right_index = T_right_index;
                %             ground.STATVAR.T_prime_right = T_prime_right;
                
                T_prime = T_prime_left + scale_factor_gamma .* (T_prime_right - T_prime_left);
                
                %special fit for low water/high matric water potential grid cells, this is accomplished
                %by linearizing the dimensionless equation in T_prime around the
                %matric water potential intercept and then solving for
                %T_prime
                low_water_range = real(log(mwp0 .* alpha)) > -0.364.* real(log(-gamma))+6.65;
                %[real(log(mwp0 .* alpha)) -0.34.* real(log(-gamma))+7.0 -0.364.* real(log(-gamma))+6.65]
%                 low_water_range = gamma<0
                T_prime0 = alpha(low_water_range) .*mwp0(low_water_range);
                T_prime(low_water_range) = (E_prime(low_water_range) - (1 + T_prime0.^n(low_water_range)).^(-m(low_water_range)) - T_prime0 .* (n(low_water_range).*m(low_water_range)).* T_prime0.^(n(low_water_range)-1) ...
                    .* (1+T_prime0.^n(low_water_range)).^(-m(low_water_range)-1)) ./ (gamma(low_water_range) -  n(low_water_range).*m(low_water_range).* T_prime0.^(n(low_water_range)-1) .* (1+T_prime0.^n(low_water_range)).^(-m(low_water_range)-1));     


                T_nonlinear = (T_prime ./ alpha - mwp0) ./ X;
                satWater_nonlinear = E_prime - gamma .* T_prime;
                matric_pot_nonlinear = T_prime ./ alpha;
                
                
                %             ground.STATVAR.waterPotential2  = -matric_pot_nonlinear ./ rho_w ./ g;
                %             ground.STATVAR.water2 = satWater_nonlinear .* porosity ;
                %             ground.STATVAR.T2 = T_nonlinear;
                
                %ground.STATVAR.waterPotential(~linear_range)  = -matric_pot_nonlinear ./ rho_w ./ g;
                ground.STATVAR.waterPotential(~linear_range) = (L_sl .* T_nonlinear ./ T0 - mwp0) ./ rho_w ./ g;
                ground.STATVAR.water(~linear_range) = satWater_nonlinear .* porosity;
                ground.STATVAR.T(~linear_range) = T_nonlinear;
                
            end

            ground.STATVAR.water(ground.STATVAR.water > waterIce) = waterIce(ground.STATVAR.water > waterIce); %elminates small rounding errors and avoids negative ice contents
            
            ground.STATVAR.water = ground.STATVAR.water .* ground.STATVAR.layerThick .* ground.STATVAR.area;

            ground.STATVAR.ice = max(0,ground.STATVAR.waterIce - ground.STATVAR.water);

            
%             for i=1:size(T_prime,1)
%                 LUT_E_prime = ground.LUT.lut_E_prime((soil_type(i,1)-1).*LUT_size_T +(1:LUT_size_T),:);
%                 LUT_T_prime = ground.LUT.lut_T_prime((soil_type(i,1)-1).*LUT_size_T +(1:LUT_size_T),:);
%                 
%                 left_gamma_index = find(gamma(i,1) - gamma_list(soil_type(i,1),:) >= 0, 1, 'last');
%                 scale_factor_gamma = (gamma(i,1) - gamma_list(soil_type(i,1), left_gamma_index) ) ./ (gamma_list(soil_type(i,1), left_gamma_index+1) - gamma_list(soil_type(i,1), left_gamma_index));
% 
%                 T_prime_left = interp1(LUT_E_prime(:,left_gamma_index), LUT_T_prime(:,left_gamma_index), E_prime(i,1), 'linear');
%                 T_prime_right = interp1(LUT_E_prime(:,left_gamma_index+1), LUT_T_prime(:,left_gamma_index+1), E_prime(i,1), 'linear');
%                 T_prime_final = T_prime_left + scale_factor_gamma .*( T_prime_right - T_prime_left);
%                 T_nonlinear(i,1) = (T_prime_final ./ alpha(i,1) - mwp0(i,1)) ./ X;
%                 satWater_nonlinear(i,1) = E_prime(i,1) - gamma(i,1) .* T_prime_final;
%                 matric_pot_nonlinear(i,1) = T_prime_final ./ alpha(i,1);
%             end
%   
%             
%             ground.STATVAR.waterPotential(~linear_range)  = -matric_pot_nonlinear ./ rho_w ./ g;
%             ground.STATVAR.water(~linear_range) = satWater_nonlinear .* porosity;
%             ground.STATVAR.T(~linear_range) = T_nonlinear;
%             
%             ground.STATVAR.water(ground.STATVAR.water > waterIce) = waterIce(ground.STATVAR.water > waterIce); %elminates small rounding errors and avoids negative ice contents
%             
%             ground.STATVAR.water = ground.STATVAR.water .* ground.STATVAR.layerThick .* ground.STATVAR.area;
%             ground.STATVAR.ice = ground.STATVAR.waterIce - ground.STATVAR.water;

        end
        
        function ground = get_T_water_freezeC_Xice(ground)
            
            sat_waterIce_min = 0.005;
            LUT_size_gamma = 2.^9;
            LUT_size_T = 2.^10;
            porosity_max = 0.95;
            porosity_min = 0.05;
            
            beta_interface = 2.2;  %make this a constant!
            
            L_sl = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            T0 = ground.CONST.Tmfw;
            
            energy = ground.STATVAR.energy ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            waterIce = ground.STATVAR.waterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            XwaterIce = ground.STATVAR.XwaterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            soil_type = ground.STATVAR.soil_type;

            %new Sebastian
            ground.STATVAR.waterPotential  = energy .*0;
            ground.STATVAR.water = energy .*0;
            ground.STATVAR.Xwater = energy .*0;
            ground.STATVAR.T = energy .*0;
            %end new Sebastian
            
            n = ground.STATVAR.n;
            alpha = ground.STATVAR.alpha ./ rho_w ./ g;
            
            m=1-1./n;
            porosity = 1 - mineral - organic - XwaterIce;
            sat_waterIce = waterIce ./ porosity;
          %  porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic) ./(ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
          %  sat_waterIce = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce) ./ porosity;
            
            %matric potential (Pa)
            mwp0 = real(1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n));
            
            linear_range = energy >= 0;
            
            %1. linear range
            ground.STATVAR.waterPotential(linear_range)  = -mwp0(linear_range) ./ rho_w ./ g;
            ground.STATVAR.water(linear_range) = waterIce(linear_range);
            ground.STATVAR.Xwater(linear_range) = XwaterIce(linear_range);
            
            ground.STATVAR.T(linear_range) = energy(linear_range) ./ ((waterIce(linear_range) + XwaterIce(linear_range)).* c_w + mineral(linear_range) .* c_m + organic(linear_range) .* c_o);
            
            %done!
            
            %2. Xice melt regime, T=0, calculate Xice and Xwater, "normal soil" unfrozen
            Xice_melt = (energy <0 & energy >= - ground.CONST.L_f .* XwaterIce);
            ground.STATVAR.waterPotential(Xice_melt)  = -mwp0(Xice_melt) ./ rho_w ./ g;
            ground.STATVAR.water(Xice_melt) = waterIce(Xice_melt);
            ground.STATVAR.Xwater(Xice_melt) = XwaterIce(Xice_melt) .*(1 + energy(Xice_melt)./ (ground.CONST.L_f .* XwaterIce(Xice_melt)));
            ground.STATVAR.T(Xice_melt) = 0;

            %3. freeze curve regime
            nonlinear_range = energy < - ground.CONST.L_f .* XwaterIce;
            
            energy = energy + ground.CONST.L_f .* XwaterIce; %subtract the latent part of the Xice energy (-Lf * XwaterIce) from energy
            %and treat the Xice part in exactly the same fashion as mineral and organic
            
            if sum(double(nonlinear_range)) > 0
                %non-linear range
                alpha = alpha(nonlinear_range);
                n= n(nonlinear_range);
                m=m(nonlinear_range);
                mwp0 = mwp0(nonlinear_range);
                porosity = porosity(nonlinear_range);
                sat_waterIce = sat_waterIce(nonlinear_range);
                
                C0 = mineral(nonlinear_range) .* c_m + organic(nonlinear_range) .* c_o + porosity .* sat_waterIce .* c_i + XwaterIce(nonlinear_range) .* c_i; %vector
                %C0 = mineral(nonlinear_range) .* c_m + organic(nonlinear_range) .* c_o + waterIce(nonlinear_range) .* c_i + XwaterIce(nonlinear_range) .* c_i; %vector
                X = -L_sl./ T0 .* beta_interface;
                L0 = porosity.* L_sl; %vector
                gamma = C0 ./ alpha ./ X ./ L0 ; %vector
                
                soil_type = ground.STATVAR.soil_type(nonlinear_range); %-1;
                
                E_prime = energy(nonlinear_range) ./ L0 + sat_waterIce + C0 .* mwp0 ./ X./ L0; %takes soil type into account, alpha is assigned correctly
                
                T_nonlinear = E_prime.*0;
                satWater_nonlinear = E_prime.*0;
                matric_pot_nonlinear = E_prime.*0;
                
                %find correct gamma
                depth_gamma_LUT = log(LUT_size_gamma)./log(2);
                gamma_left_index = (soil_type-1) .* (LUT_size_gamma+1) + 1 ; %start with first index that belongs to correct soil code
                for i = depth_gamma_LUT-1:-1:0
                    gamma_left_index = gamma_left_index + 2.^i .* double(gamma >= ground.LUT.gamma(gamma_left_index + 2.^i,1));
                end
                scale_factor_gamma = (gamma - ground.LUT.gamma(gamma_left_index, 1)) ./ (ground.LUT.gamma(gamma_left_index+1, 1) - ground.LUT.gamma(gamma_left_index, 1));
                
                %left side fit E_prime
                depth_T_LUT = log(LUT_size_T)./log(2);
                %T_left_index = (soil_type-1) .* (LUT_size_gamma+1) .* (LUT_size_T+1) + (gamma_left_index-1) .* (LUT_size_T+1) + 1; %soil_type offset + gamma offset
                T_left_index = (gamma_left_index-1) .* (LUT_size_T+1) + 1; %soil_type offset already contained in gamma offset
                for i = depth_T_LUT-1:-1:0
                    T_left_index = T_left_index + 2.^i .* double(E_prime >= ground.LUT.lut_E_prime(T_left_index + 2.^i,1));
                end
                scale_factor_T = (E_prime - ground.LUT.lut_E_prime(T_left_index, 1)) ./ (ground.LUT.lut_E_prime(T_left_index+1, 1) - ground.LUT.lut_E_prime(T_left_index, 1));
                T_prime_left = ground.LUT.lut_T_prime(T_left_index, 1) + scale_factor_T .* (ground.LUT.lut_T_prime(T_left_index+1, 1) - ground.LUT.lut_T_prime(T_left_index, 1));
                
                %right side fit E_prime
                %T_right_index = (soil_type-1) .* (LUT_size_gamma+1) .* (LUT_size_T+1) + (gamma_left_index+1-1) .* (LUT_size_T+1) + 1; %soil_type offset + gamma offset for gamma right side
                T_right_index = (gamma_left_index + 1 - 1) .* (LUT_size_T+1) + 1; %soil_type offset already contained in gamma offset
                for i = depth_T_LUT-1:-1:0
                    T_right_index = T_right_index + 2.^i .* double(E_prime >= ground.LUT.lut_E_prime(T_right_index + 2.^i,1));
                end
                scale_factor_T = (E_prime - ground.LUT.lut_E_prime(T_right_index, 1)) ./ (ground.LUT.lut_E_prime(T_right_index+1, 1) - ground.LUT.lut_E_prime(T_right_index, 1));
                T_prime_right = ground.LUT.lut_T_prime(T_right_index, 1) + scale_factor_T .* (ground.LUT.lut_T_prime(T_right_index+1, 1) - ground.LUT.lut_T_prime(T_right_index, 1));
                
                
                T_prime = T_prime_left + scale_factor_gamma .* (T_prime_right - T_prime_left);
                
                %special fit for low water/high matric water potential grid cells, this is accomplished
                %by linearizing the dimensionless equation in T_prime around the
                %matric water potential intercept and then solving for
                %T_prime
                low_water_range = real(log(mwp0 .* alpha)) > -0.364.* real(log(-gamma))+6.65;
                %[real(log(mwp0 .* alpha)) -0.34.* real(log(-gamma))+7.0 -0.364.* real(log(-gamma))+6.65]
                %                 low_water_range = gamma<0
                T_prime0 = alpha(low_water_range) .*mwp0(low_water_range);
                T_prime(low_water_range) = (E_prime(low_water_range) - (1 + T_prime0.^n(low_water_range)).^(-m(low_water_range)) - T_prime0 .* (n(low_water_range).*m(low_water_range)).* T_prime0.^(n(low_water_range)-1) ...
                    .* (1+T_prime0.^n(low_water_range)).^(-m(low_water_range)-1)) ./ (gamma(low_water_range) -  n(low_water_range).*m(low_water_range).* T_prime0.^(n(low_water_range)-1) .* (1+T_prime0.^n(low_water_range)).^(-m(low_water_range)-1));
                
                
                T_nonlinear = (T_prime ./ alpha - mwp0) ./ X;
                satWater_nonlinear = E_prime - gamma .* T_prime;
                matric_pot_nonlinear = T_prime ./ alpha;
                
                
                %ground.STATVAR.waterPotential(nonlinear_range)  = -matric_pot_nonlinear ./ rho_w ./ g;
                ground.STATVAR.waterPotential(nonlinear_range) = (L_sl .* T_nonlinear ./ T0 - mwp0) ./ rho_w ./ g;
                ground.STATVAR.water(nonlinear_range) = satWater_nonlinear .* porosity;
                ground.STATVAR.Xwater(nonlinear_range) = 0;
                ground.STATVAR.T(nonlinear_range) = T_nonlinear;
                
            end
            
            ground.STATVAR.water(ground.STATVAR.water > waterIce) = waterIce(ground.STATVAR.water > waterIce); %elminates small rounding errors and avoids negative ice contents
            
            ground.STATVAR.water = ground.STATVAR.water .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.Xwater = ground.STATVAR.Xwater .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.Xwater(ground.STATVAR.T<0) = 0;
            ground.STATVAR.ice = max(0, ground.STATVAR.waterIce - ground.STATVAR.water);
            ground.STATVAR.ice(ground.STATVAR.T>0) = 0;
            ground.STATVAR.Xice = max(0, ground.STATVAR.XwaterIce - ground.STATVAR.Xwater);
            ground.STATVAR.Xice(ground.STATVAR.T>0) = 0;
        end
         
        function ground = get_T_water_freezeC_Xice_pressure(ground)
            
            sat_waterIce_min = 0.005;
            LUT_size_gamma = 2.^9;
            LUT_size_T = 2.^10;
            %porosity_max = 0.95;
            %porosity_min = 0.05;
            
            beta_interface = 2.2;  %make this a constant!
            
            L_sl = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            T0 = ground.CONST.Tmfw;
            
            energy = ground.STATVAR.energy ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            waterIce = ground.STATVAR.waterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            XwaterIce = ground.STATVAR.XwaterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            soil_type = ground.STATVAR.soil_type;
            
            n = ground.STATVAR.n;
            alpha = ground.STATVAR.alpha ./ rho_w ./ g;
            
            m=1-1./n;
            porosity = ground.STATVAR.porosity;
            sat_waterIce = (waterIce + XwaterIce) ./ porosity;
          %  porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic) ./(ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
          %  sat_waterIce = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce) ./ porosity;
            
            %matric potential (Pa)
            mwp0 = real(1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n));
            
            linear_range = energy >= 0;
            
            %1. linear range
            ground.STATVAR.waterPotential(linear_range)  = -mwp0(linear_range) ./ rho_w ./ g;
            ground.STATVAR.water(linear_range) = waterIce(linear_range);
            ground.STATVAR.Xwater(linear_range) = XwaterIce(linear_range);
            
            ground.STATVAR.T(linear_range) = energy(linear_range) ./ ((waterIce(linear_range) + XwaterIce(linear_range)).* c_w + mineral(linear_range) .* c_m + organic(linear_range) .* c_o);
            
            %done!
            
            %2. Xice melt regime, T=0, calculate Xice and Xwater, "normal soil" unfrozen
            Xice_melt = (energy <0 & energy >= - ground.CONST.L_f .* XwaterIce);
            ground.STATVAR.waterPotential(Xice_melt)  = -mwp0(Xice_melt) ./ rho_w ./ g;
            ground.STATVAR.water(Xice_melt) = waterIce(Xice_melt);
            ground.STATVAR.Xwater(Xice_melt) = XwaterIce(Xice_melt) .*(1 + energy(Xice_melt)./ (ground.CONST.L_f .* XwaterIce(Xice_melt)));
            ground.STATVAR.T(Xice_melt) = 0;

            %3. freeze curve regime
            nonlinear_range = energy < - ground.CONST.L_f .* XwaterIce;
            
            energy = energy + ground.CONST.L_f .* XwaterIce; %subtract the latent part of the Xice energy (-Lf * XwaterIce) from energy
            %and treat the Xice part in exactly the same fashion as mineral and organic
            
            if sum(double(nonlinear_range)) > 0
                %non-linear range
                alpha = alpha(nonlinear_range);
                n= n(nonlinear_range);
                m=m(nonlinear_range);
                mwp0 = mwp0(nonlinear_range);
                porosity = porosity(nonlinear_range);
                sat_waterIce = sat_waterIce(nonlinear_range);
                
                C0 = mineral(nonlinear_range) .* c_m + organic(nonlinear_range) .* c_o + porosity .* sat_waterIce .* c_i + XwaterIce(nonlinear_range) .* c_i; %vector
                %C0 = mineral(nonlinear_range) .* c_m + organic(nonlinear_range) .* c_o + waterIce(nonlinear_range) .* c_i + XwaterIce(nonlinear_range) .* c_i; %vector
                X = -L_sl./ T0 .* beta_interface;
                L0 = porosity.* L_sl; %vector
                gamma = C0 ./ alpha ./ X ./ L0 ; %vector
                
                soil_type = ground.STATVAR.soil_type(nonlinear_range); %-1;
                
                E_prime = energy(nonlinear_range) ./ L0 + sat_waterIce + C0 .* mwp0 ./ X./ L0; %takes soil type into account, alpha is assigned correctly
                
                T_nonlinear = E_prime.*0;
                satWater_nonlinear = E_prime.*0;
                matric_pot_nonlinear = E_prime.*0;
                
                %find correct gamma
                depth_gamma_LUT = log(LUT_size_gamma)./log(2);
                gamma_left_index = (soil_type-1) .* (LUT_size_gamma+1) + 1 ; %start with first index that belongs to correct soil code
                for i = depth_gamma_LUT-1:-1:0
                    gamma_left_index = gamma_left_index + 2.^i .* double(gamma >= ground.LUT.gamma(gamma_left_index + 2.^i,1));
                end
                scale_factor_gamma = (gamma - ground.LUT.gamma(gamma_left_index, 1)) ./ (ground.LUT.gamma(gamma_left_index+1, 1) - ground.LUT.gamma(gamma_left_index, 1));
                
                %left side fit E_prime
                depth_T_LUT = log(LUT_size_T)./log(2);
                %T_left_index = (soil_type-1) .* (LUT_size_gamma+1) .* (LUT_size_T+1) + (gamma_left_index-1) .* (LUT_size_T+1) + 1; %soil_type offset + gamma offset
                T_left_index = (gamma_left_index-1) .* (LUT_size_T+1) + 1; %soil_type offset already contained in gamma offset
                for i = depth_T_LUT-1:-1:0
                    T_left_index = T_left_index + 2.^i .* double(E_prime >= ground.LUT.lut_E_prime(T_left_index + 2.^i,1));
                end
                scale_factor_T = (E_prime - ground.LUT.lut_E_prime(T_left_index, 1)) ./ (ground.LUT.lut_E_prime(T_left_index+1, 1) - ground.LUT.lut_E_prime(T_left_index, 1));
                T_prime_left = ground.LUT.lut_T_prime(T_left_index, 1) + scale_factor_T .* (ground.LUT.lut_T_prime(T_left_index+1, 1) - ground.LUT.lut_T_prime(T_left_index, 1));
                
                %right side fit E_prime
                %T_right_index = (soil_type-1) .* (LUT_size_gamma+1) .* (LUT_size_T+1) + (gamma_left_index+1-1) .* (LUT_size_T+1) + 1; %soil_type offset + gamma offset for gamma right side
                T_right_index = (gamma_left_index + 1 - 1) .* (LUT_size_T+1) + 1; %soil_type offset already contained in gamma offset
                for i = depth_T_LUT-1:-1:0
                    T_right_index = T_right_index + 2.^i .* double(E_prime >= ground.LUT.lut_E_prime(T_right_index + 2.^i,1));
                end
                scale_factor_T = (E_prime - ground.LUT.lut_E_prime(T_right_index, 1)) ./ (ground.LUT.lut_E_prime(T_right_index+1, 1) - ground.LUT.lut_E_prime(T_right_index, 1));
                T_prime_right = ground.LUT.lut_T_prime(T_right_index, 1) + scale_factor_T .* (ground.LUT.lut_T_prime(T_right_index+1, 1) - ground.LUT.lut_T_prime(T_right_index, 1));
                
                
                T_prime = T_prime_left + scale_factor_gamma .* (T_prime_right - T_prime_left);
                
                %special fit for low water/high matric water potential grid cells, this is accomplished
                %by linearizing the dimensionless equation in T_prime around the
                %matric water potential intercept and then solving for
                %T_prime
                low_water_range = real(log(mwp0 .* alpha)) > -0.364.* real(log(-gamma))+6.65;
                %[real(log(mwp0 .* alpha)) -0.34.* real(log(-gamma))+7.0 -0.364.* real(log(-gamma))+6.65]
                %                 low_water_range = gamma<0
                T_prime0 = alpha(low_water_range) .*mwp0(low_water_range);
                T_prime(low_water_range) = (E_prime(low_water_range) - (1 + T_prime0.^n(low_water_range)).^(-m(low_water_range)) - T_prime0 .* (n(low_water_range).*m(low_water_range)).* T_prime0.^(n(low_water_range)-1) ...
                    .* (1+T_prime0.^n(low_water_range)).^(-m(low_water_range)-1)) ./ (gamma(low_water_range) -  n(low_water_range).*m(low_water_range).* T_prime0.^(n(low_water_range)-1) .* (1+T_prime0.^n(low_water_range)).^(-m(low_water_range)-1));
                
                
                T_nonlinear = (T_prime ./ alpha - mwp0) ./ X;
                satWater_nonlinear = E_prime - gamma .* T_prime;
                matric_pot_nonlinear = T_prime ./ alpha;
                
                
                %ground.STATVAR.waterPotential(nonlinear_range)  = -matric_pot_nonlinear ./ rho_w ./ g;
                ground.STATVAR.waterPotential(nonlinear_range) = (L_sl .* T_nonlinear ./ T0 - mwp0) ./ rho_w ./ g;
                ground.STATVAR.water(nonlinear_range) = satWater_nonlinear .* porosity;
                ground.STATVAR.Xwater(nonlinear_range) = 0;
                ground.STATVAR.T(nonlinear_range) = T_nonlinear;
                
            end
            
            ground.STATVAR.water(ground.STATVAR.water > waterIce) = waterIce(ground.STATVAR.water > waterIce); %elminates small rounding errors and avoids negative ice contents
            
            ground.STATVAR.water = ground.STATVAR.water .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.Xwater = ground.STATVAR.Xwater .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.Xwater(ground.STATVAR.T<0) = 0;
            ground.STATVAR.ice = max(0, ground.STATVAR.waterIce - ground.STATVAR.water);
            ground.STATVAR.ice(ground.STATVAR.T>0) = 0;
            ground.STATVAR.Xice = max(0, ground.STATVAR.XwaterIce - ground.STATVAR.Xwater);
            ground.STATVAR.Xice(ground.STATVAR.T>0) = 0;
         end
        
        %         function ground = get_T_water_freezeC_Xice(ground)
        %
        %             L_sl = ground.CONST.L_f ./ ground.CONST.rho_w;
        %             c_w = ground.CONST.c_w;
        %             c_i = ground.CONST.c_i;
        %             c_o = ground.CONST.c_o;
        %             c_m = ground.CONST.c_m;
        %             rho_w = ground.CONST.rho_w;
        %             g = ground.CONST.g;
        %             Tmfw = ground.CONST.Tmfw;
        %
        %             %waterPotZero and T_star must be known here
        %             energy = ground.STATVAR.energy;
        %             mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;  %volumetric fraction of the total volume incl Xice
        %             organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
        %             waterIce = ground.STATVAR.waterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
        %             XwaterIce = ground.STATVAR.XwaterIce./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
        %
        %             n = ground.STATVAR.n;
        %             alpha = ground.STATVAR.alpha;
        %             thetaRes = ground.STATVAR.thetaRes;
        %
        %             m=1-1./n;
        %             porosity = 1 - mineral - organic - XwaterIce;
        %             waterPotZero = real(-1./alpha .*(((waterIce - thetaRes)./(porosity - thetaRes)).^(-1./m)-1).^(1./n));
        %             Tstar =  g .* Tmfw ./ L_sl .* waterPotZero;
        %
        %             %distinguish four different regimes
        %             %1. linear positive unfrozen
        %             unfrozen_positive = (energy >= 0);
        %             ground.STATVAR.waterPotential(unfrozen_positive)  = waterPotZero(unfrozen_positive);
        %             ground.STATVAR.water(unfrozen_positive) = ground.STATVAR.waterIce(unfrozen_positive);
        %             ground.STATVAR.Xwater(unfrozen_positive) = ground.STATVAR.XwaterIce(unfrozen_positive);
        %             ground.STATVAR.T(unfrozen_positive) = energy(unfrozen_positive) ./ ((ground.STATVAR.waterIce(unfrozen_positive) + ground.STATVAR.XwaterIce(unfrozen_positive)) .* c_w + ...
        %                 ground.STATVAR.mineral(unfrozen_positive) .* c_m + ground.STATVAR.organic(unfrozen_positive) .* c_o);
        %
        %             %2. Xice melt regime, T=0, calculate Xice and Xwater, "normal soil" unfrozen
        %             Xice_melt = (energy <0 & energy >= - ground.CONST.L_f .* ground.STATVAR.XwaterIce);
        %             ground.STATVAR.waterPotential(Xice_melt)  = waterPotZero(Xice_melt);
        %             ground.STATVAR.water(Xice_melt) = ground.STATVAR.waterIce(Xice_melt);
        %             ground.STATVAR.Xwater(Xice_melt) = ground.STATVAR.XwaterIce(Xice_melt) .*(1 + energy(Xice_melt)./ (ground.CONST.L_f .* ground.STATVAR.XwaterIce(Xice_melt)));
        %             ground.STATVAR.T(Xice_melt) = 0;
        %
        %             energy = energy + ground.CONST.L_f .* ground.STATVAR.XwaterIce; %subtract the latent part of the Xice energy (-Lf * XwaterIce) from energy
        %             %and treat the Xice part in exactly the same fashion as mineral and organic
        %
        %             %3. linear negative,  XwaterIce fully frozen
        %             unfrozen_negative = (energy < 0 & energy >= Tstar .* ((ground.STATVAR.waterIce + ground.STATVAR.XwaterIce) .* c_i + ...
        %                 ground.STATVAR.mineral .* c_m + ground.STATVAR.organic .* c_o));
        %             ground.STATVAR.waterPotential(unfrozen_negative)  = waterPotZero(unfrozen_negative);
        %             ground.STATVAR.water(unfrozen_negative) = ground.STATVAR.waterIce(unfrozen_negative);
        %             ground.STATVAR.Xwater(unfrozen_negative) = 0;
        %             ground.STATVAR.T(unfrozen_negative) = energy(unfrozen_negative) ./ ((ground.STATVAR.waterIce(unfrozen_negative) + ground.STATVAR.XwaterIce(unfrozen_negative)) .* c_i + ...
        %                 ground.STATVAR.mineral(unfrozen_negative) .* c_m + ground.STATVAR.organic(unfrozen_negative) .* c_o);
        %
        %             %4. freeze curve region, XwaterIce fully frozen
        %             freeze_curve = energy < Tstar .* ((ground.STATVAR.waterIce + ground.STATVAR.XwaterIce) .* c_i + ground.STATVAR.mineral .* c_m + ground.STATVAR.organic .* c_o);
        %
        %             %transform to [J/m3], at this point XwaterIce is treated like mineral/organic, latent heat of XwaterIce is already subtracted
        %             energy = energy(freeze_curve) ./ ground.STATVAR.layerThick(freeze_curve) ./ ground.STATVAR.area(freeze_curve);
        %
        %             %interpolate in look-up tables
        %             Tstar = Tstar(freeze_curve) + Tmfw;
        %             a = g .* Tstar ./ L_sl;
        %             A = Tstar - Tmfw - a .*waterPotZero(freeze_curve);
        %             b = (waterIce(freeze_curve) + XwaterIce(freeze_curve)) .* c_i + mineral(freeze_curve) .* c_m + organic(freeze_curve) .* c_o;
        %             B = rho_w .* L_sl .* (waterIce(freeze_curve) - thetaRes(freeze_curve));
        %             C = rho_w .* L_sl .* (porosity(freeze_curve) - thetaRes(freeze_curve));
        %             c = -1 ./ alpha(freeze_curve);
        %             soil_type = ground.STATVAR.soil_type(freeze_curve) - 1;
        %             size_LUT = ground.PARA.LUT_size_waterIce .* ground.PARA.LUT_size_T;
        %
        %             energy_shifted_and_scaled = real((energy - b.*A + B) ./ (a.*b.*c));
        %             scale_factor = real(C ./ (a.*b.*c));
        %
        %             sf_matrix = ground.LUT.sf_matrix;
        %             LUT = ground.LUT.lut_energy_T_water;
        %
        %             pos_sf = (scale_factor - ground.LUT.min_sf(soil_type + 1,1)) ./ (ground.LUT.max_sf(soil_type + 1,1) - ground.LUT.min_sf(soil_type + 1,1)) .* ground.PARA.LUT_size_waterIce;
        %             pos_sf(pos_sf<1) = 1;
        %             pos_sf(pos_sf > ground.PARA.LUT_size_waterIce -1) = ground.PARA.LUT_size_waterIce - 1;
        %             fraction_sf = pos_sf - floor(pos_sf);
        %             pos_sf = floor(pos_sf);
        %
        %             min_sf_interp = sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,2) + ...
        %                 fraction_sf .* (sf_matrix(pos_sf + 1 + soil_type.*ground.PARA.LUT_size_waterIce,2) - sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,2));
        %             max_sf_interp = sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,3) + ...
        %                 fraction_sf .* (sf_matrix(pos_sf + 1 + soil_type.*ground.PARA.LUT_size_waterIce,3) - sf_matrix(pos_sf + soil_type.*ground.PARA.LUT_size_waterIce,3));
        %
        %             pos_energy = (energy_shifted_and_scaled - min_sf_interp) ./ (max_sf_interp - min_sf_interp) .* ground.PARA.LUT_size_T;
        %             pos_energy(pos_energy<1) = 1;
        %             pos_energy(pos_energy > ground.PARA.LUT_size_T -1) = ground.PARA.LUT_size_T-1;
        %             fraction_energy = pos_energy - floor(pos_energy);
        %             pos_energy = floor(pos_energy);
        %
        %             left = LUT(pos_sf + (pos_energy-1).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) + ...
        %                 fraction_sf .* (LUT(pos_sf+1 + (pos_energy-1).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) - ...
        %                 LUT(pos_sf + (pos_energy-1).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT));
        %             right = LUT(pos_sf + (pos_energy).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) + ...
        %                 fraction_sf .* (LUT(pos_sf+1 + (pos_energy).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT) - ...
        %                 LUT(pos_sf + (pos_energy).*ground.PARA.LUT_size_waterIce + soil_type.*size_LUT));
        %
        %             X_interp = left + fraction_energy .*(right - left);
        %
        %             waterPot_interp = X_interp./-alpha(freeze_curve);
        %             T_interp = (waterPot_interp-waterPotZero(freeze_curve)).* Tstar.*g ./L_sl + Tstar - Tmfw;
        %             water_interp  = real(thetaRes(freeze_curve) + (porosity(freeze_curve) - thetaRes(freeze_curve)) .* (1+(X_interp).^n(freeze_curve)).^-m(freeze_curve));
        %
        %             ground.STATVAR.waterPotential(freeze_curve)  = waterPot_interp;
        %             ground.STATVAR.water(freeze_curve) = water_interp .* ground.STATVAR.layerThick(freeze_curve) .* ground.STATVAR.area(freeze_curve);
        %             ground.STATVAR.water(ground.STATVAR.water > ground.STATVAR.waterIce) = ground.STATVAR.waterIce(ground.STATVAR.water > ground.STATVAR.waterIce); %elminates small rounding errors and avoids negative ice contents
        %
        %             ground.STATVAR.Xwater(freeze_curve) = 0;
        %             ground.STATVAR.T(freeze_curve) = T_interp;
        %
        %             ground.STATVAR.ice = ground.STATVAR.waterIce - ground.STATVAR.water;
        %             ground.STATVAR.Xice = ground.STATVAR.XwaterIce - ground.STATVAR.Xwater;
        %         end
        %
        %
        %
        
        % get energy from temeprature and water contents, normally part of initializiation
        
        function ground = get_E_freezeC(ground) %required for initialization
            
            L_sl = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            T0 = ground.CONST.Tmfw;
            
            T = ground.STATVAR.T;
            mineral= ground.STATVAR.mineral ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            organic = ground.STATVAR.organic ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            waterIce = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
%             mineral= ground.STATVAR.mineral;
%             organic = ground.STATVAR.organic;
%             waterIce = ground.STATVAR.waterIce;
            layerThick = ground.STATVAR.layerThick;
            area = ground.STATVAR.area;
            soil_type = ground.STATVAR.soil_type;
            
            n = double(soil_type == 1) .* ground.CONST.vanGen_n(1) + double(soil_type == 2) .* ground.CONST.vanGen_n(2) + double(soil_type == 3) .* ground.CONST.vanGen_n(3) + double(soil_type == 4) .* ground.CONST.vanGen_n(4) + double(soil_type == 5) .* ground.CONST.vanGen_n(5);
            ground.STATVAR.n = n;
            alpha = double(soil_type == 1) .* ground.CONST.vanGen_alpha(1) + double(soil_type == 2) .* ground.CONST.vanGen_alpha(2) + double(soil_type == 3) .* ground.CONST.vanGen_alpha(3) + double(soil_type == 4) .* ground.CONST.vanGen_alpha(4) + + double(soil_type == 5) .* ground.CONST.vanGen_alpha(5);
            ground.STATVAR.alpha = alpha;
%             thetaRes = double(soil_type == 1) .* ground.CONST.vanGen_residual_wc(1) + double(soil_type == 2) .* ground.CONST.vanGen_residual_wc(2) + double(soil_type == 3) .* ground.CONST.vanGen_residual_wc(3) + double(soil_type == 4) .* ground.CONST.vanGen_residual_wc(4) +  double(soil_type == 5) .* ground.CONST.vanGen_residual_wc(5);
%             ground.STATVAR.thetaRes = thetaRes;
            
            porosity = 1-mineral-organic;
            sat_waterIce = waterIce ./ porosity;
            alpha = alpha ./ g ./ rho_w; %convert to Pa^-1;
            beta_interface = 2.2;  %make this a constant!
            m=1-1./n;
            
            mwp0 = real(1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n)); %Van genuchten
            mwp = -L_sl .* T ./ T0 .* beta_interface .* double(T<0) + mwp0; %universal step to get mwp from mwp0
            sat_water = double(mwp > 0) .* (1+(alpha.*mwp).^n).^(-m) + double(mwp <= 0);
            sat_ice = sat_waterIce - sat_water;
            energy = T.* (mineral .* c_m + organic .* c_o + porosity .* sat_waterIce .* (c_w .* double(T >= 0)+ c_i.* double(T < 0)));
            energy = energy - double(T < 0) .* L_sl .* porosity .* (sat_waterIce - sat_water);
            
            ground.STATVAR.water = sat_water .* porosity .* layerThick .* area;  % [m3]
            ground.STATVAR.ice = sat_ice .* porosity .*  layerThick .* area; %[m3]
            
            ground.STATVAR.waterIce = waterIce .* layerThick .* area; % [m3]
            ground.STATVAR.mineral = mineral .* layerThick .* area; % [m3]
            ground.STATVAR.organic = organic .* layerThick .* area; % [m3]
            ground.STATVAR.energy = energy .* layerThick .* area;  % [J]
            
            ground.STATVAR.air = (1-mineral-organic-waterIce) .* layerThick .* area;  % [m3]
            
            %ground.STATVAR.waterPotential = -mwp ./ rho_w ./ g;
            ground.STATVAR.waterPotential = (L_sl .* T ./ T0 .* double(T<0) - mwp0) ./ rho_w ./ g;
        end
        
        function ground = get_E_freezeC_pressure(ground) %required for initialization
            
            L_sl = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            T0 = ground.CONST.Tmfw;
            
            T = ground.STATVAR.T;
            mineral= (ground.STATVAR.mineral .* ground.STATVAR.initial_layerThick) ./ ground.STATVAR.layerThick; %Percentage of mineral in compressed soil
            organic= (ground.STATVAR.organic .* ground.STATVAR.initial_layerThick) ./ ground.STATVAR.layerThick; %Percentage of organic in compressed soil
            waterIce = ground.STATVAR.waterIce;
            layerThick = ground.STATVAR.layerThick;
            area = ground.STATVAR.area;
            soil_type = ground.STATVAR.soil_type;
            
            n = double(soil_type == 1) .* ground.CONST.vanGen_n(1) + double(soil_type == 2) .* ground.CONST.vanGen_n(2) + double(soil_type == 3) .* ground.CONST.vanGen_n(3) + double(soil_type == 4) .* ground.CONST.vanGen_n(4) + double(soil_type == 5) .* ground.CONST.vanGen_n(5);
            ground.STATVAR.n = n;
            alpha = double(soil_type == 1) .* ground.CONST.vanGen_alpha(1) + double(soil_type == 2) .* ground.CONST.vanGen_alpha(2) + double(soil_type == 3) .* ground.CONST.vanGen_alpha(3) + double(soil_type == 4) .* ground.CONST.vanGen_alpha(4) + + double(soil_type == 5) .* ground.CONST.vanGen_alpha(5);
            ground.STATVAR.alpha = alpha;
%             thetaRes = double(soil_type == 1) .* ground.CONST.vanGen_residual_wc(1) + double(soil_type == 2) .* ground.CONST.vanGen_residual_wc(2) + double(soil_type == 3) .* ground.CONST.vanGen_residual_wc(3) + double(soil_type == 4) .* ground.CONST.vanGen_residual_wc(4) +  double(soil_type == 5) .* ground.CONST.vanGen_residual_wc(5);
%             ground.STATVAR.thetaRes = thetaRes;
            
            porosity = ground.STATVAR.porosity;
            sat_waterIce = waterIce ./ porosity;
            alpha = alpha ./ g ./ rho_w; %convert to Pa^-1;
            beta_interface = 2.2;  %make this a constant!
            m=1-1./n;
            
            mwp0 = real(1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n)); %Van genuchten
            mwp = -L_sl .* T ./ T0 .* beta_interface .* double(T<0) + mwp0; %universal step to get mwp from mwp0
            sat_water = double(mwp > 0) .* (1+(alpha.*mwp).^n).^(-m) + double(mwp <= 0);
            sat_ice = sat_waterIce - sat_water;
            energy = T.* (mineral .* c_m + organic .* c_o + porosity .* sat_waterIce .* (c_w .* double(T >= 0)+ c_i.* double(T < 0)));
            energy = energy - double(T < 0) .* L_sl .* porosity .* (sat_waterIce - sat_water);
            
            ground.STATVAR.water = sat_water .* porosity .* layerThick .* area;  % [m3]
            ground.STATVAR.ice = sat_ice .* porosity .*  layerThick .* area; %[m3]
            
            ground.STATVAR.waterIce = waterIce .* layerThick .* area; % [m3]
            ground.STATVAR.mineral = mineral .* layerThick .* area; % [m3] Initial layerThick because mineral content does not change with changing layerThick
            ground.STATVAR.organic = organic .* layerThick .* area; % [m3]Initial layerThick because organic content does not change with changing layerThick
            ground.STATVAR.energy = energy .* layerThick .* area;  % [J]
            
            ground.STATVAR.air = (1-mineral-organic-waterIce) .* layerThick .* area;  % [m3]
            
            %ground.STATVAR.waterPotential = -mwp ./ rho_w ./ g;
            ground.STATVAR.waterPotential = (L_sl .* T ./ T0 .* double(T<0) - mwp0) ./ rho_w ./ g;
        end
        
        
        function ground = get_E_freezeC_Xice(ground) %required for initialization
            
            L_sl = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            T0 = ground.CONST.Tmfw;
            
            T = ground.STATVAR.T;
            mineral= ground.STATVAR.mineral ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
            organic = ground.STATVAR.organic ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
            waterIce = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
%             mineral= ground.STATVAR.mineral;
%             organic = ground.STATVAR.organic;
%             waterIce = ground.STATVAR.waterIce;
            layerThick = ground.STATVAR.layerThick;
            area = ground.STATVAR.area;
            soil_type = ground.STATVAR.soil_type;
            
%            Xice = ground.STATVAR.Xice .* double(T <= 0); % Xice initially only possible when frozen, provided in multiples of the "matrix", 1 means 50 vol% Xice, 50 % normal soil
            Xice = ground.STATVAR.XwaterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);        

            n = double(soil_type == 1) .* ground.CONST.vanGen_n(1) + double(soil_type == 2) .* ground.CONST.vanGen_n(2) + double(soil_type == 3) .* ground.CONST.vanGen_n(3) + double(soil_type == 4) .* ground.CONST.vanGen_n(4) + double(soil_type == 5) .* ground.CONST.vanGen_n(5);
            ground.STATVAR.n = n;
            alpha = double(soil_type == 1) .* ground.CONST.vanGen_alpha(1) + double(soil_type == 2) .* ground.CONST.vanGen_alpha(2) + double(soil_type == 3) .* ground.CONST.vanGen_alpha(3) + double(soil_type == 4) .* ground.CONST.vanGen_alpha(4) + + double(soil_type == 5) .* ground.CONST.vanGen_alpha(5);
            ground.STATVAR.alpha = alpha;
            
            porosity = 1-mineral-organic;
            sat_waterIce = waterIce ./ porosity;
            alpha = alpha ./ g ./ rho_w; %convert to Pa^-1;
            beta_interface = 2.2;  %make this a constant!
            m=1-1./n;
            
            mwp0 = real(1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n)); 
            mwp = -L_sl .* T ./ T0 .* beta_interface .* double(T<0) + mwp0;
            sat_water = double(mwp > 0) .* (1+(alpha.*mwp).^n).^(-m) + double(mwp <= 0);
            sat_ice = sat_waterIce - sat_water;
            energy = T.* (mineral .* c_m + organic .* c_o + porosity .* sat_waterIce .* (c_w .* double(T > 0)+ c_i.* double(T <= 0))); %add XwaterIce here!!!
            energy = energy - double(T <= 0) .* L_sl .* porosity .* (sat_waterIce - sat_water);
            
            ground.STATVAR.water = sat_water .* porosity .* layerThick ./ (1 + Xice) .* area;  % [m3]
            ground.STATVAR.ice = sat_ice .* porosity .*  layerThick ./ (1 + Xice) .* area; %[m3]
            
            ground.STATVAR.waterIce = waterIce .* layerThick ./ (1 + Xice) .* area; % [m3]
            ground.STATVAR.mineral = mineral .* layerThick ./ (1 + Xice) .* area; % [m3]
            ground.STATVAR.organic = organic .* layerThick ./ (1 + Xice) .* area; % [m3]
            %ground.STATVAR.energy = energy .* layerThick ./ (1 + Xice) .* area;  % [J]
            ground.STATVAR.energy = (energy + Xice .* (ground.STATVAR.T .* ground.CONST.c_i - ground.CONST.L_f)) ./ (1 + Xice) .* layerThick .* area;  % [J]
            ground.STATVAR.XwaterIce = Xice ./ (1 + Xice) .* layerThick .* area;
            ground.STATVAR.Xice = ground.STATVAR.XwaterIce .*double(T<=0);
            ground.STATVAR.Xwater = ground.STATVAR.XwaterIce .* double(T>0);
            
            ground.STATVAR.air = (1-mineral-organic-waterIce) .* layerThick ./ (1 + Xice).* area;  % [m3]
            
            %ground.STATVAR.waterPotential = -mwp ./ rho_w ./ g;
            ground.STATVAR.waterPotential = (L_sl .* T ./ T0 .* double(T<0) - mwp0) ./ rho_w ./ g;
            
            
        end
        
        
        %         function ground = get_E_freezeC_Xice(ground) %required for initialization
        %
        %             L_sl = ground.CONST.L_f ./ ground.CONST.rho_w;
        %             c_w = ground.CONST.c_w;
        %             c_i = ground.CONST.c_i;
        %             c_o = ground.CONST.c_o;
        %             c_m = ground.CONST.c_m;
        %             rho_w = ground.CONST.rho_w;
        %             g = ground.CONST.g;
        %             Tmfw = ground.CONST.Tmfw;
        %
        %             T = ground.STATVAR.T;
        %             mineral= ground.STATVAR.mineral; %properties of the matrix when Xice is removed!!!
        %             organic = ground.STATVAR.organic;
        %             waterIce = ground.STATVAR.waterIce;
        %             layerThick = ground.STATVAR.layerThick;
        %             area = ground.STATVAR.area;
        %             soil_type = ground.STATVAR.soil_type;
        %
        %             Xice = ground.STATVAR.Xice .* double(T <= 0); % Xice initially only possible when frozen, provided in multiples of the "matrix", 1 means 50 vol% Xice, 50 % normal soil
        %
        %             n = double(soil_type == 1) .* ground.CONST.vanGen_n(1) + double(soil_type == 2) .* ground.CONST.vanGen_n(2) + double(soil_type == 3) .* ground.CONST.vanGen_n(3) + double(soil_type == 4) .* ground.CONST.vanGen_n(4) + double(soil_type == 5) .* ground.CONST.vanGen_n(5);
        %             ground.STATVAR.n = n;
        %             alpha = double(soil_type == 1) .* ground.CONST.vanGen_alpha(1) + double(soil_type == 2) .* ground.CONST.vanGen_alpha(2) + double(soil_type == 3) .* ground.CONST.vanGen_alpha(3) + double(soil_type == 4) .* ground.CONST.vanGen_alpha(4) + + double(soil_type == 5) .* ground.CONST.vanGen_alpha(5);
        %             ground.STATVAR.alpha = alpha;
        %             thetaRes = double(soil_type == 1) .* ground.CONST.vanGen_residual_wc(1) + double(soil_type == 2) .* ground.CONST.vanGen_residual_wc(2) + double(soil_type == 3) .* ground.CONST.vanGen_residual_wc(3) + double(soil_type == 4) .* ground.CONST.vanGen_residual_wc(4) +  double(soil_type == 5) .* ground.CONST.vanGen_residual_wc(5);
        %             ground.STATVAR.thetaRes = thetaRes;
        %
        %             porosity = 1-mineral-organic;
        %             m=1-1./n;
        %
        %             waterPotZero = real(-1./alpha .*(((waterIce - thetaRes)./(porosity - thetaRes)).^(-1./m)-1).^(1./n));
        %             Tstar = Tmfw + g .* Tmfw ./ L_sl .* waterPotZero;
        %
        %             waterPot = waterPotZero + (L_sl./g./Tstar .* (T - Tstar + Tmfw)).*double(T < Tstar - Tmfw);
        %             waterC  = double(T < Tstar - Tmfw) .* real(thetaRes + (porosity - thetaRes) .* (1+(-alpha.*waterPot).^n).^-m) + double(T >= Tstar-Tmfw) .* waterIce;
        %
        %             energy = double (T>=0) .* (waterIce .* c_w + mineral .* c_m + organic .* c_o) .*T + double(T<0) .* ((waterIce .* c_i + mineral .* c_m + organic .* c_o) .*T - rho_w .*L_sl .* (waterIce - waterC));
        %
        %             ground.STATVAR.water = waterC .* layerThick ./ (1 + Xice) .* area;  % [m3]
        %             ground.STATVAR.ice = (waterIce - waterC) .*  layerThick ./ (1 + Xice) .* area; %[m3]
        %
        %             ground.STATVAR.waterIce = waterIce .* layerThick ./ (1 + Xice) .* area; % [m3]
        %             ground.STATVAR.mineral = mineral .* layerThick ./ (1 + Xice) .* area; % [m3]
        %             ground.STATVAR.organic = organic .* layerThick ./ (1 + Xice) .* area; % [m3]
        %             ground.STATVAR.energy = (energy + Xice .* (ground.STATVAR.T .* ground.CONST.c_i - ground.CONST.L_f)) ./ (1 + Xice) .* layerThick .* area;  % [J]
        %             ground.STATVAR.XwaterIce = Xice ./ (1 + Xice) .* layerThick .* area;
        %             ground.STATVAR.Xice = ground.STATVAR.XwaterIce .*double(T<=0);
        %             ground.STATVAR.Xwater = ground.STATVAR.XwaterIce .* double(T>0);
        %
        %             ground.STATVAR.air = (1-mineral-organic-waterIce) .* layerThick ./ (1 + Xice) .* area;  % [m3]
        %
        %             ground.STATVAR.waterPotential = waterPot;
        %         end
        %
        
        
        %---look-up tables initialization-----------
        function ground = create_LUT_freezeC(ground) % creates lookup table LUT
            
            disp('creating look-up tables')
            
            T_min = -60;
            sat_waterIce_min = 0.005;
            LUT_size_gamma = 2.^9;
            LUT_size_T = 2.^10;
            porosity_max = 0.95;
            porosity_min = 0.05;
            
            beta_interface = 2.2;

            L_sl = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            rho_w = ground.CONST.rho_w;
            g = ground.CONST.g;
            T0 = ground.CONST.Tmfw;
            
            ground.LUT.gamma = [];
            ground.LUT.lut_E_prime = [];
            ground.LUT.lut_T_prime = [];
            
            for soil_type = 1:size(ground.CONST.vanGen_alpha,2)
                
                alpha = ground.CONST.vanGen_alpha(1, soil_type) ./ rho_w ./ g; %convert to Pa^-1
                n = ground.CONST.vanGen_n(1, soil_type);
                m=1-1./n;
                
                X = -L_sl./ T0 .* beta_interface;
                
                gamma_max =  ((1-porosity_max).* c_o ) ./ alpha ./ X ./ (porosity_max.* L_sl);
                gamma_min =  ((1-porosity_min).* c_m  + c_i) ./ alpha ./ X ./ (porosity_min.* L_sl);
                
                gamma_list = linspace(log(-gamma_min), log(-gamma_max), LUT_size_gamma+1);
                gamma_list = -exp(gamma_list);
                
                %dimensionless quantities computed from the dimensionless equation
                
                LUT_T_prime = [];
                LUT_E_prime = [];
                
                for gamma_index=1:size(gamma_list,2)
                    gamma= gamma_list(1,gamma_index);
                    
                    T_prime_max = real(alpha .* X .* T_min + ((sat_waterIce_min.^(-1./m)-1)).^(1./n));
                    E_prime_min = real(gamma .* T_prime_max + (1 + T_prime_max.^n).^(-m));
                    
                    T_prime_max = T_prime_max .*1.1;
                    E_prime_min = E_prime_min .*1.1;
                    
                    T_prime2 = linspace(T_prime_max, 0, LUT_size_T+1)';
                    E_prime2 = real(gamma .* T_prime2 + (1 + T_prime2.^n).^(-m));
                    
                    for count = 1:15
                        %             plot(E_prime2)
                        %             hold on
                        %add points
                        add_points_T =[(T_prime2(1:end-1,1) + T_prime2(2:end,1))./2  (3.*T_prime2(1:end-1,1) + T_prime2(2:end,1))./4 (T_prime2(1:end-1,1) + 3.*T_prime2(2:end,1))./4];%midpoints in T
                        add_points_E = real(gamma .* add_points_T + (1 + add_points_T.^n).^(-m));
                        
                        delta_E_prime2 = repmat(E_prime2(2:end,1) - E_prime2(1:end-1,1), 1, 3) ;
                        delta_T_prime2 = repmat(T_prime2(2:end,1) - T_prime2(1:end-1,1), 1, 3);
                        add_points_T_interp = repmat(T_prime2(1:end-1),1,3) + (add_points_E - repmat(E_prime2(1:end-1),1,3)) ./ delta_E_prime2 .* delta_T_prime2;
                        
                        add_points_T = add_points_T(:);
                        add_points_E = add_points_E(:);
                        mismatch = abs(add_points_T_interp(:) - add_points_T);
                        [~,sort_index]=sort(mismatch);
                        sort_index=sort_index(end-size(sort_index,1)./3+1:end,1);
                        T_prime2 = [T_prime2; add_points_T(sort_index,1)];
                        E_prime2 = [E_prime2; add_points_E(sort_index,1)];
                        [T_prime2, sort_index] = sort(T_prime2,'descend');
                        E_prime2=E_prime2(sort_index,1);
                        
                        %reduce points
                        reduce_points_T = [T_prime2(2:4:end-3,1) T_prime2(3:4:end-2,1) T_prime2(4:4:end-1,1)];
                        reduce_points_E = [E_prime2(2:4:end-3,1) E_prime2(3:4:end-2,1) E_prime2(4:4:end-1,1)];
                        T_prime2 = T_prime2(1:4:end,1);
                        E_prime2 = E_prime2(1:4:end,1);
                        
                        delta_E_prime2 = repmat(E_prime2(2:end,1) - E_prime2(1:end-1,1), 1, 3) ;
                        delta_T_prime2 = repmat(T_prime2(2:end,1) - T_prime2(1:end-1,1), 1, 3);
                        reduce_points_T_interp = repmat(T_prime2(1:end-1),1,3) + (reduce_points_E - repmat(E_prime2(1:end-1),1,3)) ./ delta_E_prime2 .* delta_T_prime2;
                        
                        reduce_points_T = reduce_points_T(:);
                        reduce_points_E = reduce_points_E(:);
                        mismatch = abs(reduce_points_T_interp(:) - reduce_points_T);
                        [~,sort_index]=sort(mismatch);
                        sort_index=sort_index(end-size(sort_index,1)./3+1:end,1);
                        T_prime2 = [T_prime2; reduce_points_T(sort_index,1)];
                        E_prime2 = [E_prime2; reduce_points_E(sort_index,1)];
                        [T_prime2, sort_index] = sort(T_prime2,'descend');
                        E_prime2 = E_prime2(sort_index,1);
                        
                    end
                    
                    %LUT_T_prime = [LUT_T_prime T_prime2];
                    %LUT_E_prime = [LUT_E_prime E_prime2];
                    
                    LUT_T_prime = [LUT_T_prime; T_prime2];
                    LUT_E_prime = [LUT_E_prime; E_prime2];
                    
                end
                
            ground.LUT.gamma = [ground.LUT.gamma; gamma_list'];
            ground.LUT.lut_E_prime = [ground.LUT.lut_E_prime; LUT_E_prime];
            ground.LUT.lut_T_prime = [ground.LUT.lut_T_prime; LUT_T_prime];
                
            end
            
        end
        
    end
end

