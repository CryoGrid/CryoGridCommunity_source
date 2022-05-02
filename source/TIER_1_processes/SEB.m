%========================================================================
% CryoGrid TIER1 library class, functions related to surface energy balance 
% NOTE: also contains functions for penetration of short-wave radiation
% S. Westermann, October 2020
%========================================================================


classdef SEB < BASE
    
    methods
        
        %Obukhov length
        function seb = L_star(seb, forcing)

            uz = forcing.TEMP.wind;
            z =  seb.PARA.airT_height;
            z0 = seb.PARA.z0;
            Tz = forcing.TEMP.Tair+273.15;
            Lstar = seb.STATVAR.Lstar;
            p = forcing.TEMP.p;
            Qh = seb.STATVAR.Qh;
            Qe = seb.STATVAR.Qe;
            
            rho = rho_air(seb, p, Tz);
            cp = seb.CONST.cp;
            kappa = seb.CONST.kappa; %0.4;
            g = seb.CONST.g; %9.81;

            if Tz >=273.15
                L = latent_heat_evaporation(seb, Tz); %1e3.*(2500.8 - 2.36.*(Tz-273.15));  %latent heat of evaporation of water
            else
                L = latent_heat_sublimation(seb, Tz); %1e3.*2834.1; %latent heat of sublimation
            end

            u_star = real(uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)));
            L_star = real(-rho.*cp.*Tz./kappa./g.*u_star.^3./(Qh + 0.61.*cp./L.*Tz.*Qe));
            L_star=(abs(L_star)<1e-7).*L_star./abs(L_star).*1e-7 + (abs(L_star)>=1e-7).*L_star;  %limits Lstar
            
            seb.STATVAR.Lstar = L_star;
            seb.STATVAR.u_star = u_star;
            
        end
        
        %atmospheric stability functions
        function res = psi_H(seb, zeta1, zeta2) % atmospheric stability function heat/water 
            
            if zeta1<=0
                res= 1.9.*atanh((1 - 11.6.*zeta1).^0.5) + log(zeta1) - (1.9.*atanh((1 - 11.6.*zeta2).^0.5) + log(zeta2));
            else
                res=((-5 + 5^0.5).*log(-3 + 5^0.5- 2.*zeta1) - (5 + 5^0.5).*log(3 + 5^0.5 + 2.*zeta1))/2  - (((-5 + 5^0.5).*log(-3 + 5^0.5- 2.*zeta2) - (5 + 5^0.5).*log(3 + 5^0.5 + 2.*zeta2))/2);
            end
        end
        
        function res = psi_M(seb, zeta1, zeta2) %atmospheric stability function momentum 
            
            if zeta1<=0
                res=-2.*atan((1 - 19.*zeta1).^(1/4)) + 2.* log(1 + (1 - 19.*zeta1).^(1/4)) + log(1 + (1 - 19.*zeta1).^0.5) - ...
                    (-2.*atan((1 - 19.*zeta2).^(1/4)) + 2.* log(1 + (1 - 19.*zeta2).^(1/4)) + log(1 + (1 - 19.*zeta2).^0.5));
            else
                res=-19.5.*(1 + zeta1).^(1/3) - 7.5367.*atan(0.57735 - 1.72489.*(1 + zeta1).^(1/3)) + 4.35131.*log(3+4.4814.*(1+zeta1).^(1/3)) - 2.17566.*log(3 - 4.4814.*(1 + zeta1).^(1/3) + 6.69433.*(1 + zeta1).^(2/3)) - ...
                    (-19.5.*(1 + zeta2).^(1/3) - 7.5367.*atan(0.57735 - 1.72489.*(1 + zeta2).^(1/3)) + 4.35131.*log(3+4.4814.*(1+zeta2).^(1/3)) - 2.17566.*log(3 - 4.4814.*(1 + zeta2).^(1/3) + 6.69433.*(1 + zeta2).^(2/3))) ;
            end
        end
        
        %latent heat flux that can be adjusted using the empirical factor "resistance to evapotranspiration" rs
        function Q_e = Q_eq(seb, forcing) 
            
            uz = forcing.TEMP.wind;
            p = forcing.TEMP.p;
            q = forcing.TEMP.q;
            Tz = forcing.TEMP.Tair;

            z =  seb.PARA.airT_height;
            z0 = seb.PARA.z0;
            rs = seb.PARA.rs;
            kappa = seb.CONST.kappa;
            
            TForcing = seb.STATVAR.T(1);
            Lstar = seb.STATVAR.Lstar;

            Tz=Tz+273.15;
            TForcing=TForcing+273.15;

            rho = rho_air(seb, p, Tz);
            L_w = latent_heat_evaporation(seb, TForcing); % 1e3.*(2500.8 - 2.36.*(TForcing-273.15));  %latent heat of evaporation of water
            L_i = latent_heat_sublimation(seb, TForcing); % 1e3.*2834.1; %latent heat of sublimation
            
            if TForcing<=273.15
                Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q-satPresIce(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            else
                Q_e = -rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb,z./Lstar, z0./Lstar)).*(q-satPresWater(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar)  ...
                    + rs.*uz.*kappa.^2./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)));
            end
        end
        
        %latent heat flux for potential evapotranspiration
        function Q_e = Q_eq_potET(seb, forcing) 
            
            uz = forcing.TEMP.wind;
            p = forcing.TEMP.p;
            q = forcing.TEMP.q;
            Tz = forcing.TEMP.Tair;
            
            z =  seb.PARA.airT_height;
            z0 = seb.PARA.z0;
            kappa = seb.CONST.kappa;
                        
            TForcing = seb.STATVAR.T(1);
            Lstar = seb.STATVAR.Lstar;
            
            Tz=Tz+273.15;
            TForcing=TForcing+273.15;
            rho = rho_air(seb, p, Tz);
            L_w = latent_heat_evaporation(seb, TForcing); % 1e3.*(2500.8 - 2.36.*(TForcing-273.15));  %latent heat of evaporation of water
            L_i = latent_heat_sublimation(seb, TForcing); % 1e3.*2834.1; %latent heat of sublimation

            if TForcing<=273.15
                Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q-satPresIce(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            else
                Q_e = -rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q-satPresWater(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            end
        end
        
        %latent heat flux for potential evapotranspiration
        function Q_e = Q_eq_potET_snow(seb, forcing) 
            
            uz = forcing.TEMP.wind;
            p = forcing.TEMP.p;
            q = forcing.TEMP.q;
            Tz = forcing.TEMP.Tair;
            
            z =  seb.PARA.airT_height;
            z0 = seb.PARA.z0;
            kappa = seb.CONST.kappa;
            
            water_fraction = seb.STATVAR.water(1,1) ./ seb.STATVAR.waterIce(1,1);
            ice_fraction = seb.STATVAR.ice(1,1) ./ seb.STATVAR.waterIce(1,1);
                        
            TForcing = seb.STATVAR.T(1);
            Lstar = seb.STATVAR.Lstar;
            
            Tz=Tz+273.15;
            TForcing=TForcing+273.15;
            rho = rho_air(seb, p, Tz);
            L_w = latent_heat_evaporation(seb, TForcing); % 1e3.*(2500.8 - 2.36.*(TForcing-273.15));  %latent heat of evaporation of water
            L_i = latent_heat_sublimation(seb, TForcing); % 1e3.*2834.1; %latent heat of sublimation

            if TForcing<273.15 || water_fraction <= 0
                Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q-satPresIce(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
                seb.STATVAR.sublimation = -Q_e ./(seb.CONST.rho_w .* seb.CONST.L_s) .* seb.STATVAR.area(1);
                seb.TEMP.sublimation_energy = seb.STATVAR.sublimation .* (seb.STATVAR.T(1) .* seb.CONST.c_i - seb.CONST.L_f);
                seb.STATVAR.evap = 0;
                seb.STATVAR.evap_energy = 0;
            else
                evap_fraction = max(0, min(water_fraction./0.1,1));
                sublim_fraction = 1 - evap_fraction;
                Qe_sublim = -sublim_fraction .* rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q-satPresIce(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
                Qe_evap = -evap_fraction .* rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q-satPresWater(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
                
                seb.STATVAR.sublimation = - Qe_sublim ./(seb.CONST.rho_w .* L_i) .* seb.STATVAR.area(1);
                seb.TEMP.sublimation_energy = seb.STATVAR.sublimation .* (seb.STATVAR.T(1) .* seb.CONST.c_i - seb.CONST.L_f);
                seb.STATVAR.evap = - Qe_evap ./(seb.CONST.rho_w .* L_w) .* seb.STATVAR.area(1);
                seb.TEMP.evap_energy = seb.STATVAR.evap .* seb.STATVAR.T(1) .* seb.CONST.c_w;
                
                Q_e = evap_fraction .* Qe_evap + sublim_fraction .* Qe_sublim;
            end
        end
        
        
        %latent heat flux from evaporation (not transpiration!), CLM 4.5,
        %to be used with Richards Equation scheme
        
        function seb = Q_evap_CLM4_5(seb, forcing)
            uz = forcing.TEMP.wind;
            p = forcing.TEMP.p;
            q = forcing.TEMP.q;
            Tz = forcing.TEMP.Tair;
            Lstar = seb.STATVAR.Lstar;
            
            z =  seb.PARA.airT_height;
            z0 = seb.PARA.z0;
            kappa = seb.CONST.kappa; 
            
            water_fraction = seb.STATVAR.water(1,1) ./ seb.STATVAR.waterIce(1,1);
            ice_fraction = seb.STATVAR.ice(1,1) ./ seb.STATVAR.waterIce(1,1);
            
            rho = rho_air(seb, p, seb.STATVAR.T(1)+273.15);
            sat_pressure_first_cell = water_fraction .* satPresWater(seb, seb.STATVAR.T(1)+273.15) + ice_fraction .* satPresWater(seb, seb.STATVAR.T(1)+273.15);
            latent_heat = water_fraction .* latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15) + ice_fraction .* latent_heat_sublimation(seb, seb.STATVAR.T(1)+273.15);

%             sat_pressure_first_cell = double(seb.STATVAR.T(1)>=0) .* satPresWater(seb, seb.STATVAR.T(1)+273.15) + double(seb.STATVAR.T(1)<0) .* satPresWater(seb, seb.STATVAR.T(1)+273.15);
%             latent_heat = double(seb.STATVAR.T(1)>=0) .* latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15) + double(seb.STATVAR.T(1)<0) .* latent_heat_sublimation(seb, seb.STATVAR.T(1)+273.15);

            
            %saturation_fraction_air_first_cell = exp(seb.STATVAR.matric_potential(1,1) .* 9.81 ./ ((8.3145e+03 ./ 18.016) .*vegetation.mlcanopyinst.tg(p)));
            saturation_fraction_air_first_cell = exp(seb.STATVAR.waterPotential(1,1) .* seb.CONST.g ./ ((seb.CONST.R./ seb.CONST.molar_mass_w) .*(seb.STATVAR.T(1)+273.15)));
            
%             saturation_fraction_air_first_cell
            %this might be wrong if the ground is frozen?
            q_first_cell = sat_pressure_first_cell .* saturation_fraction_air_first_cell ./ p;
            
            vol_water_first_cell = seb.STATVAR.water(1,1) ./ (seb.STATVAR.layerThick(1,1) .* seb.STATVAR.area(1,1)); 
            reduce_yes_no = vol_water_first_cell < seb.STATVAR.field_capacity(1,1) && forcing.TEMP.q < q_first_cell;
            betaCLM4_5 = 1 +  double(reduce_yes_no) .* (-1 +  0.25 .* (1-(cos(pi() .* vol_water_first_cell ./ seb.STATVAR.field_capacity(1,1)))).^2);
            
%             vol_water_first_cell
%             betaCLM4_5

            seb.STATVAR.Qe = -rho.*latent_heat.*betaCLM4_5.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q - q_first_cell)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            seb.STATVAR.evap = water_fraction .* seb.STATVAR.Qe ./ (latent_heat .* seb.CONST.rho_w);
            seb.STATVAR.sublim = ice_fraction .* seb.STATVAR.Qe ./ (latent_heat .* seb.CONST.rho_w);
            seb.STATVAR.evap_energy =  seb.STATVAR.evap.*  (double(seb.STATVAR.T(1,1)>=0) .* seb.CONST.c_w .* seb.STATVAR.T(1,1) + ...
                double(seb.STATVAR.T(1,1)<0) .* seb.CONST.c_i .* seb.STATVAR.T(1,1)); 
            seb.STATVAR.sublim_energy =  seb.STATVAR.sublim .* (seb.CONST.c_i .* seb.STATVAR.T(1,1) - seb.CONST.L_f); 
        end
        
        
        %sensible heat flux
        function Q_h = Q_h(seb, forcing)
          
            cp = seb.CONST.cp; 
            kappa = seb.CONST.kappa; 
            g = seb.CONST.g; 
            sigma = seb.CONST.sigma; 
            
            uz = forcing.TEMP.wind;
            z =  seb.PARA.airT_height;
            z0 = seb.PARA.z0;
            Tz = forcing.TEMP.Tair;
            TForcing = seb.STATVAR.T(1);
            Lstar = seb.STATVAR.Lstar;
            p = forcing.TEMP.p;
                        
            Tz=Tz+273.15;
            TForcing=TForcing+273.15;
            rho = rho_air(seb, p, Tz);

            Q_h  = -rho.*cp.*kappa.* uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)) .* (Tz-TForcing)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            
        end
        
        function rho = rho_air(seb, p, T) %air density [kg m^(-3)]
            rho = p./(287.058.*T); 
        end
        
        function L_w = latent_heat_evaporation(seb, T) %specific latent heat of evaporation of water [J/kg K]
            L_w = 1e3.*(2500.8 - 2.36.*(T - 273.15));  
        end
        
        function L_i = latent_heat_sublimation(seb, T_forcing) %1e3.*2834.1; %latent heat of sublimation, constant [J/kg K]
            L_i = seb.CONST.L_s; 
        end
        
        function p = satPresWater(seb, T) %saturation pressure water, Magnus formula
            p=0.622.* 6.112 .* 100 .* exp(17.62.*(T-273.15)./(243.12-273.15+T));
        end
        
        function p = satPresIce(seb, T) %saturation pressure ice, Magnus formula
            p= 0.622.*6.112.* 100.* exp(22.46.*(T-273.15)./(272.61-273.15+T));
        end
        
        %---penetration of short-wave radiation---
        function [seb, S_up] = penetrate_SW_no_transmission(seb, S_down) % called recursively by surfaceEnergyBalance-function of TOP_CLASS, "hard" surface absorbing all SW-radiation in 1st cell
            seb.TEMP.d_energy(1,1) = seb.TEMP.d_energy(1,1) + (1 - seb.PARA.albedo) .* sum(S_down); % .* seb.STATVAR.area(1,1);  %in [W]
            S_up = seb.PARA.albedo .* S_down;  % in [W/m2]
        end
        
        function [seb, S_up] = penetrate_SW_transmission_bulk(seb, S_down)  %used with fixed albedo and SW extinction coefficient
            %S_up and S_down can in principle be spectrally resolved when provided as
            %row array, using e.g. spectral_ranges = [0.71 0.21 0.08]; 
            %SW extinction is assumed constant throughout class 
            cut_off = 0.1.* seb.STATVAR.area(1,1); %[W/m2], radiation is not penetrated further if cutoff is reached
            
            S_up = seb.PARA.albedo .* S_down;
            S_down = (1 - seb.PARA.albedo) .* S_down;
            
            i=1;
            while i<=size(seb.STATVAR.layerThick,1) && sum(S_down) >= cut_off
                reduction = exp(-seb.PARA.SW_extinction .* seb.STATVAR.layerThick(i,1));  
                seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_down .* (1-reduction)); % .* seb.STATVAR.area(i,1);
                S_down = reduction .* S_down;
                i = i+1;
            end
            
            if sum(S_down) < cut_off  %all radiation absorbed, only S_up goes out
                i = min(i, size(seb.STATVAR.layerThick,1));
                seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_down);% .* seb.STATVAR.area(i,1);
            else %end of class reached, radiation penetrated on to next class
                [seb.NEXT, S_up2] = penetrate_SW(seb.NEXT, S_down); %call mandatory function recursively for following class
                %NOTE: no check performed that NEXT is not Bottom, must always be used with class below!
                i = size(seb.STATVAR.layerThick,1); %penetrate bottom up
                while i>=1 && sum(S_up2) >= cut_off
                    reduction = exp(-seb.PARA.SW_extinction .* seb.STATVAR.layerThick(i,1));
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_up2 .* (1-reduction));% .* seb.STATVAR.area(i,1);
                    S_up2 = reduction .* S_up2;
                    i = i-1;
                end
                if sum(S_up2) < cut_off  %all radiation absorbed, only S_up goes out
                    i = max(i, 1);
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_up2); % .* seb.STATVAR.area(i,1);
                else
                    S_up = S_up + S_up2;  %add the uppwelling SW radiation to reflected, multiple reflections not accounted for!! 
                end

            end         
        end
        
        
        function [seb, S_up] = penetrate_SW_transmission_spectral(seb, S_down)  %used with variable albedo and SW exticntion coefficient
            %S_up and S_down are spectrally resolved when provided as
            %row array, using e.g. spectral_ranges = [0.71 0.21 0.08]; 
            %SW extinction is assumed constant throughout class 
            cut_off = 0.1.* seb.STATVAR.area(1,1); %[W/m2], radiation is not penetrated further if cutoff is reached
            
            S_up = seb.TEMP.spectral_albedo .* S_down;
            S_down = (1 - seb.TEMP.spectral_albedo) .* S_down;
            
            i=1;
            while i<=size(seb.STATVAR.layerThick,1) && sum(S_down) >= cut_off
                reduction = exp(-seb.TEMP.SW_extinction(i,:) .* seb.STATVAR.layerThick(i,1));  
                seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_down .* (1-reduction)); % .* seb.STATVAR.area(i,1);
                S_down = reduction .* S_down;
                i = i+1;
            end
            
            if sum(S_down) < cut_off  %all radiation absorbed, only S_up goes out
                i = min(i, size(seb.STATVAR.layerThick,1));
                seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_down); % .* seb.STATVAR.area(i,1);
            else %end of class reached, radiation penetrated on to next class
                [seb.NEXT, S_up2] = penetrate_SW(seb.NEXT, S_down); %call mandatory function recursively for following class
                %NOTE: no check performed that NEXT is not Bottom, must always be used with class below!
                i = size(seb.STATVAR.layerThick,1); %penetrate bottom up
                while i>=1 && sum(S_up2) >= cut_off
                    reduction = exp(-seb.TEMP.SW_extinction(i,:) .* seb.STATVAR.layerThick(i,1));
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_up2 .* (1-reduction));% .* seb.STATVAR.area(i,1);
                    S_up2 = reduction .* S_up2;
                    i = i-1;
                end
                if sum(S_up2) < cut_off  %all radiation absorbed, only S_up goes out
                    i = max(i, 1);
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_up2); % .* seb.STATVAR.area(i,1);
                else
                    S_up = S_up + S_up2;  %add the uppwelling SW radiation to reflected, multiple reflections not accounted for!! 
                end
            end         
        end
        
        
        %calculates derivatives of water (i.e. water fluxes) due to evapotranspiration
        function seb = calculateET(seb)
            L_v = latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15);
            
            if seb.STATVAR.Qe_pot > 0
                fraction_T = getET_fraction(seb);
                fraction_E = getET_fraction(seb);

                depth_weighting_E = exp(-1./seb.PARA.evaporationDepth .* cumsum(seb.STATVAR.layerThick));  %exponential decrease with depth
                depth_weighting_E(depth_weighting_E<0.05) = 0;
                depth_weighting_E(1,1) = depth_weighting_E(1,1) + double(depth_weighting_E(1,1) == 0); %ADDED Sebastian, prevents depth_weighting to become all zero and NaN in the next line when first grid cell is very thick
                depth_weighting_E = depth_weighting_E .* seb.STATVAR.layerThick ./ sum(depth_weighting_E .* seb.STATVAR.layerThick,1); %normalize
                               
                depth_weighting_T = exp(-1./seb.PARA.rootDepth .* cumsum(seb.STATVAR.layerThick));
                depth_weighting_T(depth_weighting_T<0.05) = 0;
                depth_weighting_T(1,1) = depth_weighting_T(1,1) + double(depth_weighting_T(1,1) == 0); %ADDED Sebastian, prevents depth_weighting to become all zero and NaN in the next line when first grid cell is very thick
                depth_weighting_T = depth_weighting_T .* seb.STATVAR.layerThick ./ sum(depth_weighting_T .* seb.STATVAR.layerThick,1);

                fraction_ET = fraction_T .* depth_weighting_T .* seb.PARA.ratioET + fraction_E .* depth_weighting_E .* (1-seb.PARA.ratioET); %make this dependent on area?
                
                seb.STATVAR.Qe = sum(fraction_ET, 1) .* seb.STATVAR.Qe_pot;
                
                fraction_ET = fraction_ET./max(1e-12, sum(fraction_ET, 1));
                
                seb.TEMP.d_water_ET = seb.TEMP.d_water_ET - seb.STATVAR.Qe ./ (L_v.*seb.CONST.rho_w) .* fraction_ET .* seb.STATVAR.area;    %in m3 water per sec
                
            else  %condensation
                seb.STATVAR.Qe = seb.STATVAR.Qe_pot .*double(seb.STATVAR.T(1)>0);
                seb.TEMP.d_water_ET(1,1) = seb.TEMP.d_water_ET(1,1) - seb.STATVAR.Qe ./ (L_v.*seb.CONST.rho_w) .* seb.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
            end
            seb.TEMP.d_water_ET_energy = seb.TEMP.d_water_ET_energy + seb.TEMP.d_water_ET .* (double(seb.STATVAR.T>=0) .* seb.CONST.c_w + double(seb.STATVAR.T<0) .* seb.CONST.c_i) .* seb.STATVAR.T; %[J/sec]
        end
        
        
        function seb = calculateET_Xice(seb)
            L_v = latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15);
            
            if seb.STATVAR.Qe_pot > 0
                fraction_T = getET_fraction_Xice(seb);
                fraction_E = getET_fraction_Xice(seb);

                depth_weighting_E = exp(-1./seb.PARA.evaporationDepth .* cumsum(seb.STATVAR.layerThick));  %exponential decrease with depth
                depth_weighting_E(depth_weighting_E<0.05) = 0;
                depth_weighting_E(1,1) = depth_weighting_E(1,1) + double(depth_weighting_E(1,1) == 0); %ADDED Sebastian, prevents depth_weighting to become all zero and NaN in the next line when first grid cell is very thick
                depth_weighting_E = depth_weighting_E .* seb.STATVAR.layerThick ./ sum(depth_weighting_E .* seb.STATVAR.layerThick,1); %normalize
                               
                depth_weighting_T = exp(-1./seb.PARA.rootDepth .* cumsum(seb.STATVAR.layerThick));
                depth_weighting_T(depth_weighting_T<0.05) = 0;
                depth_weighting_T(1,1) = depth_weighting_T(1,1) + double(depth_weighting_T(1,1) == 0); %ADDED Sebastian, prevents depth_weighting to become all zero and NaN in the next line when first grid cell is very thick
                depth_weighting_T = depth_weighting_T .* seb.STATVAR.layerThick ./ sum(depth_weighting_T .* seb.STATVAR.layerThick,1);

                fraction_ET = fraction_T .* depth_weighting_T .* seb.PARA.ratioET + fraction_E .* depth_weighting_E .* (1-seb.PARA.ratioET); %make this dependent on area?
                
                seb.STATVAR.Qe = sum(fraction_ET, 1) .* seb.STATVAR.Qe_pot;
                
                fraction_ET = fraction_ET./max(1e-12, sum(fraction_ET, 1));
                
                seb.TEMP.d_water_ET = seb.TEMP.d_water_ET - seb.STATVAR.Qe ./ (L_v.*seb.CONST.rho_w) .* fraction_ET .* seb.STATVAR.area;    %in m3 water per sec
                
            else  %condensation
                seb.STATVAR.Qe = seb.STATVAR.Qe_pot.*double(seb.STATVAR.T(1)>0);
                seb.TEMP.d_water_ET(1,1) = seb.TEMP.d_water_ET(1,1) - seb.STATVAR.Qe ./ (L_v.*seb.CONST.rho_w) .* seb.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
            end
            seb.TEMP.d_water_ET_energy = seb.TEMP.d_water_ET_energy + seb.TEMP.d_water_ET .* (double(seb.STATVAR.T>=0) .* seb.CONST.c_w + double(seb.STATVAR.T<0) .* seb.CONST.c_i) .* seb.STATVAR.T; %[J/sec]
        end
        
        
        %calculate fraction of potential evapotranspiration for each grid cell, this reduces evapotranspiration if water content is lower than field capacity 
        function fraction = getET_fraction(seb)
            %waterC = seb.STATVAR.waterIce ./ seb.STATVAR.layerThick ./ max(1e-20, seb.STATVAR.area); %area can get zero if the area of SNOW CHILD is 100%
            waterC = seb.STATVAR.waterIce ./ seb.STATVAR.layerThick ./ seb.STATVAR.area;
            waterC(isnan(waterC)) = 0;
            fraction=double(seb.STATVAR.T>0).*double(seb.STATVAR.T(1)>0).*double(waterC > 0.03) .* (double(waterC >= seb.STATVAR.field_capacity) + double(waterC < seb.STATVAR.field_capacity).*0.25.*(1-cos(pi().*waterC./seb.STATVAR.field_capacity)).^2);

        end

        function fraction = getET_fraction_Xice(seb)
            
            waterC = seb.STATVAR.waterIce ./ (seb.STATVAR.layerThick .* seb.STATVAR.area - seb.STATVAR.XwaterIce);
            %waterC = seb.STATVAR.waterIce ./ max(1e-20,seb.STATVAR.layerThick .* seb.STATVAR.area - seb.STATVAR.XwaterIce);
            waterC(isnan(waterC)) = 0;
            fraction=double(seb.STATVAR.T>0).*double(seb.STATVAR.T(1)>0).*double(waterC > 0.03).*(double(waterC >= seb.STATVAR.field_capacity) + double(waterC < seb.STATVAR.field_capacity).*0.25.*(1-cos(pi().*waterC./seb.STATVAR.field_capacity)).^2);
        end
        
        %evaporation for water body
        function seb = calculateET_LAKE(seb)
            L_v = latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15);
            F_ub_water = -seb.STATVAR.Qe ./ (L_v.*seb.CONST.rho_w) .* seb.STATVAR.area(1,1);
            F_ub_water_energy = F_ub_water .* seb.CONST.c_w .* seb.STATVAR.T(1,1); %[J/sec]
            
            seb.TEMP.d_water(1) = seb.TEMP.d_water(1)  + F_ub_water;
            seb.TEMP.d_water_energy(1) = seb.TEMP.d_water_energy(1)  + F_ub_water_energy;
        end
        
    end
end

