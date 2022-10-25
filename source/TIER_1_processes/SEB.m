 %========================================================================
% CryoGrid TIER1 library class, functions related to surface energy balance
% NOTE: also contains functions for penetration of short-wave radiation
% S. Westermann, October 2020
%========================================================================


classdef SEB < BASE
    
    methods
        
% -------------- Monin - Obukhov similarity ------------------------------        
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
        
        function seb = L_star_canopy(seb, tile) % Same as L_star(..), but allowing z0 to be a variable
            forcing = tile.FORCING;
            uz = forcing.TEMP.wind;
            z =  seb.PARA.airT_height + seb.STATVAR.layerThick;
            z0 = seb.STATVAR.z0;
            d = seb.STATVAR.d;
            Tz = forcing.TEMP.Tair+273.15;
            Lstar = seb.STATVAR.Lstar;
            p = forcing.TEMP.p;
            Qh = seb.STATVAR.Qh+seb.NEXT.STATVAR.Qh; % Qh of whole system (canopy + ground)
            Qe = seb.STATVAR.Qe+seb.NEXT.STATVAR.Qe;
            
            rho = rho_air_moist(seb, tile);
            cp = seb.CONST.cp;
            kappa = seb.CONST.kappa; %0.4;
            g = seb.CONST.g; %9.81;
            
            if Tz >=273.15
                L = latent_heat_evaporation(seb, Tz); %1e3.*(2500.8 - 2.36.*(Tz-273.15));  %latent heat of evaporation of water
            else
                L = latent_heat_sublimation(seb, Tz); %1e3.*2834.1; %latent heat of sublimation
            end
            
            u_star = real(uz.*kappa./(log(z./z0)- psi_M_CLM5(seb, z./Lstar, z0./Lstar)));
            L_star = real(-rho.*cp.*Tz./kappa./g.*u_star.^3./(Qh + 0.61.*cp./L.*Tz.*Qe));
            L_star=(abs(L_star)<1e-7).*L_star./abs(L_star).*1e-7 + (abs(L_star)>=1e-7).*L_star;  %limits Lstar
            
            L_star(Qh == 0 && Qe == 0) = 1e10; % avoid L_star = NaN for zero turbulent fluxes
            
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
        
        function result = psi_H_CLM5(seb, zeta1, zeta2) % heat/vapor stability function
            
            zeta_h = seb.PARA.zeta_h;
            if zeta1 <= 0 % Unstable
                if zeta1 < zeta_h % very unstable
                    result = -log(-zeta_h) + log(-zeta1) - 0.8*( (-zeta_h)^(-1/3) - (-zeta1)^(-1/3) ) ...
                        + 2*log( 1/2 + (1-16*zeta_h)^.5 ) - 2*log( 1/2 + (1-16*zeta2)^.5 );
                else
                    result = 2*log( 1/2 + (1-16*zeta1)^.5 ) - 2*log( 1/2 + (1-16*zeta2)^.5 );
                end
            else % stable
                if zeta1 > 1 % very stable
                    result = -4*log(zeta1) - 5 - zeta1 + 1 + 5*zeta2;
                else
                    result = -5*zeta1 + 5*zeta2 ;
                end
            end
        end
        
        function result = psi_M_CLM5(seb, zeta1, zeta2) % momentum stability function
            zeta_m = seb.PARA.zeta_m;
            
            if zeta1 <= 0 % Unstable
                
                if zeta1 < zeta_m % very unstable
                    result = -log(-zeta_m) + log(-zeta1) - 1.14*( (-zeta1)^(1/3) - (-zeta_m)^(1/3) ) ...
                        + (2*log(1/2 + (1-16*zeta_m)^.25/2) + log(1/2 + (1-16*zeta_m)^.5/2) -2*atan((1-16*zeta_m)^.25) + pi/2) ...
                        - (2*log(1/2 + (1-16*zeta2)^.25/2) + log(1/2 + (1-16*zeta2)^.5/2) -2*atan((1-16*zeta2)^.25) + pi/2);
                else
                    result = (2*log(1/2 + (1-16*zeta1)^.25/2) + log(1/2 + (1-16*zeta1)^.5/2) -2*atan((1-16*zeta1)^.25) + pi/2) ...
                        - (2*log(1/2 + (1-16*zeta2)^.25/2) + log(1/2 + (1-16*zeta2)^.5/2) -2*atan((1-16*zeta2)^.25) + pi/2);
                end
            else % stable
                if zeta1 > 1 % very stable
                    result = -4*log(zeta1) - 5 - zeta1 + 1 + 5*zeta2;
                else
                    result = 5*zeta2 - 5*zeta1;
                end
            end
        end
        
% ------------- Latent heat fluxes ---------------------------------------        
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
                Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q-.622.*satPresIce(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            else
                Q_e = -rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb,z./Lstar, z0./Lstar)).*(q-.622.*satPresWater(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar)  ...
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
                Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q-.622.*satPresIce(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            else
                Q_e = -rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q-.622.*satPresWater(seb, TForcing)./p)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
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
            sat_pressure_first_cell = water_fraction .* satPresWater(seb, seb.STATVAR.T(1)+273.15) + ice_fraction .* satPresIce(seb, seb.STATVAR.T(1)+273.15);
            latent_heat = water_fraction .* latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15) + ice_fraction .* latent_heat_sublimation(seb, seb.STATVAR.T(1)+273.15);
            
            %             sat_pressure_first_cell = double(seb.STATVAR.T(1)>=0) .* satPresWater(seb, seb.STATVAR.T(1)+273.15) + double(seb.STATVAR.T(1)<0) .* satPresWater(seb, seb.STATVAR.T(1)+273.15);
            %             latent_heat = double(seb.STATVAR.T(1)>=0) .* latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15) + double(seb.STATVAR.T(1)<0) .* latent_heat_sublimation(seb, seb.STATVAR.T(1)+273.15);
            
            
            %saturation_fraction_air_first_cell = exp(seb.STATVAR.matric_potential(1,1) .* 9.81 ./ ((8.3145e+03 ./ 18.016) .*vegetation.mlcanopyinst.tg(p)));
            saturation_fraction_air_first_cell = exp(seb.STATVAR.waterPotential(1,1) .* seb.CONST.g ./ ((seb.CONST.R./ seb.CONST.molar_mass_w) .*(seb.STATVAR.T(1)+273.15)));
            
            %             saturation_fraction_air_first_cell
            %this might be wrong if the ground is frozen?
            q_first_cell = sat_pressure_first_cell .* saturation_fraction_air_first_cell ./ p;
            
            %             betaCLM4_5 = soil_resistance_beta(seb.NEXT,forcing); Does not
            %             work both with and without canopy
            betaCLM4_5 = soil_resistance_beta(seb,forcing);
            %             vol_water_first_cell = seb.STATVAR.water(1,1) ./ (seb.STATVAR.layerThick(1,1) .* seb.STATVAR.area(1,1));
            %             reduce_yes_no = vol_water_first_cell < seb.STATVAR.field_capacity(1,1) && forcing.TEMP.q < q_first_cell;
            %             betaCLM4_5 = 1 +  double(reduce_yes_no) .* (-1 +  0.25 .* (1-(cos(pi() .* vol_water_first_cell ./ seb.STATVAR.field_capacity(1,1)))).^2);
            %
            %             vol_water_first_cell
            %             betaCLM4_5
            
            seb.STATVAR.Qe = -rho.*latent_heat.*betaCLM4_5.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q - q_first_cell)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            seb.STATVAR.evap = water_fraction .* seb.STATVAR.Qe ./ (latent_heat .* seb.CONST.rho_w);
            seb.STATVAR.sublim = ice_fraction .* seb.STATVAR.Qe ./ (latent_heat .* seb.CONST.rho_w);
            seb.STATVAR.evap_energy =  seb.STATVAR.evap.*  (double(seb.STATVAR.T(1,1)>=0) .* seb.CONST.c_w .* seb.STATVAR.T(1,1) + ...
                double(seb.STATVAR.T(1,1)<0) .* seb.CONST.c_i .* seb.STATVAR.T(1,1));
            seb.STATVAR.sublim_energy =  seb.STATVAR.sublim .* (seb.CONST.c_i .* seb.STATVAR.T(1,1) - seb.CONST.L_f);
        end
        
        function seb = Q_evap_CLM4_5_Xice(seb, forcing)
            uz = forcing.TEMP.wind;
            p = forcing.TEMP.p;
            q = forcing.TEMP.q;
            Tz = forcing.TEMP.Tair;
            Lstar = seb.STATVAR.Lstar;
            
            z =  seb.PARA.airT_height;
            z0 = seb.PARA.z0;
            kappa = seb.CONST.kappa;
            
            water_fraction = (seb.STATVAR.water(1,1)+seb.STATVAR.Xwater(1,1) ) ./ (seb.STATVAR.waterIce(1,1) + seb.STATVAR.XwaterIce(1,1));
            ice_fraction = (seb.STATVAR.ice(1,1) + seb.STATVAR.Xice(1,1)) ./ (seb.STATVAR.waterIce(1,1) + seb.STATVAR.XwaterIce(1,1));
            
            rho = rho_air(seb, p, seb.STATVAR.T(1)+273.15);
            sat_pressure_first_cell = water_fraction .* satPresWater(seb, seb.STATVAR.T(1)+273.15) + ice_fraction .* satPresWater(seb, seb.STATVAR.T(1)+273.15);
            latent_heat = water_fraction .* latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15) + ice_fraction .* latent_heat_sublimation(seb, seb.STATVAR.T(1)+273.15);
            
            saturation_fraction_air_first_cell = exp(seb.STATVAR.waterPotential(1,1) .* seb.CONST.g ./ ((seb.CONST.R./ seb.CONST.molar_mass_w) .*(seb.STATVAR.T(1)+273.15)));
            
            %             saturation_fraction_air_first_cell
            %this might be wrong if the ground is frozen?
            q_first_cell = sat_pressure_first_cell .* saturation_fraction_air_first_cell ./ p;
            
            vol_water_first_cell = seb.STATVAR.water(1,1) ./ (seb.STATVAR.layerThick(1,1) .* seb.STATVAR.area(1,1) - seb.STATVAR.XwaterIce(1,1));
            reduce_yes_no = vol_water_first_cell < seb.STATVAR.field_capacity(1,1) && forcing.TEMP.q < q_first_cell && ~(seb.STATVAR.XwaterIce(1,1) > 1e-5.*seb.STATVAR.area(1,1));
            betaCLM4_5 = 1 +  double(reduce_yes_no) .* (-1 +  0.25 .* (1-(cos(pi() .* vol_water_first_cell ./ seb.STATVAR.field_capacity(1,1)))).^2);
            
            seb.STATVAR.Qe = -rho.*latent_heat.*betaCLM4_5.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q - q_first_cell)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            seb.STATVAR.evap = water_fraction .* seb.STATVAR.Qe ./ (latent_heat .* seb.CONST.rho_w);
            seb.STATVAR.sublim = ice_fraction .* seb.STATVAR.Qe ./ (latent_heat .* seb.CONST.rho_w);
            seb.STATVAR.evap_energy =  seb.STATVAR.evap.*  (double(seb.STATVAR.T(1,1)>=0) .* seb.CONST.c_w .* seb.STATVAR.T(1,1) + ...
                double(seb.STATVAR.T(1,1)<0) .* seb.CONST.c_i .* seb.STATVAR.T(1,1));
            seb.STATVAR.sublim_energy =  seb.STATVAR.sublim .* (seb.CONST.c_i .* seb.STATVAR.T(1,1) - seb.CONST.L_f);
        end
        
        function seb = Q_evap_CLM4_5_Xice(seb, forcing)
            uz = forcing.TEMP.wind;
            p = forcing.TEMP.p;
            q = forcing.TEMP.q;
            Tz = forcing.TEMP.Tair;
            Lstar = seb.STATVAR.Lstar;
            
            z =  seb.PARA.airT_height;
            z0 = seb.PARA.z0;
            kappa = seb.CONST.kappa; 
            
            water_fraction = (seb.STATVAR.water(1,1)+seb.STATVAR.Xwater(1,1) ) ./ (seb.STATVAR.waterIce(1,1) + seb.STATVAR.XwaterIce(1,1));
            ice_fraction = (seb.STATVAR.ice(1,1) + seb.STATVAR.Xice(1,1)) ./ (seb.STATVAR.waterIce(1,1) + seb.STATVAR.XwaterIce(1,1));
            
            rho = rho_air(seb, p, seb.STATVAR.T(1)+273.15);
            sat_pressure_first_cell = water_fraction .* satPresWater(seb, seb.STATVAR.T(1)+273.15) + ice_fraction .* satPresWater(seb, seb.STATVAR.T(1)+273.15);
            latent_heat = water_fraction .* latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15) + ice_fraction .* latent_heat_sublimation(seb, seb.STATVAR.T(1)+273.15);

            saturation_fraction_air_first_cell = exp(seb.STATVAR.waterPotential(1,1) .* seb.CONST.g ./ ((seb.CONST.R./ seb.CONST.molar_mass_w) .*(seb.STATVAR.T(1)+273.15)));
            
%             saturation_fraction_air_first_cell
            %this might be wrong if the ground is frozen?
            q_first_cell = sat_pressure_first_cell .* saturation_fraction_air_first_cell ./ p;
            
            vol_water_first_cell = seb.STATVAR.water(1,1) ./ (seb.STATVAR.layerThick(1,1) .* seb.STATVAR.area(1,1) - seb.STATVAR.XwaterIce(1,1)); 
            reduce_yes_no = vol_water_first_cell < seb.STATVAR.field_capacity(1,1) && forcing.TEMP.q < q_first_cell && ~(seb.STATVAR.XwaterIce(1,1) > 1e-5.*seb.STATVAR.area(1,1));
            betaCLM4_5 = 1 +  double(reduce_yes_no) .* (-1 +  0.25 .* (1-(cos(pi() .* vol_water_first_cell ./ seb.STATVAR.field_capacity(1,1)))).^2);

            seb.STATVAR.Qe = -rho.*latent_heat.*betaCLM4_5.*kappa.*uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)).*(q - q_first_cell)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            seb.STATVAR.evap = water_fraction .* seb.STATVAR.Qe ./ (latent_heat .* seb.CONST.rho_w);
            seb.STATVAR.sublim = ice_fraction .* seb.STATVAR.Qe ./ (latent_heat .* seb.CONST.rho_w);
            seb.STATVAR.evap_energy =  seb.STATVAR.evap.*  (double(seb.STATVAR.T(1,1)>=0) .* seb.CONST.c_w .* seb.STATVAR.T(1,1) + ...
                double(seb.STATVAR.T(1,1)<0) .* seb.CONST.c_i .* seb.STATVAR.T(1,1)); 
            seb.STATVAR.sublim_energy =  seb.STATVAR.sublim .* (seb.CONST.c_i .* seb.STATVAR.T(1,1) - seb.CONST.L_f); 
        end
        % Latent heat flux from evaporation + transpiration as in CLM5
        function seb = Q_e_CLM5(seb, tile)
            forcing = tile.FORCING;
            L = seb.STATVAR.LAI; % Leaf area index
            S = seb.STATVAR.SAI; % Stem area index
            Tmfw = forcing.CONST.Tmfw;
            Tv = seb.STATVAR.T(1); % leaf temperature
            q_atm = forcing.TEMP.q; % atm. specific humidity
            p = forcing.TEMP.p; % surface pressure
            rho_atm = seb.TEMP.rho_atm; % air density
            f_dry = seb.STATVAR.f_dry;
            f_wet = seb.STATVAR.f_wet;
            
            q_g = get_humidity_ground(seb.IA_NEXT, tile); % ground specific humidity
            e_v = double(Tv>=0).*satPresWater(seb,Tv+Tmfw) + double(Tv<0).*satPresIce(seb,Tv+Tmfw); % saturation water pressure of leafs
            qs_Tv = .622.*e_v./p; % saturation water vapor specific humidity at leaf temperature
            
            r_aw = seb.TEMP.r_a; % aerodynamic resistance against water vapor transport canopy air - atmosphere
            r_aw_prime = seb.TEMP.r_a_prime; % aerodynamic resistance against water vapor transport soil - canopy air
            r_b = seb.TEMP.r_b; % leaf boudnary layer resitance against vapor transport
            r_total = seb.TEMP.r_total;
            r_soil = seb.TEMP.r_soil;
            
            ca = 1/r_aw;
            cv = 1/r_total;
            cg = 1./(r_aw_prime + r_soil);
            
            q_s = (q_atm.*ca + q_g.*cg + qs_Tv.*cv)./(ca + cv + cg); % Eq. 5.108, canopy specific humidity
            
            seb.STATVAR.q_s = q_s;
            seb.STATVAR.qs_Tv = qs_Tv;
            
            water_fraction = seb.STATVAR.water(1) ./ seb.STATVAR.waterIce(1);
            ice_fraction = seb.STATVAR.ice(1) ./ seb.STATVAR.waterIce(1);
            water_fraction(seb.STATVAR.waterIce == 0) = double(Tv>=0);
            ice_fraction(seb.STATVAR.waterIce == 0) = double(Tv<0);
            
            if seb.TEMP.Ev_pot > 0 % evaporation/transpiration in m3 water /(m2*s), not kg/(m2*s) as in CLM5
                seb.TEMP.evap = -rho_atm.*(q_s - qs_Tv).* f_wet.*water_fraction.*(L+S) ./r_b ./seb.CONST.rho_w;
                seb.TEMP.sublim = -rho_atm.*(q_s - qs_Tv).* f_wet.*ice_fraction.*(L+S) ./r_b ./seb.CONST.rho_w; % Is this correct? or use f_snow?
                seb.TEMP.transp_sun = -rho_atm.*(q_s - qs_Tv).* f_dry.*L_sun ./(r_b+r_sun) ./seb.CONST.rho_w;
                seb.TEMP.transp_sha = -rho_atm.*(q_s - qs_Tv).* f_dry.*L_sha ./(r_b+r_sha) ./seb.CONST.rho_w;
            else % condensation/deposition
                seb.TEMP.evap = -rho_atm.*(q_s - qs_Tv).* water_fraction.*(L+S) ./r_b ./seb.CONST.rho_w;
                seb.TEMP.sublim = -rho_atm.*(q_s - qs_Tv).* ice_fraction.*(L+S) ./r_b ./seb.CONST.rho_w;
                seb.TEMP.transp_sun = 0;
                seb.TEMP.transp_sha = 0;
            end
            seb.TEMP.Qe_canopy = latent_heat_evaporation(seb, Tv+Tmfw).*seb.TEMP.evap.*seb.CONST.rho_w + latent_heat_sublimation(seb, Tv+Tmfw).*seb.TEMP.sublim.*seb.CONST.rho_w;
            seb.TEMP.Qt_sun = latent_heat_evaporation(seb, Tv+Tmfw).*seb.TEMP.transp_sun.*seb.CONST.rho_w;
            seb.TEMP.Qt_sha = latent_heat_evaporation(seb, Tv+Tmfw).*seb.TEMP.transp_sha.*seb.CONST.rho_w;
            
            seb.STATVAR.Qe = seb.TEMP.Qe_canopy + seb.TEMP.Qt_sun + seb.TEMP.Qt_sha;
            seb.TEMP.evap_energy = seb.TEMP.evap.*( double(Tv>=0).*seb.CONST.c_w.*Tv + double(Tv<0).*seb.CONST.c_i.*Tv);
            seb.TEMP.sublim_energy = seb.TEMP.sublim .* (seb.CONST.c_i.*Tv - seb.CONST.L_f);
            seb.TEMP.transp_energy = (seb.TEMP.transp_sun+seb.TEMP.transp_sha).*seb.CONST.c_w.*Tv; % assumes water within the leaf is always unfrozen when transpiring
            
            if seb.STATVAR.waterIce > seb.PARA.Wmax.*(L+S).*seb.STATVAR.area && seb.TEMP.evap < 0
                seb.TEMP.evap = 0;
                seb.TEMP.evap_energy = 0;
            elseif seb.STATVAR.waterIce > seb.PARA.Wmax.*(L+S).*seb.STATVAR.area && seb.TEMP.sublim < 0
                seb.TEMP.sublim = 0;
                seb.TEMP.sublim_energy = 0;
            end
        end
        
        % Latent heat from evaporation + transpiration as in CLM5, but with
        % stomatal resistance from Stewart (1988)
        function seb = Q_e_CLM5_Stewart(seb, tile)
            % evaporation and transpiration following CLM5, but with
            % simpler stomata/leaf resistance after Stewart (1988)
            % R. B. Zweigel, December 2021
            forcing = tile.FORCING;
            L = seb.STATVAR.LAI; % Leaf area index
            S = seb.STATVAR.SAI; % Stem area index
            Tmfw = forcing.CONST.Tmfw;
            Tv = seb.STATVAR.T(1); % leaf temperature
            q_atm = forcing.TEMP.q; % atm. specific humidity
            p = forcing.TEMP.p; % surface pressure
            rho_atm = seb.TEMP.rho_atm; % air density 
            f_wet = seb.STATVAR.f_wet; % wetted fraction of canopy
            f_dry = seb.STATVAR.f_dry; % dry and transpiring fraction of leaves
            
            % Resistance terms
            r_aw = seb.TEMP.r_a; % aerodynamic resistance against water vapor transport canopy air - atmosphere
            r_aw_prime = seb.TEMP.r_a_prime; % aerodynamic resistance against water vapor transport soil - canopy air
            r_b = seb.TEMP.r_b; % leaf boudnary layer resitance against vapor transport
            r_total = seb.TEMP.r_total;
            r_soil = seb.TEMP.r_soil;
            r_canopy = seb.TEMP.r_canopy;
            
            % Humidity
            q_g = get_humidity_surface(seb.IA_NEXT, tile); % ground/snow surface specific humidity
            e_v = double(Tv>=0).*satPresWater(seb,Tv+Tmfw) + double(Tv<0).*satPresIce(seb,Tv+Tmfw); % saturation water pressure of leafs
            qs_Tv = .622.*e_v./p; % saturation water vapor specific humidity at leaf temperature
            ca = 1/r_aw;
            cv = 1/r_total;
            cg = 1./(r_aw_prime + r_soil);
            
            q_s = (q_atm.*ca + q_g.*cg + qs_Tv.*cv)./(ca + cv + cg); % Eq. 5.108, canopy air specific humidity
            
            seb.STATVAR.q_s = q_s;
            seb.STATVAR.qs_TV = qs_Tv;
            seb.TEMP.q_atm = q_atm;
            
            water_fraction = seb.STATVAR.water(1) ./ seb.STATVAR.waterIce(1);
            ice_fraction = seb.STATVAR.ice(1) ./ seb.STATVAR.waterIce(1);
            water_fraction(seb.STATVAR.waterIce == 0) = double(Tv>=0);
            ice_fraction(seb.STATVAR.waterIce == 0) = double(Tv<0);
            
            if q_s - qs_Tv < 0 % evaporation/transpiration in m3 water /(m2*s), not kg/(m2*s) as in CLM5
                seb.TEMP.evap = -rho_atm.*(q_s - qs_Tv).* f_wet.*water_fraction.*(L+S) ./r_b ./seb.CONST.rho_w;
                seb.TEMP.sublim = -rho_atm.*(q_s - qs_Tv).* f_wet.*ice_fraction.*(L+S) ./r_b ./seb.CONST.rho_w; % Is this correct? or use f_snow?
                seb.TEMP.transp = -rho_atm.*(q_s - qs_Tv).* f_dry.*(L+S) ./ (r_b + r_canopy) ./ seb.CONST.rho_w;
            else % condensation/deposition
                seb.TEMP.evap = -rho_atm.*(q_s - qs_Tv).* water_fraction.*(L+S) ./r_b ./seb.CONST.rho_w;
                seb.TEMP.sublim = -rho_atm.*(q_s - qs_Tv).* ice_fraction.*(L+S) ./r_b ./seb.CONST.rho_w;
                seb.TEMP.transp = 0;
            end
            
            seb.TEMP.Qe_evap = latent_heat_evaporation(seb, Tv+Tmfw).*seb.TEMP.evap.*seb.CONST.rho_w + latent_heat_sublimation(seb, Tv+Tmfw).*seb.TEMP.sublim.*seb.CONST.rho_w;
            seb.TEMP.Qe_transp = latent_heat_evaporation(seb, Tv+Tmfw).*seb.TEMP.transp.*seb.CONST.rho_w;
            
            seb.STATVAR.Qe = seb.TEMP.Qe_evap + seb.TEMP.Qe_transp;
            seb.TEMP.evap_energy = seb.TEMP.evap.*( double(Tv>=0).*seb.CONST.c_w.*Tv + double(Tv<0).*seb.CONST.c_i.*Tv);
            seb.TEMP.sublim_energy = seb.TEMP.sublim .* (seb.CONST.c_i.*Tv - seb.CONST.L_f);
            seb.TEMP.transp_energy = seb.TEMP.transp.*seb.CONST.c_w.*Tv; 
        end
        
% -------------- Sensible heat fluxes ------------------------------------        
        %sensible heat flux for ground surface
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
            
            Tz=Tz+seb.CONST.Tmfw;
            TForcing=TForcing+seb.CONST.Tmfw;
            rho = rho_air(seb, p, Tz);
            
            Q_h  = -rho.*cp.*kappa.* uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)) .* (Tz-TForcing)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            
        end
        
        % Sensible heat flux for canopy as in CLM5
        function seb = Q_h_CLM5(seb, tile)
            forcing = tile.FORCING;
            cp = seb.CONST.cp;
            L = seb.STATVAR.LAI;
            S = seb.PARA.SAI;
            Tv = seb.STATVAR.T(1)+forcing.CONST.Tmfw;
%             Tg = seb.NEXT.STATVAR.T(1)+forcing.CONST.Tmfw; % ground temperature
            Tz = forcing.TEMP.Tair+forcing.CONST.Tmfw; % air temperature in Kelvin
            z = forcing.PARA.airT_height; % height above ground for measurements
            z = z + sum(seb.STATVAR.layerThick); % adjust so forcing height is above canopy
            Gamma = seb.CONST.Gamma_dry; % negative of dry adiabatic lapse rate
            r_b = seb.TEMP.r_b;
            r_ah = seb.TEMP.r_a; % r_ah = r_aw, resistance to sensible heat transfer canopy air - atmosphere
            r_ah_prime = seb.TEMP.r_a_prime;
            rho_atm = seb.TEMP.rho_atm; % air density
            Tz_pot = Tz+Gamma.*z; % Eq. 5.7 potential temperature
            Tg = get_surface_T(seb.NEXT, tile); % ground/surface temperature
            Tg = Tg + forcing.CONST.Tmfw;
            
            Ts = (Tz_pot./r_ah + Tg./r_ah_prime + Tv.*2*(L+S)./r_b)./(1./r_ah + 1/r_ah_prime + 2*(L+S)./r_b); % Eq. 5.93 canopy air temperature
            
            seb.STATVAR.Ts = Ts - forcing.CONST.Tmfw;
            seb.STATVAR.Qh = -rho_atm.*cp.*(Ts-Tv).*2*(L+S)./r_b; % Eq. 5.88   
            seb.TEMP.Tair = forcing.TEMP.Tair;
        end
        
% -------------- Evaporation/transpiration resistances -------------------
        function seb = canopy_resistances_CLM5(seb, tile)
            Throw error
            forcing = tile.FORCING;
            p = forcing.TEMP.p; % surface pressure
            L = seb.STATVAR.LAI; % Leaf area index
            S = seb.STATVAR.SAI; % Stem area index
            q_s = seb.STATVAR.q_s; % canopy air specific humidity (last timestep)
            L_sun = seb.TEMP.L_sun; % sunlit LAI
            L_sha = seb.TEMP.L_sha; % shaded LAI
            f_dry = seb.STATVAR.f_dry; % dry fraction of canopy
            f_wet = seb.STATVAR.f_wet; % wet fraction of canopy
            Tv = seb.STATVAR.T(1)+forcing.CONST.Tmfw; % leaf temperature
            kappa = seb.CONST.kappa; % van Karman constant
            u_star = seb.STATVAR.u_star; % from Monin-Ubukhov
            Lstar = seb.STATVAR.Lstar; % Monin-Ubukhov length
            Cv = seb.PARA.Cv; % Turbulent transfer coefficient canopy surface - canopy air
            d_leaf = seb.PARA.d_leaf; % characteristic dimension of the leaves in the direction of wind flow [m]
            ypsilon = seb.CONST.ypsilon; % kinematic viscosity of air
            Cs_dense = seb.PARA.Cs_dense; %  dense canopy turbulent transfer coefficient
            z = forcing.PARA.airT_height; % height above ground for measurements
            z0v = seb.STATVAR.z0; % roughness length of vegetation
            z0g = get_z0_ground(seb.IA_NEXT); % Roughness lenght of ground
            z = z + sum(seb.STATVAR.layerThick); % adjust so forcing height is above canopy
            beta_t = seb.TEMP.beta_t; % soil moisture stress function
            
            e_v = double(Tv>=forcing.CONST.Tmfw).*satPresWater(seb,Tv) + double(Tv<forcing.CONST.Tmfw).*satPresIce(seb,Tv); % saturation water pressure of leafs
            qs_Tv = .622.*e_v./p;
            seb.TEMP.rho_atm = rho_air_moist(seb,tile);
            Cs_bare = kappa./.13.*(z0g.*u_star./ypsilon).^(-.45); % Eq.5.121 bare soil turbulent transfer coefficient
            W = exp(-(L+S));% Eq. 5.119
            Cs = Cs_bare.*W + Cs_dense.*(1-W); % Eq. 5.118 turbulent transfer coefficient between soil and canopy air
            
            seb.TEMP.r_a = 1./(kappa^2.*u_star).*(log(z./z0v)- psi_M(seb, z./Lstar, z0v./Lstar)).*(log(z./z0v)- psi_H(seb, z./Lstar, z0v./Lstar)); % aerodynamic resistance to heat/water vapor fransper canopy air - atmosphere. From original CG3 publication
            seb.TEMP.r_b = 1./Cv*(u_star./d_leaf).^(-.5); % Eq. 5.122 leaf boundary layer resistance
            seb.TEMP.r_a_prime = 1/(Cs.*u_star); % Eq. 5.116 aerodynamic resistance to heat/water vapor transfer soil - canopy air
            seb.TEMP.r_soil = ground_resistance_evap(seb.IA_NEXT, tile); % resistance to water vapor flux within the soil matrix
            seb.TEMP.rr_dry = f_dry.*seb.TEMP.r_b./L .* (L_sun./(seb.TEMP.r_b+seb.TEMP.r_sun) + L_sha./(seb.TEMP.r_b+seb.TEMP.r_sha));
            seb.TEMP.Ev_pot = -seb.TEMP.rho_atm.*(q_s - qs_Tv)./seb.TEMP.r_b;
            seb.TEMP.rr = double(seb.TEMP.Ev_pot<=0) + double(seb.TEMP.Ev_pot>0).*(f_wet + double(beta_t>0).*seb.TEMP.rr_dry);
            if L > 0
                seb.TEMP.r_total = seb.TEMP.r_b./(seb.TEMP.rr.*(L+S));
            else % No leaves - no transpiration
                seb.TEMP.r_total = seb.TEMP.r_b./S;
            end
        end
        
        function seb = canopy_resistances_CLM5_Stewart(seb, tile)
            % similar to canopy_resistances_CLM5(), but with resistances 
            % against transpiration from Stewart (1988) as in Dingman (2015)     
            forcing = tile.FORCING;
            uz = forcing.TEMP.wind; % atm. wind speed
            p = forcing.TEMP.p; % surface pressure
            L = seb.STATVAR.LAI; % Leaf area index
            S = seb.STATVAR.SAI; % Stem area index
            q_s = seb.STATVAR.q_s; % canopy air specific humidity (last timestep)
            f_dry = seb.STATVAR.f_dry; % dry fraction of canopy
            f_wet = seb.STATVAR.f_wet; % wet fraction of canopy
            Tv = seb.STATVAR.T(1)+forcing.CONST.Tmfw; % leaf temperature
            kappa = seb.CONST.kappa; % van Karman constant
            Lstar = seb.STATVAR.Lstar; % Monin-Ubukhov length
            Cv = seb.PARA.Cv; % Turbulent transfer coefficient canopy surface - canopy air
            d_leaf = seb.PARA.d_leaf; % characteristic dimension of the leaves in the direction of wind flow [m]
            ypsilon = seb.CONST.ypsilon; % kinematic viscosity of air
            Cs_dense = seb.PARA.Cs_dense; %  dense canopy turbulent transfer coefficient
            z = forcing.PARA.airT_height + sum(seb.STATVAR.layerThick); % height above ground for measurements
            z0v = seb.STATVAR.z0; % roughness length of vegetation
            z0g = get_z0_surface(seb.NEXT); % Roughness lenght of ground
            d = seb.STATVAR.d; % displacement heigh
            k_s = seb.PARA.k_shelter; % canopy sheltering coefficient
            
            u_star = real(uz.*kappa./(log((z-d)./z0v)- psi_M_CLM5(seb, (z-d)./Lstar, z0v./Lstar)));
            
            e_v = double(Tv>=forcing.CONST.Tmfw).*satPresWater(seb,Tv) + double(Tv<forcing.CONST.Tmfw).*satPresIce(seb,Tv); % saturation water pressure of leafs
            qs_Tv = .622.*e_v./p;
            seb.TEMP.rho_atm = rho_air_moist(seb,tile);
            Cs_bare = kappa./.13.*(z0g.*u_star./ypsilon).^(-.45); % Eq.5.121 bare soil turbulent transfer coefficient
            W = exp(-(L+S));% Eq. 5.119
            Cs = Cs_bare.*W + Cs_dense.*(1-W); % Eq. 5.118 turbulent transfer coefficient between soil and canopy air
            seb.TEMP.C_leaf = leaf_conductance_Stewart(seb, tile);
            
            seb.TEMP.r_a = 1./(kappa^2.*u_star).*(log((z-d)./z0v)- psi_M_CLM5(seb, (z-d)./Lstar, z0v./Lstar)).*(log((z-d)./z0v)- psi_H_CLM5(seb, (z-d)./Lstar, z0v./Lstar)); % aerodynamic resistance to heat/water vapor fransper canopy air - atmosphere. From original CG3 publication
            seb.TEMP.r_b = 1./Cv*(u_star./d_leaf).^(-.5); % Eq. 5.122 leaf boundary layer resistance
            seb.TEMP.r_a_prime = 1/(Cs.*u_star); % Eq. 5.116 aerodynamic resistance to heat/water vapor transfer soil - canopy air
            seb.TEMP.r_soil = ground_resistance_evap(seb.IA_NEXT, tile); % resistance to water vapor flux within the soil matrix
            if seb.STATVAR.LAI > 0 % Transpiring leaves present
                seb.TEMP.r_canopy = 1./(seb.TEMP.C_leaf.*k_s);
            else % No leaves -> no transpiration
                seb.TEMP.r_canopy = Inf;
            end
            seb.TEMP.Ev_pot = -seb.TEMP.rho_atm.*(q_s - qs_Tv)./seb.TEMP.r_b;
            if seb.TEMP.Ev_pot <= 0 % Condensation -> on whole leaf + stem area, NO transpiration
                seb.TEMP.r_total = seb.TEMP.r_b./(2*(L+S));
            elseif seb.TEMP.Ev_pot > 0 % evaporation + transpiration
                seb.TEMP.r_total = seb.TEMP.r_b./(2*(L+S)) ./ ( f_wet + f_dry.*seb.TEMP.r_b./(seb.TEMP.r_b + seb.TEMP.r_canopy) ) ;  
            end
        end
        
        function beta = soil_resistance_beta(seb, forcing)
            error
            uz = forcing.TEMP.wind;
            kappa = seb.CONST.kappa;
            field_capacity = seb.STATVAR.field_capacity(1);
            q = forcing.TEMP.q;
            p = forcing.TEMP.p;
            T_surf = seb.STATVAR.T(1)+seb.CONST.Tmfw;
            
            water_fraction = seb.STATVAR.water(1) ./ seb.STATVAR.waterIce(1);
            ice_fraction = seb.STATVAR.ice(1) ./ seb.STATVAR.waterIce(1);
            saturation_factor_surf= exp(seb.STATVAR.waterPotential(1) .* seb.CONST.g ./ ((seb.CONST.R./ seb.CONST.molar_mass_w) .*T_surf));
            e_sat_surface = water_fraction .* satPresWater(seb, T_surf) + ice_fraction .* satPresIce(seb, T_surf);
            q_surf = .622.*e_sat_surface.*saturation_factor_surf./p;
            
            water_cont_surf = seb.STATVAR.water(1) ./ (seb.STATVAR.layerThick(1) .* seb.STATVAR.area(1));
            
            reduce =  water_cont_surf >= field_capacity | q > q_surf;
            beta = 1 + double(reduce) .* (-1 +  0.25 .* (1-(cos(pi .* water_cont_surf ./ field_capacity))).^2);
        end
        
        function C_leaf = leaf_conductance_Stewart(seb, tile)
            forcing = tile.FORCING;
            C_leaf_max = seb.PARA.C_leaf_max;
            Sin = forcing.TEMP.Sin;
            p = forcing.TEMP.p;
            q_s = seb.STATVAR.q_s;
            Tmfw = seb.CONST.Tmfw;
            rho_air = seb.TEMP.rho_atm;
            T_leaf = max(0,min(40,seb.STATVAR.T(1))); % only valid for 0 <= T_leaf <= 40
            
            e_sat = satPresWater(seb,seb.STATVAR.Ts+Tmfw); 
            rho_v_saturated = e_sat./p.*rho_air;
            rho_v_actual = q_s.*rho_air;
            Delta_rho_v = rho_v_saturated - rho_v_actual;
            
            f_Sin = 1.1046.*Sin./(Sin+104.4);
            f_Sin = max(0,min(f_Sin,1)); % limit to [0 1]
            f_rho_v = 1 - 66.6.*Delta_rho_v;  
            f_rho_v = max(0.233,min(f_rho_v,1));
            f_T_leaf = T_leaf.*(40-T_leaf).^1.18 ./ 691;
            f_T_leaf = max(0,min(f_T_leaf,1));
            f_psi_soil = get_soil_moisture_stress(seb.IA_GROUND);
            f_psi_soil = max(0,min(f_psi_soil,1));
            
            C_leaf = C_leaf_max.*f_Sin.*f_rho_v.*f_T_leaf.*f_psi_soil;
            
            seb.TEMP.Delta_rho_v = Delta_rho_v;
            seb.TEMP.f_Sin = f_Sin;
            seb.TEMP.f_rho_v = f_rho_v;
            seb.TEMP.f_T_leaf = f_T_leaf;
            seb.TEMP.f_psi_soil = f_psi_soil;
        end
        
% -------------- Support functions ---------------------------------------        

        
        function z0 = get_z0_surface(seb)
            z0 = seb.PARA.z0;
        end
        
        function rho = rho_air(seb, p, T) %air density [kg m^(-3)]
            rho = p./(287.058.*T);
        end
        
        function rho = rho_air_moist(seb, tile)
            forcing = tile.FORCING;
            R_da = seb.CONST.R_spec; % gas constant for dry air = 287.0423 (CLM5 documentation)
            p = forcing.TEMP.p;
            Tz= forcing.TEMP.Tair+forcing.CONST.Tmfw;
            q = forcing.TEMP.q;
            
            e = q.*p./(0.622+0.378*q); % atmospheric vapor pressure
            rho = (p - 0.378*e)./(R_da * Tz);
        end
        
        function L_w = latent_heat_evaporation(seb, T) %specific latent heat of evaporation of water [J/kg]
            L_w = 1e3.*(2500.8 - 2.36.*(T - 273.15));
        end
        
        function L_i = latent_heat_sublimation(seb, T_forcing) %1e3.*2834.1; %latent heat of sublimation, constant [J/kg]
            L_i = seb.CONST.L_s;
        end
        
        function p = satPresWater(seb, T) %saturation pressure water, Magnus formula
            p= 6.112 .* 100 .* exp(17.62.*(T-273.15)./(243.12-273.15+T)); % Removed the factor .622 to give actual vapor pressure, RBZ Nov. 2021
        end
        
        function p = satPresIce(seb, T) %saturation pressure ice, Magnus formula
            p= 6.112.* 100.* exp(22.46.*(T-273.15)./(272.61-273.15+T)); % Removed the factor .622 to give actual vapor pressure, RBZ Nov. 2021
        end
        
%--------------- penetration of LW --------------------------------------
        function [seb, Lout] = penetrate_LW_no_transmission(seb, Lin)
            % Lin is in W, not W/m2!
            Lout = (1-seb.PARA.epsilon) .* Lin + seb.PARA.epsilon .* seb.CONST.sigma .* (seb.STATVAR.T(1)+ seb.CONST.Tmfw).^4 .*seb.STATVAR.area(1);
            seb.TEMP.L_abs = Lin - Lout;
            seb.STATVAR.Lin = Lin./seb.STATVAR.area(1);
            seb.STATVAR.Lout = Lout./seb.STATVAR.area(1);
            seb.TEMP.d_energy(1,1) = seb.TEMP.d_energy(1,1) + seb.TEMP.L_abs;
        end
        
        function [seb, Lout] = penetrate_LW_no_transmission_GROUND_snow(seb, Lin)
            % For GROUND that is below a snow in CHILD phase
            % Lin is in W, not W/m2!
            ground_area = seb.STATVAR.area(1) - seb.CHILD.STATVAR.area;
            Lout = (1-seb.PARA.epsilon) .* Lin + seb.PARA.epsilon .* seb.CONST.sigma .* (seb.STATVAR.T(1)+ seb.CONST.Tmfw).^4 .*ground_area;
            seb.TEMP.L_abs = Lin - Lout;
            seb.STATVAR.Lin = Lin./seb.STATVAR.area(1);
            seb.STATVAR.Lout = Lout./seb.STATVAR.area(1);
            seb.TEMP.d_energy(1,1) = seb.TEMP.d_energy(1,1) + seb.TEMP.L_abs;
        end
        
        function [seb, Lout] = penetrate_LW_simpleShading(seb, Lin)
            % Lin is in W, not W/m2!!!!
            fractional_canopy_cover = seb.PARA.fractional_canopy_cover;
            canopy_transmissivity = seb.PARA.canopy_transmissivity;
            canopy_emissivity = seb.PARA.canopy_emissivity;
            Tmfw = seb.CONST.Tmfw;
            sigma = seb.CONST.sigma;
            
            % transmitted (directly and through canopy) and emitted LW downwards
            L_down = fractional_canopy_cover * (canopy_transmissivity * Lin + ...
                (1 - canopy_transmissivity) * canopy_emissivity * sigma * (seb.STATVAR.T + Tmfw)^4.*seb.STATVAR.area) + ...
                (1 - fractional_canopy_cover) * Lin;
            
            % penetrate LW to class below
            [seb.NEXT, L_up] = penetrate_LW(seb.NEXT, L_down);
            
            % transmitted/emitted/reflected upward
            Lout = L_up * (fractional_canopy_cover * canopy_transmissivity + (1 - fractional_canopy_cover)) ...
                + fractional_canopy_cover * (1 - canopy_transmissivity) * canopy_emissivity * sigma * (seb.STATVAR.T + Tmfw)^4.*seb.STATVAR.area ...
                + Lin * (1 - canopy_emissivity) * fractional_canopy_cover + L_up.*(1-canopy_emissivity)*fractional_canopy_cover;
            % Last term is upward LW reflected by the base of the canopy, which is added to Lout because multiple backstattering is not considered.
            
            seb.TEMP.L_abs = (Lin + L_up.*canopy_emissivity*fractional_canopy_cover - L_down - Lout); % L_up which is reflected by the canopy is discarded
            seb.TEMP.d_energy(1) = seb.TEMP.d_energy(1) + seb.TEMP.L_abs;
            seb.TEMP.L_down = L_down./seb.STATVAR.area(1);
            seb.TEMP.L_up = L_up./seb.STATVAR.area(1);
            seb.STATVAR.Lin = Lin./seb.STATVAR.area(1);
            seb.STATVAR.Lout = Lout./seb.STATVAR.area(1);
        end
        
        function [seb, Lout] = penetrate_LW_CLM5(seb,Lin)
            % Lin is in W, not W/m2!!!!
            canopy_emissivity = seb.STATVAR.emissivity;
            Tmfw = seb.CONST.Tmfw;
            sigma = seb.CONST.sigma;
            
            % Transmitted and emitted downward below canopy
            L_down = (1 - canopy_emissivity).*Lin + canopy_emissivity.*sigma.*(seb.STATVAR.T + Tmfw)^4.*seb.STATVAR.area;
            
            % penetrate LW to class below
            [seb.NEXT, L_up] = penetrate_LW(seb.NEXT, L_down);
            
            % transmitted/emitted/reflected upward
            Lout = (1- canopy_emissivity).*L_up + canopy_emissivity.*sigma.*(seb.STATVAR.T + Tmfw)^4.*seb.STATVAR.area;
            
            seb.TEMP.L_abs = Lin + L_up - Lout - L_down;
            seb.TEMP.d_energy(1) = seb.TEMP.d_energy(1) + seb.TEMP.L_abs;
            seb.TEMP.L_down = L_down./seb.STATVAR.area(1);
            seb.TEMP.L_up = L_up./seb.STATVAR.area(1); 
            seb.STATVAR.Lin = Lin./seb.STATVAR.area(1);
            seb.STATVAR.Lout = Lout./seb.STATVAR.area(1);
        end
        
% -------------- Penetration of short-wave radiation ---------------------
        function [seb, Sout] = penetrate_SW_no_transmission(seb, Sin) % called recursively by surfaceEnergyBalance-function of TOP_CLASS, "hard" surface absorbing all SW-radiation in 1st cell
            seb.TEMP.d_energy(1,1) = seb.TEMP.d_energy(1,1) + (1 - seb.PARA.albedo) .* sum(Sin); %in [W], not [W/m2]
            Sout = seb.PARA.albedo .* Sin;
            
            if seb.TEMP.SW_split == 0 % First time this function is called
                seb.STATVAR.Sin = sum(Sin)./seb.STATVAR.area(1);
                seb.STATVAR.Sout = sum(Sout)./seb.STATVAR.area(1);
                seb.TEMP.S_abs = (1 - seb.PARA.albedo) .* sum(Sin);
            else % SW is penetrated below child, and these variables are already populated
                seb.STATVAR.Sin = seb.STATVAR.Sin + sum(Sin)./seb.STATVAR.area(1);
                seb.STATVAR.Sout = seb.STATVAR.Sout + sum(Sout)./seb.STATVAR.area(1);
                seb.TEMP.S_abs = seb.TEMP.S_abs + (1 - seb.PARA.albedo) .* sum(Sin);
            end
            
        end
        
        function [seb, Sout] = penetrate_SW_CLM5(seb, Sin)
            % Documentation: https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/tech_note/Surface_Albedos/CLM50_Tech_Note_Surface_Albedos.html
            spectral_weights = [seb.PARA.nir_fraction 1-seb.PARA.nir_fraction ];
            sun_angle = seb.TEMP.sun_angle;
            zenith = min(89,90 - sun_angle); % sun angles close to 90 produce NaNs
            Sin(zenith == 89) = 0;
            Sin_dir = Sin.*seb.TEMP.Sin_dir_fraction;
            Sin_dif = Sin - Sin_dir;
            
            L = seb.STATVAR.LAI; % Leaf area index
            S = seb.STATVAR.SAI; % Stem area index
            alpha_leaf = [seb.PARA.alpha_leaf_nir seb.PARA.alpha_leaf_vis]; % leaf reflectances
            alpha_stem = [seb.PARA.alpha_stem_nir seb.PARA.alpha_stem_vis]; % stem reflectances
            tau_leaf = [seb.PARA.tau_leaf_nir seb.PARA.tau_leaf_vis]; % leaf transmittances
            tau_stem = [seb.PARA.tau_stem_nir seb.PARA.tau_stem_vis]; % stem transmittances
            Khi_L = seb.PARA.Khi_L; % departure of leaf angles from a random distribution
            alpha_g = get_albedo(seb.NEXT); % ground albedo, direct and diffuse
            
            w_leaf = L/(L+S); % leaf weighting
            w_stem = S/(L+S); % stem weighting
            alpha = alpha_leaf.*w_leaf + alpha_stem.*w_stem; % canopy reflectance, Eq. 3.11
            tau = tau_leaf.*w_leaf + tau_stem.*w_stem; % canopy transmittance, Eq.
            phi1 = .5 - .633.*Khi_L - .33.*Khi_L.^2; % for -.4 <= Khi_L <=.6
            phi2 = .877 .* (1-2.*phi1);
            my = cosd(zenith);
            G = phi1 + phi2.*my; % Relative projecter area of canopy (leaf and stem), Eq. 3.3
            K = G./my; % optical depth
            my_bar = 1./phi2*( 1- phi1./phi2 * log( (phi1+phi2)./phi1) ); % average inverse diffuse optical depth per unit leaf and stem area, Eq. 3.4
            omega = alpha + tau; % scattering coefficient
            cos_theta = (1+Khi_L)./2; % theta = mean leaf inclination angle relative to the horizontal plane, Eq. 3.14
            omega_beta = .5.*(alpha + tau + (alpha - tau).*cos_theta.^2); % upscatter for diffuse radiation, Eq. 3.13
            as = omega./2 .* G./max(my.*phi2+G,1e-6) .* (1 - my.*phi1./max(my.*phi2+G,1e-6) .* log((my.*phi1+max(my.*phi2+G,1e-6))./(my.*phi1)) ); % single scatter albedo, Eq. 3.16
            omega_beta_0 = (1+my_bar.*K)./(my_bar.*K).*as; % upscatter for direct beam radiation, Eq. 3.15
            beta_0 = omega_beta_0./omega;
            
            b = 1 - omega + omega_beta; % Eq. 3.31
            c = omega_beta; % Eq. 3.32
            d = omega.*my_bar.*K.*beta_0; % Eq.3.33
            f = omega.*my_bar.*K.*(1-beta_0); % Eq. 3.34
            h = sqrt(b.^2-c.^2)./my_bar; % Eq. 3.35
            sigma = (my_bar.*K).^2 + c.^2 - b.^2; % Eq. 3.36
            u1 = b - c./alpha_g; % Eq. 3.37
            u2 = b - c.*alpha_g; % Eq. 3.38
            u3 = f + c.*alpha_g; % Eq. 3.39
            s1 = exp(-min(h.*(L+S),40)); % Eq. 3.40
            s2 = exp(-min(K*(L+S),40)); % Eq. 3.41
            p1 = b + my_bar.*h; % Eq. 3.42
            p2 = b - my_bar.*h; % Eq. 3.43
            p3 = b + my_bar*K; % Eq. 3.44
            p4 = b - my_bar*K; % Eq. 3.45
            d1 = p1.*(u1-my_bar*h)./s1 - p2.*(u1+my_bar*h).*s1; % Eq. 3.46
            d2 = (u2 + my_bar.*h)./s1 - (u2 - my_bar*h).*s1; % Eq. 3.47
            h1 = -d.*p4 - c.*f; % Eq. 3.48
            h2 = 1./d1.*( (d-h1./sigma.*p3).*(u1-my_bar.*h)./s1 - p2.*(d-c-h1./sigma.*(u1+my_bar.*K)).*s2 ); % Eq. 3.50
            h3 = -1./d1.*( (d-h1./sigma.*p3).*(u1+my_bar.*h).*s1 - p1.*(d-c-h1./sigma.*(u1+my_bar.*K)).*s2 ); % Eq. 3.51
            h4 = -f.*p3 - c.*d; % Eq. 3.51
            h5 = -1./d2.*( h4.*(u2+my_bar.*h)./(sigma.*s1) + (u3-h4./sigma.*(u2-my_bar*K)).*s2 ); % Eq. 3.52
            h6 = 1./d2.*( h4./sigma.*(u2-my_bar.*h).*s1 + (u3-h4./sigma.*(u2 - my_bar.*K)).*s2 ); % Eq. 3.53
            h7 = c.*(u1-my_bar.*h)./(d1.*s1); % Eq. 3.54
            h8 = -c.*(u1+my_bar*h).*s1./d1; % Eq. 3.55
            h9 = (u2 + my_bar.*h)./(d2.*s1); % Eq. 3.56
            h10 = -s1.*(u2 - my_bar*h)./d2; % Eq. 3.57
            
            % Downwelling shortwave fluxes
            I_down_from_dir = (h4./sigma.*exp(-K*(L+S)) + h5.*s1 + h6./s1)*spectral_weights'; % Eq. 3.19
            I_down_from_dif = (h9.*s1 + h10./s1)*spectral_weights'; % Eq. 3.20
            I_transmitted = exp(-K*(L+S));
            
            S_down =  Sin_dir*(I_down_from_dir + I_transmitted) + Sin_dif*I_down_from_dif;
            
            [seb.NEXT, S_up] = penetrate_SW(seb.NEXT, S_down);
            S_up = sum(S_up); % snow crocus resolves SW spectraly
            
            % Upwelling shortwave fluxes
            I_out_from_dir = (h1./sigma + h2 + h3);
            I_out_from_dir = I_out_from_dir*spectral_weights'; % Eq. 3.17
            I_out_from_dif = (h7 + h8)*spectral_weights'; % Eq. 3.18
            
            Sout = Sin_dif.*I_out_from_dif + Sin_dir.*I_out_from_dir;
            
            seb.TEMP.S_abs = Sin + S_up - S_down - Sout;
            seb.TEMP.d_energy(1) = seb.TEMP.d_energy(1) + seb.TEMP.S_abs;
            seb.TEMP.S_down = S_down./seb.STATVAR.area(1);
            seb.TEMP.S_up = S_up./seb.STATVAR.area(1);
            seb.STATVAR.Sin = Sin./seb.STATVAR.area(1);
            seb.STATVAR.Sout = Sout./seb.STATVAR.area(1);
        end
        
        function [seb, Sout] = penetrate_SW_simpleShading(seb, Sin)
            albedo = seb.PARA.canopy_albedo;
            fractional_canopy_cover = seb.PARA.fractional_canopy_cover;
            canopy_height = seb.PARA.canopy_height;
            canopy_extinction_coefficient = seb.PARA.canopy_extinction_coefficient;
            canopy_transmissivity = seb.PARA.canopy_transmissivity;
            sun_angle = seb.TEMP.sun_angle;
            
            % split Sin into direct and diffuse radiation
            S_down_dir = Sin .* seb.TEMP.Sin_dir_fraction;
            S_down_dif = Sin .* (1 - seb.TEMP.Sin_dir_fraction);
            
            % Transmitted SW down
            S_down_dir = fractional_canopy_cover * (S_down_dir * ...
                exp(-canopy_extinction_coefficient * canopy_height /max(1e-12,sind(sun_angle)))) ...
                + (1 - fractional_canopy_cover) * S_down_dir;
            S_down_dif = fractional_canopy_cover * (S_down_dif * canopy_transmissivity) + ...
                (1 - fractional_canopy_cover) * S_down_dif;
            S_down = S_down_dir + S_down_dif;
            
            % penetrate SW to class below
            [seb.NEXT, S_up] = penetrate_SW(seb.NEXT, S_down);
            
            % transmitted & reflected upward
            % original parameterization only includes reduction of Sin, and has no albedo /Sout calculation!!
            Sout = fractional_canopy_cover .* (S_up .* canopy_transmissivity) + (1 - fractional_canopy_cover) .* S_up ...
                + albedo .* fractional_canopy_cover .* Sin + S_up*albedo*fractional_canopy_cover;
            % Last term is backstattered SW reflected by the canopy, and is added to Sout because multiple backstatter is not considered!
            
            seb.TEMP.S_abs = Sin + S_up - S_down - Sout;
            seb.TEMP.d_energy(1) = seb.TEMP.d_energy(1) + seb.TEMP.S_abs;
            seb.TEMP.S_backscat_discarded = albedo*fractional_canopy_cover*S_up;
            seb.TEMP.S_down = S_down./seb.STATVAR.area(1);
            seb.TEMP.S_up = S_up./seb.STATVAR.area(1);
            seb.STATVAR.Sin = Sin./seb.STATVAR.area(1);
            seb.STATVAR.Sout = Sout./seb.STATVAR.area(1);
        end
        
        function [seb, Sout] = penetrate_SW_transmission_bulk(seb, Sin)  %used with fixed albedo and SW extinction coefficient
            %S_up and S_down can in principle be spectrally resolved when provided as
            %row array, using e.g. spectral_ranges = [0.71 0.21 0.08]; 
            %SW extinction is assumed constant throughout class 
            cut_off = 0.1.* seb.STATVAR.area(1,1); %[W/m2], radiation is not penetrated further if cutoff is reached
            
            Sout = seb.STATVAR.albedo .* Sin;
            S_down = (1 - seb.STATVAR.albedo) .* Sin;
            
            i=1;
            while i<=size(seb.STATVAR.layerThick,1) && sum(S_down) >= cut_off
                reduction = exp(-seb.PARA.SW_extinction .* seb.STATVAR.layerThick(i,1));
                seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_down .* (1-reduction));
                S_down = reduction .* S_down;
                i = i+1;
            end
            
            if sum(S_down) < cut_off  %all radiation absorbed, only Sout goes out
                i = min(i, size(seb.STATVAR.layerThick,1));
                seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_down);
            else %end of class reached, radiation penetrated on to next class
                [seb.NEXT, S_up] = penetrate_SW(seb.NEXT, S_down); %call mandatory function recursively for following class
                S_up2 = S_up;
                %NOTE: no check performed that NEXT is not Bottom, must always be used with class below!
                i = size(seb.STATVAR.layerThick,1); %penetrate bottom up
                while i>=1 && sum(S_up2) >= cut_off
                    reduction = exp(-seb.PARA.SW_extinction .* seb.STATVAR.layerThick(i,1));
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_up2 .* (1-reduction));
                    S_up2 = reduction .* S_up2;
                    i = i-1;
                end
                if sum(S_up2) < cut_off  %all radiation absorbed, only S_up goes out
                    i = max(i, 1);
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_up2); % .* seb.STATVAR.area(i,1);
                else
                    Sout = Sout + S_up2;  %add the uppwelling SW radiation to reflected, multiple reflections not accounted for!!
                end
                seb.TEMP.S_down = S_down./seb.STATVAR.area(1);
                seb.TEMP.S_up = S_up./seb.STATVAR.area(1);
            end
            seb.STATVAR.Sin = Sin./seb.STATVAR.area(1);
            seb.STATVAR.Sout = Sout./seb.STATVAR.area(1);
        end
        
        function [seb, Sout] = penetrate_SW_transmission_spectral(seb, Sin)  %used with variable albedo and SW exticntion coefficient
            %S_up and S_down are spectrally resolved when provided as
            %row array, using e.g. spectral_ranges = [0.71 0.21 0.08]; 
            %SW extinction is assumed constant throughout class
            
            Sin = sum(Sin) .* seb.PARA.spectral_ranges; % Moved this here to allow for use with abovelying classes, ie. vegetation. RBZ 25022022
            
            cut_off = 0.1*seb.STATVAR.area(1); % = 0.1 W/m2, radiation is not penetrated further if cutoff is reached
            
            Sout = seb.TEMP.spectral_albedo .* Sin;
            S_down = (1 - seb.TEMP.spectral_albedo) .* Sin;
            
            i=1;
            while i<=size(seb.STATVAR.layerThick,1) && sum(S_down) >= cut_off
                reduction = exp(-seb.TEMP.SW_extinction(i,:) .* seb.STATVAR.layerThick(i,1));
                seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_down .* (1-reduction));
                S_down = reduction .* S_down;
                i = i+1;
            end
            
            if sum(S_down) < cut_off  %all radiation absorbed, only S_up goes out
                i = min(i, size(seb.STATVAR.layerThick,1));
                seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_down);
            else %end of class reached, radiation penetrated on to next class
                [seb.NEXT, S_up] = penetrate_SW(seb.NEXT, S_down); %call mandatory function recursively for following class
                S_up2 = S_up;
                %NOTE: no check performed that NEXT is not Bottom, must always be used with class below!
                i = size(seb.STATVAR.layerThick,1); %penetrate bottom up
                while i>=1 && sum(S_up2) >= cut_off
                    reduction = exp(-seb.TEMP.SW_extinction(i,:) .* seb.STATVAR.layerThick(i,1));
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_up2 .* (1-reduction));
                    S_up2 = reduction .* S_up2;
                    i = i-1;
                end
                if sum(S_up2) < cut_off  %all radiation absorbed, only S_up goes out
                    i = max(i, 1);
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_up2); % .* seb.STATVAR.area(i,1);
                else
                    Sout = Sout + S_up2;  %add the uppwelling SW radiation to reflected, multiple reflections not accounted for!!
                end
                seb.TEMP.S_down = sum(S_down)./seb.STATVAR.area(1);
                seb.TEMP.S_up = sum(S_up)./seb.STATVAR.area(1);
            end
            seb.STATVAR.Sin = sum(Sin)./seb.STATVAR.area(1);
            seb.STATVAR.Sout = sum(Sout)./seb.STATVAR.area(1);
        end
        
% -------------- Water fluxes from evapotranspiration --------------------
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

