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
            
            Tz=Tz+forcing.CONST.Tmfw;
            TForcing=TForcing+forcing.CONST.Tmfw;
            rho = rho_air(seb, p, Tz);
            
            Q_h  = -rho.*cp.*kappa.* uz.*kappa./(log(z./z0)- psi_M(seb, z./Lstar, z0./Lstar)) .* (Tz-TForcing)./(log(z./z0)- psi_H(seb, z./Lstar, z0./Lstar));
            
        end
        
        function Q_h = Q_h_CLM5(seb, forcing)
            Cp = seb.CONST.cp; % Specific heat capacity of air
            L = seb.PARA.LAI; % Leaf area index
            S = seb.PARA.SAI; % Stem area index
            Tv = seb.STATVAR.T(1)+forcing.CONST.Tmfw; % leaf temperature
            Tg = seb.NEXT.STATVAR.T(1)+forcing.CONST.Tmfw; % ground temperature
            z = forcing.PARA.airT_height;
            z0v = seb.PARA.z0; % roughness length of vegetation
            z0g = seb.NEXT.PARA.z0; % Roughness lenght of ground
            %R_z0 = seb.PARA.R_z0; % not needed when zo
            Tz = forcing.TEMP.Tair+forcing.CONST.Tmfw; % air temperature in Kelvin
            Gamma = seb.CONST.Gamma_dry; % negative of dry adiabatic lapse rate
            kappa = seb.CONST.kappa;
            u_star = seb.STATVAR.u_star;
            Lstar = seb.STATVAR.Lstar;
            Cv = seb.PARA.Cv; % Turbulent transfer coefficient canopy surface - canopy air
            d_leaf = seb.PARA.d_leaf;
            z_top = sum(seb.STATVAR.layerThick);
            ypsilon = seb.CONST.ypsilon; % kinematic viscosity of air
            Cs_dense = seb.PARA.Cs_dense;
            z = z + z_top;
            
            Tz_pot = Tz+Gamma.*z; % Eq. 5.7 potential temperature
            Cs_bare = kappa./.13.*(z0g.*u_star./ypsilon).^(-.45); % Eq.5.121 bare soil turbulent transfer coefficient
            W = exp(-(L+S));% Eq. 5.119
            Cs = Cs_bare.*W + Cs_dense.*(1-W); % Eq. 5.118 turbulent transfer coefficient between soil and canopy air
            r_ah = 1./(kappa^2.*u_star).*(log(z./z0v)- psi_M(seb, z./Lstar, z0v./Lstar)).*(log(z./z0v)- psi_H(seb, z./Lstar, z0v./Lstar)); % From original CG3 publication
            r_b = 1./Cv*(u_star./d_leaf).^(-.5); % Eq. 5.122 leaf boundary layer resistance
            r_ah2 = 1/(Cs.*u_star); % Eq. 5.116 aerodynamic resistance to heat transfer soil - canopy air
            Ts = (Tz_pot./r_ah + Tg./r_ah2 + Tv.*(L+S)./r_b)./(1./r_ah + 1/r_ah2 + (L+S)./r_b); % Eq. 5.93 canopy air temperature
            rho = rho_air_moist(seb,forcing); % Eq. 5.8-9, moist air density
            Q_h = -rho.*Cp.*(Ts-Tv).*(L+S)./r_b; % Eq. 5.88
            
            Q_h_ground  = -rho.*Cp.*(Ts-Tg)./r_ah2;
            seb.NEXT.TEMP.d_energy(1) = seb.NEXT.TEMP.d_energy(1) - Q_h_ground;

        end
        
        function rho = rho_air(seb, p, T) %air density [kg m^(-3)]
            rho = p./(287.058.*T);
        end
        
        function rho = rho_air_moist(seb, forcing)
            R_da = seb.CONST.R_spec; % gas constant for dry air = 287.0423 (CLM5 documentation)
            p = forcing.TEMP.p;
            Tz= forcing.TEMP.Tair+forcing.CONST.Tmfw;
            q = forcing.TEMP.q;
            
            e = q.*p./(0.622+0.378*q); % atmospheric vapor pressure
            rho = (p - 0.378*e)./(R_da * Tz);
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
        
        % ---------- penetration of LW ----------------------------------
        function [seb, Lout] = penetrate_LW_no_transmission(seb, Lin)
            % Lin is in W, not W/m2!
            Lout = (1-seb.PARA.epsilon) .* Lin + seb.PARA.epsilon .* seb.CONST.sigma .* (seb.STATVAR.T(1)+ seb.CONST.Tmfw).^4 .*seb.STATVAR.area(1);
            seb.TEMP.L_abs = Lin - Lout;
            seb.STATVAR.Lin = Lin;
            seb.STATVAR.Lout = Lout;
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
            
            seb.TEMP.L_abs = (Lin + L_up.*canopy_emissivity*fractional_canopy_cover - L_down - Lout)./seb.STATVAR.area; % L_up which is reflected by the canopy is discarded
            seb.TEMP.d_energy(1) = seb.TEMP.d_energy(1) + seb.TEMP.L_abs;
            seb.TEMP.L_down = L_down./seb.STATVAR.area;
            seb.TEMP.L_up = L_up./seb.STATVAR.area;
        end
        
        function [seb, Lout] = penetrate_LW_CLM5(seb,Lin)
            % Lin is in W, not W/m2!!!!
            canopy_emissivity = seb.PARA.canopy_emissivity;
            Tmfw = seb.CONST.Tmfw;
            sigma = seb.CONST.sigma;
            
            % Transmitted and emitted downward below canopy
            L_down = (1 - canopy_emissivity).*Lin + canopy_emissivity.*sigma.*(seb.STATVAR.T + Tmfw)^4.*seb.STATVAR.area;
            
            % penetrate LW to class below
            [seb.NEXT, L_up] = penetrate_LW(seb.NEXT, L_down);
            
            % transmitted/emitted/reflected upward
            Lout = (1- canopy_emissivity).*L_up + canopy_emissivity.*sigma.*(seb.STATVAR.T + Tmfw)^4.*seb.STATVAR.area;
            
            seb.TEMP.L_abs = (Lin + L_up - Lout - L_down)./seb.STATVAR.area;
            seb.TEMP.d_energy(1) = seb.TEMP.d_energy(1) + seb.TEMP.L_abs;
            seb.TEMP.L_down = L_down./seb.STATVAR.area;
            seb.TEMP.L_up = L_up./seb.STATVAR.area; % STATVARS are assigned in canopy_energy_balance
        end
        
        % --------- penetration of short-wave radiation ----------------
        function [seb, Sout] = penetrate_SW_no_transmission(seb, Sin) % called recursively by surfaceEnergyBalance-function of TOP_CLASS, "hard" surface absorbing all SW-radiation in 1st cell
            seb.TEMP.d_energy(1,1) = seb.TEMP.d_energy(1,1) + (1 - seb.PARA.albedo) .* sum(Sin); %in [W], not [W/m2]
            Sout = seb.PARA.albedo .* Sin;
            seb.STATVAR.Sin = Sin;
            seb.STATVAR.Sout = Sout;
            seb.TEMP.S_abs = (1 - seb.PARA.albedo) .* sum(Sin);
        end
        
        function [seb, Sout] = penetrate_SW_CLM5(seb, Sin)
            % Documentation: https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/tech_note/Surface_Albedos/CLM50_Tech_Note_Surface_Albedos.html
            spectral_weights = [seb.PARA.nir_fraction 1-seb.PARA.nir_fraction ];
            zenith = min(89,90 - seb.TEMP.sun_angle); % sun angles close to 90 produce NaNs
            Sin(zenith == 89) = 0;
            Sin_dir = Sin.*seb.TEMP.Sin_dir_fraction;
            Sin_dif = Sin - Sin_dir;
            
            L = seb.PARA.LAI; % Leaf area index
            S = seb.PARA.SAI; % Stem area index
            alpha_leaf = [seb.PARA.alpha_leaf_nir seb.PARA.alpha_leaf_vis]; % leaf reflectances
            alpha_stem = [seb.PARA.alpha_stem_nir seb.PARA.alpha_leaf_vis]; % stem reflectances
            tau_leaf = [seb.PARA.tau_leaf_nir seb.PARA.tau_leaf_vis]; % leaf transmittances
            tau_stem = [seb.PARA.tau_stem_nir seb.PARA.tau_stem_vis]; % stem transmittances
            Khi_L = seb.PARA.Khi_L; % departure of leaf angles from a random distribution
            alpha_g = seb.NEXT.PARA.albedo; % ground albedo, direct and diffuse
            
            w_leaf = L/(L+S); % leaf weighting
            w_stem = S/(L+S); % stem weighting
            alpha = alpha_leaf.*w_leaf + alpha_stem.*w_stem; % canopy reflectance
            tau = tau_leaf.*w_leaf + tau_stem.*w_stem; % canopy transmittance
            phi1 = .5 - .633.*Khi_L - .33.*Khi_L.^2;
            phi2 = .877 .* (1-2.*phi1);
            my = cosd(zenith);
            G = phi1 + phi2.*my; % Relative projecter area of canopy (leaf and stem)
            K = G./my; % optical depth
            my_bar = 1./phi2*( 1- phi1./phi2 * log( (phi1+phi2)./phi1) ); % average inverse diffuse optical depth per unit leaf and stem area
            omega = alpha + tau; % scattering coefficient
            cos_theta = (1+Khi_L)./2; % theta = mean leaf inclination angle relative to the horizontal plane
            omega_beta = .5.*(alpha + tau + (alpha - tau).*cos_theta.^2); % upscatter for diffuse radiation
            as = omega./2 .* G./max(my.*phi2+G,1e-6) .* (1 - my.*phi1./max(my.*phi2+G,1e-6) .* log((my.*phi1+max(my.*phi2+G,1e-6))./(my.*phi1)) );
            omega_beta_0 = (1+my_bar.*K)./(my_bar.*K).*as;
            beta_0 = omega_beta_0./omega;
            
            b = 1 - omega + omega_beta;
            c = omega_beta;
            d = omega.*my_bar.*K.*beta_0;
            f = omega.*my_bar.*K.*(1-beta_0);
            h = sqrt(b.^2+c.^2)./my_bar;
            sigma = (my_bar.*K).^2 + c.^2 + b.^2;
            u1 = b - c./alpha_g;
            u2 = b - c.*alpha_g;
            u3 = f + c.*alpha_g;
            s1 = exp(-min(h.*(L+S),40));
            s2 = exp(-min(K*(L+S),40));
            p1 = b + my_bar*h;
            p2 = b - my_bar*h;
            p3 = b + my_bar*K;
            p4 = b - my_bar*K;
            d1 = p1.*(u1-my_bar*h)./s1 - p2.*(u1+my_bar*h).*s1; %
            d2 = (u2 + my_bar.*h)./s1 - (u2 - my_bar*h).*s1;%
            h1 = -d.*p4 - c.*d;
            h2 = 1./d1.*( (d-h1./sigma.*p3).*(u1-my_bar.*h)./s1 - p2.*(d-c-h1/sigma.*(u1+my_bar.*K)).*s2 );%
            h3 = -1./d1.*( (d-h1./sigma.*p3).*(u1+my_bar.*h).*s1 - p1.*(d-c-h1/sigma.*(u1+my_bar.*K)).*s2 );%
            h4 = -f.*p3 - c.*d;
            h5 = -1./d2.*( h4.*(u2+my_bar.*h)./(sigma.*s1) + (u3-h4./sigma.*(u2-my_bar*K)).*s2 );%
            h6 = 1./d2.*( h4./sigma.*(u2-my_bar.*h).*s1 + (u3-h4./sigma.*(u2 - my_bar.*K)).*s2 );%
            h7 = c.*(u1-my_bar.*h)./(d1.*s1);%
            h8 = -c.*(u1+my_bar*h).*s1./d1;%
            h9 = (u2 + my_bar.*h)/(d2.*s1);%
            h10 = -s1.*(u2 - my_bar*h)./d2;%
            
            % Downwelling shortwave fluxes
            I_down_from_dir = (h4./sigma.*exp(-K*(L+S)) + h5.*s1 + h6./s1)*spectral_weights'; %
            I_down_from_dif = (h9.*s1 + h10./s1)*spectral_weights'; %
            I_transmitted = exp(-K*(L+S));
            
            S_down =  Sin_dir*(I_down_from_dir + I_transmitted) + Sin_dif*I_down_from_dif;
            
            [seb.NEXT, S_up] = penetrate_SW(seb.NEXT, S_down);
            
            % Upwelling shortwave fluxes
            I_out_from_dir = (h1./sigma + h2 + h3)*spectral_weights';
            I_out_from_dif = (h7 + h8)*spectral_weights';
            
            Sout = Sin_dif.*I_out_from_dif + Sin_dir.*I_out_from_dir;
            
            seb.TEMP.S_abs = Sin + S_up - S_down - Sout;
            seb.TEMP.d_energy(1) = seb.TEMP.d_energy(1) + seb.TEMP.S_abs;
            seb.TEMP.S_down = S_down;
            seb.TEMP.S_up = S_up;
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
            seb.TEMP.S_down = S_down./seb.STATVAR.area;
            seb.TEMP.S_up = S_up./seb.STATVAR.area;
        end
        
        function [seb, Sout] = penetrate_SW_transmission_bulk(seb, Sin)  %used with fixed albedo and SW extinction coefficient
            %S_up and S_down can in principle be spectrally resolved when provided as
            %row array, using e.g. spectral_ranges = [0.71 0.21 0.08];
            %SW extinction is assumed constant throughout class
            cut_off = 0.1*seb.STATVAR.area(1); % = 0.1 W/m2, radiation is not penetrated further if cutoff is reached
            
            Sout = seb.PARA.albedo .* Sin;
            S_down = (1 - seb.PARA.albedo) .* Sin;
            
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
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_up2);
                else
                    Sout = Sout + S_up2;  %add the uppwelling SW radiation to reflected, multiple reflections not accounted for!!
                end
                seb.TEMP.S_down = S_down;
                seb.TEMP.S_up = S_up;
            end
            seb.STATVAR.Sin = Sin;
            seb.STATVAR.Sout = Sout;
        end
        
        function [seb, Sout] = penetrate_SW_transmission_spectral(seb, Sin)  %used with variable albedo and SW exticntion coefficient
            %S_up and S_down are spectrally resolved when provided as
            %row array, using e.g. spectral_ranges = [0.71 0.21 0.08];
            %SW extinction is assumed constant throughout class
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
                    seb.TEMP.d_energy(i,1) = seb.TEMP.d_energy(i,1) + sum(S_down);
                else
                    Sout = Sout + S_up2;  %add the uppwelling SW radiation to reflected, multiple reflections not accounted for!!
                end
                seb.TEMP.S_down = S_down;
                seb.TEMP.S_up = seb.STATVAR.S_up;
            end
            seb.STATVAR.Sin = Sin;
            seb.STATVAR.Sout = Sout;
        end
        
        %calculates derivatives of water (i.e. water fluxes) due to evapotranspiration
        function seb = calculateET(seb)
            L_v = latent_heat_evaporation(seb, seb.STATVAR.T(1)+273.15);
            
            if seb.STATVAR.Qe_pot > 0
                fraction_T = getET_fraction(seb);
                fraction_E = getET_fraction(seb);
                
                depth_weighting_E = exp(-1./seb.PARA.evaporationDepth .* cumsum(seb.STATVAR.layerThick));  %exponential decrease with depth
                depth_weighting_E(depth_weighting_E<0.05) = 0;
                depth_weighting_E = depth_weighting_E .* seb.STATVAR.layerThick ./ sum(depth_weighting_E .* seb.STATVAR.layerThick,1); %normalize
                
                depth_weighting_T = exp(-1./seb.PARA.rootDepth .* cumsum(seb.STATVAR.layerThick));
                depth_weighting_T(depth_weighting_T<0.05) = 0;
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
                depth_weighting_E = depth_weighting_E .* seb.STATVAR.layerThick ./ sum(depth_weighting_E .* seb.STATVAR.layerThick,1); %normalize
                
                depth_weighting_T = exp(-1./seb.PARA.rootDepth .* cumsum(seb.STATVAR.layerThick));
                depth_weighting_T(depth_weighting_T<0.05) = 0;
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
            fraction=double(seb.STATVAR.T>0).*double(seb.STATVAR.T(1)>0).*(double(waterC >= seb.STATVAR.field_capacity) + double(waterC < seb.STATVAR.field_capacity).*0.25.*(1-cos(pi().*waterC./seb.STATVAR.field_capacity)).^2);
            
        end
        
        function fraction = getET_fraction_Xice(seb)
            
            waterC = seb.STATVAR.waterIce ./ (seb.STATVAR.layerThick .* seb.STATVAR.area - seb.STATVAR.XwaterIce);
            %waterC = seb.STATVAR.waterIce ./ max(1e-20,seb.STATVAR.layerThick .* seb.STATVAR.area - seb.STATVAR.XwaterIce);
            waterC(isnan(waterC)) = 0;
            fraction=double(seb.STATVAR.T>0).*double(seb.STATVAR.T(1)>0).*(double(waterC >= seb.STATVAR.field_capacity) + double(waterC < seb.STATVAR.field_capacity).*0.25.*(1-cos(pi().*waterC./seb.STATVAR.field_capacity)).^2);
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

