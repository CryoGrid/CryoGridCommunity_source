%========================================================================
% CryoGrid TIER1 INTERACTION (IA) class for functions related to surface
% energy balance extending through the top class (e.g. canopy/vegetation)
% R. B. Zweigel, August 2021
%========================================================================

classdef IA_SEB < IA_BASE
    
    methods
        
%         function ia_seb_water = get_boundary_condition_Qh_m(ia_seb_water, tile)
%             % equivalent to function Q_h(seb, forcing) in SEB
%             forcing = tile.FORCING;
%             stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
%             stratigraphy2 = ia_seb_water.NEXT; % ground
%             
%             cp = stratigraphy2.CONST.cp;
%             kappa = stratigraphy2.CONST.kappa;
%             
%             uz = forcing.TEMP.wind;
%             z =  stratigraphy2.PARA.airT_height;
%             z0 = stratigraphy2.PARA.z0;
%             Tair = forcing.TEMP.Tair;
%             Tsurf = stratigraphy2.STATVAR.T(1);
%             Lstar = stratigraphy2.STATVAR.Lstar;
%             p = forcing.TEMP.p;
%             
%             Tair=Tair+forcing.CONST.Tmfw;
%             Tsurf=Tsurf+forcing.CONST.Tmfw;
%             rho = rho_air(stratigraphy2,p, Tair);
%             
%             Q_h  = -rho.*cp.*kappa.* uz.*kappa./(log(z./z0)- psi_M(stratigraphy2, z./Lstar, z0./Lstar)) .* (Tair-Tsurf)./(log(z./z0)- psi_H(stratigraphy2, z./Lstar, z0./Lstar));
%             stratigraphy2.STATVAR.Qh = Q_h;
%         end
        
%         function ia_seb_water = get_boundary_condition_Qe_m(ia_seb_water, tile)
%             forcing = tile.FORCING;
%             stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
%             stratigraphy2 = ia_seb_water.NEXT; % ground
%             
%             stratigraphy2 = Q_evap_CLM4_5(stratigraphy2, forcing);
%             
%             stratigraphy2.TEMP.d_water_ET(1,1) = stratigraphy2.TEMP.d_water_ET(1,1) + stratigraphy2.STATVAR.evaporation.* stratigraphy2.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
%             stratigraphy2.TEMP.d_water_ET_energy(1,1) = stratigraphy2.TEMP.d_water_ET_energy(1,1) + stratigraphy2.TEMP.evaporation_energy.* stratigraphy2.STATVAR.area(1,1);
%         end
        
                
        function ia_seb_water = get_boundary_condition_Qh_CLM5_m(ia_seb_water, tile)
            % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % ground/snow
            Tg = stratigraphy2.STATVAR.T(1)+tile.FORCING.CONST.Tmfw;
            Cp = stratigraphy1.CONST.cp;
            Ts = stratigraphy1.STATVAR.Ts+tile.FORCING.CONST.Tmfw; % canopy air temperature (K)
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density
            r_ah_prime = stratigraphy1.TEMP.r_a_prime; % aerodynamic resistance to heat transfer soil - canopy air
            
            stratigraphy2.STATVAR.Qh  = -rho_atm.*Cp.*(Ts-Tg)./r_ah_prime; % Eq. 5.90
        end
        
        function ia_seb_water = get_boundary_condition_Qe_CLM5_m(ia_seb_water, tile)
            % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % ground
            q_g = stratigraphy2.STATVAR.q_g; % ground specific humidity
            q_s = stratigraphy1.STATVAR.q_s; % canopy air specific humidity
            Tg = stratigraphy2.STATVAR.T(1);
            Tmfw = tile.FORCING.CONST.Tmfw;
            r_aw_prime = stratigraphy1.TEMP.r_a_prime;
            r_soil = stratigraphy1.TEMP.r_soil;
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density
            
            water_fraction = stratigraphy2.STATVAR.water(1) ./ stratigraphy2.STATVAR.waterIce(1);
            ice_fraction = stratigraphy2.STATVAR.ice(1) ./ stratigraphy2.STATVAR.waterIce(1);
            latent_heat = water_fraction.*latent_heat_evaporation(stratigraphy2, Tg+Tmfw) + ice_fraction.*latent_heat_sublimation(stratigraphy2, Tg+Tmfw);
            
            stratigraphy2.STATVAR.Qe = -rho_atm.*latent_heat.*(q_s-q_g)./(r_aw_prime+r_soil);
            
            stratigraphy2.STATVAR.evaporation = -water_fraction .* stratigraphy2.STATVAR.Qe ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.STATVAR.sublimation = -ice_fraction .* stratigraphy2.STATVAR.Qe ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.TEMP.evaporation_energy =  stratigraphy2.STATVAR.evaporation.*  (double(Tg>=0) .* stratigraphy2.CONST.c_w .* Tg + ...
                double(Tg<0) .* stratigraphy2.CONST.c_i .* Tg);
            stratigraphy2.TEMP.sublimation_energy =  stratigraphy2.STATVAR.sublimation .* (stratigraphy2.CONST.c_i .* stratigraphy2.STATVAR.T(1,1) - stratigraphy2.CONST.L_f);
            
            stratigraphy2.TEMP.d_water_ET(1) = stratigraphy2.TEMP.d_water_ET(1) + stratigraphy2.STATVAR.evaporation.*stratigraphy2.STATVAR.area(1);
            stratigraphy2.TEMP.d_water_ET_energy(1) = stratigraphy2.TEMP.d_water_ET_energy(1) + stratigraphy2.TEMP.evaporation_energy.*stratigraphy2.STATVAR.area(1);
                        
        end
        
        function ia_seb_water = get_boundary_condition_Qe_CLM5_m_Xice(ia_seb_water, tile)
                        % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % ground
            q_g = stratigraphy2.STATVAR.q_g; % ground specific humidity
            q_s = stratigraphy1.STATVAR.q_s; % canopy air specific humidity
            Tg = stratigraphy2.STATVAR.T(1);
            Tmfw = tile.FORCING.CONST.Tmfw;
            r_aw_prime = stratigraphy1.TEMP.r_a_prime;
            r_soil = stratigraphy1.TEMP.r_soil;
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density
            
            water_fraction = (stratigraphy2.STATVAR.water(1,1)+stratigraphy2.STATVAR.Xwater(1,1) ) ./ (stratigraphy2.STATVAR.waterIce(1,1) + stratigraphy2.STATVAR.XwaterIce(1,1));
            ice_fraction = (stratigraphy2.STATVAR.ice(1,1) + stratigraphy2.STATVAR.Xice(1,1)) ./ (stratigraphy2.STATVAR.waterIce(1,1) + stratigraphy2.STATVAR.XwaterIce(1,1));
            water_fraction(stratigraphy2.STATVAR.waterIce == 0) = double(Tg>=0);
            ice_fraction(stratigraphy2.STATVAR.waterIce == 0) = double(Tg<0);
            latent_heat = water_fraction.*latent_heat_evaporation(stratigraphy2, Tg+Tmfw) + ice_fraction.*latent_heat_sublimation(stratigraphy2, Tg+Tmfw);
            
            stratigraphy2.STATVAR.Qe = -rho_atm.*latent_heat.*(q_s-q_g)./(r_aw_prime+r_soil);
            
            stratigraphy2.STATVAR.evaporation = -water_fraction .* stratigraphy2.STATVAR.Qe ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.STATVAR.sublimation = -ice_fraction .* stratigraphy2.STATVAR.Qe ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.TEMP.evaporation_energy =  stratigraphy2.STATVAR.evaporation.*  (double(Tg>=0) .* stratigraphy2.CONST.c_w .* Tg + ...
                double(Tg<0) .* stratigraphy2.CONST.c_i .* Tg);
            stratigraphy2.TEMP.sublimation_energy =  stratigraphy2.STATVAR.sublimation .* (stratigraphy2.CONST.c_i .* Tg - stratigraphy2.CONST.L_f);
            
            stratigraphy2.TEMP.d_water_ET(1) = stratigraphy2.TEMP.d_water_ET(1) + stratigraphy2.STATVAR.evaporation.*stratigraphy2.STATVAR.area(1);
            stratigraphy2.TEMP.d_water_ET_energy(1) = stratigraphy2.TEMP.d_water_ET_energy(1) + stratigraphy2.TEMP.evaporation_energy.*stratigraphy2.STATVAR.area(1);
        end
        
        function ia_seb_water = get_boundary_condition_Qe_CLM5_SNOW_m(ia_seb_water, tile)
            % equivalent to Q_eq_potET_snow, but for snow below canopy
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % snow
            q_g = stratigraphy2.STATVAR.q_g; % snow specific humidity
            q_s = stratigraphy1.STATVAR.q_s; % canopy air specific humidity
            Tsurf = stratigraphy2.STATVAR.T(1);
            Tmfw = tile.FORCING.CONST.Tmfw;
            r_aw_prime = stratigraphy1.TEMP.r_a_prime;
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density
            
%             water_fraction = stratigraphy2.STATVAR.water(1) ./ stratigraphy2.STATVAR.waterIce(1);
%             ice_fraction = stratigraphy2.STATVAR.ice(1) ./ stratigraphy2.STATVAR.waterIce(1);
            %saturation = stratigraphy2.STATVAR.water(1)./(stratigraphy2.STATVAR.layerThick(1).*stratigraphy2.STATVAR.area(1) - stratigraphy2.STATVAR.ice(1));
            porespace = stratigraphy2.STATVAR.layerThick(1).*stratigraphy2.STATVAR.area(1) - stratigraphy2.STATVAR.ice(1);
            saturation = stratigraphy2.STATVAR.water(1)./porespace;
            saturation(porespace == 0) = 0;
            L_w = latent_heat_evaporation(stratigraphy2, Tsurf+Tmfw); 
            L_i = latent_heat_sublimation(stratigraphy2, Tsurf+Tmfw);
            
            if Tsurf < 0 || saturation <= 0 % All sublimation
                stratigraphy2.STATVAR.Qe = -rho_atm.*L_i.*(q_s-q_g)./r_aw_prime;
                stratigraphy2.STATVAR.sublimation = - stratigraphy2.STATVAR.Qe ./ (L_i .* stratigraphy2.CONST.rho_w) .*stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.sublimation_energy =  stratigraphy2.STATVAR.sublimation .* (stratigraphy2.CONST.c_i .* stratigraphy2.STATVAR.T(1,1) - stratigraphy2.CONST.L_f);
                stratigraphy2.STATVAR.evaporation = 0;
                stratigraphy2.TEMP.evaporation_energy = 0;
            else
                evap_fraction = max(0, min(saturation./(2.*stratigraphy2.PARA.field_capacity),1)); % all evap when saturation is more than twice of field capacity (in fraction of porespace)
                sublim_fraction = 1 - evap_fraction;
                Qe_sublim = -sublim_fraction.*rho_atm.*L_i.*(q_s-q_g)./r_aw_prime;
                Qe_evap = -evap_fraction.*rho_atm.*L_w.*(q_s-q_g)./r_aw_prime;
                stratigraphy2.STATVAR.sublimation = - Qe_sublim ./(stratigraphy2.CONST.rho_w .* L_i) .* stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.sublimation_energy = stratigraphy2.STATVAR.sublimation .* (Tsurf .* stratigraphy2.CONST.c_i - stratigraphy2.CONST.L_f);
                stratigraphy2.STATVAR.evaporation = - Qe_evap ./(stratigraphy2.CONST.rho_w .* L_w) .* stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.evaporation_energy = stratigraphy2.STATVAR.evaporation .* stratigraphy2.CONST.c_w;
                
                stratigraphy2.STATVAR.Qe = Qe_evap + Qe_sublim; % Note this is not like in Q_eq_potET_snow! RBZ
            end
            % Note: no d_water_ET for SNOW, evap/sublimation is added to water fluxes in advance_prognostics(..)
        end
        
        function ia_seb_water = get_boundary_condition_m_CHILD(ia_seb_water, tile)
            % SEB and water balance of child below a canopy, excluding radiation
            forcing = tile.FORCING;
            snow = ia_seb_water.NEXT.CHILD;
            snow = get_boundary_condition_allSNOW_rain_canopy_m(snow, tile); %add full snow, but rain only for snow-covered part
            snow = get_boundary_condition_water_SNOW_canopy_m(snow, tile);
            
            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing);
            
            ia_seb_water = get_boundary_condition_Qh_CLM5_CHILD_m(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_Qe_CLM5_CHILD_m(ia_seb_water, tile);
            
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + (-snow.STATVAR.Qh - snow.STATVAR.Qe).*snow.STATVAR.area(1);
        end
        
        function ia_seb_water = get_boundary_condition_m_vegetation_snow_CHILD_below(ia_seb_water, tile)
            % SEB and water balance of child below a canopy with snow, excluding radiation
            forcing = tile.FORCING;
            snow = ia_seb_water.NEXT.CHILD;
            vegetation = ia_seb_water.NEXT.PREVIOUS;
            
            %replaces
            %snow = get_boundary_condition_allSNOW_rain_canopy_m(snow, tile); %add full snow, but rain only for snow-covered part
            snow.TEMP.snowfall = vegetation.TEMP.snow_thru .* (snow.PARENT.STATVAR.area(1,1) + snow.STATVAR.area); %snowfall is in mm/day -> [m3/sec]
            %snow.TEMP.rainfall = vegetation.TEMP.rain_thru .* snow.STATVAR.area;
            snow.TEMP.snow_energy = vegetation.TEMP.snow_thru .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);  %[J/sec]
            %snow.TEMP.rain_energy = vegetation.TEMP.rain_thru .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;
        
            %rainfall handled here
            snow = get_boundary_condition_water_SNOW_canopy_m(snow, tile);
            
            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing);
            
            ia_seb_water = get_boundary_condition_Qh_CLM5_CHILD_m(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_Qe_CLM5_CHILD_m(ia_seb_water, tile);
            
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + (-snow.STATVAR.Qh - snow.STATVAR.Qe).*snow.STATVAR.area(1);
        end
        
        
        function ia_seb_water = get_boundary_condition_Qh_CLM5_CHILD_m(ia_seb_water, tile)
            % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT.CHILD; % snow CHILD
            Tg = stratigraphy2.STATVAR.T(1)+tile.FORCING.CONST.Tmfw;
            Cp = stratigraphy1.CONST.cp;
            Ts = stratigraphy1.STATVAR.Ts+tile.FORCING.CONST.Tmfw; % canopy air temperature (K)
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density
            r_ah_prime = stratigraphy1.TEMP.r_a_prime; % aerodynamic resistance to heat transfer soil - canopy air
            
            stratigraphy2.STATVAR.Qh  = -rho_atm.*Cp.*(Ts-Tg)./r_ah_prime; % Eq. 5.90
        end
        
        function ia_seb_water = get_boundary_condition_Qe_CLM5_CHILD_m(ia_seb_water, tile)
            % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT.CHILD; % snow CHILD
            q_g = stratigraphy2.STATVAR.q_g; % ground specific humidity
            q_s = stratigraphy1.STATVAR.q_s; % canopy air specific humidity
            Tsurf = stratigraphy2.STATVAR.T(1);
            Tmfw = tile.FORCING.CONST.Tmfw;
            r_aw_prime = stratigraphy1.TEMP.r_a_prime;
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density
            
            porespace = stratigraphy2.STATVAR.layerThick(1).*stratigraphy2.STATVAR.area(1) - stratigraphy2.STATVAR.ice(1);
            saturation = stratigraphy2.STATVAR.water(1)./porespace;
            saturation(porespace == 0) = 0;
            L_w = latent_heat_evaporation(stratigraphy2, Tsurf+Tmfw); 
            L_i = latent_heat_sublimation(stratigraphy2, Tsurf+Tmfw);
            
            if Tsurf < 0 || saturation <= 0 % All sublimation
                stratigraphy2.STATVAR.Qe = -rho_atm.*L_i.*(q_s-q_g)./r_aw_prime;
                stratigraphy2.STATVAR.sublimation = - stratigraphy2.STATVAR.Qe ./ (L_i .* stratigraphy2.CONST.rho_w) .*stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.sublimation_energy =  stratigraphy2.STATVAR.sublimation .* (stratigraphy2.CONST.c_i .* stratigraphy2.STATVAR.T(1,1) - stratigraphy2.CONST.L_f);
                stratigraphy2.STATVAR.evaporation = 0;
                stratigraphy2.TEMP.evaporation_energy = 0;
            else
                evap_fraction = max(0, min(saturation./(2.*stratigraphy2.PARA.field_capacity),1)); % all evap when saturation is more than twice of field capacity (in fraction of porespace)
                sublim_fraction = 1 - evap_fraction;
                Qe_sublim = -sublim_fraction.*rho_atm.*L_i.*(q_s-q_g)./r_aw_prime;
                Qe_evap = -evap_fraction.*rho_atm.*L_w.*(q_s-q_g)./r_aw_prime;
                stratigraphy2.STATVAR.sublimation = - Qe_sublim ./(stratigraphy2.CONST.rho_w .* L_i) .* stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.sublimation_energy = stratigraphy2.STATVAR.sublimation .* (Tsurf .* stratigraphy2.CONST.c_i - stratigraphy2.CONST.L_f);
                stratigraphy2.STATVAR.evaporation = - Qe_evap ./(stratigraphy2.CONST.rho_w .* L_w) .* stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.evaporation_energy = stratigraphy2.STATVAR.evaporation .* stratigraphy2.CONST.c_w;
                
                stratigraphy2.STATVAR.Qe = Qe_evap + Qe_sublim; % Note this is not like in Q_eq_potET_snow! RBZ
            end
        end
        
        %SERVICE FUNCTIONS,MOVE TO SEB.m and remove everything with canopy,
        %called in Qe_CLM5, function
        function q_g = get_humidity_surface_GROUND(ia_seb_water, tile)
            % Specific humidity of ground surface as in CLM5 documentation
            % (sect. 5.2)
            forcing = tile.FORCING;
            canopy = ia_seb_water.PREVIOUS;
            ground = ia_seb_water.NEXT;
            Tmfw = canopy.CONST.Tmfw; % Freezing T of water (K)
            R = ground.CONST.R; % Univ. gas constant [J/K/mol]
            M = ground.CONST.molar_mass_w; % [g/mol] SW,Dec 2022, must be in kg/mol, otherwise units in alpha_soil do not match! 
            g = ground.CONST.g;
            p = forcing.TEMP.p;
            Tg = ground.STATVAR.T(1) + Tmfw; % Tg in [K]
            psi_surf = ground.STATVAR.waterPotential(1);
            f_water = ground.STATVAR.water(1)./ground.STATVAR.waterIce(1);
            f_ice = ground.STATVAR.ice(1)./ground.STATVAR.waterIce(1);

            e_sat = f_water.*satPresWater(ground,Tg) + f_ice.*satPresIce(ground,Tg);
            alpha_soil = exp(psi_surf.*g./(R./M.*Tg)); % Eq. 5.73
            q_g_sat = .622.*e_sat./p;
            q_g = q_g_sat.*alpha_soil; % Eq. 5.72
            ground.STATVAR.q_g = q_g; % save for calculation of ground Qe
        end
        
        function q_g = get_humidity_surface_SNOW(ia_seb_water, tile)
            forcing = tile.FORCING;
            p = forcing.TEMP.p;
            Tmfw = ia_seb_water.NEXT.CONST.Tmfw;
            T_surf = ia_seb_water.NEXT.STATVAR.T(1);
            
            e = double(T_surf>=0).*satPresWater(ia_seb_water.NEXT,T_surf+Tmfw) + double(T_surf<0).*satPresIce(ia_seb_water.NEXT,T_surf+Tmfw); % saturation water pressure of snow
            q_g = .622.*e./p; % saturation water vapor specific humidity at snow temperature
            ia_seb_water.NEXT.STATVAR.q_g = q_g;
        end
        
        function q_g = get_humidity_surface_SNOW_CHILD(ia_seb_water, tile)
            forcing = tile.FORCING;
            snow = ia_seb_water.NEXT.CHILD;
            p = forcing.TEMP.p;
            Tmfw = snow.CONST.Tmfw;
            T_surf = snow.STATVAR.T(1);
            
            e = double(T_surf>=Tmfw).*satPresWater(snow,T_surf+Tmfw) + double(T_surf<Tmfw).*satPresIce(snow,T_surf+Tmfw); % saturation water pressure of snow
            q_g = .622.*e./p; % saturation water vapor specific humidity at snow temperature
            snow.STATVAR.q_g = q_g;
        end
        
    end
end
