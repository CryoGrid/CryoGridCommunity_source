%========================================================================
% CryoGrid TIER1 INTERACTION (IA) class for functions related to surface
% energy balance extending through the top class (e.g. canopy/vegetation)
% R. B. Zweigel, August 2021
%========================================================================

classdef IA_SEB < IA_BASE
    
    methods
        
        function ia_seb_water = get_boundary_condition_Qh_m(ia_seb_water, tile)
            % equivalent to function Q_h(seb, forcing) in SEB
            forcing = tile.FORCING;
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % ground
            
            cp = stratigraphy2.CONST.cp;
            kappa = stratigraphy2.CONST.kappa;
            
            uz = forcing.TEMP.wind;
            z =  stratigraphy2.PARA.airT_height;
            z0 = stratigraphy2.PARA.z0;
            Tz = forcing.TEMP.Tair;
            TForcing = stratigraphy2.STATVAR.T(1);
            Lstar = stratigraphy2.STATVAR.Lstar;
            p = forcing.TEMP.p;
            
            Tz=Tz+forcing.CONST.Tmfw;
            TForcing=TForcing+forcing.CONST.Tmfw;
            rho = rho_air(stratigraphy2,p, Tz);
            
            Q_h  = -rho.*cp.*kappa.* uz.*kappa./(log(z./z0)- psi_M(stratigraphy2, z./Lstar, z0./Lstar)) .* (Tz-TForcing)./(log(z./z0)- psi_H(stratigraphy2, z./Lstar, z0./Lstar));
            stratigraphy2.STATVAR.Qh = Q_h;
        end
        
        function ia_seb_water = get_boundary_condition_Qe_m(ia_seb_water, tile)
            forcing = tile.FORCING;
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % ground
            
            stratigraphy2 = Q_evap_CLM4_5(stratigraphy2, forcing);
            
            stratigraphy2.TEMP.d_water_ET(1,1) = stratigraphy2.TEMP.d_water_ET(1,1) -  stratigraphy2.STATVAR.evap.* stratigraphy2.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
            stratigraphy2.TEMP.d_water_ET_energy(1,1) = stratigraphy2.TEMP.d_water_ET_energy(1,1) -  stratigraphy2.STATVAR.evap_energy.* stratigraphy2.STATVAR.area(1,1);
        end
        
        function ia_seb_water = get_boundary_condition_Qh_CLM5_m(ia_seb_water, tile)
            % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % ground
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
            
            stratigraphy2.TEMP.evap = water_fraction .* stratigraphy2.STATVAR.Qe ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.TEMP.sublim = ice_fraction .* stratigraphy2.STATVAR.Qe ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.TEMP.evap_energy =  stratigraphy2.TEMP.evap.*  (double(Tg>=0) .* stratigraphy2.CONST.c_w .* Tg + ...
                double(Tg<0) .* stratigraphy2.CONST.c_i .* Tg);
            stratigraphy2.TEMP.sublim_energy =  stratigraphy2.TEMP.sublim .* (stratigraphy2.CONST.c_i .* stratigraphy2.STATVAR.T(1,1) - stratigraphy2.CONST.L_f);
            
            stratigraphy2.TEMP.d_water_ET(1) = stratigraphy2.TEMP.d_water_ET(1) - (stratigraphy2.TEMP.evap + stratigraphy2.TEMP.sublim).*stratigraphy2.STATVAR.area(1);
            stratigraphy2.TEMP.d_water_ET_energy(1) = stratigraphy2.TEMP.d_water_ET_energy(1) - (stratigraphy2.TEMP.evap_energy + stratigraphy2.TEMP.sublim_energy).*stratigraphy2.STATVAR.area(1);
                        
        end
        
    end
end
