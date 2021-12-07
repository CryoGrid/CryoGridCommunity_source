%========================================================================
% CryoGrid INTERACTION (IA) class for surface energy balance of a GROUND
% class below a shading VEGETATION class
% R. B. Zwegiel, October 2021
%========================================================================

classdef IA_HEAT01_SEB11_vegetation_CLM5 < IA_SEB & IA_WATER % & IA_HEAT 
    
    methods
        
        function ia_seb_water = get_boundary_condition_m(ia_seb_water, tile)
            ia_seb_water = get_boundary_condition_Qh_CLM5_m(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_Qe_CLM5_m(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_RichardsEq_canopy_m(ia_seb_water, tile);
            ia_seb_water = get_water_transpiration(ia_seb_water);
            
            % add fluxes to uppermost cell (ratiative fluxes are added by penetration)
            ia_seb_water.NEXT.TEMP.d_energy(1) = ia_seb_water.NEXT.TEMP.d_energy(1) + (-ia_seb_water.NEXT.STATVAR.Qh - ia_seb_water.NEXT.STATVAR.Qe) .* ia_seb_water.NEXT.STATVAR.area(1);
            
            ia_seb_water.NEXT.TEMP.d_water_ET(1) = ia_seb_water.NEXT.TEMP.d_water_ET(1) -  ia_seb_water.NEXT.STATVAR.evap.* ia_seb_water.NEXT.STATVAR.area(1); %in m3 water per sec, put everything in uppermost grid cell
            ia_seb_water.NEXT.TEMP.d_water_ET_energy(1) = ia_seb_water.NEXT.TEMP.d_water_ET_energy(1) -  ia_seb_water.NEXT.STATVAR.evap_energy.* ia_seb_water.NEXT.STATVAR.area(1);
        end
        
        function z0 = get_z0_ground(ia_seb_water)
            z0 = ia_seb_water.NEXT.PARA.z0;
        end
        
        function q_g = get_humidity_ground(ia_seb_water, forcing)
            canopy = ia_seb_water.PREVIOUS;
            ground = ia_seb_water.NEXT;
            Tmfw = canopy.CONST.Tmfw; % Freezing T of water (K)
            R = ground.CONST.R; % Univ. gas constant [J/K/mol]
            M = ground.CONST.molar_mass_w; % [g/mol]
            g = ground.CONST.g;
            p = forcing.TEMP.p;
            Tg = ground.STATVAR.T(1) + Tmfw; % Tg in [K]
            psi_surf = ground.STATVAR.waterPotential(1);
            f_water = ground.STATVAR.water(1)./ground.STATVAR.waterIce(1);
            f_ice = ground.STATVAR.ice(1)./ground.STATVAR.waterIce(1);

            e_sat = f_water.*satPresWater(ground,Tg) + f_ice.*satPresIce(ground,Tg);
            alpha_soil = exp(psi_surf.*g./(R./M.*Tg));
            q_g_sat = .622.*e_sat./p;
            q_g = q_g_sat.*alpha_soil;
            ground.STATVAR.q_g = q_g; % save for calculation of ground Qe
        end
        
        function r_soil = ground_resistance_evap(ia_seb_water, forcing)
            % soil resistance to evapotranspiration - based on Dry Surface Layer parameterization in CLM5 / Swenson and Lawrence (2014)
            canopy = ia_seb_water.PREVIOUS;
            ground = ia_seb_water.NEXT;
            Tmfw = forcing.CONST.Tmfw;
            Dmax = canopy.PARA.Dmax; % maximum thickness of DSL
            tau = canopy.CONST.tau; %  tortuosity of the vapor flow paths through the soil matrix 
            Tg = ground.STATVAR.T(1) + Tmfw; % Tg in [K]
            theta_init = ground.STATVAR.field_capacity(1);
            theta_surf = ground.STATVAR.water(1)./ground.STATVAR.layerThick(1);
            
            Dv = 2.12.* 10^(-5).*(Tg./Tmfw)^1.75; % molecular diffusivity of water vapor in air
            DSL = double(theta_surf<theta_init).*Dmax.*(theta_init-theta_surf)./theta_init;
            r_soil = DSL./(Dv*tau);
        end
        
        function beta_t = get_soil_water_stress(ia_seb_water)
            % adopted from section 8.4 in CLM4.5! 
            stratigraphy2 = ia_seb_water.NEXT; % ground
            r = stratigraphy2.STATVAR.f_root; % layer root fraction
            
            w = double(stratigraphy2.STATVAR.water>0).*double(stratigraphy2.STATVAR.T>-2); % Plant wilting factor
            
            beta_t = sum(r.*w);
            
            stratigraphy2.TEMP.w = w;
        end
            
    end
end