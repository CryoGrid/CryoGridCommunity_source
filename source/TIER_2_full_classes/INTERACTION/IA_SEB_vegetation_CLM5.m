%========================================================================
% CryoGrid INTERACTION (IA) class for surface energy balance of a GROUND
% class below a shading VEGETATION class
% R. B. Zwegiel, October 2021
%========================================================================

classdef IA_SEB_vegetation_CLM5 < IA_SEB & IA_WATER % & IA_HEAT
    
    methods
        
        function ia_seb_water = get_boundary_condition_m(ia_seb_water, tile)
            ia_seb_water = get_boundary_condition_Qh_CLM5_m(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_Qe_CLM5_m(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_RichardsEq_canopy_m(ia_seb_water, tile);
            
            % add fluxes to uppermost cell (ratiative fluxes are added by penetration)
            ia_seb_water.NEXT.TEMP.F_ub = (-ia_seb_water.NEXT.STATVAR.Qh - ia_seb_water.NEXT.STATVAR.Qe) .* ia_seb_water.NEXT.STATVAR.area(1);
            ia_seb_water.NEXT.TEMP.d_energy(1) = ia_seb_water.NEXT.TEMP.d_energy(1) + ia_seb_water.NEXT.TEMP.F_ub; % cant skip F_ub, it is required for the ground-snow IA_class
        end
        
        function q_g = get_humidity_surface(ia_seb_water, tile)
            q_g = get_humidity_surface_GROUND(ia_seb_water, tile);
        end
        
        function r_soil = ground_resistance_evap(ia_seb_water, tile)
            % soil resistance to evapotranspiration - based on Dry Surface Layer parameterization in CLM5 / Swenson and Lawrence (2014)
            forcing = tile.FORCING;
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
        
        function ia_seb_water = distribute_roots(ia_seb_water)
            %             biomass_root = ia_seb_water.PARA.biomass_root;
            %             density_root = ia_seb_water.PARA.density_root;
            %             r_root = ia_seb_water.PARA.r_root;
            beta = ia_seb_water.PARA.beta_root;
            dz = ia_seb_water.NEXT.STATVAR.layerThick;
            z = cumsum(dz);
            
            % Root fraction per soil layer
            f_root = beta.^([0; z(1:end-1)].*100) - beta.^(z*100); % Eq. 11.1
            
            %             % Root spacing
            %             CA_root = pi.*r_root.^2; % Eq. 11.5, fine root cross sectional area
            %             B_root = biomass_root.*f_root./dz; % Eq. 11.4, root biomass density
            %             L_root = B_root./(density_root.*CA_root); % Eq. 11.3, root length density
            %             dx_root = (pi.*L_root).^(-1/2);
            
            %             ia_seb_water.NEXT.STATVAR.dx_root = dx_root;
            ia_seb_water.NEXT.STATVAR.f_root = f_root;
            %             ia_seb_water.NEXT.TEMP.w = f_root.*0;
            
            tkjerht % Throw error!
            
        end
        
        function ia_seb_water = canopy_drip(ia_seb_water, tile)
            ia_seb_water = canopy_drip@IA_WATER(ia_seb_water, tile);
        end
        
    end
end