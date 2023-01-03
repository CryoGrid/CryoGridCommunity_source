%========================================================================
% CryoGrid vegetation class VEGETATION_CLM5_seb_snow
%
% as VEGETATION_CLM5_seb, but with canopy interception of snowfall and snow
% processes within the canopy
%
% S. Westermann, R. Zweigel, Dec 2022
%========================================================================

classdef VEGETATION_CLM5_seb_snow < VEGETATION_CLM5_seb
    properties
        CHILD
        IA_CHILD
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

       function vegetation = provide_PARA(vegetation)  
            vegetation = provide_PARA@VEGETATION_CLM5_seb(vegetation);
            vegetation.PARA.max_snow_fraction = 0.5; %set to zero when no interception shall occur
            vegetation.PARA.SWE_max = 0.01; % in [m]
       end
       
       function vegetation = provide_CONST(vegetation)  
           vegetation = provide_CONST@VEGETATION_CLM5_seb(vegetation);
       end
       
       function vegetation = provide_STATVAR(vegetation)  
           vegetation = provide_STATVAR@VEGETATION_CLM5_seb(vegetation);
       end
       
       function vegetation = finalize_init(vegetation, tile)
           vegetation = finalize_init@VEGETATION_CLM5_seb(vegetation, tile);
           vegetation.CHILD = 0; % no snow
           vegetation.IA_CHILD = 0;
           vegetation.STATVAR.max_vegetation_height = sum(vegetation.STATVAR.layerThick); %use to recompute layerThick when snow accumulates
           %change also LAI and SAI
           %recompute emissivity, z0 and d and the energy
       end

       
        %---time integration------
        
        function vegetation = get_boundary_condition_u(vegetation, tile)
            forcing = tile.FORCING;
            vegetation.TEMP.snow_thru = 0;
            vegetation.TEMP.rain_thru = 0;
            if vegetation.CHILD == 0  %CHILD does not exist
                vegetation = get_boundary_condition_u@VEGETATION_CLM5_seb(vegetation, tile); %call the native function for the vegetation class
                
                if forcing.TEMP.snowfall > 0 && vegetation.PARA.max_snow_fraction > 0 %create CHILD 
                    
                    vegetation.CHILD = copy(tile.STORE.SNOW);
                    vegetation.CHILD.PARENT = vegetation;
                    vegetation.CHILD.NEXT = vegetation; 
                    vegetation.IA_CHILD = get_IA_class(class(vegetation.CHILD), class(vegetation));
                    vegetation.IA_CHILD.PREVIOUS = vegetation.CHILD;
                    vegetation.IA_CHILD.NEXT = vegetation; %SNOW CHILD created
                    vegetation.CHILD.PARA.swe_per_cell = vegetation.PARA.SWE_max;
                    
                    %modify forcing and then make back
                    f_intr_snow = tanh(vegetation.STATVAR.LAI + vegetation.STATVAR.SAI); %same as in "get_boundary_condition_u_water_canopy(canopy, tile)"
                    tile.FORCING.TEMP.snowfall = tile.FORCING.TEMP.snowfall .* f_intr_snow;
                    vegetation.CHILD = get_boundary_condition_u_create_CHILD(vegetation.CHILD, tile);  %initialize with fresh snowfall
                    tile.FORCING.TEMP.snowfall = tile.FORCING.TEMP.snowfall ./ f_intr_snow;
                    
                    vegetation.TEMP.snow_thru = tile.FORCING.TEMP.snowfall ./1000 ./(24.*3600) .* (1 - f_intr_snow);
                    %make IA functions for using snow and rain throughfall 
                end
            else %CHILD exists
                snow_fraction = vegetation.CHILD.STATVAR.area./vegetation.STATVAR.area;
                 
                % 1. Longwave penetration
                [vegetation, L_up] = penetrate_LW(vegetation, forcing.TEMP.Lin .* vegetation.STATVAR.area(1));
                
                % 2. Shortwave penetration
                vegetation.TEMP.sun_angle = forcing.TEMP.sunElevation * double(forcing.TEMP.sunElevation > 0);
                vegetation.TEMP.Sin_dir_fraction = forcing.TEMP.Sin_dir ./ forcing.TEMP.Sin;
                vegetation.TEMP.Sin_dir_fraction(isnan(vegetation.TEMP.Sin_dir_fraction)) = 0;
                [vegetation, S_up] = penetrate_SW(vegetation, forcing.TEMP.Sin .* vegetation.STATVAR.area(1));
                
                % 3. Sensible and latent heat
                vegetation = canopy_resistances_CLM5_Stewart_snow(vegetation, tile);
                
                vegetation = Q_h_CLM5_snow(vegetation, tile);
                vegetation = Q_e_CLM5_Stewart_snow(vegetation, tile);
                
                vegetation.TEMP.d_energy = vegetation.TEMP.d_energy - (vegetation.STATVAR.Qh + vegetation.STATVAR.Qe).*(1-snow_fraction) .* vegetation.STATVAR.area;
                vegetation.TEMP.d_water_ET = vegetation.TEMP.d_water_ET - (vegetation.STATVAR.evap + vegetation.STATVAR.sublim).*(1-snow_fraction).*vegetation.STATVAR.area; % transpired water is removed from soil, not canopy
                vegetation.TEMP.d_water_ET_energy = vegetation.TEMP.d_water_ET_energy - (vegetation.STATVAR.evap_energy + vegetation.STATVAR.sublim_energy).*(1-snow_fraction).*vegetation.STATVAR.area;
                
                vegetation.CHILD.TEMP.d_energy = vegetation.CHILD.TEMP.d_energy - (vegetation.CHILD.STATVAR.Qh + vegetation.CHILD.STATVAR.Qe).*vegetation.CHILD.STATVAR.area;
                %sublimation is removed in advance_prognostic
%                 vegetation.CHILD.TEMP.d_water = vegetation.CHILD.TEMP.d_water + vegetation.CHILD.STATVAR.sublimation .*vegetation.CHILD.STATVAR.area; % transpired water is removed from soil, not vegetation
%                 vegetation.CHILD.TEMP.d_water_energy = vegetation.CHILD.TEMP.d_water_ET_energy + vegetation.CHILD.TEMP.sublimation_energy.*vegetation.CHILD.STATVAR.area;

                 vegetation.STATVAR.Qh = (1-snow_fraction) .* vegetation.STATVAR.Qh + snow_fraction .* vegetation.CHILD.STATVAR.Qh;
                 vegetation.STATVAR.Qe = (1-snow_fraction) .* vegetation.STATVAR.Qe + snow_fraction .* vegetation.CHILD.STATVAR.Qe; 

                
                %rain and snowfall
                f_intr = tanh(vegetation.STATVAR.LAI + vegetation.STATVAR.SAI); %same as in "get_boundary_condition_u_water_canopy(canopy, tile)"
                snow_fraction = vegetation.CHILD.STATVAR.area(1)./vegetation.STATVAR.area(1);

                %snow
                snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600);
                vegetation.TEMP.snow_thru = snowfall .* (1 - f_intr); %in per m2
                vegetation.CHILD.TEMP.snowfall =  snowfall  .* f_intr .* vegetation.STATVAR.area; %-> all goes to the snow CHILD
                vegetation.CHILD.TEMP.snow_energy = vegetation.CHILD.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* vegetation.CHILD.CONST.c_i - vegetation.CHILD.CONST.L_f);  %[J/sec]
                snow_property_function = str2func(vegetation.CHILD.PARA.snow_property_function);
                vegetation.CHILD = snow_property_function(vegetation.CHILD, forcing); 
                
                rainfall = forcing.TEMP.rainfall ./1000 ./(24.*3600);
                vegetation.TEMP.rain_thru = rainfall .* (1 - f_intr); 
                vegetation.CHILD.TEMP.rainfall =  rainfall  .* f_intr .* snow_fraction .* vegetation.STATVAR.area;            
                vegetation.CHILD.TEMP.rain_energy = vegetation.CHILD.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* vegetation.CHILD.CONST.c_w;
                vegetation.CHILD.TEMP.d_water = vegetation.CHILD.TEMP.d_water + vegetation.CHILD.TEMP.rainfall;
                vegetation.CHILD.TEMP.d_water_energy =  vegetation.CHILD.TEMP.d_water_energy + vegetation.CHILD.TEMP.rain_energy;
                
                %vegetation
                rain_vegetation = rainfall .* f_intr .* (1 - snow_fraction) .* vegetation.STATVAR.area;
                vegetation.TEMP.d_water = vegetation.TEMP.d_water + rain_vegetation;
                vegetation.TEMP.d_water_energy = vegetation.TEMP.d_water_energy(1) + rain_vegetation .* max(0, forcing.TEMP.Tair) .* vegetation.CONST.c_w;
            end
        end
        
        function [vegetation, S_up] = penetrate_SW(vegetation, S_down)  %mandatory function when used with class that features SW penetration
            if vegetation.CHILD ~= 0
                snow_fraction = vegetation.CHILD.STATVAR.area(1)./vegetation.STATVAR.area(1);
                vegetation_fraction = 1 - snow_fraction;
                [vegetation, S_up] = penetrate_SW@VEGETATION_CLM5_seb(vegetation, S_down.*vegetation_fraction);
                [vegetation.CHILD, S_up2] = penetrate_SW(vegetation.CHILD, S_down.*snow_fraction); %should first pass the snow, then vegetation, and then snow again 
                %alternative: just use reflection and omit transmission
%                 vegetation.CHILD.d_energy = vegetation.CHILD.d_energy + S_down.*snow_fraction .*(1-vegetation.CHILD.STATVAR.albedo);
%                 S_up2 = S_down.*snow_fraction .* vegetation.CHILD.STATVAR.albedo;
                
                S_up = S_up + sum(S_up2); % snow_crocus splits SW into spectral bands
            else
                [vegetation, S_up] = penetrate_SW@VEGETATION_CLM5_seb(vegetation, S_down);
            end
        end
        
        function [vegetation, S_up] = penetrate_SW_PARENT(vegetation, S_down) %called below the snow cover
            [vegetation, S_up] = penetrate_SW_PARENT@VEGETATION_CLM5_seb(vegetation, S_down);
        end
        
        function [vegetation, Lout] = penetrate_LW(vegetation, Lin)  %mandatory function when used with class that features SW penetration
            % Lin is in W, not W/m2!
            
             if vegetation.CHILD ~= 0
                 % Lin is in W, not W/m2!!!!
                 canopy_emissivity = vegetation.STATVAR.emissivity; %possibly change emssivity
                 Tmfw = vegetation.CONST.Tmfw;
                 sigma = vegetation.CONST.sigma;
                 snow_fraction = vegetation.CHILD.STATVAR.area(1)./vegetation.STATVAR.area(1);
                 
                 % Transmitted and emitted downward below canopy
                 L_down = (1 - canopy_emissivity).*Lin + canopy_emissivity.*sigma.*((1-snow_fraction).*vegetation.STATVAR.T + snow_fraction.*vegetation.CHILD.STATVAR.T + Tmfw)^4.*vegetation.STATVAR.area;
                 
                 % penetrate LW to class below
                 [seb.NEXT, L_up] = penetrate_LW(vegetation.NEXT, L_down);
                 
                 % transmitted/emitted/reflected upward
                 Lout = (1- canopy_emissivity).*L_up + canopy_emissivity.*sigma.*((1-snow_fraction).*vegetation.STATVAR.T + snow_fraction.*vegetation.CHILD.STATVAR.T + Tmfw)^4.*vegetation.STATVAR.area;
                 
                 vegetation.TEMP.L_abs = Lin + L_up - Lout - L_down;
                 vegetation.TEMP.d_energy(1) = vegetation.TEMP.d_energy(1) + (1-snow_fraction).*vegetation.TEMP.L_abs;
                 vegetation.CHILD.TEMP.d_energy = vegetation.CHILD.TEMP.d_energy + snow_fraction.*vegetation.TEMP.L_abs;
                 vegetation.TEMP.L_down = L_down./vegetation.STATVAR.area(1);
                 vegetation.TEMP.L_up = L_up./vegetation.STATVAR.area(1);
                 vegetation.STATVAR.Lin = Lin./vegetation.STATVAR.area(1);
                 vegetation.STATVAR.Lout = Lout./vegetation.STATVAR.area(1);
                 
                 %test lateral
                 exchange_area = min(snow_fraction, 1-snow_fraction).* vegetation.STATVAR.area;
                 snow_L_out = sigma.*(vegetation.CHILD.STATVAR.T + Tmfw).^4.* exchange_area;
                 vegetation_L_out = sigma.*(vegetation.STATVAR.T + Tmfw).^4.* exchange_area;
                 coupling_strength = 1;
                 snow_L_abs = coupling_strength .* (vegetation_L_out - snow_L_out);
                 vegetation_L_abs = coupling_strength .*(-vegetation_L_out + snow_L_out);
                 vegetation.TEMP.d_energy = vegetation.TEMP.d_energy + vegetation_L_abs;
                 vegetation.CHILD.TEMP.d_energy = vegetation.CHILD.TEMP.d_energy + snow_L_abs ;
             else
                 [vegetation, Lout] = penetrate_LW@VEGETATION_CLM5_seb(vegetation, Lin);
             end
        end
        
        function vegetation = canopy_resistances_CLM5_Stewart_snow(vegetation, tile)
            % similar to canopy_resistances_CLM5(), but with resistances 
            % against transpiration from Stewart (1988) as in Dingman (2015)     
            forcing = tile.FORCING;
            
            snow_fraction = vegetation.CHILD.STATVAR.area./vegetation.STATVAR.area;
            uz = forcing.TEMP.wind; % atm. wind speed
            p = forcing.TEMP.p; % surface pressure
            L = vegetation.STATVAR.LAI; % Leaf area index
            S = vegetation.STATVAR.SAI; % Stem area index

            q_s = vegetation.STATVAR.q_s; % canopy air specific humidity (last timestep)
            f_dry = vegetation.STATVAR.f_dry; % dry fraction of canopy
            f_wet = vegetation.STATVAR.f_wet; % wet fraction of canopy
            Tv = vegetation.STATVAR.T + forcing.CONST.Tmfw; % leaf temperature
            kappa = vegetation.CONST.kappa; % van Karman constant
            Lstar = vegetation.STATVAR.Lstar; % Monin-Ubukhov length
            Cv = vegetation.PARA.Cv; % Turbulent transfer coefficient canopy surface - canopy air
            d_leaf = vegetation.PARA.d_leaf; % characteristic dimension of the leaves in the direction of wind flow [m]
            ypsilon = vegetation.CONST.ypsilon; % kinematic viscosity of air
            Cs_dense = vegetation.PARA.Cs_dense; %  dense canopy turbulent transfer coefficient
            z = forcing.PARA.airT_height + sum(vegetation.STATVAR.layerThick); % height above ground for measurements
            z0v = vegetation.STATVAR.z0; % roughness length of vegetation ---------------- consider makeing a get_z0_surface_vegetation function!
            z0g = get_z0_surface(vegetation.NEXT); % Roughness lenght of ground
            d = vegetation.STATVAR.d; % displacement heigh
            k_s = vegetation.PARA.k_shelter; % canopy sheltering coefficient
            
            u_star = real(uz.*kappa./(log((z-d)./z0v)- psi_M_CLM5(vegetation, (z-d)./Lstar, z0v./Lstar)));
            
            e_v = double(Tv>=forcing.CONST.Tmfw).*satPresWater(vegetation,Tv) + double(Tv<forcing.CONST.Tmfw).*satPresIce(vegetation,Tv); % saturation water pressure of leafs
            qs_Tv = .622.*e_v./p;
            vegetation.TEMP.rho_atm = rho_air_moist(vegetation,tile);
            Cs_bare = kappa./.13.*(z0g.*u_star./ypsilon).^(-.45); % Eq.5.121 bare soil turbulent transfer coefficient
            W = exp(-(L+S));% Eq. 5.119
            Cs = Cs_bare.*W + Cs_dense.*(1-W); % Eq. 5.118 turbulent transfer coefficient between soil and canopy air
            vegetation.TEMP.C_leaf = leaf_conductance_Stewart(vegetation, tile);
            
            vegetation.TEMP.r_a = 1./(kappa^2.*uz).*(log((z-d)./z0v)- psi_M_CLM5(vegetation, (z-d)./Lstar, z0v./Lstar)).*(log((z-d)./z0v)- psi_H_CLM5(vegetation, (z-d)./Lstar, z0v./Lstar)); % aerodynamic resistance to heat/water vapor fransper canopy air - atmosphere. From original CG3 publication
            vegetation.TEMP.r_b = 1./Cv*(u_star./d_leaf).^(-.5); % Eq. 5.122 leaf boundary layer resistance
            vegetation.TEMP.r_a_prime = 1/(Cs.*u_star); % Eq. 5.116 aerodynamic resistance to heat/water vapor transfer soil - canopy air

            vegetation.TEMP.r_soil = double(get_humidity_surface(vegetation.IA_NEXT, tile) > q_s) .* ground_resistance_evap(vegetation.IA_NEXT, tile); % resistance to water vapor flux within the soil matrix
            if vegetation.STATVAR.LAI > 0 % Transpiring leaves present
                vegetation.TEMP.r_canopy = 1./(vegetation.TEMP.C_leaf.*k_s);   
            else % No leaves -> no transpiration
                vegetation.TEMP.r_canopy = Inf;
            end
            vegetation.TEMP.Ev_pot = -vegetation.TEMP.rho_atm.*(q_s - qs_Tv)./vegetation.TEMP.r_b;
            if vegetation.TEMP.Ev_pot <= 0 % Condensation -> on whole leaf + stem area, NO transpiration
                vegetation.TEMP.r_total = vegetation.TEMP.r_b./(2*(L+S));
            elseif vegetation.TEMP.Ev_pot > 0 % evaporation + transpiration
                vegetation.TEMP.r_total = vegetation.TEMP.r_b./(2*(L+S)) ./ ( f_wet + f_dry.*vegetation.TEMP.r_b./(vegetation.TEMP.r_b + vegetation.TEMP.r_canopy) ) ;  
            end
            vegetation.TEMP.r_total_snow = vegetation.TEMP.r_b./(2.*(L+S)); %snow only one-sided ?
        end
        
        % Sensible heat flux for canopy as in CLM5
        function vegetation = Q_h_CLM5_snow(vegetation, tile)
            forcing = tile.FORCING;
            
            snow_fraction = vegetation.CHILD.STATVAR.area./vegetation.STATVAR.area;
            Tv = (1-snow_fraction) .* vegetation.STATVAR.T + snow_fraction .* vegetation.CHILD.STATVAR.T + forcing.CONST.Tmfw; % leaf temperature

            cp = vegetation.CONST.cp;
            L = vegetation.STATVAR.LAI;
            S = vegetation.PARA.SAI;
%             Tv = vegetation.STATVAR.T(1)+forcing.CONST.Tmfw;
            Tz = forcing.TEMP.Tair+forcing.CONST.Tmfw; % air temperature in Kelvin
            z = forcing.PARA.airT_height; % height above ground for measurements
            z = z + sum(vegetation.STATVAR.layerThick); % adjust so forcing height is above canopy
            Gamma = vegetation.CONST.Gamma_dry; % negative of dry adiabatic lapse rate
            r_b = vegetation.TEMP.r_b;
            r_ah = vegetation.TEMP.r_a; % r_ah = r_aw, resistance to sensible heat transfer canopy air - atmosphere
            r_ah_prime = vegetation.TEMP.r_a_prime;
            rho_atm = vegetation.TEMP.rho_atm; % air density
            Tz_pot = Tz+Gamma.*z; % Eq. 5.7 potential temperature
            Tg = get_surface_T(vegetation.NEXT, tile); % ground/surface temperature
            Tg = Tg + forcing.CONST.Tmfw;
            
            Ts = (Tz_pot./r_ah + Tg./r_ah_prime + Tv.*2*(L+S)./r_b)./(1./r_ah + 1/r_ah_prime + 2*(L+S)./r_b); % Eq. 5.93 canopy air temperature
            
            vegetation.STATVAR.Ts = Ts - forcing.CONST.Tmfw;
            vegetation.STATVAR.Qh = - rho_atm.*cp.*(vegetation.STATVAR.Ts - vegetation.STATVAR.T).*2*(L+S)./r_b; % Eq. 5.88
            vegetation.CHILD.STATVAR.Qh = - rho_atm.*cp.*(vegetation.STATVAR.Ts - vegetation.CHILD.STATVAR.T).*2*(L+S)./r_b; % Eq. 5.88
            vegetation.TEMP.Tair = forcing.TEMP.Tair;
        end
        
        function vegetation = Q_e_CLM5_Stewart_snow(vegetation, tile)
            % evaporation and transpiration following CLM5, but with
            % simpler stomata/leaf resistance after Stewart (1988)
            % R. B. Zweigel, December 2021
            forcing = tile.FORCING;
            
            snow_fraction = vegetation.CHILD.STATVAR.area./vegetation.STATVAR.area;
            
            Tmfw = forcing.CONST.Tmfw;
            p = forcing.TEMP.p; % surface pressure
            Tsnow = vegetation.CHILD.STATVAR.T +forcing.CONST.Tmfw;
            e_v = double(Tsnow>=Tmfw).*satPresWater(vegetation,Tsnow) + double(Tsnow<Tmfw).*satPresIce(vegetation,Tsnow); % saturation water pressure of leafs
            qs_snow = .622.*e_v./p;
            
            L = vegetation.STATVAR.LAI; % Leaf area index
            S = vegetation.STATVAR.SAI; % Stem area index
            Tv = vegetation.STATVAR.T(1); % leaf temperature
            q_atm = forcing.TEMP.q; % atm. specific humidity
            rho_atm = vegetation.TEMP.rho_atm; % air density 
            f_wet = vegetation.STATVAR.f_wet; % wetted fraction of canopy
            f_dry = vegetation.STATVAR.f_dry; % dry and transpiring fraction of leaves
            
            % Resistance terms
            r_aw = vegetation.TEMP.r_a; % aerodynamic resistance against water vapor transport canopy air - atmosphere
            r_aw_prime = vegetation.TEMP.r_a_prime; % aerodynamic resistance against water vapor transport soil - canopy air
            r_b = vegetation.TEMP.r_b; % leaf boudnary layer resitance against vapor transport
            r_total = vegetation.TEMP.r_total;
            r_total_snow = vegetation.TEMP.r_total_snow;            
            r_soil = vegetation.TEMP.r_soil;
            r_canopy = vegetation.TEMP.r_canopy;
            
            % Humidity
            q_g = get_humidity_surface(vegetation.IA_NEXT, tile); % ground/snow surface specific humidity
            e_v = double(Tv>=0).*satPresWater(vegetation,Tv+Tmfw) + double(Tv<0).*satPresIce(vegetation,Tv+Tmfw); % saturation water pressure of leafs
            qs_Tv = .622.*e_v./p; % saturation water vapor specific humidity at leaf temperature
            ca = 1./r_aw;
            cv = (1-snow_fraction)./r_total;
            cv_snow = snow_fraction./r_total_snow;
            cg = 1./(r_aw_prime + r_soil);
            
            q_s = (q_atm.*ca + q_g.*cg + qs_Tv.*cv + qs_snow.*cv_snow)./(ca + cv + cv_snow + cg); % Eq. 5.108, canopy air specific humidity
            
            vegetation.STATVAR.q_s = q_s;
            vegetation.STATVAR.qs_TV = qs_Tv;
            vegetation.STATVAR.qs_snow = qs_snow;
            vegetation.TEMP.q_atm = q_atm;
            
            water_fraction = vegetation.STATVAR.water(1) ./ vegetation.STATVAR.waterIce(1);
            ice_fraction = vegetation.STATVAR.ice(1) ./ vegetation.STATVAR.waterIce(1);
            water_fraction(vegetation.STATVAR.waterIce == 0) = double(Tv>=0);
            ice_fraction(vegetation.STATVAR.waterIce == 0) = double(Tv<0);
            
            if q_s - qs_Tv < 0 % evaporation/transpiration in m3 water /(m2*s), not kg/(m2*s) as in CLM5
                vegetation.STATVAR.evap = -rho_atm.*(q_s - qs_Tv).* f_wet.*water_fraction.*2.*(L+S) ./r_b ./vegetation.CONST.rho_w ;
                vegetation.STATVAR.sublim = -rho_atm.*(q_s - qs_Tv).* f_wet.*ice_fraction.*2.*(L+S) ./r_b ./vegetation.CONST.rho_w ; 
                vegetation.STATVAR.transp = -rho_atm.*(q_s - qs_Tv).* f_dry.*2.*(L+S) ./ (r_b + r_canopy) ./ vegetation.CONST.rho_w ; 
            else % condensation/deposition
                vegetation.STATVAR.evap = -rho_atm.*(q_s - qs_Tv).* water_fraction.*2.*(L+S) ./r_b ./vegetation.CONST.rho_w;
                vegetation.STATVAR.sublim = -rho_atm.*(q_s - qs_Tv).* ice_fraction.*2.*(L+S) ./r_b ./vegetation.CONST.rho_w;
                vegetation.STATVAR.transp = 0;
            end
            
            vegetation.TEMP.Qe_evap = latent_heat_evaporation(vegetation, Tv+Tmfw).*vegetation.STATVAR.evap.*vegetation.CONST.rho_w + latent_heat_sublimation(vegetation, Tv+Tmfw).*vegetation.STATVAR.sublim.*vegetation.CONST.rho_w;
            vegetation.TEMP.Qe_transp = latent_heat_evaporation(vegetation, Tv+Tmfw).*vegetation.STATVAR.transp .* vegetation.CONST.rho_w;
            
            vegetation.STATVAR.Qe = vegetation.TEMP.Qe_evap + vegetation.TEMP.Qe_transp;
            vegetation.STATVAR.evap_energy = vegetation.STATVAR.evap.*( double(Tv>=0).*vegetation.CONST.c_w.*Tv + double(Tv<0).*vegetation.CONST.c_i.*Tv);
            vegetation.STATVAR.sublim_energy = vegetation.STATVAR.sublim .* (vegetation.CONST.c_i.*Tv - vegetation.CONST.L_f);
            vegetation.STATVAR.transp_energy = vegetation.STATVAR.transp.*vegetation.CONST.c_w.*Tv; 
            
            %snow CHILD, only sublimation for now, sign convention as for
            %snow, sublimation proportional to d_water
            vegetation.CHILD.STATVAR.sublimation = rho_atm.*(q_s - qs_snow) .* 2.* (L+S) ./r_b ./vegetation.CONST.rho_w .* vegetation.CHILD.STATVAR.area;
            vegetation.CHILD.TEMP.sublimation_energy = vegetation.CHILD.STATVAR.sublimation .* (vegetation.CONST.c_i.*vegetation.CHILD.STATVAR.T - vegetation.CONST.L_f);
            vegetation.CHILD.STATVAR.Qe = -latent_heat_sublimation(vegetation, vegetation.CHILD.STATVAR.T + Tmfw).*vegetation.CHILD.STATVAR.sublimation.*vegetation.CONST.rho_w ./ vegetation.CHILD.STATVAR.area;
        end
        
        
        function vegetation = get_boundary_condition_l(vegetation, tile)
              vegetation = get_boundary_condition_l@VEGETATION_CLM5_seb(vegetation, tile);
        end
        
        function vegetation = get_derivatives_prognostic(vegetation, tile)
            if vegetation.CHILD == 0  
                vegetation = get_derivatives_prognostic@VEGETATION_CLM5_seb(vegetation, tile); %call normal function
            else
                vegetation.CHILD = get_derivatives_prognostic_CHILD(vegetation.CHILD, tile);
                vegetation = get_derivatives_prognostic@VEGETATION_CLM5_seb(vegetation, tile); 
            end
        end
        
        function timestep = get_timestep(vegetation, tile) 
            if vegetation.CHILD == 0
                timestep = get_timestep@VEGETATION_CLM5_seb(vegetation, tile);
            else 
                timestep_snow = get_timestep_CHILD(vegetation.CHILD, tile);
                timestep_vegetation =  get_timestep@VEGETATION_CLM5_seb(vegetation, tile);
                timestep = timestep_vegetation + double(timestep_snow > 0 && timestep_snow < timestep_vegetation) .* (timestep_snow - timestep_vegetation);
            end
        end
        
        function vegetation = advance_prognostic(vegetation, tile) 
            if vegetation.CHILD == 0
                vegetation =  advance_prognostic@VEGETATION_CLM5_seb(vegetation, tile);
            else                
                vegetation.CHILD = advance_prognostic_CHILD(vegetation.CHILD, tile);
                vegetation =  advance_prognostic@VEGETATION_CLM5_seb(vegetation, tile);
            end
        end
        
        function vegetation = compute_diagnostic_first_cell(vegetation, tile)
            forcing = tile.FORCING;
            vegetation = compute_diagnostic_first_cell@VEGETATION_CLM5_seb(vegetation, tile);
        end
        
        function vegetation = compute_diagnostic(vegetation, tile)
            if vegetation.CHILD == 0
                vegetation = compute_diagnostic@VEGETATION_CLM5_seb(vegetation, tile);
            else
                vegetation.CHILD = compute_diagnostic_CHILD(vegetation.CHILD, tile); %requires additional function for snow "drip" and water drio from the snow CHILD
                %in IA class between snow CHILD and canopy, there is
                %funtion remove_excess_water, this routes excess
                %meltwatre to canopy water and then drip from there will go
                %to the class below - check how this is handled!
                if vegetation.CHILD ~= 0
                    remove_excessWater_CHILD(vegetation.IA_CHILD); 
                    vegetation.IA_NEXT = canopy_snow_drip(vegetation.IA_NEXT, tile);
                end
                vegetation = compute_diagnostic@VEGETATION_CLM5_seb(vegetation, tile);  %this function is fine, handles drip from canopy itself
                
                %handle snow drip
%                 if vegetation.CHILD ~= 0
%                 end
            end
        end
        
        function vegetation = check_trigger(vegetation, tile)
            if vegetation.CHILD ~= 0
                %delete CHILD
                if vegetation.CHILD.STATVAR.area ./ vegetation.STATVAR.area(1,1) < 1e-6 %cutoff to get rid of remaining snow
                   vegetation.CHILD = 0;
                   vegetation.IA_CHILD = 0;
                %make SNOW CHILD full class   
                elseif vegetation.CHILD.STATVAR.area ./ vegetation.STATVAR.area(1,1) > 1 
                    %transforms dimensions and STAVAR
                    snow_volume = vegetation.CHILD.STATVAR.area .* vegetation.CHILD.STATVAR.layerThick;
                    vegetation.CHILD.STATVAR.area = vegetation.STATVAR.area(1,1);
                    vegetation.CHILD.STATVAR.layerThick = snow_volume ./ vegetation.CHILD.STATVAR.area;
                   
                    %make snow a real class
                    vegetation.CHILD.PARENT = 0;
                    vegetation.CHILD.PREVIOUS = vegetation.PREVIOUS;
                    vegetation.CHILD.NEXT = vegetation;
                    vegetation.PREVIOUS.NEXT = vegetation.CHILD;
                    ia_class = get_IA_class(class(vegetation.PREVIOUS), class(vegetation.CHILD));
                    vegetation.PREVIOUS.IA_NEXT = ia_class;
                    vegetation.CHILD.IA_PREVIOUS = ia_class;
                    vegetation.CHILD.IA_PREVIOUS.NEXT = vegetation.CHILD;
                    vegetation.CHILD.IA_PREVIOUS.PREVIOUS = vegetation.PREVIOUS;
                    finalize_init(vegetation.CHILD.IA_PREVIOUS, tile);
                    
                    vegetation.PREVIOUS = vegetation.CHILD;
                    vegetation.CHILD = 0;
                    vegetation.IA_PREVIOUS = vegetation.IA_CHILD; 
                    vegetation.PREVIOUS.IA_NEXT = vegetation.IA_CHILD;
                    vegetation.IA_CHILD = 0;
                    
                end
            end
        end
        

        %----------
        %reset timestamp when changing TILES
        function vegetation = reset_timestamps(vegetation, tile)
            if vegetation.CHILD ~= 0
                vegetation.CHILD = reset_timestamps(vegetation.CHILD, tile);
            end
        end
        
        
        
        
        %-------------param file generation-----
        function canopy = param_file_info(canopy)
            canopy = provide_PARA(canopy);

            canopy.PARA.STATVAR = [];
            canopy.PARA.class_category = 'VEGETATION';
            canopy.PARA.options = [];

            canopy.PARA.comment.max_snow_fraction = {'maximum aereal fraction of canopy that can be covered by snow; no interception occurs if set to zero'};
            canopy.PARA.default_value.max_snow_fraction = {0.5};
            
            canopy.PARA.comment.SWE_max = {'maximum SWE in [m] within canopy; 0.01m corresponds to 5-10cm of snow'};
            canopy.PARA.default_value.SWE_max = {0.01};
            
            canopy.PARA.comment.t_leafsprout = {'DayOfYear when leaves emerge (135 = 15. May)'};
            canopy.PARA.default_value.t_leafsprout = {'152'};

            canopy.PARA.comment.t_leaffall = {'DayOfYear when leaves fall (274 = 1. Oct.)'};
            canopy.PARA.default_value.t_leaffall = {'274'};

            canopy.PARA.comment.kv = {'parameter to adjust for non-cylinderness of trees'};
            canopy.PARA.default_value.kv = {0.5};

            canopy.PARA.comment.D_bh = {'mean breast-height diameter of trees'};
            canopy.PARA.default_value.D_bh = {0.28};

            canopy.PARA.comment.N_tree = {'stand density (number of trees per area)'};
            canopy.PARA.default_value.N_tree = {0.15};

            canopy.PARA.comment.rho_wood = {'Density of dry wood'};
            canopy.PARA.default_value.rho_wood = {500};

            canopy.PARA.comment.SLA = {'Specific leaf area (leaf area per unit mass of carbon)'};
            canopy.PARA.default_value.SLA = {23};
            
            canopy.PARA.comment.f_carbon = {'mass of carbon per leaf dry mass'};
            canopy.PARA.default_value.f_carbon = {0.5};

            canopy.PARA.comment.f_water = {'mass of water per leaf total mass'};
            canopy.PARA.default_value.f_water = {0.45};

            canopy.PARA.comment.nir_fraction = {'fraction of Sin that is in the NIR spectrum'};
            canopy.PARA.default_value.nir_fraction = {0.3};

            canopy.PARA.comment.LAI = {'Leaf area index'};
            canopy.PARA.default_value.LAI = {1.5};

            canopy.PARA.comment.SAI = {'Stem area index'};
            canopy.PARA.default_value.SAI = {1};

            canopy.PARA.comment.Khi_L = {'departure of leaf angles from a random distribution'};
            canopy.PARA.default_value.Khi_L = {0.01};
           
            canopy.PARA.comment.alpha_leaf_vis = {'leaf reflectence in the VIS'};
            canopy.PARA.default_value.alpha_leaf_vis = {0.07};

            canopy.PARA.comment.alpha_stem_vis = {'stem reflectence in the VIS'};
            canopy.PARA.default_value.alpha_stem_vis = {0.16};

            canopy.PARA.comment.tau_leaf_vis = {'leaf transmittance in the VIS'};
            canopy.PARA.default_value.tau_leaf_vis = {0.05};

            canopy.PARA.comment.tau_stem_vis = {'stem transmittance in the VIS'};
            canopy.PARA.default_value.tau_stem_vis = {0.001};

            canopy.PARA.comment.alpha_leaf_nir = {'leaf reflectence in the NIR'};
            canopy.PARA.default_value.alpha_leaf_nir = {0.35};

            canopy.PARA.comment.alpha_stem_nir = {'stem reflectence in the NIR'};
            canopy.PARA.default_value.alpha_stem_nir = {0.39};

            canopy.PARA.comment.tau_leaf_nir = {'leaf transmittance in the NIR'};
            canopy.PARA.default_value.tau_leaf_nir = {0.1};

            canopy.PARA.comment.tau_stem_nir = {'stem transmittance in the NIR'};
            canopy.PARA.default_value.tau_stem_nir = {0.001};

            canopy.PARA.comment.dT_max = {'max temperature change per timestep [K]'};
            canopy.PARA.default_value.dT_max = {1};

            canopy.PARA.comment.dt_max = {'maximum timestep [sec]'};
            canopy.PARA.default_value.dt_max = {3600};

            canopy.PARA.comment.Cv = {'turbulent transfer coefficient between canopy surface and canopy air'};
            canopy.PARA.default_value.Cv = {0.01};

            canopy.PARA.comment.d_leaf = {'characteristic dimension of leaves in direction of wind flow'};
            canopy.PARA.default_value.d_leaf = {0.04};
            
            canopy.PARA.comment.Cs_dense = {'dense canopy turbulent transfer coefficient (Dickinson et al. 1993)'};
            canopy.PARA.default_value.Cs_dense = {0.004};

            canopy.PARA.comment.R_z0 = {'ratio of momentum roughness length'};
            canopy.PARA.default_value.R_z0 = {0.055};

            canopy.PARA.comment.R_d = {'ratio of displacement height to canopy height'};
            canopy.PARA.default_value.R_d = {0.67};
            
            canopy.PARA.comment.Dmax = {'Maximum dry soil layer thickness, 15 mm'};
            canopy.PARA.default_value.Dmax = {0.015};

            canopy.PARA.comment.Wmax = {'water holding capacity of canopy per unit area (L+S)'}; 
            canopy.PARA.default_value.Wmax = {0.0001}; 
            
            canopy.PARA.comment.beta_root = {'root distribution parameter (CLM5)'};
            canopy.PARA.default_value.beta_root = {0.943};

            canopy.PARA.comment.C_leaf_max = {'maximum leaf conductance'};
            canopy.PARA.default_value.C_leaf_max = {0.01};

            canopy.PARA.comment.k_shelter = {'shelter factor, between 0.5 - 1 (Carlson, 1991)'};
            canopy.PARA.default_value.k_shelter = {0.5};

            canopy.PARA.comment.psi_wilt = {'water potential at permanent wilting point (-15 bar)'};
            canopy.PARA.default_value.psi_wilt = {-150};

            canopy.PARA.comment.zeta_m = {'threshold for very unstable atm. Conditions wrt. Heat/vapor fluxes'};
            canopy.PARA.default_value.zeta_m = {-0.465};

            canopy.PARA.comment.zeta_h = {'threshold for very unstable atm. Conditions wrt. mass fluxes'};
            canopy.PARA.default_value.zeta_h = {-1.574};

        end
    end
end