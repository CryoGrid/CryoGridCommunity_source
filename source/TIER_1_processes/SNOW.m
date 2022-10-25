%========================================================================
% CryoGrid TIER1 library class for functions related to snow cover
% representation
% S. Westermann, R. Zweigel, October 2020
%========================================================================

classdef SNOW < BASE
    
    methods
        
        function snow = initialize_zero_snow_BASE(snow)
            snow.STATVAR.T = 0;
            snow.STATVAR.energy = 0;
            snow.STATVAR.waterIce = 0;
            snow.STATVAR.mineral = 0;
            snow.STATVAR.organic = 0;
            snow.STATVAR.layerThick = 0;
            snow.STATVAR.ice=0;
            snow.STATVAR.water=0;
            snow.STATVAR.area = 0;
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
        end
        
%-----------boundary conditions----------------
        function snow = get_boundary_condition_SNOW_u(snow, forcing) %snow and rain from entire grid cell area
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600) .* snow.STATVAR.area(1,1); %snowfall is in mm/day -> [m3/sec]
            snow.TEMP.rainfall = forcing.TEMP.rainfall ./1000 ./(24.*3600) .* snow.STATVAR.area(1,1);
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);  %[J/sec]
            snow.TEMP.rain_energy = snow.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;
        end
        
        %used in CHILD phase
        function snow = get_boundary_condition_allSNOW_rain_u(snow, forcing) %snow from entire area (PARENT+CHILD), rain only from snow-covered part (CHILD)
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600) .* (snow.PARENT.STATVAR.area(1,1) + snow.STATVAR.area); %snowfall is in mm/day -> [m3/sec]
            snow.TEMP.rainfall = forcing.TEMP.rainfall ./1000 ./(24.*3600) .* snow.STATVAR.area;
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);  %[J/sec]
            snow.TEMP.rain_energy = snow.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;
        end
        
        function snow = get_boundary_condition_allSNOW_rain_canopy_m(snow, tile) %snow from entire area (PARENT+CHILD), rain only from snow-covered part (CHILD)
            forcing = tile.FORCING;
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600) .* (snow.PARENT.STATVAR.area(1,1) + snow.STATVAR.area); %snowfall is in mm/day -> [m3/sec]
            snow.TEMP.rainfall = snow.PARENT.PREVIOUS.TEMP.rain_thru .* snow.STATVAR.area;
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);  %[J/sec]
            snow.TEMP.rain_energy = snow.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;
        end
        
        %when creating CHILD
        function snow = get_boundary_condition_allSNOW_u(snow, forcing) %only snow from entire area, i.e. from PARENT
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600) .* snow.PARENT.STATVAR.area(1,1); %snowfall is in mm/day -> [m3/sec]
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);  %[J/sec]
        end
        
        function snow = get_sublimation(snow, forcing)
            snow.STATVAR.sublim = -snow.STATVAR.Qe ./(snow.CONST.rho_w .* snow.CONST.L_s) .* snow.STATVAR.area(1);
            snow.TEMP.sublim_energy = snow.STATVAR.sublim .* (snow.STATVAR.T(1) .* snow.CONST.c_i - snow.CONST.L_f);
        end
        
%--------timesteps----
        function timestep = get_timestep_SNOW_mass_balance(snow) %at maximum timestep the entire cell is melted
            timestep = min((-snow.STATVAR.energy ./ snow.TEMP.d_energy) .*double(snow.TEMP.d_energy>0) + double(snow.TEMP.d_energy<=0).*snow.PARA.dt_max); %when snow is melting, do not melt more than there is in a grid cell
            timestep(isnan(timestep)) = snow.PARA.dt_max;
        end
        
        function timestep = get_timestep_SNOW(snow) %at maximum timestep maximum half the cell is melted under melting conditions, otherwise normal heat conduction timestep
%             timestep = min(double(snow.TEMP.d_energy>0 & snow.STATVAR.T >=0) .*0.5 .*  (-snow.STATVAR.energy ./ snow.TEMP.d_energy) + ...
%                 double(snow.TEMP.d_energy<=0 | snow.STATVAR.T < 0) .* snow.PARA.dE_max ./ (abs(snow.TEMP.d_energy) ./ snow.STATVAR.layerThick./ snow.STATVAR.area));
            timestep = min(double(snow.TEMP.d_energy>0 & snow.STATVAR.energy >= - 0.95.* snow.STATVAR.ice.*snow.CONST.L_f) .*0.5 .*  (-snow.STATVAR.energy ./ snow.TEMP.d_energy) + ...
                double(snow.TEMP.d_energy<0 & snow.STATVAR.energy >= - 0.95.* snow.STATVAR.ice.*snow.CONST.L_f) .*0.5 .*  ((-snow.STATVAR.ice.*snow.CONST.L_f -snow.STATVAR.energy) ./ snow.TEMP.d_energy) + ...
                double( snow.STATVAR.energy < - 0.95.* snow.STATVAR.ice.*snow.CONST.L_f) .* snow.CONST.c_i .* (snow.STATVAR.ice./snow.STATVAR.layerThick./ snow.STATVAR.area) .*0.5 ./ (abs(snow.TEMP.d_energy) ./ snow.STATVAR.layerThick./ snow.STATVAR.area));
             
            timestep(isnan(timestep)) = snow.PARA.dt_max;
        end
        
        function timestep = get_timestep_SNOW_CHILD(snow)
            %timestep = (-snow.STATVAR.energy ./ snow.TEMP.d_energy) .* double(snow.TEMP.d_energy>0) + double(snow.TEMP.d_energy<=0).*snow.PARA.dt_max;
            
%             timestep = double(snow.TEMP.d_energy>0 & snow.STATVAR.T >=0) .*0.5 .*  (-snow.STATVAR.energy ./ snow.TEMP.d_energy) + ...
%                 double(snow.TEMP.d_energy<=0 | snow.STATVAR.T < 0) .* snow.PARA.dE_max ./ (abs(snow.TEMP.d_energy) ./ snow.STATVAR.layerThick./ snow.STATVAR.area);
            timestep = double(snow.TEMP.d_energy>0 & snow.STATVAR.energy >= - 0.95) .* snow.STATVAR.ice.*snow.CONST.L_f .*0.5 .*  (-snow.STATVAR.energy ./ snow.TEMP.d_energy) + ...
                double(snow.TEMP.d_energy<0 & snow.STATVAR.energy >= - 0.95.* snow.STATVAR.ice.*snow.CONST.L_f) .*0.5 .*  ((-snow.STATVAR.ice.*snow.CONST.L_f -snow.STATVAR.energy) ./ snow.TEMP.d_energy) + ...
                double( snow.STATVAR.energy < - 0.95.* snow.STATVAR.ice.*snow.CONST.L_f) .* snow.CONST.c_i .* (snow.STATVAR.ice./snow.STATVAR.layerThick./ snow.STATVAR.area) .*0.5 ./ (abs(snow.TEMP.d_energy) ./ snow.STATVAR.layerThick./ snow.STATVAR.area);
         
         %a=double(snow.TEMP.d_energy>0 & snow.STATVAR.T >=0) .*0.5 .*  (-snow.STATVAR.energy ./ snow.TEMP.d_energy);
         %b =double(snow.TEMP.d_energy<=0 | snow.STATVAR.T < 0) .* snow.PARA.dE_max ./ (abs(snow.TEMP.d_energy) ./ snow.STATVAR.layerThick./ snow.STATVAR.area);
         
         
            timestep(isnan(timestep)) = snow.PARA.dt_max;
        end
        
        function timestep = get_timestep_SNOW_sublimation(snow) 
            timestep = double(snow.TEMP.sublim_energy > 0) .*0.25 .*  (-snow.STATVAR.energy(1,1) ./ snow.TEMP.sublim_energy) +  double(snow.TEMP.sublim_energy <= 0) .* snow.PARA.dt_max;
             
                
            timestep(isnan(timestep)) = snow.PARA.dt_max;
        end
        
        
%--------diagnostic step------------
        %remove or reroute meltwater in excess of snow matrix - check if
        %subtract_water is redundant
        function snow = subtract_water(snow)
            snow.STATVAR.layerThick = min(snow.STATVAR.layerThick, snow.STATVAR.ice ./ snow.STATVAR.target_density ./snow.STATVAR.area); %adjust so that old density is maintained; do not increase layerThick (when water refreezes)
            difference = max(0, snow.STATVAR.waterIce - snow.STATVAR.layerThick .*  snow.STATVAR.area);
            snow.STATVAR.waterIce = min(snow.STATVAR.layerThick .*  snow.STATVAR.area, snow.STATVAR.waterIce); % Remove water that is in excess of cell volume (drains water out of the system)
            snow.STATVAR.water = max(0, min(snow.STATVAR.water, snow.STATVAR.waterIce - snow.STATVAR.ice));
            snow.STATVAR.excessWater = snow.STATVAR.excessWater + sum(difference);
        end
        
        function snow = subtract_water2(snow)
            snow.STATVAR.layerThick = min(snow.STATVAR.layerThick, snow.STATVAR.ice ./ snow.STATVAR.target_density ./snow.STATVAR.area); %adjust so that old density is maintained; do not increase layerThick (when water refreezes)
            difference = max(0, snow.STATVAR.waterIce - snow.STATVAR.layerThick .*  snow.STATVAR.area);
            snow.STATVAR.waterIce = snow.STATVAR.waterIce - difference;
            snow.STATVAR.water = max(0, snow.STATVAR.waterIce - snow.STATVAR.ice);
            snow.STATVAR.excessWater = snow.STATVAR.excessWater + sum(difference);
        end
        
        
        function snow = subtract_water_CHILD(snow)
            snow.STATVAR.area = min(snow.STATVAR.area, snow.STATVAR.ice ./ snow.STATVAR.target_density ./snow.STATVAR.layerThick); %adjust so that old density is maintained; do not increase layerThick (when water refreezes)
            difference = snow.STATVAR.waterIce;  
            snow.STATVAR.waterIce = min(snow.STATVAR.layerThick .* snow.STATVAR.area, snow.STATVAR.waterIce); % Remove water that is in excess of cell volume (drains water out of the system)
            snow.STATVAR.water = max(0, min(snow.STATVAR.water, snow.STATVAR.waterIce - snow.STATVAR.ice));
            difference = difference - snow.STATVAR.waterIce;
            snow.STATVAR.excessWater = snow.STATVAR.excessWater + difference;
        end
        
        function snow = subtract_water_CHILD2(snow)  %unclear if and where needed
            snow.STATVAR.area = min(snow.STATVAR.area, snow.STATVAR.ice ./ snow.STATVAR.target_density ./snow.STATVAR.layerThick); %adjust so that old density is maintained; do not increase layerThick (when water refreezes)
            difference = max(0, snow.STATVAR.waterIce - snow.STATVAR.layerThick .*  snow.STATVAR.area);
            snow.STATVAR.waterIce = max(0, snow.STATVAR.waterIce - difference);
            snow.STATVAR.water = max(0, min(snow.STATVAR.water, snow.STATVAR.waterIce - snow.STATVAR.ice));
            snow.STATVAR.excessWater = snow.STATVAR.excessWater + sum(difference);
            snow.STATVAR.excessWater_energy = double(snow.STATVAR.energy > 0) .* snow.STATVAR.energy; %rare case that snow is all liquid, only possible through lateral flow.
            snow.STATVAR.energy =  snow.STATVAR.energy - snow.STATVAR.excessWater_energy;
            snow.STATVAR.excessWater_energy = sum(snow.STATVAR.excessWater_energy,1);
        end
        
        
        %albedo parametrizations
        function snow = calculate_albedo_simple(snow, timestep)  %albedo formulation used in the "old" CryoGrid 3
            if snow.STATVAR.T(1) == 0    % melting conditions
                snow.STATVAR.albedo = snow.PARA.min_albedo + (snow.STATVAR.albedo - snow.PARA.min_albedo) .* exp(-snow.PARA.tau_f .* timestep ./ snow.PARA.tau_1);
            else
                snow.STATVAR.albedo = max(snow.STATVAR.albedo - snow.PARA.tau_a .* timestep ./ snow.PARA.tau_1, snow.PARA.min_albedo);
            end
            
            %reset albedo after snowfall
            snow.TEMP.newSnow =  snow.TEMP.newSnow - snow.TEMP.newSnow.*0.1./24./3600.*timestep + snow.TEMP.snowfall .*timestep ./ snow.STATVAR.area(1);
            if snow.TEMP.newSnow >= 0.001
                snow.STATVAR.albedo = snow.PARA.max_albedo;
                snow.TEMP.newSnow = 0;
            end
        end
        
        
        function snow = calculate_albedo_crocus(snow, forcing)  %find better way to include the pressure? Maybe include forcing? Or altitude
            
            
            rho=snow.STATVAR.waterIce./snow.STATVAR.layerThick ./snow.STATVAR.area .*1000; % WET DENSITY ?? Check Vionet 2012/Sebastian
            d=snow.STATVAR.d;
            s=snow.STATVAR.s;
            gs=snow.STATVAR.gs;
            %snow_age = forcing.TEMP.t - snow.STATVAR.time_snowfall(1);
            snow_age = forcing.TEMP.t - snow.STATVAR.top_snow_date(1,1);
            %spectral_ranges = [0.71; 0.21; 0.08];
            
            d_opt = double(d>0).*1e-4.*(d+(1-d).*(4-s)) + double(d==0).*(gs.*s+(1-s).*max(gs.*0+4e-4, gs/2));
            
            %albedo1 = max(0.6, min(0.92,0.96-1.58.*(d_opt(1).^0.5)) - min(1,max(forcing.TEMP.p./870e2,0.5)).*0.2./60.*snow_age);
            %albedo1 = max(0.6, min(0.92,0.96-1.58.*(d_opt(1).^0.5)) - 0.0.*0.2./60.*snow_age);
            albedo1 = max(0.6, min(0.92,0.96-1.58.*(d_opt(1).^0.5)) - snow.PARA.albedo_age_factor.*snow_age);
            albedo2 = max(0.3, 0.9-15.4.*(d_opt(1).^(0.5)));
            albedo3 = 346.3.*min(d_opt(1), 0.0023) - 32.1.*((min(d_opt(1), 0.0023)).^0.5) + 0.88;
            
            spectral_albedo = [albedo1 albedo2 albedo3];
            
            beta1 = max(40+d_opt.*0, 0.00192.*rho./(d_opt.^0.5));
            beta2 = max(100, 0.01098.*rho./(d_opt.^0.5));
            
            SW_extinction = [beta1 beta2 beta2.*Inf];
            
            total_albedo = sum(spectral_albedo.*snow.PARA.spectral_ranges);
            
            snow.STATVAR.albedo = total_albedo;
            snow.TEMP.spectral_albedo = spectral_albedo;
            snow.TEMP.SW_extinction = SW_extinction;
        end
        
        %properties of new snow CROCUS, Vionnet et al., 2012
        function snow = get_snow_properties_crocus(snow, forcing)
            if snow.TEMP.snowfall >0
                T_air=min(forcing.TEMP.Tair, 0);
                windspeed = forcing.TEMP.wind;
                
                T_fus=0;  %degree C
                a_rho=109;  %kg/m3
                b_rho=6; % kg/m3K
                %c_rho=26; % kg m7/2s1/2
                c_rho = snow.PARA.wind_factor_fresh_snow; 
                min_snowDensity=50; %kg/m3
                
                snow.TEMP.newSnow.STATVAR.density = max(min_snowDensity, a_rho+b_rho.*(T_air-T_fus)+ c_rho.*windspeed.^0.5);  %initial snow density
                %test Sebastian
                snow.TEMP.newSnow.STATVAR.density = snow.TEMP.newSnow.STATVAR.density .*1000 ./920;
                
                snow.TEMP.newSnow.STATVAR.d = min(max(1.29-0.17.*windspeed,0.2),1);
                snow.TEMP.newSnow.STATVAR.s = min(max(0.08.*windspeed + 0.38,0.5),0.9);
                snow.TEMP.newSnow.STATVAR.gs = 0.1e-3+(1-snow.TEMP.newSnow.STATVAR.d).*(0.3e-3-0.1e-3.*snow.TEMP.newSnow.STATVAR.s);
                snow.TEMP.newSnow.STATVAR.time_snowfall = forcing.TEMP.t;
                
                %new Sebastian
                snow.TEMP.newSnow.STATVAR.top_snow_date = forcing.TEMP.t;
                snow.TEMP.newSnow.STATVAR.bottom_snow_date = forcing.TEMP.t; % a few minites earluer
            end
        end
        
        function snow = get_snow_properties_vanKampenhout(snow, forcing)
            % Snow properties function similar to
            % get_snow_properties_crocus(...), but with new snow density
            % according to van Kampenhout et al. (2017) https://doi.org/10.1002/2017MS000988
            % R.B. Zweigel, July 2022
            
            if snow.TEMP.snowfall >0
                windspeed = forcing.TEMP.wind;
                T_air = forcing.TEMP.Tair;
                
                T_fus=0;  %degree C
                rho_Tair = double(T_air > T_fus+2).*( 50 + 1.7*17^(3/2) ) ...
                    + double(T_air > T_fus-15 & T_air <= T_fus+2).*( 50 + 1.7*(T_air-T_fus+15)^(3/2) ) ...
                    + double(T_air <= T_fus-15).*( -3.8328*(T_air-T_fus) - 0.0333*(T_air-T_fus)^2);
                rho_wind = 266.861*(.5*(1+tanh(windspeed/5)))^8.8;
                snow.TEMP.newSnow.STATVAR.density = rho_Tair + rho_wind;
                snow.TEMP.newSnow.STATVAR.density = snow.TEMP.newSnow.STATVAR.density .*1000 ./920;
                
                snow.TEMP.newSnow.STATVAR.d = min(max(1.29-0.17.*windspeed,0.2),1);
                snow.TEMP.newSnow.STATVAR.s = min(max(0.08.*windspeed + 0.38,0.5),0.9);
                snow.TEMP.newSnow.STATVAR.gs = 0.1e-3+(1-snow.TEMP.newSnow.STATVAR.d).*(0.3e-3-0.1e-3.*snow.TEMP.newSnow.STATVAR.s);
                snow.TEMP.newSnow.STATVAR.time_snowfall = forcing.TEMP.t;
            end
        end
        
        function snow = get_snow_properties_none(snow, forcing)
            % Function for snow classes that do not feature any dynamic new
            % snow properties (density, dendricity etc.)
        end
        
        %crocus snow microphycics, Vionnet et al., 2012
        function snow = get_T_gradient_snow(snow)
            
            delta=[snow.STATVAR.layerThick; snow.NEXT.STATVAR.layerThick(1,1)];
            T = [snow.STATVAR.T; snow.NEXT.STATVAR.T(1,1)];
            dT = (T(1:end-2) - T(3:end))./(0.5.*delta(1:end-2) +  delta(2:end-1) + 0.5.*delta(3:end));
            dT =abs([ (T(1) - T(2))./(0.5.*delta(1) + 0.5.*delta(2)) ; dT]);
            snow.TEMP.dT = dT;
        end
        
        function snow = get_T_gradient_snow_single_cell(snow) %when SNOW is a sigle cell
            
            delta=[snow.STATVAR.layerThick; snow.NEXT.STATVAR.layerThick(1,1)];
            T = [snow.STATVAR.T; snow.NEXT.STATVAR.T(1,1)];
            dT = abs((T(1) - T(2))./(0.5.*delta(1) + 0.5.*delta(2)));
            snow.TEMP.dT = dT;
        end
        
        function snow = prog_metamorphism(snow)
            T = snow.STATVAR.T;
            D = snow.STATVAR.layerThick;
            D_water = snow.STATVAR.water ./ snow.STATVAR.area;
            dT = snow.TEMP.dT;
            d = snow.STATVAR.d;
            s = snow.STATVAR.s;
            gs = snow.STATVAR.gs;
            
            daysec = 60*60*24;
            rho = max(50, snow.STATVAR.waterIce ./ (snow.STATVAR.layerThick .* snow.STATVAR.area) .*920);
            %rho = max(50, snow.STATVAR.waterIce ./ (snow.STATVAR.layerThick .* snow.STATVAR.area) .*1000);
            W_liq=D_water.*920; %W_liq is mass density, and the densities cancel out in theta_cubed if ice density is assumed here
            %W_liq=D_water.*1000;
            
            small_gradient=(dT<=5);
            gradient_term=double(small_gradient)+double(~small_gradient).*dT.^0.4;
            a=2e8.*exp(-6000./(T+273.15)).*gradient_term;
            
            %new Sebastian
            %D=max(D, D.*0+5e-3);
            theta_cubed=(min(10, 100.*W_liq./rho./D)).^3; %limit to 10%, Brun 1989, Annals Glac. Fig. 6;
            
            % Wet snow metamporphism
            d_d_wet = -1/16.*theta_cubed;
            d_s_wet = +1/16.*theta_cubed;
            d_gs_wet= double(s>=1 & d<=0) .*  2./3.14./gs.^2 .* (1.28e-17 + 4.22e-19.*theta_cubed);
            
            % Dry snow metamporhism
            d_d_dry = -a;
            d_s_dry = 5.*a.*double(small_gradient)- a.*double(~small_gradient);
            
            f=double(T<-40) + double(T>=-40 & T<-22).*0.011.*(T+40) +  double(T>=-22 & T<-6).*(0.2+0.05.*(T+22)) + double(T>=-6).*(1-0.05.*T);
            h=max(rho.*0, min(rho.*0+1, 1-0.004.*(rho-150)));
            g=min(dT.*0+1, double(dT>=15 & dT<25).*0.01.*(dT-15) + double(dT>=25 & dT<40).*(0.1+0.037.*(dT-25)) + double(dT>=40 & dT<50).*(0.65+0.02.*(dT-40)) + double(dT>=50).*(0.85+0.0075.*(dT-50)));
            d_gs_dry = double(s<=0 & d<=0) .* 1.0417e-9.*f.*h.*g;
            
            snow.TEMP.metam_d_d = (double(W_liq>0).*d_d_wet + double(W_liq==0).*d_d_dry)./daysec;
            snow.TEMP.metam_d_s = (double(W_liq>0).*d_s_wet + double(W_liq==0).*d_s_dry)./daysec;
            snow.TEMP.metam_d_gs = d_gs_wet + d_gs_dry; % No need to divide by daysec!
        end
        
        %compaction due to wind drift
        function snow = prog_wind_drift(snow)
            timescale = snow.PARA.timescale_winddrift;
            rho = max(50, snow.STATVAR.waterIce ./ (snow.STATVAR.layerThick .* snow.STATVAR.area) .*920);
            %rho = max(50, snow.STATVAR.waterIce ./ (snow.STATVAR.layerThick .* snow.STATVAR.area) .*1000);
            d = snow.STATVAR.d;
            s = snow.STATVAR.s;
            gs = snow.STATVAR.gs;
            D = snow.STATVAR.layerThick;
            
            F = 1.25-0.0042.*(max(50,rho)-50);
            Mo = double(d>0).*(0.34.*(0.75.*d-0.5.*s+0.5)+0.66.*F) + double(d<=0).*(0.34.*(-0.583.*gs.*1000-0.833.*s+0.833)+0.66.*F);
            Si = (-2.868 .* exp(-0.085.*snow.TEMP.wind_surface)+1+Mo) .* double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f);
            
            
            zi=Si.*0;  %Problem fixed Sebastian Westermann June 2020!
            i=2;
            while i<=size(Si,1) && sum(Si(1:i)<=0,1)==0
                zi(i)=sum(D(1:i-1).*(3.25-Si(1:i-1)));
                i=i+1;
            end
            one_over_tau=max(0,Si.*exp(-zi./0.1))./(timescale.*3600);
            one_over_tau(i:end,1) = 0;
            
            snow.TEMP.wind_d_d = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* -d./2.*one_over_tau;
            snow.TEMP.wind_d_s = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* (1-s).*one_over_tau;
            snow.TEMP.wind_d_gs = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* double(d==0).*5e-4.*one_over_tau;
            %d_rho = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* double(rho<350).* (350-rho).*one_over_tau;
            d_rho = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* double(rho<snow.PARA.max_wind_slab_density).* (snow.PARA.max_wind_slab_density - rho).*one_over_tau; %new value Barrere et al., 2017, GMD
            snow.TEMP.wind_d_D = -D./rho .*d_rho;
            snow.TEMP.Si=Si;
            snow.TEMP.one_over_tau = one_over_tau;
        end
        
        %snow compaction due to overlying snow, Vionnet et al., 2012
        function snow = compaction(snow)
            rho = max(50, snow.STATVAR.waterIce ./ (snow.STATVAR.layerThick .* snow.STATVAR.area) .*1000);
            gs = snow.STATVAR.gs;
            D = snow.STATVAR.layerThick;
            W_liq = snow.STATVAR.water*1000;
            T = snow.STATVAR.T;
            
            %sigma=9.81.*cosd(snow.PARA.slope).*rho.*D;
            sigma=9.81.*cosd(snow.PARA.slope).*snow.STATVAR.waterIce ./ snow.STATVAR.area .*1000;
            sigma(end)=sum(sigma(1:end-1));
            for i=size(sigma,1)-1:-1:2
                sigma(i)=sigma(i+1)-sigma(i);
            end
            sigma(1)=0.5.*sigma(1);
            
            f1=(1+60.*W_liq./1000./D).^(-1);
            f2= min(4, exp(min(0.4, gs.*1000-0.2)/0.1));
            eta=f1.*f2.*7.62237e6./250.*rho.*exp(-T.*0.1+ 0.023.*rho);
            
            snow.TEMP.compact_d_D = - sigma./eta.*D;
        end
        
        %advance prognostic variable for new snow falling during a timestep
        function  snow = advance_prognostic_new_snow_crocus(snow, timestep)
            snow.TEMP.newSnow.STATVAR.waterIce = timestep .* snow.TEMP.snowfall;
            snow.TEMP.newSnow.STATVAR.ice = timestep .* snow.TEMP.snowfall;
            snow.TEMP.newSnow.STATVAR.layerThick = timestep .* snow.TEMP.snowfall ./snow.STATVAR.area(1) ./ (snow.TEMP.newSnow.STATVAR.density ./1000);
            snow.TEMP.newSnow.STATVAR.energy = timestep .* snow.TEMP.snow_energy;
        end
        
        function  snow = advance_prognostic_new_snow_CHILD_crocus(snow, timestep)
            snow.TEMP.newSnow.STATVAR.waterIce = timestep .* snow.TEMP.snowfall;
            snow.TEMP.newSnow.STATVAR.ice = timestep .* snow.TEMP.snowfall;
            snow.TEMP.newSnow.STATVAR.volume = timestep .* snow.TEMP.snowfall ./ (snow.TEMP.newSnow.STATVAR.density ./1000);
            snow.TEMP.newSnow.STATVAR.energy = timestep .* snow.TEMP.snow_energy;
        end
                
        function z0 = get_z0_surface(snow)
            z0 = snow.PARA.z0;
        end
        
        %-------triggers-----------------
        %make SNOW a CHILD again
        function snow = make_SNOW_CHILD(snow)
            if size(snow.STATVAR.layerThick,1) == 1 && snow.STATVAR.ice(1,1) ./ snow.STATVAR.area(1,1) < 0.5 .* snow.PARA.swe_per_cell
                
                ground = snow.NEXT;
                
                ground.PREVIOUS = snow.PREVIOUS; %reassign ground
                ground.PREVIOUS.NEXT = ground;
                ground.CHILD = snow;
                snow.PARENT = ground;
                ground.IA_CHILD = snow.IA_NEXT;
                ground.IA_CHILD.NEXT = ground;
                ground.IA_CHILD.PREVIOUS = snow;
                
                if ~isempty(snow.PREVIOUS.IA_NEXT) % Class above snow is a "real" stratigraphy class
                    ia_class = get_IA_class(class(snow.PREVIOUS),class(ground));
                    ground.IA_PREVIOUS = ia_class();
                    ground.PREVIOUS.IA_NEXT = ia_class();
                    ground.IA_PREVIOUS.NEXT = ground;
                    ground.IA_PREVIOUS.PREVIOUS = snow.PREVIOUS;
                else
                ground.IA_PREVIOUS=[]; %change to get_ia_class, if there is a possibility for another class on top of the snow cover
                end
                %snow.NEXT =[];  %cut all dependencies, except for snow.NEXT which keeps being pointed to snow.PARENT, so that SW radiation can be transmitted
                snow.PREVIOUS =[];
                snow.IA_NEXT =[];
                snow.IA_PREVIOUS =[];
                
                %change to constant layerThick, variable area
                volume = snow.STATVAR.layerThick .* snow.STATVAR.area;
                snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ snow.STATVAR.target_density; %[m] constant layerThick
                snow.STATVAR.area = volume ./ snow.STATVAR.layerThick;
                
                snow = ground; %assign snow pointer to ground to return to regular stratigraphy
            end
        end
        
        function snow = make_SNOW_CHILD2(snow, tile)
            if size(snow.STATVAR.layerThick,1) == 1 && snow.STATVAR.ice(1,1) ./ snow.STATVAR.area(1,1) < 0.5 .* snow.PARA.swe_per_cell
                
                ground = snow.NEXT;
                
                ground.PREVIOUS = snow.PREVIOUS; %reassign ground
                ground.PREVIOUS.NEXT = ground;
                ground.CHILD = snow;
                snow.PARENT = ground;
                ground.IA_CHILD = snow.IA_NEXT;
                ground.IA_CHILD.NEXT = ground;
                ground.IA_CHILD.PREVIOUS = snow;
                
                ground.IA_PREVIOUS=[]; %change to get_ia_class, if there is a possibility for another class on top of the snow cover
                
                %snow.NEXT =[];  %cut all dependencies, except for snow.NEXT which keeps being pointed to snow.PARENT, so that SW radiation can be transmitted
                snow.PREVIOUS =[];
                snow.IA_NEXT =[];
                snow.IA_PREVIOUS =[];
                
                %change to constant layerThick, variable area
                volume = snow.STATVAR.layerThick .* snow.STATVAR.area;
                snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ snow.STATVAR.target_density; %[m] constant layerThick
                snow.STATVAR.area = volume ./ snow.STATVAR.layerThick;
                
                ground.CHILD = compute_diagnostic_CHILD(ground.CHILD, tile);
                %ground = check_trigger(ground, tile);
                
                snow = ground; %assign snow pointer to ground to return to regular stratigraphy
            end
        end
        
        function snow = make_SNOW_CHILD_ubT(snow)
            if size(snow.STATVAR.layerThick,1) == 1 && snow.STATVAR.ice(1,1) ./ snow.STATVAR.area(1,1) < 0.5 .* snow.PARA.swe_per_cell
                
                ground = snow.NEXT;
                
                ground.PREVIOUS = snow.PREVIOUS; %reassign ground
                ground.PREVIOUS.NEXT = ground;
                ground.CHILD = snow;
                snow.PARENT = ground;
                ground.IA_CHILD = snow.IA_NEXT;
                ground.IA_CHILD.NEXT = ground;
                ground.IA_CHILD.PREVIOUS = snow;
                
                ground.IA_PREVIOUS=[]; %change to get_ia_class, if there is a possibility for another class on top of the snow cover
                
                %snow.NEXT =[];  %cut all dependencies, except for snow.NEXT which keeps being pointed to snow.PARENT, so that SW radiation can be transmitted
                snow.PREVIOUS =[];
                snow.IA_NEXT =[];
                snow.IA_PREVIOUS =[];
                
                %change to constant layerThick, variable area
                volume = snow.STATVAR.layerThick .* snow.STATVAR.area;
                snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ snow.PARA.density; %[m] constant layerThick
                snow.STATVAR.area = volume ./ snow.STATVAR.layerThick;
                
                snow = ground; %assign snow pointer to ground to return to regular stratigraphy
            end

            if isempty(snow.STATVAR.layerThick) %snow depth has become zero in a ssigle time step, jumping over the CHILD phase 
                
                ground = snow.NEXT;
                
                ground.PREVIOUS = snow.PREVIOUS; %reassign ground
                ground.PREVIOUS.NEXT = ground;
                ground.CHILD = 0;
               
                ground.IA_CHILD = [];
                
                ground.IA_PREVIOUS=[]; %change to get_ia_class, if there is a possibility for another class on top of the snow cover
                
                %snow.NEXT =[];  %cut all dependencies, except for snow.NEXT which keeps being pointed to snow.PARENT, so that SW radiation can be transmitted
                snow.PREVIOUS =[];
                snow.IA_NEXT =[];
                snow.IA_PREVIOUS =[];

                snow = ground; %assign snow pointer to ground to return to regular stratigraphy
            end
        end
        
    end
end

