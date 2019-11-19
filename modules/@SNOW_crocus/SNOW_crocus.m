


classdef SNOW_crocus < SNOW_base_class
    
    methods
        
        %------ Mandatory Functions ------% 
        
         function xls_out = write_excel(snow)
            xls_out = {} %Fill out later
         end
        
         function snow = provide_variables(snow)  %initializes the subvariables as empty arrays 
           snow = provide_variables@SNOW_base_class(snow); 
           snow = provide_PARA(snow);
           snow = provide_CONST(snow);
           snow = provide_STATVAR(snow);
         end
        
         function variable = initialize_from_file(snow, variable, section)
            variable = initialize_from_file@SNOW_base_class(snow, variable, section);
         end
        
         function snow = assign_global_variables(snow, forcing)
            snow.PARA.airT_height = forcing.PARA.airT_height;
         end
         
         function snow = initialize_zero_snow(snow, parentGround)
            snow = initialize_zero_snow@SNOW_base_class(snow, parentGround);
         end
        
         function snow = get_boundary_condition_u(snow, forcing)
             [snow.TEMP.albedo, snow.TEMP.spectral_ranges, snow.TEMP.spectral_albedo, snow.TEMP.damping_factor] = get_albedo(snow, FORCING.TEMP.t, forcing.TEMP.p)
             
             snow.STATVAR.Sout = snow.PARA.albedo .*  forcing.TEMP.Sin;
             snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
             snow.STATVAR.Qh = Q_h(snow, forcing);
             snow.STATVAR.Qe = Q_eq_potET(snow, forcing);
             
             snow.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe;
             snow.TEMP.U = forcing.TEMP.wind;
             snow.TEMP.rain_energy = snow.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;

%             snow.TEMP.T_rainWater = forcing.TEMP.Tair;
%             snow.TEMP.F_ub_water = snow.TEMP.rainfall;
            snow.TEMP.d_ice_sublim = -snow.STATVAR.Qe ./ snow.CONST.L_s; % snow.CONST.L_s ~= L_i IN CRYOGRID 3!
%             snow.TEMP.d_E_sublim =  snow.TEMP.d_ice_sublim .* (snow.CONST.c_w .*snow.STATVAR.T(1) - snow.CONST.L_f);
%           Copied from SNOW_simple_seb_bucketW
         end
         
         function snow = get_boundary_condition_l(snow, forcing) 
            snow = get_boundary_condition_l@SNOW_base_class(snow, forcing);
         end
         
         function snow = get_derivatives_prognostic(snow)
            % Radiative(solar) and conductive (T gradient) energy transfer
            snow = get_derivatives_prognostic@SNOW_base_class(snow);
            snow.TEMP.d_energy = snow.TEMP.d_energy + snow.TEMP.d_energy_seb;
            
            % Snow layer properties
            snow = get_T_water(snow);
            snow = get_T_gradient_snow(snow);
            snow = prog_metamorphism(snow);
            snow = prog_wind_drift(snow); % Add blowing snow sublimation according to Gordon etAl 2006
            snow = compaction(snow);
         end
         
         function timestep = get_timestep(snow)  %could involve check for several state variables
             timestep1 = get_timestep@SNOW_base_class(snow);
             timestep2 = min((-snow.STATVAR.energy ./ snow.TEMP.d_energy) .*double(snow.TEMP.d_energy>0) + double(snow.TEMP.d_energy<=0).*1e5); %when snow is melting, do not melt more than there is in a grid cell
             timestep = min(timestep1, timestep2);
         end
         
         function snow = advance_prognostic(snow, timestep)
            snow = advance_prognostic@SNOW_base_class(snow, timestep); % OR MAYBE NOT?
            
            snow.STATVAR.d = max(snow.STATVAR.d.*0, snow.STATVAR.d + timestep .*(snow.TEMP.metam_d_d + snow.TEMP.wind_d_d));
            snow.STATVAR.s = max(snow.STATVAR.s.*0, min(snow.STATVAR.s.*0+1, snow.STATVAR.s + timestep .*(snow.TEMP.metam_d_s + snow.TEMP.wind_d_s))); 
        	snow.STATVAR.gs = max(snow.STATVAR.gs, snow.STATVAR.gs + timestep .*(snow.TEMP.metam_d_gs + snow.TEMP.wind_d_gs));
            % CAN'T gs DECREASE WHEN DRIFTING/ BREAKING APPART? 
            snow.STATVAR.layerThick = min(snow.STATVAR.layerThick, max(snow.STATVAR.ice, snow.STATVAR.layerThick + timestep .*(snow.TEMP.wind_d_D + snow.TEMP.compact_d_D)));
            % CANT'T D BECOME SMALLER THAN snow.STATVAR.ice WHEN DRIFTING?
            snow.TEMP.ice_fraction_old_snow = snow.STATVAR.ice./snow.STATVAR.layerThick;
         end
         
        function snow = compute_diagnostic_first_cell(snow, forcing)
            snow = L_star(snow, forcing); % From SNOW_simple_seb
            % WHY HERE % WHAT IS L_star?  
        end
        
        function snow = compute_diagnostic(snow, forcing) % number refers to # in Cryogrid3
            snow = infiltrate_rain_energy(snow, timestep); % #1
            snow = sublimation_loss(snow, timestep); % #1a timestep IS NOT AVAILABLE HERE, MOVE TO advance_prognistics??
            
            snow = get_T_water(snow); % #2
            snow = adjust_layerThick(snow); % #3
            
            snow = melt_infiltrate(snow, timestep); % #4
            % HOW TO INFILTRATE INTO GROUND BELOW?
            if snow.TEMP.meltwater>0 %pool up from bottom
                snow = infiltrateBottom2Top(snow);
            end
            
            snow = get_T_water(snow);
            snow = reduce_snow_grid(snow); % #5         
            
            snow = sublimation_gain(snow, timestep); % #6
            
            snow = snowfall(snow, forcing, timestep); % #7&8
            
            snow = get_T_water(snow); % #9
            
            snow = conductivity_snow_Yen(snow);
            snow.STATVAR.upperPos = snow.STATVAR.lowerPos + cumsum(snow.STATVAR.layerThick);
        end
        
%---------------- Non-Mandatory Functions ----------------% 
        
        function [total_albedo, spectral_ranges, spectral_albedo, damping_factor] = get_albedo_snow(snow, t, pressure)

            rho=snow.STATVAR.waterIce./snow.STATVAR.layerThick.*1000;
            d=snow.STATVAR.d;
            s=snow.STATVAR.s;
            gs=snow.STATVAR.gs;
            snow_age=t-snow.STATVAR.time_snowFall(1);

            spectral_ranges = [0.71; 0.21; 0.08];

            d_opt=double(d>0).*1e-4.*(d+(1-d).*(4-s)) + double(d==0).*(gs.*s+(1-s).*max(gs.*0+4e-4, gs/2));

            albedo1=max(0.6, min(0.92,0.96-1.58.*(d_opt(1).^0.5)) - min(1,max(pressure./870e2,0.5)).*0.2./60.*snow_age);
            albedo2=max(0.3, 0.9-15.4.*(d_opt(1).^(0.5)));
            albedo3=346.3.*min(d_opt(1), 0.0023) - 32.1.*((min(d_opt(1), 0.0023)).^0.5) + 0.88;

            spectral_albedo=[albedo1; albedo2; albedo3];

            beta1=max(40+d_opt.*0, 0.00192.*rho./(d_opt.^0.5));
            beta2=max(100, 0.01098.*rho./(d_opt.^0.5));

            damping_factor=[beta1 beta2 beta2.*Inf];

            total_albedo = sum(spectral_albedo.*spectral_ranges);
        end
            
        function snow = energy_infiltration(snow,forcing)
            Q_solar=(1-snow.TEMP.spectral_albedo).*snow.TEMP.spectral_ranges.*forcing.TEMP.Sin;
            residual=0.1;
            i=1;
            d_energy = snow.STATVAR.energy.*0;
            
            while i<=size(snow.STATVAR.layerThick,1) && sum(Q_solar)>=residual
                reduction=exp(-snow.TEMP.damping_factor(i,:)'.*snow.STATVAR.layerThick(i));  %vector
                d_energy(i) = sum(Q_solar .* (1-reduction));
                Q_solar=reduction.*Q_solar;  %vector

                i=i+1;
            end
            if sum(Q_solar) < residual
                d_energy(i) = sum(Q_solar); 
            else
                % HOW TO TRANSFER REST OF ENERGY TO TOP OF SOIL???
                snow.TEMP.Qsolar_trans_lb = sum(Q_solar); % Solar energy transmitted through snowpack
            end
            
            snow.TEMP.d_energy_seb = d_energy;
        end
        
        function snow = get_T_water(snow)
            E_frozen = - snow.STATVAR.waterIce .* snow.CONST.L_f;
            
            % WHAT IF ENERGY = E_frozen (T = 0)?
            snow.STATVAR.T = double(snow.STATVAR.energy<=E_frozen) .* (snow.STATVAR.energy-E_frozen) ./ (snow.STATVAR.waterIce .* snow.CONST.c_i);
            snow.STATVAR.water = double(snow.STATVAR.energy > E_frozen) .*  snow.STATVAR.waterIce .* (E_frozen - snow.STATVAR.energy) ./E_frozen;
            snow.STATVAR.ice = double(snow.STATVAR.energy > E_frozen) .*  snow.STATVAR.waterIce .* (snow.STATVAR.energy) ./E_frozen + double(snow.STATVAR.energy <= E_frozen) .* snow.STATVAR.waterIce;
        
        end
        
        function snow = get_T_gradient_snow(snow)

            delta=[snow.STATVAR.layerThick; snow.NEXT.STATVAR.layerThick(1,1)];
            T = [snow.STATVAR.T; snow.NEXT.STATVAR.T(1,1)];
            dT = (T(1:end-2) - T(2:end))./(0.5.*delta(1:end-2) +  delta(2:end-1) + 0.5.*delta(3:end));
            dT =abs([ (T(1) - T(2))./(0.5.*delta(1) + 0.5.*delta(2)) ; dT]);
            % WHY FORWARD DIFFERENCE FOR UPPERMOST CELL? 
            snow.TEMP.dT = dT;
        end
        
        function snow = prog_metamorphism(snow)
            T = snow.STATVAR.T;
            D = snow.STATVAR.layerThick;
            D_water = snow.STATVAR.water;
            dT = snow.TEMP.dT;
            d = snow.STATVAR.d;
            s = snow.STATVAR.s;
            gs = snow.STATVAR.gs;
            
            rho = snow.STATVAR.waterIce./snow.STATVAR.layerThick .*1000;
            W_liq=D_water.*1000;
            
            small_gradient=(dT<=5);
            gradient_term=double(small_gradient)+double(~small_gradient).*dT.^0.4;
            a=2e8.*exp(-6000./(T+273.15)).*gradient_term;
            
            D=max(D, D.*0+5e-3);
            theta_cubed=(100.*W_liq./rho./D).^3;
            
            % Wet snow metamporphism
            d_d_wet = -1/16.*theta_cubed;
            d_s_wet = +1/16.*theta_cubed;
            d_gs_wet= double(s>=1 & d<=0) .*  2./3.14./gs.^2 .* (1.28e-17 + 4.22e-19.*theta_cubed).*day_sec;
            % WHY TIMES day_sec IN d_gs_wet?
            
            % Dry snow metamporhism
            d_d_dry = -a;
            d_s_dry = 5.*a.*double(small_gradient)- a.*double(~small_gradient);
            
            f=double(T<-40) + double(T>=-40 & T<-22).*0.011.*(T+40) +  double(T>=-22 & T<-6).*(0.2+0.05.*(T+22)) + double(T>=-6).*(1-0.05.*T);
            h=max(rho.*0, min(rho.*0+1, 1-0.004.*(rho-150)));
            g=min(dT.*0+1, double(dT>=15 & dT<25).*0.01.*(dT-15) + double(dT>=25 & dT<40).*(0.1+0.037.*(dT-25)) + double(dT>=40 & dT<50).*(0.65+0.02.*(dT-40)) + double(dT>=50).*(0.85+0.0075.*(dT-50)));
            d_gs_dry = double(s<=0 & d<=0) .* 1.0417e-9.*f.*h.*g.*day_sec;

            snow.TEMP.metam_d_d = double(W_liq>0).*d_d_wet + double(W_liq==0).*d_d_dry;
            snow.TEMP.metam_d_s = double(W_liq>0).*d_s_wet + double(W_liq==0).*d_s_dry;
            snow.TEMP.metam_d_gs = d_gs_wet + d_gs_dry;
        end
        
        function snow = prog_wind_drift(snow)
            timescale = snow.PARA.timescale_Winddrift;
            rho = snow.STATVAR.waterIce./snow.STATVAR.layerThick .*1000;
            d = snow.STATVAR.d;
            s = snow.STATVAR.s;
            gs = snow.STATVAR.gs;
            D = snow.STATVAR.layerThick;
            
            F = 1.25-0.0042.*(max(50,rho)-50);
            Mo = double(d>0).*(0.34.*(0.75.*d-0.5.*s+0.5)+0.66.*F) + double(d<=0).*(0.34.*(-0.583.*gs.*1000-0.833.*s+0.833)+0.66.*F);
            Si=-2.868 .* exp(-0.085.*snow.TEMP.U)+1+Mo;
            
            zi=Si.*0;
            i=2;
            while i<=size(Si,1) && Si(i)>0  
                zi(i)=sum(D(1:i-1).*(3.25-Si(1:i-1)));
                i=i+1;
            end
            one_over_tau=max(0,Si.*exp(-zi./0.1))./(timescale.*3600);
            
            snow.TEMP.wind_d_d = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* -d./2.*one_over_tau.*day_sec;
            snow.TEMP.wind_d_s = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* (1-s).*one_over_tau.*day_sec;
            snow.TEMP.wind_d_gs = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* double(d==0).*5e-4.*one_over_tau.*day_sec;
            d_rho = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* double(rho<350).* (350-rho).*one_over_tau.*day_sec;
            snow.TEMP.wind_d_D = -D./rho .*d_rho;
            snow.TEMP.Si=Si;
        end
        
        function snow = compaction(snow)
            rho = snow.STATVAR.waterIce./snow.STATVAR.layerThick .*1000;
            d = snow.STATVAR.d;
            s = snow.STATVAR.s;
            gs = snow.STATVAR.gs;
            D = snow.STATVAR.layerThick;
            W_liq = snow.STATVAR.water*1000;
            
            sigma=9.81.*cosd(snow.PARA.slope).*rho.*D;
            sigma(end)=sum(sigma(1:end-1));
            for i=size(sigma,1)-1:-1:2
                sigma(i)=sigma(i+1)-sigma(i);
            end
            sigma(1)=0.5.*sigma(1);
            
            f1=(1+60.*W_liq./1000./D).^(-1);
            f2= min(4, exp(min(0.4, gs.*1000-0.2)/0.1));
            eta=f1.*f2.*7.62237e6./250.*rho.*exp(-T.*0.1+ 0.023.*rho);
            
            snow.TEMP.compact_d_D = - day_sec.*sigma./eta.*D;
        end
       
            
        function snow = infiltrate_rain_energy(snow,timestep)
            energy = snow.STATVAR.energy;
            
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.rain_energy.*timestep;
            
            while sum(energy(1:end-1)>0)~=0  %all ends up in last snow grid cell
                change = double(energy>=0) .* energy;
                energy(1:end-1) = energy(1:end-1) - change(1:end-1);
                energy(2:end) = energy(2:end) + change(1:end-1);
            end

            snow.STATVAR.energy = energy;
        end
        
        function snow = sublimation_loss(snow,timestep)
            Lf = snow.CONST.L_f;
            i = 1; 
            sublim_ice = snow.TEMP.d_ice_sublim.*timestep;
            while sublim_ice < 0 && i<=size(snow.STATVAR.waterIce,1)
                sublim_ice = -sublim_ice;
                D_ice=min(-snow.STATVAR.energy(i)./Lf, snow.STATVAR.waterIce(i));
                change = min(sublim_ice, D_ice);
               
                snow.STATVAR.energy(i) = double(snow.STATVAR.energy(i) < -snow.STATVAR.waterIce(i).*Lf) .* snow.STATVAR.energy(i).*(snow.STATVAR.waterIce(i)-change)./snow.STATVAR.waterIce(i) + ...
                    double(snow.STATVAR.energy(i) >=  -snow.STATVAR.waterIce(i).*Lf) .* (snow.STATVAR.energy(i) - change.*Lf);
                snow.STATVAR.waterIce(i)=snow.STATVAR.waterIce(i)-change;
                snow.STATVAR.layerThick(i) = snow.STATVAR.layerThick(i).*(D_ice-change)./D_ice;
                
                sublim_ice = sublim_ice-change;
                i=i+1;
            end
        end
            
        function snow = adjust_layerThick(snow)
            snow.STATVAR.layerThick = snow.STATVAR.ice./snow.TEMP.ice_fraction_old_snow;
        end

        function snow = melt_infiltrate(snow, timestep)
            Lf = snow.CONST.L_f;
            max_water = snow.PARA.field_capacity.*(snow.STATVAR.layerThick - snow.STATVAR.ice);  % field capacity of 5%of the porosity
            meltwater = max(0, snow.STATVAR.water - max_water);
            water = min(snow.STATVAR.water, max_water);

            meltwater = [snow.TEMP.rainfall.*timestep; meltwater]; %has one element more than snow

            if sum(meltwater)>0
                infiltration_capacity = min(snow.STATVAR.layerThick-snow.STATVAR.ice, (max_water - snow.STATVAR.water) + ...
                    double(snow.STATVAR.energy < - snow.STATVAR.waterIce.*Lf) .* (-snow.STATVAR.waterIce.*Lf - snow.STATVAR.energy)./Lf );

                for i=1:size(snow.STATVAR.layerThick,1)
                    water(i)=water(i) +  min(infiltration_capacity(i), meltwater(i));   %
                    meltwater(i+1) =  meltwater(i+1) + max(0, meltwater(i)-infiltration_capacity(i));
                    meltwater(i)=0;
                end

            end
            snow.STATVAR.waterIce = snow.STATVAR.ice + water;
            snow.TEMP.meltwater = meltwater(end);
        end

        function snow = infiltrateBottom2Top(snow)
            meltwater = snow.TEMP.meltwater;
            i=size(snow.STATVAR.layerThick,1);
            while meltwater>0 && i>0
                infiltration_capacity =  snow.STATVAR.layerThick - snow.STATVAR.waterIce;
                infiltration = min(meltwater, infiltration_capacity);
                meltwater=meltwater-infiltration;
                snow.STATVAR.waterIce(i)=snow.STATVAR.waterIce(i)+infiltration;
                i=i-1; 
            end
        end
        
        function snow = reduce_snow_grid(snow)
            i=1;
            while i<size(snow.STATVAR.layerThick,1)
                if snow.STATVAR.ice(i) < .5.*snow.PARA.swe_per_cell
                    % Linearly mixed properties
                   snow.STATVAR.s(i+1) = (snow.STATVAR.ice(i).* snow.STATVAR.s(i) + snow.STATVAR.ice(i+1).* snow.STATVAR.s(i+1)) ./ (snow.STATVAR.ice(i) + snow.STATVAR.ice(i+1));  
                   snow.STATVAR.d(i+1) = (snow.STATVAR.ice(i).* snow.STATVAR.d(i) + snow.STATVAR.ice(i+1).* snow.STATVAR.d(i+1)) ./ (snow.STATVAR.ice(i) + snow.STATVAR.ice(i+1));  
                   snow.STATVAR.gs(i+1) = (snow.STATVAR.ice(i).* snow.STATVAR.gs(i) + snow.STATVAR.ice(i+1).* snow.STATVAR.gs(i+1)) ./ (snow.STATVAR.ice(i) + snow.STATVAR.ice(i+1));  

                   snow.STATVAR.s(i) = [];
                   snow.STATVAR.d(i) = [];
                   snow.STATVAR.g(i) = [];
                   snow.STATVAR.time_snowfall(i) = [];

                    % Additive properties
                   snow.STATVAR.energy(i+1) = snow.STATVAR.energy(i) + snow.STATVAR.energy(i+1);
                   snow.STATVAR.energy(i) = [];

                   snow.STATVAR.layerThick(i+1) = snow.STATVAR.layerthick(i) + snow.STATVAR.layerThick(i+1);
                   snow.STATVAR.layerThick(i) = [];

                   snow.STATVAR.waterIce(i+1) = snow.STATVAR.waterIce(i) + snow.STATVAR.waterIce(i+1);
                   snow.STATVAR.waterIce(i) = [];

                   snow.STATVAR.water(i+1) = snow.STATVAR.water(i) + snow.STATVAR.water(i+1);
                   snow.STATVAR.water(i) = [];

                   snow.STATVAR.ice(i+1) = snow.STATVAR.ice(i) + snow.STATVAR.ice(i+1);
                   snow.STATVAR.ice(i) = [];

                else 
                   i = i+1;
                end
            end
            
            if snow.STATVAR.ice(end) < .5*snow.PARA.swe_per_cell % bottom cell
                if size(snow.STATVAR.s) > 1
                   snow.STATVAR.s(end-1) = (snow.STATVAR.ice(end).* snow.STATVAR.s(end) + snow.STATVAR.ice(end-1).* snow.STATVAR.s(end-1)) ./ (snow.STATVAR.ice(end) + snow.STATVAR.ice(end-1));  
                   snow.STATVAR.d(end-1) = (snow.STATVAR.ice(end).* snow.STATVAR.d(end) + snow.STATVAR.ice(end-1).* snow.STATVAR.d(end-1)) ./ (snow.STATVAR.ice(end) + snow.STATVAR.ice(end-1));  
                   snow.STATVAR.gs(end-1) = (snow.STATVAR.ice(end).* snow.STATVAR.gs(end) + snow.STATVAR.ice(end-1).* snow.STATVAR.gs(end-1)) ./ (snow.STATVAR.ice(end) + snow.STATVAR.ice(end-1));  

                   snow.STATVAR.s(end) = [];
                   snow.STATVAR.d(end) = [];
                   snow.STATVAR.gs(end) = [];
                   snow.STATVAR.time_snowfall(end) = [];
                   
                   snow.STATVAR.energy(end-1) = snow.STATVAR.energy(end) + snow.STATVAR.energy(end-1);
                   snow.STATVAR.energy(end) = [];

                   snow.STATVAR.layerThick(end-1) = snow.STATVAR.layerthick(end) + snow.STATVAR.layerThick(end-1);
                   snow.STATVAR.layerThick(end) = [];

                   snow.STATVAR.waterIce(end-1) = snow.STATVAR.waterIce(end) + snow.STATVAR.waterIce(end-1);
                   snow.STATVAR.waterIce(end) = [];

                   snow.STATVAR.water(end-1) = snow.STATVAR.water(end) + snow.STATVAR.water(end-1);
                   snow.STATVAR.water(end) = [];

                   snow.STATVAR.ice(end-1) = snow.STATVAR.ice(end) + snow.STATVAR.ice(end-1);
                   snow.STATVAR.ice(end) = [];
                   
                else % only one cell
                   snow.STATVAR.s = [];
                   snow.STATVAR.d = [];
                   snow.STATVAR.g = [];
                   snow.STATVAR.time_snowfall = [];
                   snow.STATVAR.energy = [];
                   snow.STATVAR.layerThick = [];
                   snow.STATVAR.waterIce = [];
                   snow.STATVAR.water = [];
                   snow.STATVAR.ice = [];
                end
                % FROM SEBASTIAN:
                   %the water of the last cell is lost in the current version, should
                   %be added to the meltwater flux - on Samoilov, this can be something
                   %like 10-20% of the annual snowpack - also conservation of energy
                   %is violated to the extent of 0.5.*PARA.technical.SWEperCell .* Lf
                   %set threshold lower in the future!
            end
        end
        
        function snow = sublimation_gain(snow,timestep)
            Lf = snow.CONST.L_f;
            sublim_ice = snow.TEMP.d_ice_sublim.*timestep;
            
            if sublim_ice > 0
                ice_fraction = snow.STATVAR.ice(1)./snow.STATVAR.layerThick(1);
                
                snow.STATVAR.layerThick(1) = (snow.STATVAR.ice(1) + sublim_ice)./ice_fraction;
                snow.STATVAR.energy(1) = double(snow.STATVAR.energy(1) <  -snow.STATVAR.waterIce(1).*Lf).*snow.STATVAR.energy(1).* (snow.STATVAR.waterIce(1) + sublim_ice)./snow.STATVAR.waterIce(1) ...
                    + double(snow.STATVAR.energy(1) >=  -snow.STATVAR.waterIce(1).*Lf).*(snow.STATVAR.energy(1) + sublim_ice.*Lf);
                snow.STATVAR.waterIce = snow.STATVAR.waterIce + sublim_ice;
                % WHAT HAPPENS TO snow.STATVAR.water ./ice ??
            end
            
        end
        
        function snow = snowfall(snow,forcing,timestep)
            if isempty(snow.STATVAR.layerThick) % no previous snowcover
                snow.TEMP.initialSWE = snow.TEMP.initialSWE + snow.TEMP.snowfall.*timestep./1000 - snow.TEMP.initialSWE .* 0.1 .* timestep;
                % WHY MINUS 0.1*timestep???
                if snow.TEMP.initialSWE >= snow.PARA.swe_per_cell
                    snow = setInitialSnow(snow,forcing);
                    snow.TEMP.initialSWE = 0;
                end
            else % snowcover exists
                newSnow = setInitialSnow(snow,forcing);
                snow = combineNewSnowLayer(snow,newSnow);
            end
        end
        
        function snow = setInitialSnow(snow,forcing)
            T_air=min(forcing.TEMP.Tair, 0);
            windspeed = forcing.TEMP.wind;
            
            T_fus=0;  %degree C
            a_rho=109;  %kg/m3
            b_rho=6; % kg/m3K
            c_rho=26; % kg m7/2s1/2
            min_snowDensity=50; %kg/m3
            rho_w=1000;
            
            rho=max(min_snowDensity, a_rho+b_rho.*(T_air-T_fus)+ c_rho.*windspeed.^0.5);  %initial snow density
            
            snow.STATVAR.waterIce = snow.TEMP.initialSWE;
            snow.STATVAR.water = 0; % dry by definition
            snow.STATVAR.ice = snow.TEMP.initialSWE; 
            snow.STATVAR.layerThick = snow.TEMP.initalSWE.*rho_w./rho;
            snow.STATVAR.d=min(max(1.29-0.17.*windspeed,0.2),1);
            snow.STATVAR.s=min(max(0.08.*windspeed + 0.38,0.5),0.9);
            snow.STATVAR.gs=0.1e-3+(1-snow.STATVAR.d).*(0.3e-3-0.1e-3.*snow.STATVAR.s);
            snow.STATVAR.time_snowFall = timestamp;
            snow.STATVAR.energy = get_energy_first_snow_cell(T_air, rho).*snow.STATVAR.layerThick;
        end
        
        function energy = get_energy_first_snow_cell(T, rho)
            rho_i=1000; % WHY? RHO ICE = 917??
            E_max=-rho.*snow.CONST.L_f;

            energy=E_max + T.*c_i.*rho./rho_i;
        end
        
        function snow = combineNewSnowLayer(snow,newSnow)
                % WHAT IF snowfall*timestep >= 0.5*SWE_per_cell?

            if (snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce(1)) < 1.5 .* snow.PARA.swe_per_cell % Layers can be merged
                % Linear mixing
                snow.STATVAR.d(1) = (snow.STATVAR.d(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.d.*snow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.s(1) = (snow.STATVAR.s(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.s.*snow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.gs(1) = (snow.STATVAR.gs(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.gs.*snow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.time_snowfall(1) = (snow.STATVAR.time_snowfall(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.time_snowfall.*snow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                % Additive
                snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + newSnow.STATVAR.energy;
                snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce;
                snow.STATVAR.ice(1) = snow.STATVAR.ice(1) + newSnow.STATVAR.ice;
                snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + newSnow.STATVAR.layerThick;
            else % split layers before adding snow (as above)
                fraction = (snow.STATVAR.waterIce(1)-snow.PARA.swe_per_cell)./snow.STATVAR.waterIce(1);
                snow.STATVAR.waterIce = [snow.STATVAR.waterIce(1)*fraction; snow.STATVAR.waterIce(1)*(1-fraction); snow.STATVAR.waterIce(2:end)];
                snow.STATVAR.water = [snow.STATVAR.water(1)*fraction; snow.STATVAR.water(1)*(1-fraction); snow.STATVAR.water(2:end)];
                snow.STATVAR.ice = [snow.STATVAR.ice(1)*fraction; snow.STATVAR.ice(1)*(1-fraction); snow.STATVAR.ice(2:end)];
                snow.STATVAR.layerThick = [snow.STATVAR.layerThick(1)*fraction; snow.STATVAR.layerThick(1)*(1-fraction); snow.STATVAR.layerThick(2:end)];
                snow.STATVAR.energy = [snow.STATVAR.energy(1)*fraction; snow.STATVAR.energy(1)*(1-fraction); snow.STATVAR.energy(2:end)];
                
                snow.STATVAR.d = [snow.STATVAR.d(1); snow.STATVAR.d]; 
                snow.STATVAR.s = [snow.STATVAR.s(1); snow.STATVAR.s];
                snow.STATVAR.gs = [snow.STATVAR.gs(1); snow.STATVAR.gs];
                snow.STATVAR.time_snowfall = [snow.STATVAR.time_snowfall(1); snow.STATVAR.time_snowfall];
                
                snow.STATVAR.d(1) = (snow.STATVAR.d(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.d.*snow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.s(1) = (snow.STATVAR.s(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.s.*snow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.gs(1) = (snow.STATVAR.gs(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.gs.*snow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.time_snowfall(1) = (snow.STATVAR.time_snowfall(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.time_snowfall.*snow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                
                snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + newSnow.STATVAR.energy;
                snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce;
                snow.STATVAR.ice(1) = snow.STATVAR.ice(1) + newSnow.STATVAR.ice;
                snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + newSnow.STATVAR.layerThick;
            end
        end
        
    end
end
