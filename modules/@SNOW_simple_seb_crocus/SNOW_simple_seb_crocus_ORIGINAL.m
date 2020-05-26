% inhertits from SNOW_base_class and adds SEB as upper boundary
% designed to function as CHILD of a GROUND class that is compatible with
% SNOW classes; compatible with interaction classes IA_SNOW_GROUND and IA_SNOW_GROUND_fcSimple_salt

classdef SNOW_simple_seb_crocus< SNOW_simple_seb
    
    methods
        
        %mandatory functions for each class
        
        function xls_out = write_excel(~)
            xls_out = {{'CLASS','index',NaN,NaN,NaN;'SNOW_simple_seb_crocus',1,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN;NaN,'value','default','unit',NaN;'density',300,350,'[kg/m3]','snow density';'albedo_max',0.850000000000000,0.850000000000000,'[-]','not active';'albedo_min',0.550000000000000,0.550000000000000,'[-]','not active';'albedo',0.800000000000000,0.800000000000000,'[-]','surface albedo';'epsilon',0.990000000000000,0.990000000000000,'[-]','surface emissivity';'z0',0.000100000000000000,0.000100000000000000,'[m]','roughness length';'field_capacity',0.0500000000000000,0.0500000000000000,'[-]','%fraction of porosity that can be filled with water before draining';'hydraulicConductivity',1,'    ','[m/sec]','    ';'dt_max',3600,3600,'[sec]','longest possible timestep';'dE_max',50000,50000,'[J/m3]','maximum change of energy per timestep';'swe_per_cell',0.0100000000000000,0.0100000000000000,'[m]','target SWE regulating grid cell size, 0.01m is ca. 3cm ';'slope',0,0,'[°]',NaN;'timescale_winddrift',24,48,'[hours]',NaN;'CLASS_END',NaN,NaN,NaN,NaN}};
        end
        
        function snow = provide_variables(snow)  %initializes the subvariables as empty arrays
            snow = provide_variables@SNOW_simple_seb(snow);
            snow = provide_PARA(snow);
            snow = provide_CONST(snow);
            snow = provide_STATVAR(snow);
        end
        
        function variable = initialize_from_file(snow, variable, section)
            variable = initialize_from_file@SNOW_simple_seb(snow, variable, section);
        end
        
        function snow = assign_global_variables(snow, forcing)
            snow = assign_global_variables@SNOW_simple_seb(snow, forcing);
        end
        
        function snow = initialize_zero_snow(snow, parentGround)
            snow = initialize_zero_snow@SNOW_simple_seb(snow, parentGround);
            snow.STATVAR.d = [];
            snow.STATVAR.s = [];
            snow.STATVAR.gs = [];
            snow.STATVAR.time_snowfall = [];
            snow.STATVAR.ice_fraction = [];
            snow.STATVAR.target_density = [];
            snow.STATVAR.energy = [];
        end
        
        function snow = get_boundary_condition_u(snow, forcing) 
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600); %snowfall is in mm/day
            snow.TEMP.rainfall = forcing.TEMP.rainfall ./1000 ./(24.*3600);
            snow.TEMP.wind = forcing.TEMP.wind;
            
            if snow.TEMP.snowfall > 0
                snow = get_snow_properties(snow,forcing);
            end
            
            if ~isempty(snow.STATVAR.time_snowfall) 
                snow = get_albedo_snow(snow, forcing.TEMP.t, forcing.TEMP.p);
                snow = energy_infiltration(snow,forcing.TEMP.Sin);
            else
                snow.TEMP.albedo = snow.PARA.albedo;
                snow.TEMP.d_energy_seb = 0; % No transmission when only fractional snowcover
            end
            
            snow.STATVAR.Sout = snow.TEMP.albedo .*  forcing.TEMP.Sin;
            snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
            snow.STATVAR.Qh = Q_h(snow, forcing);
            snow.STATVAR.Qe = Q_eq_potET(snow, forcing);
            
            snow.TEMP.F_ub = forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Qh - snow.STATVAR.Qe;
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);
            snow.TEMP.rain_energy = snow.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;
            snow.TEMP.d_ice_sublim = -snow.STATVAR.Qe ./ snow.CONST.L_s;
            snow.TEMP.d_E_sublim =  snow.TEMP.d_ice_sublim .* (snow.CONST.c_w .*snow.STATVAR.T(1) - snow.CONST.L_f);
            %snow.TEMP.d_E_sublim = snow.STATVAR.Qe;
            snow.TEMP.T_rainWater = forcing.TEMP.Tair; % same as below
            snow.TEMP.F_ub_water = snow.TEMP.rainfall; % Needed when using the bucketW water infiltration scheme, change advance_prognostics
            
%             fraction_snow = snow.IA_PARENT.FRACTIONAL_SNOW_COVER;
%             snow.TEMP.d_ice_sublim = fraction_snow .* snow.TEMP.d_ice_sublim;
%             snow.TEMP.d_E_sublim = fraction_snow .* snow.TEMP.d_E_sublim;

        end
        
        function snow = get_boundary_condition_u_CHILD(snow, forcing)
             snow = get_boundary_condition_u(snow, forcing); %same function as  for normal snow class
             
             fraction_snow = snow.IA_PARENT.FRACTIONAL_SNOW_COVER; %sublimation must be scaled!
             snow.TEMP.d_ice_sublim = fraction_snow .* snow.TEMP.d_ice_sublim;
             snow.TEMP.d_E_sublim = fraction_snow .* snow.TEMP.d_E_sublim;
        end
        
        function snow = get_boundary_condition_u_create_CHILD(snow, forcing) 
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600); %snowfall is in mm/day
            snow = get_snow_properties(snow,forcing);
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);
            
        end
        
        function snow = get_boundary_condition_l(snow, forcing)
            snow = get_boundary_condition_l@SNOW_simple_seb(snow, forcing);
            snow.TEMP.F_lb_water = 0; % zero flux if used as bottom class
        end
        
        
        function snow = get_derivatives_prognostic(snow)
            snow = get_derivatives_prognostic@SNOW_simple_seb(snow);
            snow.TEMP.d_energy = snow.TEMP.d_energy + snow.TEMP.d_energy_seb;
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.rain_energy;
            snow = get_derivative_water(snow);
            
            snow = get_T_gradient_snow(snow);
            snow = prog_metamorphism(snow);
            snow = prog_wind_drift(snow); % Add blowing snow sublimation according to Gordon etAl 2006
            snow = compaction(snow);
            
        end
        
        function snow = get_derivatives_prognostic_CHILD(snow)   
           
            snow.TEMP.d_energy = snow.TEMP.F_ub + snow.TEMP.snow_energy + snow.TEMP.rain_energy + snow.TEMP.d_E_sublim + snow.TEMP.d_energy_seb;
            %snow.TEMP.d_waterIce = snow.TEMP.snowfall + snow.TEMP.rainfall + snow.TEMP.d_ice_sublim;
  
            snow.TEMP.dT = 0; %assuming zero gradient
            snow = prog_metamorphism(snow);
            snow = prog_wind_drift(snow); 
            %snow = compaction(snow); compaction does not work yet with only one cell, check below
            snow.TEMP.compact_d_D = 0;
        end
        
        
        function timestep = get_timestep(snow)  %could involve check for several state variables
            timestep = get_timestep@SNOW_simple_seb(snow);
        end
        
        
        function timestep = get_timestep_CHILD(snow)  %will be ignored if it has a negative value
            timestep = -snow.STATVAR.energy ./ snow.TEMP.d_energy;
            
        end
         
         
        function snow = advance_prognostic(snow, timestep) %real timestep derived as minimum of several classes
            % All prognistic variables are advanced BEFORE the new snow is
            % added, i.e. rain/sublimation occur at the "old" surface
            
            snow = advance_prognostic@SNOW_base_class(snow, timestep);
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* snow.TEMP.d_E_sublim;  %rainfall energy already added, snowfall energy is in get_new_snow
            snow.STATVAR.waterIce(1)  = snow.STATVAR.waterIce(1) + timestep .* (snow.TEMP.d_ice_sublim); % + snow.TEMP.rainfall) ; %rainfall needs to be removed if advance_prognostic_water is used
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* (snow.TEMP.d_ice_sublim ./ (snow.STATVAR.ice(1) ./ snow.STATVAR.layerThick(1)));
            
            snow.STATVAR.d = max(snow.STATVAR.d.*0, snow.STATVAR.d + timestep .*(snow.TEMP.metam_d_d + snow.TEMP.wind_d_d));
            snow.STATVAR.s = max(snow.STATVAR.s.*0, min(snow.STATVAR.s.*0+1, snow.STATVAR.s + timestep .*(snow.TEMP.metam_d_s + snow.TEMP.wind_d_s)));
            snow.STATVAR.gs = max(snow.STATVAR.gs, snow.STATVAR.gs + timestep .*(snow.TEMP.metam_d_gs + snow.TEMP.wind_d_gs));
            snow.STATVAR.layerThick = min(snow.STATVAR.layerThick, max(snow.STATVAR.ice, snow.STATVAR.layerThick + timestep .*(snow.TEMP.compact_d_D + snow.TEMP.wind_d_D)));
            
            snow = advance_prognostic_water(snow, timestep);
            
            % Add dry snow on top
            if snow.TEMP.snowfall > 0
                newSnow = get_new_snow(snow,timestep);
                snow = combineNewSnowLayer(snow,newSnow);
            end
            
            snow.STATVAR.target_density = min(1, snow.STATVAR.ice ./ snow.STATVAR.layerThick);
            snow.STATVAR.target_density(1) = min(1, (snow.STATVAR.ice(1) + timestep .* snow.TEMP.d_ice_sublim) ./ snow.STATVAR.layerThick(1));
            
            if sum(snow.STATVAR.layerThick<=0)~=0 || sum(snow.STATVAR.waterIce<=0)~=0
                dff
            end
            
        end
        
        function snow = advance_prognostic_create_CHILD(snow, timestep)
            newSnow = get_new_snow(snow,timestep);
            snow.STATVAR.d = newSnow.STATVAR.d;
            snow.STATVAR.s = newSnow.STATVAR.s;
            snow.STATVAR.gs = newSnow.STATVAR.gs;
            snow.STATVAR.time_snowfall = newSnow.STATVAR.time_snowfall;
            snow.STATVAR.energy = newSnow.STATVAR.energy;
            snow.STATVAR.waterIce = newSnow.STATVAR.waterIce;
            snow.STATVAR.ice =  newSnow.STATVAR.ice;
            snow.STATVAR.layerThick = newSnow.STATVAR.layerThick;
            snow.STATVAR.water = 0;
            snow.STATVAR.target_density = snow.STATVAR.ice ./ snow.STATVAR.layerThick;
            
     
            
        end
        
        function snow = advance_prognostic_CHILD(snow, timestep)
            
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* snow.TEMP.d_energy;
            
            if snow.STATVAR.waterIce < 0
                test2
            end
            
            %snow.STATVAR.energy = snow.STATVAR.energy + timestep .*(snow.TEMP.F_ub + snow.TEMP.snow_energy + snow.TEMP.rain_energy + snow.TEMP.d_E_sublim);
            snow.STATVAR.waterIce  = snow.STATVAR.waterIce + timestep .* (snow.TEMP.d_ice_sublim + snow.TEMP.rainfall); %rainfall added directly here, snowfall added below
            snow.STATVAR.layerThick = snow.STATVAR.layerThick + timestep .* (snow.TEMP.d_ice_sublim ./ (snow.STATVAR.ice ./ snow.STATVAR.layerThick));
            
            snow.STATVAR.d = max(snow.STATVAR.d.*0, snow.STATVAR.d + timestep .*(snow.TEMP.metam_d_d + snow.TEMP.wind_d_d));
            snow.STATVAR.s = max(snow.STATVAR.s.*0, min(snow.STATVAR.s.*0+1, snow.STATVAR.s + timestep .*(snow.TEMP.metam_d_s + snow.TEMP.wind_d_s)));
            snow.STATVAR.gs = max(snow.STATVAR.gs, snow.STATVAR.gs + timestep .*(snow.TEMP.metam_d_gs + snow.TEMP.wind_d_gs));
            snow.STATVAR.layerThick = min(snow.STATVAR.layerThick, max(snow.STATVAR.ice, snow.STATVAR.layerThick + timestep .*(snow.TEMP.compact_d_D + snow.TEMP.wind_d_D)));
            
           if snow.STATVAR.waterIce < 0
                test1
            end
            
            %add dry snow
            if snow.TEMP.snowfall > 0
                newSnow = get_new_snow(snow,timestep);
                
                snow.STATVAR.d(1) = (snow.STATVAR.d(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.d.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.s(1) = (snow.STATVAR.s(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.s.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.gs(1) = (snow.STATVAR.gs(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.gs.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.time_snowfall(1) = (snow.STATVAR.time_snowfall(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.time_snowfall.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                
                %snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + newSnow.STATVAR.energy;
                snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce;
                snow.STATVAR.ice(1) = snow.STATVAR.ice(1) + newSnow.STATVAR.ice;
                snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + newSnow.STATVAR.layerThick;
            end
            snow.STATVAR.target_density = min(1, (snow.STATVAR.ice + timestep .* snow.TEMP.d_ice_sublim) ./ snow.STATVAR.layerThick);

            
        end
        
        
        function snow = compute_diagnostic_first_cell(snow, forcing)
            snow = compute_diagnostic_first_cell@SNOW_simple_seb(snow, forcing);
        end
        
        function snow = compute_diagnostic(snow, forcing)
%           snow = compute_diagnostic@SNOW_simple_seb(snow, forcing);
            snow = compute_diagnostic@SNOW_base_class(snow, forcing); % changed 171019
            snow = check_trigger(snow); %checks if snow goes back to CHILD

            snow.STATVAR.upperPos = snow.STATVAR.lowerPos + sum(snow.STATVAR.layerThick);
            
        end
        
        function snow = compute_diagnostic_CHILD(snow, forcing)
            
            snow = conductivity(snow);
            snow = get_T_water(snow);
            snow.STATVAR.upperPos = snow.STATVAR.lowerPos + snow.STATVAR.layerThick;

            if snow.STATVAR.waterIce <1e-15   %resets STATUS back to zero
                snow.IA_PARENT.IA_PARENT_GROUND.IA_CHILD.STATUS = 0;
                snow = initialize_zero_snow(snow, snow.IA_PARENT.IA_PARENT_GROUND); %set all variables to zero
            end
            
            
            if snow.STATVAR.ice >= snow.PARA.swe_per_cell./2
                
                snow.IA_PARENT.IA_PARENT_GROUND.PREVIOUS.NEXT = snow; 
                snow.PREVIOUS = snow.IA_PARENT.IA_PARENT_GROUND.PREVIOUS;
                snow.NEXT = snow.IA_PARENT.IA_PARENT_GROUND;
                snow.IA_PARENT.IA_PARENT_GROUND.PREVIOUS = snow;
                snow.IA_PARENT.IA_PARENT_GROUND.IA_CHILD.STATUS = -1;
                snow.IA_PARENT.IA_PARENT_GROUND.IA_CHILD.FRACTIONAL_SNOW_COVER = 0;

                
                snow.IA_NEXT = get_IA_class(class(snow.NEXT), class(snow));
                snow.IA_PARENT.IA_PARENT_GROUND.IA_PREVIOUS = snow.IA_NEXT;
                snow.IA_NEXT.PREVIOUS = snow;
                snow.IA_NEXT.NEXT = snow.IA_PARENT.IA_PARENT_GROUND;
                
                %snow.IA_PARENT.IA_PARENT_GROUND.IA_CHILD.IA_CHILD_SNOW = [];
                snow.NEXT.IA_CHILD.IA_CHILD_SNOW = [];  %does not work yet to cut the connection between ground CHILD and snow
                snow.IA_PARENT = [];
            end
            % checks if snow CHILD needs to become full
            %snow class and rearrange the stratigraphy

        end
        
        
        function ground = troubleshoot(ground)
            ground = checkNaN(ground);
        end

        
        % -----------------   non-mandatory  ------------------
        
        function snow = conductivity(snow)   
            snow = conductivity@SNOW_base_class(snow);
        end
            
        function snow = get_T_water(snow)
            snow = get_T_water@SNOW_base_class(snow);
        end
        
        function snow = get_snow_properties(snow,forcing)
            
            T_air=min(forcing.TEMP.Tair, 0);
            windspeed = forcing.TEMP.wind;
            
            T_fus=0;  %degree C
            a_rho=109;  %kg/m3
            b_rho=6; % kg/m3K
            c_rho=26; % kg m7/2s1/2
            min_snowDensity=50; %kg/m3
            
            snow.TEMP.rho_falling_snow = max(min_snowDensity, a_rho+b_rho.*(T_air-T_fus)+ c_rho.*windspeed.^0.5);  %initial snow density
            snow.TEMP.d=min(max(1.29-0.17.*windspeed,0.2),1);
            snow.TEMP.s=min(max(0.08.*windspeed + 0.38,0.5),0.9);
            snow.TEMP.gs=0.1e-3+(1-snow.TEMP.d).*(0.3e-3-0.1e-3.*snow.TEMP.s);
            snow.TEMP.time_snowfall = forcing.TEMP.t;
            
        end
        
        function snow = get_albedo_snow(snow, t, pressure)
            
            rho=snow.STATVAR.waterIce./snow.STATVAR.layerThick.*1000; % WET DENSITY ?? Check Vionet 2012/Sebastian
            d=snow.STATVAR.d;
            s=snow.STATVAR.s;
            gs=snow.STATVAR.gs;
            snow_age=t-snow.STATVAR.time_snowfall(1);
            
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
            
            snow.TEMP.albedo = min(max(snow.PARA.albedo_min,total_albedo),snow.PARA.albedo_max);
            snow.TEMP.spectral_ranges = spectral_ranges;
            snow.TEMP.spectral_albedo = spectral_albedo;
            snow.TEMP.damping_factor = damping_factor;
            
        end
        
        function snow = energy_infiltration(snow,Sin)
            Q_solar=(1-snow.TEMP.spectral_albedo).*snow.TEMP.spectral_ranges.*Sin;
            residual=0.1;
            i=1;
            d_energy = snow.STATVAR.energy.*0;
            
            while i<=size(snow.STATVAR.layerThick,1) && sum(Q_solar)>=residual
                reduction=exp(-snow.TEMP.damping_factor(i,:)'.*snow.STATVAR.layerThick(i));  %vector
                d_energy(i) = sum(Q_solar .* (1-reduction));
                Q_solar=reduction.*Q_solar;  %vector
                
                i=i+1;
            end
            
            if sum(Q_solar) < residual && i<=size(snow.STATVAR.layerThick,1)
                d_energy(i) = sum(Q_solar);
            elseif sum(Q_solar) >= residual
                d_energy(end) = d_energy(end) + sum(Q_solar);
            end
            
            snow.TEMP.d_energy_seb = d_energy;
        end
        
        function snow = get_T_gradient_snow(snow)
            
            delta=[snow.STATVAR.layerThick; snow.NEXT.STATVAR.layerThick(1,1)];  
            T = [snow.STATVAR.T; snow.NEXT.STATVAR.T(1,1)];
            dT = (T(1:end-2) - T(3:end))./(0.5.*delta(1:end-2) +  delta(2:end-1) + 0.5.*delta(3:end));
            dT =abs([ (T(1) - T(2))./(0.5.*delta(1) + 0.5.*delta(2)) ; dT]);
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
            
            daysec = 60*60*24;
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
        
        function snow = prog_wind_drift(snow)
            timescale = snow.PARA.timescale_winddrift;
            rho = snow.STATVAR.waterIce./snow.STATVAR.layerThick .*1000;
            d = snow.STATVAR.d;
            s = snow.STATVAR.s;
            gs = snow.STATVAR.gs;
            D = snow.STATVAR.layerThick;
            
            F = 1.25-0.0042.*(max(50,rho)-50);
            Mo = double(d>0).*(0.34.*(0.75.*d-0.5.*s+0.5)+0.66.*F) + double(d<=0).*(0.34.*(-0.583.*gs.*1000-0.833.*s+0.833)+0.66.*F);
            Si=-2.868 .* exp(-0.085.*snow.TEMP.wind)+1+Mo;
            
            zi=Si.*0;
            i=2;
            while i<=size(Si,1) && Si(i)>0
                zi(i)=sum(D(1:i-1).*(3.25-Si(1:i-1)));
                i=i+1;
            end
            one_over_tau=max(0,Si.*exp(-zi./0.1))./(timescale.*3600);
            
            snow.TEMP.wind_d_d = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* -d./2.*one_over_tau;
            snow.TEMP.wind_d_s = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* (1-s).*one_over_tau;
            snow.TEMP.wind_d_gs = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* double(d==0).*5e-4.*one_over_tau;
            d_rho = double(snow.STATVAR.energy < -snow.STATVAR.waterIce.*snow.CONST.L_f) .* double(rho<350).* (350-rho).*one_over_tau;
            snow.TEMP.wind_d_D = -D./rho .*d_rho;
            snow.TEMP.Si=Si;
            snow.TEMP.one_over_tau = one_over_tau;
        end
        
        function snow = compaction(snow)
            rho = snow.STATVAR.waterIce./snow.STATVAR.layerThick .*1000;
            gs = snow.STATVAR.gs;
            D = snow.STATVAR.layerThick;
            W_liq = snow.STATVAR.water*1000;
            T = snow.STATVAR.T;
            
            sigma=9.81.*cosd(snow.PARA.slope).*rho.*D;
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
        
        function newSnow = get_new_snow(snow,timestep)
            
            newSnow.STATVAR.d = snow.TEMP.d;
            newSnow.STATVAR.s = snow.TEMP.s;
            newSnow.STATVAR.gs = snow.TEMP.gs;
            newSnow.STATVAR.time_snowfall = snow.TEMP.time_snowfall;
            newSnow.STATVAR.waterIce = timestep .* snow.TEMP.snowfall;
            newSnow.STATVAR.ice = timestep .* snow.TEMP.snowfall;
            %newSnow.STATVAR.target_density = snow.TEMP.rho_falling_snow./1000;
            newSnow.STATVAR.layerThick = timestep .* (snow.TEMP.snowfall ./ (snow.TEMP.rho_falling_snow ./1000));
            newSnow.STATVAR.energy = timestep .* snow.TEMP.snow_energy;
        end
        
        function snow = combineNewSnowLayer(snow,newSnow)
            if (snow.STATVAR.ice(1) + newSnow.STATVAR.ice) < 1.5*snow.PARA.swe_per_cell % Layers can be merged
                % Linear mixing
                snow.STATVAR.d(1) = (snow.STATVAR.d(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.d.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.s(1) = (snow.STATVAR.s(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.s.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.gs(1) = (snow.STATVAR.gs(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.gs.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.time_snowfall(1) = (snow.STATVAR.time_snowfall(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.time_snowfall.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
               % snow.STATVAR.target_density(1) = (snow.STATVAR.target_density(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.target_density.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                % Additive
                snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + newSnow.STATVAR.energy;
                snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce;
                snow.STATVAR.ice(1) = snow.STATVAR.ice(1) + newSnow.STATVAR.ice;
                snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + newSnow.STATVAR.layerThick;
            else % split layers before adding snow (as above)
                fraction = (snow.STATVAR.ice(1)-snow.PARA.swe_per_cell)./snow.STATVAR.ice(1);
                snow.STATVAR.waterIce = [snow.STATVAR.waterIce(1)*fraction; snow.STATVAR.waterIce(1)*(1-fraction); snow.STATVAR.waterIce(2:end)];
                snow.STATVAR.water = [snow.STATVAR.water(1)*fraction; snow.STATVAR.water(1)*(1-fraction); snow.STATVAR.water(2:end)];
                snow.STATVAR.ice = [snow.STATVAR.ice(1)*fraction; snow.STATVAR.ice(1)*(1-fraction); snow.STATVAR.ice(2:end)];
                snow.STATVAR.layerThick = [snow.STATVAR.layerThick(1)*fraction; snow.STATVAR.layerThick(1)*(1-fraction); snow.STATVAR.layerThick(2:end)];
                snow.STATVAR.energy = [snow.STATVAR.energy(1)*fraction; snow.STATVAR.energy(1)*(1-fraction); snow.STATVAR.energy(2:end)];
                
                snow.STATVAR.d = [snow.STATVAR.d(1); snow.STATVAR.d];
                snow.STATVAR.s = [snow.STATVAR.s(1); snow.STATVAR.s];
                snow.STATVAR.gs = [snow.STATVAR.gs(1); snow.STATVAR.gs];
                snow.STATVAR.time_snowfall = [snow.STATVAR.time_snowfall(1); snow.STATVAR.time_snowfall];
              %  snow.STATVAR.target_density = [snow.STATVAR.target_density(1); snow.STATVAR.target_density];
                
                snow.STATVAR.d(1) = (snow.STATVAR.d(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.d.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.s(1) = (snow.STATVAR.s(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.s.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.gs(1) = (snow.STATVAR.gs(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.gs.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                snow.STATVAR.time_snowfall(1) = (snow.STATVAR.time_snowfall(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.time_snowfall.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                %snow.STATVAR.target_density(1) = (snow.STATVAR.target_density(1).*snow.STATVAR.waterIce(1) + newSnow.STATVAR.target_density.*newSnow.STATVAR.waterIce)./(snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce);
                
                snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + newSnow.STATVAR.energy;
                snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + newSnow.STATVAR.waterIce;
                snow.STATVAR.ice(1) = snow.STATVAR.ice(1) + newSnow.STATVAR.ice;
                snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + newSnow.STATVAR.layerThick;
            end
        end
        
        
        function snow = get_derivative_water(snow)
            snow.TEMP.F_lb_water =0; %CHANGE LATER
            
            saturation = snow.STATVAR.water ./ max(1e-12, snow.STATVAR.layerThick - snow.STATVAR.ice);
            waterMobile = double(saturation > snow.PARA.field_capacity);
            snow.TEMP.d_water_out = waterMobile .* snow.PARA.hydraulicConductivity .* snow.STATVAR.water ./ snow.STATVAR.layerThick;
            snow.TEMP.d_water_in = snow.TEMP.d_water_out .*0;
            snow.TEMP.d_water_in(2:end,1) = snow.TEMP.d_water_out(1:end-1,1);
            snow.TEMP.d_water_in(1,1) = snow.TEMP.F_ub_water;
            if ~isempty(snow.IA_NEXT) % Must be a nicer way to do this!
                if strcmp(snow.IA_NEXT,'IA_HEAT_WATER_SNOW_GROUND')
                    get_boundary_condition_water_m(snow.IA_NEXT); 
                else
                    snow.TEMP.d_water_out(end,1) = snow.TEMP.F_lb_water; %CHANGE LATER
                end
            else
                snow.TEMP.d_water_out(end,1) = snow.TEMP.F_lb_water; %CHANGE LATER
            end
        end
        
        function snow = advance_prognostic_water(snow, timestep)
            %snow.STATVAR.water  = snow.STAVAR.water + snow.TEMP.d_water_ET .* timestep; %add this function when SEB is present
            
            snow.TEMP.d_water_in = snow.TEMP.d_water_in .* timestep;
            snow.TEMP.d_water_out = snow.TEMP.d_water_out .* timestep;
            %limit outflow to field capacity
            snow.TEMP.d_water_out  = min(snow.TEMP.d_water_out, max(0, snow.STATVAR.water - snow.PARA.field_capacity .* (snow.STATVAR.layerThick -  snow.STATVAR.ice)));
            snow.TEMP.d_water_in(2:end,1) = snow.TEMP.d_water_out(1:end-1,1);
            %limit inflow so that unity is not exceeded
            snow.TEMP.d_water_in = max(0, min(snow.TEMP.d_water_in, snow.STATVAR.layerThick  - snow.STATVAR.waterIce));
            snow.TEMP.d_water_out(1:end-1,1) = snow.TEMP.d_water_in(2:end,1);
            
            %            finalize_boundary_condition_water_m(snow.IA_NEXT);
            
            energy_out = snow.CONST.c_w .* snow.STATVAR.T .* snow.TEMP.d_water_out;
            energy_in = energy_out.*0;
            energy_in (2:end,1) = energy_out(1:end-1,1);
            energy_in(1,1) = snow.CONST.c_w .* snow.TEMP.T_rainWater .* snow.TEMP.d_water_in(1,1);
            
            snow.STATVAR.waterIce = snow.STATVAR.waterIce - snow.TEMP.d_water_out + snow.TEMP.d_water_in;
            snow.STATVAR.energy = snow.STATVAR.energy - energy_out + energy_in;
            
        end
        
        function snow = check_trigger(snow)
            if size(snow.STATVAR.energy,1) ==1 && snow.STATVAR.ice < snow.PARA.swe_per_cell/2
                
                ground = snow.NEXT;
                
                ground.PREVIOUS = snow.PREVIOUS; %reassign ground
                ground.PREVIOUS.NEXT = ground;
                ground.IA_PREVIOUS=[];
                
                snow.NEXT =[]; %reassign snow -> cut all dependencies
                snow.PREVIOUS =[];
                snow.IA_NEXT =[];
                snow.IA_PREVIOUS =[];
                
                %ground.IA_CHILD = IA_SNOW_GROUND_crocus();  %reinitialize interaction class
                ground.IA_CHILD.STATUS = 2; %snow initially active
                %ground.IA_CHILD.IA_PARENT_GROUND = ground;  %attach snow and ground to interaction class
                ground.IA_CHILD.IA_CHILD_SNOW = snow;
                ground.IA_CHILD.IA_CHILD_SNOW.IA_PARENT = ground.IA_CHILD;
                
                snow = ground; %assign snow pointer to ground to return to regular stratigraphy
            end
        end
        
    end
end