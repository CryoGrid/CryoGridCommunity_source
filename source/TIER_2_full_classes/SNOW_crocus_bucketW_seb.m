%========================================================================
% CryoGrid GROUND class SNOW_crocus_bucketW_seb
% CROCUS snow model Vionnet et al., 2012, but with simpler layer splitting and regridding scheme compared to CROCUS 
% temperature and windspeed-dependent initial snow density, snow microstructure (dendricity, sphericity, grain size), 
% compaction, sublimation, water flow, refreezing, variable albedo.
% Meltwater exceeding the available pore space within the snow cover is automatically removed.
% R. Zweigel, S. Westermann, October 2020
%========================================================================

classdef SNOW_crocus_bucketW_seb < SEB & HEAT_CONDUCTION & WATER_FLUXES & HEAT_FLUXES_LATERAL & WATER_FLUXES_LATERAL & SNOW & SNOW_FLUXES_LATERAL  & REGRID

    properties
        PARENT
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function snow = provide_PARA(snow)

            snow.PARA.epsilon = []; %surface emissivity [-]
            snow.PARA.z0 = [];      %roughness length [m]
            
            snow.PARA.SW_spectral_range1 = []; %fraction of incoming short-wave radiation in first spectral band [-], see Vionnet et al.,2012
            snow.PARA.SW_spectral_range2 = []; %fraction of incoming short-wave radiation in second spectral band [-], fraction of third spectral band calculated automatically
            
            snow.PARA.crocus_version = []; % 'normal' or 'arctic' (Royer et al., 2021) 
            snow.PARA.field_capacity = []; %snow field capacity in fraction of available pore space [-] NOTE: the definition is different for GROUND_XX classes
            snow.PARA.hydraulicConductivity = []; %hydraulic conductivity of snow [m/sec]
            snow.PARA.swe_per_cell = [];  %target SWE per grid cell [m]
            
%             snow.PARA.slope = [];  %slope angle [-]
            snow.PARA.timescale_winddrift = []; %timescale of snow compaction for wind drift [hours!!]
            snow.PARA.max_wind_slab_density = [];
            snow.PARA.wind_factor_fresh_snow = [];
            snow.PARA.albedo_age_factor = [];
            
            snow.PARA.conductivity_function = [];
            snow.PARA.dt_max = [];  %maximum possible timestep [sec]
            snow.PARA.dE_max = [];  %maximum possible energy change per timestep [J/m3]
            
            snow.PARA.snow_property_function = []; % function to be used for new snoe properties
            snow.PARA.conductivity_function = []; % function to be used for snow thermal conductivity
        end
        
        function snow = provide_STATVAR(snow)
            
            snow.STATVAR.upperPos = [];  % upper surface elevation [m]
            snow.STATVAR.lowerPos = [];  % lower surface elevation [m]
            snow.STATVAR.layerThick = [];  % thickness of grid cells [m]
            snow.STATVAR.area = [];     % grid cell area [m2]
            
            snow.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            snow.STATVAR.mineral = [];  % total volume of minerals [m3]
            snow.STATVAR.organic = [];  % total volume of organics [m3]
            snow.STATVAR.energy = [];   % total internal energy[J]
            
            snow.STATVAR.T = [];      % temperature [degree C]
            snow.STATVAR.water = [];  % total volume of water [m3]
            snow.STATVAR.ice = [];    %total volume of ice [m3]
            snow.STATVAR.air = [];    % total volume of air [m3] - NOT USED
            snow.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            snow.STATVAR.hydraulicConductivity = []; %hydraulic conductivity of snow [m/sec]
            snow.STATVAR.albedo = []; %snow albedo [-]
            
            snow.STATVAR.d = []; %dendricity [-], range from 1 (dendric, i.e. crystal-shaped snow particles) to 0 (non-dendric, i.e. round snow particles) 
            snow.STATVAR.s = []; %sphericity [-], range from 1 (round snow particles) to 0 (elongated snow particles) 
            snow.STATVAR.gs = [];  % snow grain size [m]
            snow.STATVAR.time_snowfall = []; % average time of snowfall of a layers (i.e. grid cell) [Matlab time, days]
            snow.STATVAR.target_density = []; %ice fraction prior to  melt in diagnostic step [-]
            snow.STATVAR.excessWater = []; %water volume exceeding snow pore space [m3]
            
            snow.STATVAR.Lstar = []; %Obukhov length [m]
            snow.STATVAR.Qh = [];    %sensible heat flux [W/m2]
            snow.STATVAR.Qe = [];    %latent heat flux [W/m2]
        end
    
        function snow = provide_CONST(snow)
            
            snow.CONST.L_f = []; % volumetric latent heat of fusion, freezing
            
            snow.CONST.c_w = []; % volumetric heat capacity water
            snow.CONST.c_i = []; % volumetric heat capacity ice
            snow.CONST.c_o = []; % volumetric heat capacity organic
            snow.CONST.c_m = []; % volumetric heat capacity mineral
            
            snow.CONST.k_a = [];  % thermal conductivity air
            snow.CONST.k_w = [];  % thermal conductivity water
            snow.CONST.k_i = [];  % thermal conductivity ice 
            snow.CONST.k_o = [];  % thermal conductivity organic 
            snow.CONST.k_m = [];  % thermal conductivity mineral 
            
            snow.CONST.sigma = []; %Stefan-Boltzmann constant
            snow.CONST.kappa = []; % von Karman constant
            snow.CONST.L_s = [];  %latent heat of sublimation, evaporation handled in a dedicated function
            
            snow.CONST.cp = [];  % specific heat capacity at constant pressure of air
            snow.CONST.g = [];   % gravitational acceleration Earth surface
            snow.CONST.Tmfw = []; % Freezing temperature of water in C
            
            snow.CONST.rho_w = []; % water density
            snow.CONST.rho_i = []; % ice density
        end
        
        function snow = finalize_init(snow, tile) %assign all variables, that must be calculated or assigned otherwise for initialization
            if isempty(snow.PARA.conductivity_function) || sum(isnan(snow.PARA.conductivity_function))>0
                snow.PARA.conductivity_function = 'conductivity_snow_Yen';
            end
            
            if ~isempty(snow.PARA.crocus_version) || sum(isnan(snow.PARA.crocus_version))==0
                if strcmp(snow.PARA.crocus_version, 'normal')
                    snow.PARA.conductivity_function = 'conductivity_snow_Yen';
                    snow.PARA.timescale_winddrift = 48; %timescale of snow compaction for wind drift [hours!!]
                    snow.PARA.max_wind_slab_density = 350;
                    snow.PARA.wind_factor_fresh_snow = 26;
                elseif strcmp(snow.PARA.crocus_version, 'arctic')
                    snow.PARA.conductivity_function = 'conductivity_snow_Sturm';
                    snow.PARA.timescale_winddrift = 48./3; %timescale of snow compaction for wind drift [hours!!]
                    snow.PARA.max_wind_slab_density = 600;
                    snow.PARA.wind_factor_fresh_snow = 26.*2;
                end
            end
            if isempty(snow.PARA.albedo_age_factor) || sum(isnan(snow.PARA.albedo_age_factor))>0
                snow.PARA.albedo_age_factor = min(1,max(mean(tile.FORCING.DATA.p)./870e2, 0.5)).*0.2./60; %as in Crocus, Vionnet et al., 2011
            end
            
            
            snow.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            snow.PARA.airT_height = tile.FORCING.PARA.airT_height;
            snow.PARA.slope = tile.FORCING.PARA.slope_angle;
            
            snow = initialize_zero_snow_BASE(snow);  %initialize all values to be zero
            snow.PARA.spectral_ranges = [snow.PARA.SW_spectral_range1 snow.PARA.SW_spectral_range2 1 - snow.PARA.SW_spectral_range1 - snow.PARA.SW_spectral_range2];
            
            snow.TEMP.d_energy = snow.STATVAR.energy .*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
            
        end

        %---time integration------
        %separate functions for CHILD phase of snow cover
        
        function snow = get_boundary_condition_u(snow, tile) 
            forcing = tile.FORCING;
            snow = get_boundary_condition_SNOW_u(snow, forcing);
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            
            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing);    
            
            snow = surface_energy_balance(snow, forcing);
%             snow = get_sublimation(snow, forcing);
            
            snow.TEMP.wind = forcing.TEMP.wind;
            snow.TEMP.wind_surface = forcing.TEMP.wind;
            snow.STATVAR.Sin = forcing.TEMP.Sin;
        end
        
        function snow = get_boundary_condition_u_CHILD(snow, tile)
            forcing = tile.FORCING;
            snow = get_boundary_condition_allSNOW_rain_u(snow, forcing); %add full snow, but rain only for snow-covered part
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            
            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing);
            
            snow = surface_energy_balance(snow, forcing); %this works including penetration of SW radiation through the CHILD snow
         %   snow = get_sublimation(snow, forcing);
            
            snow.TEMP.wind = forcing.TEMP.wind;
            snow.TEMP.wind_surface = forcing.TEMP.wind;
            snow.STATVAR.Sin = forcing.TEMP.Sin;
        end
        
        function snow = get_boundary_condition_u_create_CHILD(snow, tile)
            forcing = tile.FORCING;
            snow = get_boundary_condition_allSNOW_u(snow, forcing); %add all snow, no rain
            
            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing);
            
            snow.TEMP.F_ub = 0;
            snow.TEMP.F_lb = 0;
            snow.TEMP.F_ub_water = 0;
            snow.TEMP.F_lb_water = 0;
            snow.TEMP.F_ub_water_energy = 0;
            snow.TEMP.F_lb_water_energy = 0;
          
          % Figure out a consistent naming and storage scheme for latent fluxes with Sebastian
          % Currently we are have two competing schemes
            snow.STATVAR.sublimation = 0;
            snow.TEMP.sublimation_energy = 0;
            snow.STATVAR.evap = 0;
            snow.TEMP.evap_energy = 0;
            snow.TEMP.sublim = 0;
            snow.TEMP.sublim_energy = 0;

            snow.TEMP.rain_energy = 0;
            snow.TEMP.rainfall = 0;
            
            snow.TEMP.d_energy = 0;
            snow.TEMP.d_water = 0;
            snow.TEMP.d_water_energy = 0;
            
            snow.STATVAR.d = 0;
            snow.STATVAR.s = 0;
            snow.STATVAR.gs = 0;
            snow.STATVAR.time_snowfall = 0;
            %new
            snow.STATVAR.top_snow_date = forcing.TEMP.t;
            snow.STATVAR.bottom_snow_date = forcing.TEMP.t-0.1./24;
            
            
            snow.TEMP.metam_d_d = 0;
            snow.TEMP.wind_d_d = 0;
            snow.TEMP.metam_d_s = 0;
            snow.TEMP.wind_d_s = 0;
            snow.TEMP.metam_d_gs = 0;
            snow.TEMP.compact_d_D = 0;
            snow.TEMP.wind_d_D = 0;           
            snow.TEMP.wind = forcing.TEMP.wind;
            snow.TEMP.wind_surface = forcing.TEMP.wind;
            
            %start with  non-zero values for area and layerThick
            snow.STATVAR.area = 1;
            snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ (snow.TEMP.newSnow.STATVAR.density ./1000); %[m] layerThick adjusted so that always 0.5 .* snow.PARA.swe_per_cell is contained
            snow.STATVAR.energy = 0;
            snow.STATVAR.waterIce = 0;
            snow.STATVAR.ice = 0;
            snow.STATVAR.excessWater = 0;
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos;

        end
        
        function snow = get_boundary_condition_m_CHILD(snow, tile)
            forcing = tile.FORCING;
            snow = get_boundary_condition_allSNOW_rain_canopy_m(snow, tile);%add full snow, but rain only for snow-covered part
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            
            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing); %makes a TEMP variable newSnow that contains all information on the fresh snow - which is merged in the diagnostic step
            
            snow.STATVAR.Qh = get_boundary_condition_m_Qh_CLM5(snow, forcing);
            snow.STATVAR.Qe = Q_e_CLM5_snow(snow, forcing);
            snow = get_sublimation(snow, forcing);
            
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + (-snow.STATVAR.Qh - snow.STATVAR.Qe);
        end
        
        function [snow, S_up] = penetrate_SW(snow, S_down)  %mandatory function when used with class that features SW penetration
            [snow, S_up] = penetrate_SW_transmission_spectral(snow, S_down);
        end
        
        function [snow, L_up] = penetrate_LW(snow, L_down)
            [snow, L_up] = penetrate_LW_no_transmission(snow, L_down);
        end
        
        function snow = get_boundary_condition_l(snow, tile)
            forcing = tile.FORCING;
            
            snow.TEMP.F_lb = forcing.PARA.heatFlux_lb .* snow.STATVAR.area(end);
            snow.TEMP.d_energy(end) = snow.TEMP.d_energy(end) + snow.TEMP.F_lb;
            
            snow.TEMP.F_lb_water = 0;
            snow.TEMP.F_lb_water_energy = 0;
            
        end
        
        function snow = get_derivatives_prognostic(snow, tile)
            if size(snow.STATVAR.layerThick,1) > 1
                snow = get_derivative_energy(snow);
                snow = get_derivative_water_SNOW(snow);
                snow = get_T_gradient_snow(snow);
            else
                snow = get_T_gradient_snow_single_cell(snow);
            end
            snow = prog_metamorphism(snow);
            snow = prog_wind_drift(snow); % possibly add blowing snow sublimation according to Gordon et al. 2006
            snow = compaction(snow);
            
        end
        
        function snow = get_derivatives_prognostic_CHILD(snow, tile)
            
            snow = get_T_gradient_snow_single_cell(snow);
            snow = prog_metamorphism(snow);
            snow = prog_wind_drift(snow); % possibly add blowing snow sublimation according to Gordon et al. 2006
            snow = compaction(snow);
            
        end
        
        function timestep = get_timestep(snow, tile) 
            timestep = get_timestep_SNOW(snow);
            %timestep1 = get_timestep_heat_coduction(snow);
            %timestep2 = get_timestep_SNOW_mass_balance(snow);
            timestep2 = get_timestep_SNOW_sublimation(snow);
            timestep3 = get_timestep_water_SNOW(snow);

            timestep = min(timestep, timestep2);
            
            timestep = min(timestep, timestep3);
        end
        
        function timestep = get_timestep_CHILD(snow, tile)  
            timestep = get_timestep_SNOW_CHILD(snow);
            
            timestep = min(timestep, get_timestep_SNOW_sublimation(snow));
            timestep = min(timestep, get_timestep_water_SNOW(snow));
            %timestep1 = get_timestep_heat_coduction(snow);
            %timestep2 = get_timestep_SNOW_CHILD(snow);
            %timestep = min(timestep1, timestep2);
            
        end
        
        function snow = advance_prognostic(snow, tile)
            timestep = tile.timestep;
            
            snow.STATVAR.layerThick = min(snow.STATVAR.layerThick, max(snow.STATVAR.ice ./ snow.STATVAR.area, snow.STATVAR.layerThick + timestep .*(snow.TEMP.compact_d_D + snow.TEMP.wind_d_D)));
            %energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy + snow.TEMP.d_water_energy);
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* snow.TEMP.sublim_energy;  %snowfall energy added below, when new snow layer is merged
            %mass
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* snow.TEMP.d_water;
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + timestep .* snow.TEMP.sublim;
            %snow.STATVAR.ice(1) = snow.STATVAR.ice(1) + timestep .* snow.STATVAR.sublim;
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* snow.TEMP.sublim ./snow.STATVAR.area(1,1) ./ (snow.STATVAR.ice(1) ./ snow.STATVAR.layerThick(1) ./snow.STATVAR.area(1));
            
            snow.STATVAR.ice(1) = snow.STATVAR.ice(1) + timestep .* snow.TEMP.sublim;
            %microphysics
            snow.STATVAR.d = max(snow.STATVAR.d.*0, snow.STATVAR.d + timestep .*(snow.TEMP.metam_d_d + snow.TEMP.wind_d_d));
            snow.STATVAR.s = max(snow.STATVAR.s.*0, min(snow.STATVAR.s.*0+1, snow.STATVAR.s + timestep .*(snow.TEMP.metam_d_s + snow.TEMP.wind_d_s)));
            snow.STATVAR.gs = max(snow.STATVAR.gs, snow.STATVAR.gs + timestep .*(snow.TEMP.metam_d_gs + snow.TEMP.wind_d_gs));
            %snow.STATVAR.layerThick = min(snow.STATVAR.layerThick, max(snow.STATVAR.ice ./ snow.STATVAR.area, snow.STATVAR.layerThick + timestep .*(snow.TEMP.compact_d_D + snow.TEMP.wind_d_D)));
            
            %subtract water and energy from evaporation (not checked in
            %get_timestep()
            evap = min(snow.STATVAR.evap.* timestep, snow.STATVAR.water(1,1));
            evap_energy = evap ./ (snow.STATVAR.evap.* timestep) .* snow.TEMP.evap_energy;
            evap_energy(isnan(evap_energy)) = 0;
            snow.STATVAR.waterIce(1,1) = snow.STATVAR.waterIce(1,1) + evap;
            snow.STATVAR.water(1,1) = snow.STATVAR.water(1,1) + evap;
            snow.STATVAR.energy(1,1) = snow.STATVAR.energy(1,1) + evap_energy;
            
            %new snow
            if snow.TEMP.snowfall >0
                snow = advance_prognostic_new_snow_crocus(snow, timestep);
                %merge with uppermost layer
                snow = merge_cells_intensive2(snow, 1, snow.TEMP.newSnow, 1, {'d'; 's'; 'gs'; 'time_snowfall'}, 'ice');
                snow = merge_cells_extensive2(snow, 1, snow.TEMP.newSnow, 1, {'waterIce'; 'energy'; 'layerThick'; 'ice'});
                snow = merge_cells_snowfall_times2(snow, 1, snow.TEMP.newSnow, 1); %specific function merginging bottom and top snow dates
            end
            
            %store "old" density - ice is updated for new snowfall and sublimation losses
            snow.STATVAR.target_density = min(1, snow.STATVAR.ice ./ snow.STATVAR.layerThick ./ snow.STATVAR.area);
            
            %add surface runoff from non infiltrating surface runoff to precip 
            snow.STATVAR.excessWater = snow.STATVAR.excessWater + snow.TEMP.surface_runoff .* timestep;
        end
        
        function snow = advance_prognostic_CHILD(snow, tile)
            timestep = tile.timestep;
            %energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy  + snow.TEMP.d_water_energy + snow.TEMP.sublim_energy);
            %mass
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* (snow.TEMP.d_water + snow.TEMP.sublim);
            %snow.STATVAR.ice = snow.STATVAR.ice + timestep .* snow.STATVAR.sublim;
            snow.STATVAR.volume = snow.STATVAR.layerThick .* snow.STATVAR.area;
            snow.STATVAR.volume = min(snow.STATVAR.volume, max(snow.STATVAR.ice, snow.STATVAR.volume + timestep .* snow.STATVAR.area .* (snow.TEMP.compact_d_D + snow.TEMP.wind_d_D))); %mass is conserved, reduce layerthick
            %snow.STATVAR.volume = snow.STATVAR.volume + timestep .* snow.STATVAR.sublim ./ max(50, snow.STATVAR.ice ./ snow.STATVAR.layerThick); 
            
            %snow.STATVAR.volume = snow.STATVAR.volume + timestep .* snow.STATVAR.sublim ./ (snow.STATVAR.ice ./ snow.STATVAR.layerThick./snow.STATVAR.area); 
            snow.STATVAR.volume = snow.STATVAR.volume + timestep .* snow.TEMP.sublim ./ (snow.STATVAR.ice ./ snow.STATVAR.volume); 
%                         disp('Hallo2')
%             disp(snow.STATVAR.volume)
            
            snow.STATVAR.ice = snow.STATVAR.ice + timestep .* snow.TEMP.sublim;
            
            %microphysics
            snow.STATVAR.d = max(snow.STATVAR.d.*0, snow.STATVAR.d + timestep .*(snow.TEMP.metam_d_d + snow.TEMP.wind_d_d));
            snow.STATVAR.s = max(snow.STATVAR.s.*0, min(snow.STATVAR.s.*0+1, snow.STATVAR.s + timestep .*(snow.TEMP.metam_d_s + snow.TEMP.wind_d_s)));
            snow.STATVAR.gs = max(snow.STATVAR.gs, snow.STATVAR.gs + timestep .*(snow.TEMP.metam_d_gs + snow.TEMP.wind_d_gs));
            %snow.STATVAR.volume = min(snow.STATVAR.volume, max(snow.STATVAR.ice, snow.STATVAR.volume + timestep .* snow.STATVAR.area .* (snow.TEMP.compact_d_D + snow.TEMP.wind_d_D))); %mass is conserved, reduce layerthick
            
            %subtract water and energy from evaporation (not checked in
            %get_timestep()
            evap = min(snow.STATVAR.evap.* timestep, snow.STATVAR.water(1,1));
            evap_energy = evap ./ (snow.STATVAR.evap.* timestep) .* snow.TEMP.evap_energy;
            evap_energy(isnan(evap_energy)) = 0;
            snow.STATVAR.waterIce(1,1) = snow.STATVAR.waterIce(1,1) + evap;
            snow.STATVAR.water(1,1) = snow.STATVAR.water(1,1) + evap;
            snow.STATVAR.energy(1,1) = snow.STATVAR.energy(1,1) + evap_energy;
            
            %new snow and merge
            if snow.TEMP.snowfall >0
                snow = advance_prognostic_new_snow_CHILD_crocus(snow, timestep);  %add new snow with the new layerThick
                %merge with uppermost layer
                snow = merge_cells_intensive2(snow, 1, snow.TEMP.newSnow, 1, {'d'; 's'; 'gs'; 'time_snowfall'}, 'ice');
                snow = merge_cells_extensive2(snow, 1, snow.TEMP.newSnow, 1, {'waterIce'; 'energy'; 'volume'; 'ice'});
                snow = merge_cells_snowfall_times2(snow, 1, snow.TEMP.newSnow, 1); %specific function merginging bottom and top snow dates
            end
            
            %store "old" density - ice is updated for new snowfall and sublimation losses
            snow.STATVAR.target_density = min(1, snow.STATVAR.ice ./ snow.STATVAR.volume);

            
            %adjust layerThick, so that exactly 0.5 .* snow.PARA.swe_per_cell is contained
            snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ snow.STATVAR.target_density;
            snow.STATVAR.area = snow.STATVAR.volume ./ snow.STATVAR.layerThick;
            
%             %add surface runoff from non infiltrating surface runoff to precip
%             snow.STATVAR.excessWater = snow.STATVAR.excessWater + snow.TEMP.surface_runoff .* timestep;
        end
        
        function snow = compute_diagnostic_first_cell(snow, tile)
            forcing = tile.FORCING;
            snow = L_star(snow, forcing);
        end
       
        function snow = compute_diagnostic(snow, tile)
            forcing = tile.FORCING;
            ice_before = max(1e-12, snow.STATVAR.ice(1,1));
            T_before = snow.STATVAR.T;
            snow = get_T_water_freeW(snow);
            melt_fraction = 1- max(0, min(1, snow.STATVAR.ice(1,1)./ice_before));
            if melt_fraction > 0
                snow = split_first_cell_snowfall_times(snow, melt_fraction, 'reduce');
            end
            
%             if tile.t>datenum(2008,6,4)
%                 datestr(tile.t)
%                 'Hei';
%             end
            
            snow = subtract_water2(snow);
            
            [snow, regridded_yesNo] = regrid_snow_crocus(snow, {'waterIce'; 'energy'; 'layerThick'; 'mineral'; 'organic'; 'ice'}, {'area'; 'target_density'; 'd'; 's'; 'gs'; 'time_snowfall'}, 'ice');
            if regridded_yesNo
                snow = get_T_water_freeW(snow);
            end
            
            
            snow = conductivity(snow);
            snow = calculate_hydraulicConductivity_SNOW(snow);
            
            snow.STATVAR.upperPos = snow.NEXT.STATVAR.upperPos + sum(snow.STATVAR.layerThick);
            
            snow = calculate_albedo_crocus(snow, forcing); %albedo calculation is a diagnostic operation
                        
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
        end
        
        function snow = compute_diagnostic_CHILD(snow, tile)
            forcing = tile.FORCING;
            
            ice_before = max(1e-12, snow.STATVAR.ice);
            snow = get_T_water_freeW(snow);
            melt_fraction = 1- max(0, min(1, snow.STATVAR.ice./ice_before));
            if melt_fraction >0
                snow = split_first_cell_snowfall_times(snow, melt_fraction, 'reduce');
            end
            
            snow = subtract_water_CHILD(snow);

            snow = conductivity(snow);
            snow = calculate_hydraulicConductivity_SNOW(snow);
            
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos + (snow.STATVAR.layerThick .* snow.STATVAR.area ./ snow.PARENT.STATVAR.area(1,1));
            
            snow = calculate_albedo_crocus(snow, forcing);
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
            
        end
        
        function snow = check_trigger(snow, tile)
            forcing = tile.FORCING;
            
            snow = make_SNOW_CHILD(snow); 
        end
        
        %-----non-mandatory functions-------
        
        function snow = surface_energy_balance(snow, forcing)
            snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
            
            [snow, S_up] = penetrate_SW(snow, forcing.TEMP.Sin .* snow.STATVAR.area(1) ); % moved splitting of SW to subfunction, to allow for use below other classes. RBZ 25022022
            snow.STATVAR.Sout = sum(S_up) ./ snow.STATVAR.area(1);
            
            snow.STATVAR.Qh = Q_h(snow, forcing);
            snow.STATVAR.Qe = Q_eq_potET_snow(snow, forcing);
            
            snow.TEMP.F_ub = (forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Qh - snow.STATVAR.Qe) .* snow.STATVAR.area(1);
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.F_ub;
        end
        
%         
%         function snow = conductivity(snow)
%             snow = conductivity_snow_Yen(snow);
%         end
        
        function snow = conductivity(snow)
            conductivity_function = str2func(snow.PARA.conductivity_function);
            snow = conductivity_function(snow);
            
        end
        
        %reset timestamp when changing TILES
        function snow = reset_timestamps(snow, tile)
            snow.STATVAR.time_snowfall = snow.STATVAR.time_snowfall - tile.TEMP.time_difference;
            snow.STATVAR.top_snow_date = snow.STATVAR.top_snow_date - tile.TEMP.time_difference;
            snow.STATVAR.bottom_snow_date = snow.STATVAR.bottom_snow_date - tile.TEMP.time_difference;
        end
        
        %-----LATERAL-------------------
        
        %-----LAT_REMOVE_SURFACE_WATER-----
        function snow = lateral_push_remove_surfaceWater(snow, lateral)
            snow = lateral_push_remove_surfaceWater_simple(snow, lateral);
        end
        
        %----LAT_SEEPAGE_FACE----------
        function snow = lateral_push_remove_water_seepage(snow, lateral)
            snow = lateral_push_remove_water_seepage_snow(snow, lateral);
        end

        %----LAT_WATER_RESERVOIR------------        
        function snow = lateral_push_water_reservoir(snow, lateral)
            snow = lateral_push_water_reservoir_snow(snow, lateral);
        end
        
        %----LAT3D_WATER_UNCONFINED_AQUIFER------------        
        function snow = lateral3D_pull_water_unconfined_aquifer(snow, lateral)
            snow = lateral3D_pull_water_unconfined_aquifer_snow(snow, lateral);
        end
        
        function snow = lateral3D_push_water_unconfined_aquifer(snow, lateral)
            snow = lateral3D_push_water_unconfined_aquifer_snow(snow, lateral);
        end
        
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(snow, lateral)
            [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_snow(snow, lateral);
        end
        
        %----LAT3D_WATER_UNCONFINED_AQUIFER_RICHARDS_EQ------------
        function snow = lateral3D_pull_water_unconfined_aquifer_RichardsEq(snow, lateral)
            snow = lateral3D_pull_water_unconfined_aquifer_RichardsEq_snow(snow, lateral);
        end
        
        %LAT3D_WATER_RESERVOIR and LAT3D_WATER_SEEPAGE_FACE do not require specific functions
        
        %----LAT3D_WATER_UNCONFINED_AQUIFER_RICHARDS_EQ------------
        function snow = lateral3D_pull_water_unconfined_aquifer_RichardsEq(snow, lateral)
            snow = lateral3D_pull_water_unconfined_aquifer_RichardsEq_snow(snow, lateral);
        end
        
        %----LAT3D_HEAT------------
        function snow = lateral3D_pull_heat(snow, lateral)
            snow = lateral3D_pull_heat_simple(snow, lateral);
        end
        
        function snow = lateral3D_push_heat(snow, lateral)
            snow = lateral3D_push_heat_simple(snow, lateral);
        end
        
        %----LAT3D_SNOW_CROCUS------------
        function snow = lateral3D_pull_snow(snow, lateral)
            snow = lateral3D_pull_snow_crocus(snow, lateral);
        end
        
        function snow = lateral3D_push_snow(snow, lateral)
            snow = lateral3D_push_snow_crocus(snow, lateral);
        end
        
        %----LAT3D_SNOW_CROCUS_snow_dump------------
        function snow = lateral3D_pull_snow_dump(snow, lateral)
            snow = lateral3D_pull_snow_crocus_dump(snow, lateral);
        end
        
        function snow = lateral3D_push_snow_dump(snow, lateral)
            snow = lateral3D_push_snow_crocus_dump(snow, lateral);
        end
        
        
        
        %----inherited Tier 1 functions ------------
        
        function snow = get_derivative_energy(snow)
           snow = get_derivative_energy@HEAT_CONDUCTION(snow); 
        end
        
%         function snow = conductivity_snow_Yen(snow); Do we need this when
%         class inherits from HEAT_CONDUCTION? RBZ
%             snow = conductivity_snow_Yen@HEAT_CONDUCTION(snow);
%         end
        
        function flux = Q_h(snow, forcing)
           flux = Q_h@SEB(snow, forcing);
        end
    
        function flux = Q_eq_potET(snow, forcing)
            flux = Q_eq_potET@SEB(snow, forcing);
        end
        
        function timestep = get_timestep_heat_coduction(snow)
            timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(snow);
        end
        
        function timestep = get_timestep_SNOW_mass_balance(snow)
            timestep = get_timestep_SNOW_mass_balance@SNOW(snow);
        end
        
        function timestep = get_timestep_SNOW_CHILD(snow)
            timestep = get_timestep_SNOW_CHILD@SNOW(snow);
        end
        
        function snow = L_star(snow, forcing)
           snow = L_star@SEB(snow, forcing); 
        end
        
         function [snow, S_up] = penetrate_SW_transmission_spectral(snow, S_down)
             [snow, S_up] = penetrate_SW_transmission_spectral@SEB(snow, S_down);
         end
        
        function snow = get_T_water_freeW(snow)
            snow = get_T_water_freeW@HEAT_CONDUCTION(snow);
        end
        
        function snow = subtract_water(snow)
            snow = subtract_water@SNOW(snow);
        end
        
        function snow = subtract_water_CHILD(snow)
            snow = subtract_water_CHILD@SNOW(snow);
        end
        
        function snow = get_boundary_condition_SNOW_u(snow, forcing)
            snow = get_boundary_condition_SNOW_u@SNOW(snow, forcing);
        end
        
        function snow = get_boundary_condition_allSNOW_rain_u(snow, forcing)
            snow = get_boundary_condition_allSNOW_rain_u@SNOW(snow, forcing);
        end
        
        function snow = get_boundary_condition_allSNOW_u(snow, forcing) 
            snow = get_boundary_condition_allSNOW_u@SNOW(snow, forcing);
        end
        
        function snow = initialize_zero_snow_BASE(snow)
            snow = initialize_zero_snow_BASE@SNOW(snow);
        end
        
        function snow = make_SNOW_CHILD(snow)
            snow = make_SNOW_CHILD@SNOW(snow);
        end
        
        function [snow, regridded_yesNo] = regrid_snow(snow, extensive_variables, intensive_variables, intensive_scaling_variable)
            [snow, regridded_yesNo] = regrid_snow@REGRID(snow, extensive_variables, intensive_variables, intensive_scaling_variable);
        end
        
        function z0 = get_z0_surface(snow)
            z0 = get_z0_surface@SNOW(snow);
        end
        
        function albedo = get_albedo(snow)
           albedo = snow.STATVAR.albedo;
        end
        
        function Ts = get_surface_T(snow, tile)
            Ts = snow.STATVAR.T(1);
        end
                
        %-------------param file generation-----
         function ground = param_file_info(ground)
             ground = param_file_info@BASE(ground);
             
             ground.PARA.class_category = 'SNOW';
             
             ground.PARA.STATVAR = {''};
             

             ground.PARA.default_value.epsilon = {0.99};
             ground.PARA.comment.epsilon = {'surface emissivity [-]'};
             
             ground.PARA.default_value.z0 = {0.01};
             ground.PARA.comment.z0 = {'roughness length [m]'};

             ground.PARA.default_value.SW_spectral_range1 = {0.71};
             ground.PARA.comment.SW_spectral_range1 = {'fraction of incoming short-wave radiation in first spectral band [-], see Vionnet et al.,2012'};
             
             ground.PARA.default_value.SW_spectral_range2 = {0.21};
             ground.PARA.comment.SW_spectral_range2 = {'fraction of incoming short-wave radiation in second spectral band [-], fraction of third spectral band calculated automatically'};
            
             ground.PARA.default_value.crocus_version = {''};
             ground.PARA.comment.crocus_version = {'normal or arctic (Royer et al., 2021), automatically selects timescale_winddrift, max_wind_slab_density,wind_factor_fresh_snow, and conductivity_fucntion'}; 
             
             ground.PARA.default_value.field_capacity = {0.05};
             ground.PARA.comment.field_capacity = {'snow field capacity in fraction of available pore space [-] NOTE: the definition is different for GROUND_XX classes'};
            
             ground.PARA.default_value.hydraulicConductivity = {1e-4};
             ground.PARA.comment.hydraulicConductivity = {'hydraulic conductivity of snow [m/sec]'};
             
             ground.PARA.default_value.swe_per_cell = {0.02};
             ground.PARA.comment.swe_per_cell = {'target SWE per grid cell [m]'};
             
             ground.PARA.default_value.slope = {0};
             ground.PARA.comment.slope = {'slope angle [-]'};
             
             ground.PARA.default_value.timescale_winddrift = {48};
             ground.PARA.comment.timescale_winddrift = {'timescale of snow compaction for wind drift [hours!!]'};
             
             ground.PARA.default_value.max_wind_slab_density = {350};
             ground.PARA.comment.max_wind_slab_density = {'maximum density achievable by wind compaction'};
             
             ground.PARA.default_value.wind_factor_fresh_snow = {26};
             ground.PARA.comment.wind_factor_fresh_snow = {'factor for wind speed dependency in fresh snow density equation'};
            
             ground.PARA.default_value.conductivity_function = {'conductivity_snow_Yen'};
             ground.PARA.comment.conductivity_function = {'function employed to calculate thermal conductivity, leave empty for default'};
             
             ground.PARA.snow.PARA.albedo_age_factor = {''};
             ground.PARA.comment.albedo_age_factor = {'if empty, compute automatically as function of mean pressure, i.e. 0.2./60 day^-1 for altitude 870hPa'};
            
             ground.PARA.default_value.dt_max = {3600};
             ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
             
             ground.PARA.default_value.dE_max = {50000};
             ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end

    end
    
end
