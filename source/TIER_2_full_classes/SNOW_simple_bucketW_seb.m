%========================================================================
% CryoGrid GROUND class SNOW_simple_bucketW_seb
% largely corresponds to snow scheme in the old CryoGrid 3, but without short-wave radiation penetration
% constant initial density, sublimation, water flow, refreezing, variable albedo
% S. Westermann, October 2020
%========================================================================

classdef SNOW_simple_bucketW_seb < SEB & HEAT_CONDUCTION & WATER_FLUXES & WATER_FLUXES_LATERAL & SNOW & REGRID 

    properties
        PARENT
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function snow = provide_PARA(snow)
            
            snow.PARA.max_albedo = []; %maximum surface albedo (fresh snow) [-]
            snow.PARA.min_albedo = []; %minimum surface albedo (old snow) [-]
            snow.PARA.tau_1 = []; % time constants for transient albedo
            snow.PARA.tau_a = [];
            snow.PARA.tau_f = [];
            snow.PARA.SW_extinction = [];
            
            snow.PARA.epsilon = []; %surface emissivity [-]
            snow.PARA.z0 = []; %roughness length [m]
            
            snow.PARA.density = []; %(initial) snow density [kg/m3]
            snow.PARA.field_capacity = []; %snow field capacity in fraction of available pore space [-] NOTE: the definition is different for GROUND_XX classes
            snow.PARA.hydraulicConductivity = []; %hydraulic conductivity of snow [m/sec]
            snow.PARA.swe_per_cell = []; %target SWE per grid cell [m]
            
            snow.PARA.dt_max = []; %maximum possible timestep [sec]
            snow.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
            
            snow.PARA.snow_property_function = [];
            snow.PARA.conductivity_function = [];
        end
        
        function snow = provide_STATVAR(snow)
            
            snow.STATVAR.upperPos = []; % upper surface elevation [m]
            snow.STATVAR.lowerPos = []; % lower surface elevation [m]
            snow.STATVAR.layerThick = [];  % thickness of grid cells [m]
            
            snow.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            snow.STATVAR.mineral = []; % total volume of minerals [m3]
            snow.STATVAR.organic = []; % total volume of organics [m3]
            snow.STATVAR.energy = [];  % total internal energy[J]
            
            snow.STATVAR.T = [];      % temperature [degree C]
            snow.STATVAR.water = [];  % total volume of water [m3]
            snow.STATVAR.ice = [];    %total volume of ice [m3]
            snow.STATVAR.air = [];    % total volume of air [m3] - NOT USED
            snow.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            snow.STATVAR.hydraulicConductivity = []; %hydraulic conductivity of snow [m/sec]
            snow.STATVAR.albedo = []; %snow albedo [-]
            
            snow.STATVAR.Lstar = []; %Obukhov length [m]
            snow.STATVAR.Qh = [];    %sensible heat flux [W/m2]
            snow.STATVAR.Qe = [];    % latent heat flux [W/m2]
        end
    
        function snow = provide_CONST(snow)
            
            snow.CONST.L_f = []; % volumetric latent heat of fusion, freezing
            snow.CONST.c_w = []; % volumetric heat capacity water
            snow.CONST.c_i = []; % volumetric heat capacity ice
            snow.CONST.c_o = []; % volumetric heat capacity organic
            snow.CONST.c_m = []; % volumetric heat capacity mineral
            
            snow.CONST.k_a = [];   % thermal conductivity air
            snow.CONST.k_w = [];   % thermal conductivity water
            snow.CONST.k_i = [];   % thermal conductivity ice 
            snow.CONST.k_o = [];   % thermal conductivity organic 
            snow.CONST.k_m = [];   % thermal conductivity mineral 
            
            snow.CONST.sigma = []; %Stefan-Boltzmann constant
            snow.CONST.kappa = []; % von Karman constant
            snow.CONST.L_s = [];   %latent heat of sublimation, evaporation handled in a dedicated function
            
            snow.CONST.cp = []; % specific heat capacity at constant pressure of air
            snow.CONST.g = []; % gravitational acceleration Earth surface
            
            snow.CONST.rho_w = []; % water density
            snow.CONST.rho_i = []; %ice density
            snow.CONST.Tmfw = [];
        end
        
        function snow = finalize_init(snow, tile)
            snow.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            snow.PARA.airT_height = tile.FORCING.PARA.airT_height;
            
            snow = initialize_zero_snow_BASE(snow); 
            snow.TEMP.d_energy = snow.STATVAR.energy .*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
            snow.TEMP.newSnow = 0;
            snow.STATVAR.albedo = snow.PARA.max_albedo;
            snow.STATVAR.excessWater = 0;
        end
        
        %---time integration------
        %separate functions for CHILD pphase of snow cover
        
        function snow = get_boundary_condition_u(snow, tile) 
            forcing = tile.FORCING;
            snow = surface_energy_balance(snow, forcing);
            snow = get_boundary_condition_SNOW_u(snow, forcing);
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            snow = get_sublimation(snow, forcing);
        end
        
        function snow = get_boundary_condition_u_CHILD(snow, tile)
            forcing = tile.FORCING;
            snow = surface_energy_balance(snow, forcing);
            snow = get_boundary_condition_allSNOW_rain_u(snow, forcing); %add full snow, but rain only for snow-covered part
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            snow = get_sublimation(snow, forcing);
        end
        
        function snow = get_boundary_condition_u_create_CHILD(snow, tile)
            forcing = tile.FORCING;
            snow = get_boundary_condition_allSNOW_u(snow, forcing); %add full snow, no rain
            
            snow.TEMP.d_energy = 0;
            snow.TEMP.d_water = 0;
            snow.TEMP.d_water_energy = 0;
            
            snow.TEMP.F_ub = 0;
            snow.TEMP.F_lb = 0;
            snow.TEMP.F_ub_water = 0;
            snow.TEMP.F_lb_water = 0;
            snow.TEMP.F_ub_water_energy = 0;
            snow.TEMP.F_lb_water_energy = 0;
            snow.STATVAR.sublimation = 0;
            snow.TEMP.sublimation_energy = 0;
            snow.TEMP.rain_energy = 0;
            snow.TEMP.rainfall = 0;
            snow.STATVAR.area = 0;
            snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ (snow.PARA.density ./1000); %[m] constant layerThick during CHILD phase
            snow.STATVAR.ice = 0;
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos;
        end
        
        function [snow, S_up] = penetrate_SW(snow, S_down)  %mandatory function when used with class that features SW penetration
            [snow, S_up] = penetrate_SW_transmission_bulk(snow, S_down);
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
            else
                %handled in advance_prognostic()
            end
        end
        
        function snow = get_derivatives_prognostic_CHILD(snow, tile)
            %do nothing
        end
        
        function timestep = get_timestep(snow, tile) 
            timestep = get_timestep_SNOW(snow);
            %timestep1 = get_timestep_heat_coduction(snow);
            %timestep2 = get_timestep_SNOW_mass_balance(snow);
            timestep3 = get_timestep_water_SNOW(snow);

            %timestep = min(timestep1, timestep2);
            timestep = min(timestep, timestep3);
        end
        
        function timestep = get_timestep_CHILD(snow, tile)  
            timestep = get_timestep_SNOW_CHILD(snow);
        end
        
        function snow = advance_prognostic(snow, tile)
            timestep = tile.timestep;
            %energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy + snow.TEMP.d_water_energy);
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* (snow.TEMP.snow_energy + snow.TEMP.sublimation_energy);  %rainfall energy already added
            %mass
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* snow.TEMP.d_water;
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + timestep .* (snow.TEMP.snowfall + snow.STATVAR.sublimation);
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* snow.STATVAR.sublimation ./snow.STATVAR.area(1,1) ./ (snow.STATVAR.ice(1) ./ snow.STATVAR.layerThick(1) ./ snow.STATVAR.area(1));
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* snow.TEMP.snowfall ./snow.STATVAR.area(1,1) ./ (snow.PARA.density ./1000);
            %store "old" density
            snow.STATVAR.target_density = snow.STATVAR.ice ./ snow.STATVAR.layerThick ./ snow.STATVAR.area;
            snow.STATVAR.target_density(1) = (snow.STATVAR.ice(1) + timestep .* (snow.TEMP.snowfall + snow.STATVAR.sublimation)) ./ snow.STATVAR.layerThick(1) ./ snow.STATVAR.area(1,1);
        
            snow = calculate_albedo_simple(snow, timestep);
        end
        
        function snow = advance_prognostic_CHILD(snow, tile)
            timestep = tile.timestep;
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy + snow.TEMP.snow_energy + snow.TEMP.d_water_energy + snow.TEMP.sublimation_energy);
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* (snow.TEMP.snowfall + snow.TEMP.d_water + snow.STATVAR.sublimation);
            snow.STATVAR.area = snow.STATVAR.area + timestep .* snow.STATVAR.sublimation ./snow.STATVAR.layerThick ./ max(50/1000, snow.STATVAR.ice ./ snow.STATVAR.layerThick ./ snow.STATVAR.area);
            snow.STATVAR.area = snow.STATVAR.area + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./1000) ./  snow.STATVAR.layerThick ; %[m2]
            snow.STATVAR.target_density = min(1,(snow.STATVAR.ice + timestep .* (snow.TEMP.snowfall + snow.STATVAR.sublimation)) ./ snow.STATVAR.layerThick ./ snow.STATVAR.area);
            
            snow = calculate_albedo_simple(snow, timestep);
        end
        
        function snow = compute_diagnostic_first_cell(snow, tile)
            forcing = tile.FORCING;
            snow = L_star(snow, forcing);
        end
       
        function snow = compute_diagnostic(snow, tile)
            snow = get_T_water_freeW(snow);
            snow = subtract_water2(snow);
            
            [snow, regridded_yesNo] = regrid_snow(snow, {'waterIce'; 'energy'; 'layerThick'; 'mineral'; 'organic'}, {'area'; 'target_density'}, 'ice');
            if regridded_yesNo
                snow = get_T_water_freeW(snow);
            end
            
            snow = conductivity(snow);
            snow = calculate_hydraulicConductivity_SNOW(snow);
            
            snow.STATVAR.upperPos = snow.NEXT.STATVAR.upperPos + sum(snow.STATVAR.layerThick);
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
        end
        
        function snow = compute_diagnostic_CHILD(snow, tile)

            snow = get_T_water_freeW(snow);
            snow = subtract_water_CHILD(snow);

            snow = conductivity(snow);
            snow = calculate_hydraulicConductivity_SNOW(snow);
            
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos + (snow.STATVAR.layerThick .* snow.STATVAR.area ./ snow.PARENT.STATVAR.area(1,1));
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
            
        end
        
        function snow = check_trigger(snow, tile)
            snow = make_SNOW_CHILD(snow); 
        end
        
        %-----non-mandatory functions-------
        function snow = surface_energy_balance(snow, forcing)
            Lin = forcing.TEMP.Lin .* snow.STATVAR.area(1);
            Sin = forcing.TEMP.Sin .* snow.STATVAR.area(1); 
            [snow, Lout] = penetrate_LW(snow, Lin);
            [snow, Sout] = penetrate_SW(snow, Sin);
            snow.STATVAR.Qh = Q_h(snow, forcing);
            snow.STATVAR.Qe = Q_eq_potET(snow, forcing);
            
            snow.TEMP.F_ub = (- snow.STATVAR.Qh - snow.STATVAR.Qe) .* snow.STATVAR.area(1);
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.F_ub;
        end
        
        function snow = conductivity(snow)
            conductivity_function = str2func(snow.PARA.conductivity_function);
            snow = conductivity_function(snow);
        end
        
        
        %-----LATERAL-------------------
        
        %-----LAT_REMOVE_SURFACE_WATER-----
        function ground = lateral_push_remove_surfaceWater(ground, lateral)
            ground = lateral_push_remove_surfaceWater_simple(ground, lateral);
        end
                
        %----LAT_SEEPAGE_FACE----------
        function snow = lateral_push_remove_water_seepage(snow, lateral)
            snow = lateral_push_remove_water_seepage_snow(snow, lateral);
        end
        
        %----LAT_WATER_RESERVOIR------------
        function snow = lateral_push_water_reservoir(snow, lateral)
            snow = lateral_push_water_reservoir_snow(snow, lateral);
        end
        
        %-----------------------------
        
        %----inherited Tier 1 functions ------------
        
        function snow = get_derivative_energy(snow)
           snow = get_derivative_energy@HEAT_CONDUCTION(snow); 
        end
        
        function snow = conductivity_snow_Yen(snow)
            snow = conductivity_snow_Yen@HEAT_CONDUCTION(snow);
        end
        
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
        
        function [snow, S_up] = penetrate_SW_no_transmission(snow, S_down)
             [snow, S_up] = penetrate_SW_no_transmission@SEB(snow, S_down);
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
             
             
             ground.PARA.default_value.max_albedo = {0.85};
             ground.PARA.comment.max_albedo = {'maximum surface albedo (fresh snow) [-]'};
             
             ground.PARA.default_value.min_albedo = {0.5};
             ground.PARA.comment.min_albedo = {'minimum surface albedo (old snow) [-]'};
             
             ground.PARA.default_value.tau_1 = {86400};
             ground.PARA.comment.tau_1 = {'time constants for transient albedo'};
             
             ground.PARA.default_value.tau_a = {0.008};
             ground.PARA.default_value.tau_f = {0.24};

             ground.PARA.default_value.epsilon = {0.99};
             ground.PARA.comment.epsilon = {'surface emissivity [-]'};
             
             ground.PARA.default_value.z0 = {0.01};
             ground.PARA.comment.z0 = {'roughness length [m]'};
             
             ground.PARA.default_value.density = {350};
             ground.PARA.comment.density = {'(initial) snow density [kg/m3]'};
             
             ground.PARA.default_value.field_capacity = {0.05};
             ground.PARA.comment.field_capacity = {'snow field capacity in fraction of available pore space [-] NOTE: the definition is different for GROUND_XX classes'};
            
             ground.PARA.default_value.hydraulicConductivity = {1e-4};
             ground.PARA.comment.hydraulicConductivity = {'hydraulic conductivity of snow [m/sec]'};
             
             ground.PARA.default_value.swe_per_cell = {0.02};
             ground.PARA.comment.swe_per_cell = {'target SWE per grid cell [m]'};
            
             ground.PARA.default_value.dt_max = {3600};
             ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
             
             ground.PARA.default_value.dE_max = {50000};
             ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end
        
    end
    
end
