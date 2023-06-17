%========================================================================
% CryoGrid GROUND class SNOW_simple_seb
% simple snow scheme with little value for science
% constant initial density, no sublimation, no water flow, refreezing of melt/rain water within a cell, constant albedo
% S. Westermann, October 2020
%========================================================================

classdef SNOW_simple_seb < SEB & HEAT_CONDUCTION & SNOW & WATER_FLUXES_LATERAL & REGRID 

    properties
        PARENT
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

        
        function snow = provide_PARA(snow)
            
            snow.PARA.albedo = []; %surface albedo [-]
            snow.PARA.epsilon = [];  % surface emissivity [-]
            snow.PARA.z0 = []; %roughness length [m]
            snow.PARA.density = []; %(initial) snow density [kg/m3]
            snow.PARA.swe_per_cell = []; %target SWE per grid cell [m]
            snow.PARA.dt_max = []; %maximum possible timestep [sec]
            snow.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
        end
        
        function snow = provide_STATVAR(snow)
            
            snow.STATVAR.upperPos = []; % upper surface elevation [m]
            snow.STATVAR.lowerPos = []; % lower surface elevation [m]
            snow.STATVAR.layerThick = []; % thickness of grid cells [m]
            
            snow.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            snow.STATVAR.mineral = []; % total volume of minerals [m3]
            snow.STATVAR.organic = []; % total volume of organics [m3]
            snow.STATVAR.energy = []; % total internal energy[J]
            
            snow.STATVAR.T = []; % temperature [degree C]
            snow.STATVAR.water = [];  % total volume of water [m3]
            snow.STATVAR.ice = []; %total volume of ice [m3]
            snow.STATVAR.air = [];   % total volume of air [m3] - NOT USED
            snow.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            
            snow.STATVAR.Lstar = []; %Obukhov length [m]
            snow.STATVAR.Qh = []; %sensible heat flux [W/m2]
            snow.STATVAR.Qe = []; % latent heat flux [W/m2]
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
            
            snow.CONST.cp = [];   % specific heat capacity at constant pressure of air
            snow.CONST.g = [];    % gravitational acceleration Earth surface
            
            snow.CONST.rho_w = [];  % water density
            snow.CONST.rho_i = [];  %ice density
        end
        
        
        function snow = finalize_init(snow, tile)
            snow.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            snow.PARA.airT_height = tile.FORCING.PARA.airT_height;
            
            snow = initialize_zero_snow_BASE(snow); 
            snow.STATVAR.excessWater = 0;
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
        end
        
        
        %---time integration------
        %separate functions for CHILD pphase of snow cover
        
        function snow = get_boundary_condition_u(snow, tile) 
            forcing = tile.FORCING;
            snow = surface_energy_balance(snow, forcing);
            snow = get_boundary_condition_SNOW_u(snow, forcing);
        end
        
        function snow = get_boundary_condition_u_CHILD(snow, tile)
            forcing = tile.FORCING;
            snow = surface_energy_balance(snow, forcing);
            snow = get_boundary_condition_allSNOW_rain_u(snow, forcing); %add full snow, but rain only for snow-covered part
        end
        
        function snow = get_boundary_condition_u_create_CHILD(snow, tile)
            forcing = tile.FORCING;
            snow = get_boundary_condition_allSNOW_u(snow, forcing); %add full snow, no rain
            snow.TEMP.d_energy = 0;
            snow.TEMP.F_ub = 0;
            snow.TEMP.F_lb = 0;
            snow.TEMP.rain_energy = 0;
            snow.TEMP.rainfall = 0;
            snow.STATVAR.area = 0;
            snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ (snow.PARA.density ./1000); %[m] constant layerThick during CHILD phase
            snow.STATVAR.ice = 0;
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos;
        end
        
        function snow = get_boundary_condition_l(snow, tile)
            forcing = tile.FORCING;
            snow.TEMP.F_lb = forcing.PARA.heatFlux_lb .* snow.STATVAR.area(end);
            snow.TEMP.d_energy(end) = snow.TEMP.d_energy(end) + snow.TEMP.F_lb;
        end
        
        function snow = get_derivatives_prognostic(snow, tile)
            if size(snow.STATVAR.layerThick,1) > 1
                snow = get_derivative_energy(snow);
            else
                %handled in advance_prognostic()
            end
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.rain_energy;  %add rain energy here, since it can melt snow and does not change layerThick - must be taken into account for timestep calculation
        end
        
        function snow = get_derivatives_prognostic_CHILD(snow, tile)
            snow.TEMP.d_energy = snow.TEMP.d_energy + snow.TEMP.rain_energy;
        end
        
        function timestep = get_timestep(snow, tile) 
%             timestep1 = get_timestep_heat_coduction(snow);
%             timestep2 = get_timestep_SNOW_mass_balance(snow);
%             timestep = min(timestep1, timestep2);
            timestep = get_timestep_SNOW(snow);
        end
        
        function timestep = get_timestep_CHILD(snow, tile)  
            timestep = get_timestep_SNOW_CHILD(snow);
        end
        
        function snow = advance_prognostic(snow, tile)
            timestep = tile.timestep;
            %energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* snow.TEMP.d_energy;
            %mass
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* snow.TEMP.snow_energy;  %rainfall energy already added
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + timestep .* (snow.TEMP.snowfall + snow.TEMP.rainfall);
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* snow.TEMP.snowfall ./snow.STATVAR.area(1,1) ./ (snow.PARA.density ./1000);
            %store "old" density
            snow.STATVAR.target_density = snow.STATVAR.ice ./ snow.STATVAR.layerThick ./ snow.STATVAR.area;
            snow.STATVAR.target_density(1) = (snow.STATVAR.ice(1) + timestep .* snow.TEMP.snowfall) ./ snow.STATVAR.layerThick(1) ./ snow.STATVAR.area(1,1);
        end
        
        function snow = advance_prognostic_CHILD(snow, tile)
            timestep = tile.timestep;
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy + snow.TEMP.snow_energy);
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* (snow.TEMP.snowfall + snow.TEMP.rainfall);
            snow.STATVAR.area = snow.STATVAR.area + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./1000) ./  snow.STATVAR.layerThick ; 
            snow.STATVAR.target_density = min(1,(snow.STATVAR.ice + timestep .* snow.TEMP.snowfall) ./ snow.STATVAR.layerThick ./ snow.STATVAR.area);
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
            snow.STATVAR.upperPos = snow.NEXT.STATVAR.upperPos + sum(snow.STATVAR.layerThick);
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
        end
        
        function snow = compute_diagnostic_CHILD(snow, tile)
            snow = get_T_water_freeW(snow);
            snow = subtract_water_CHILD(snow);

            snow = conductivity(snow);
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos + (snow.STATVAR.layerThick .* snow.STATVAR.area ./ snow.PARENT.STATVAR.area(1,1));
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            
        end
        
        function snow = check_trigger(snow, tile)
            snow = make_SNOW_CHILD(snow); 
        end
        
        
        %-----non-mandatory functions-------
        function snow = surface_energy_balance(snow, forcing)
            snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
            snow.STATVAR.Sout = snow.PARA.albedo .*  forcing.TEMP.Sin;
            snow.STATVAR.Qh = Q_h(snow, forcing);
            snow.STATVAR.Qe = Q_eq_potET(snow, forcing);
            
            snow.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe) .* snow.STATVAR.area(1);
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.F_ub;
        end
        
        function snow = conductivity(snow)
            snow = conductivity_snow_Yen(snow);
        end
        
        function yesNo = is_ground_surface(snow)
            yesNo = 0;
        end
                
        %-----LATERAL-------------------
        
        function gse = get_groundSurfaceElevation(ground)
           gse = []; 
        end
        
        
        %-----LAT_REMOVE_SURFACE_WATER-----
        function snow = lateral_push_remove_surfaceWater(snow, lateral)
            snow = lateral_push_remove_surfaceWater_simple(snow, lateral);
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
        
        
        %-------------param file generation-----
         function ground = param_file_info(ground)
             ground = param_file_info@BASE(ground);
             
             ground.PARA.class_category = 'SNOW';
             
             ground.PARA.STATVAR = {''};
             
             ground.PARA.default_value.albedo = {0.8};
             ground.PARA.comment.albedo = {'surface albedo [-]'};
             
             ground.PARA.default_value.epsilon = {0.99};
             ground.PARA.comment.epsilon = {'surface emissivity [-]'};
             
             ground.PARA.default_value.z0 = {0.01};
             ground.PARA.comment.z0 = {'roughness length [m]'};
             
             ground.PARA.default_value.density = {350};
             ground.PARA.comment.density = {'(initial) snow density [kg/m3]'};
             
             ground.PARA.default_value.swe_per_cell = {0.02};
             ground.PARA.comment.swe_per_cell = {'target SWE per grid cell [m]'};
            
             ground.PARA.default_value.dt_max = {3600};
             ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
             
             ground.PARA.default_value.dE_max = {50000};
             ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end
    end
    
end
