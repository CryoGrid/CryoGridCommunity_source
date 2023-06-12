%========================================================================
% CryoGrid GROUND class SNOW_simple_ubtf_mf
% Simple snow scheme with upper boundary temperature forcing (Dirichlet
% boundary condition) and snow melt based on melt factors.
% See TIER1 class SNOW_MELTFACTOR for description of snow melt scheme.
%
% T. Ingeman-Nielsen, S. Westermann, December 2021
%========================================================================


classdef SNOW_simple_ubtf_mf < HEAT_CONDUCTION & UB_TEMPERATURE_FORCING & SNOW & SNOW_MELTFACTOR & REGRID


    properties
        PARENT
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function snow = provide_PARA(snow)
            
            snow.PARA.density = []; % (initial) snow density [kg/m3]
            snow.PARA.swe_per_cell = []; % target SWE per grid cell [m]
            snow.PARA.dt_max = []; % maximum possible timestep [sec]
            snow.PARA.dE_max = []; % maximum possible energy change per timestep [J/m3]
            snow.PARA.melt_threshold = []; % threshold air temperature for snow melt to occur [degC]
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
            snow.STATVAR.water = []; % total volume of water [m3]
            snow.STATVAR.ice = []; % total volume of ice [m3]
            snow.STATVAR.air = []; % total volume of air [m3] - NOT USED
            snow.STATVAR.thermCond = []; % thermal conductivity [W/mK]
            
        end
    
        
        function snow = provide_CONST(snow)
            
            snow.CONST.L_f = []; % volumetric latent heat of fusion, freezing
            snow.CONST.c_w = []; % volumetric heat capacity water
            snow.CONST.c_i = []; % volumetric heat capacity ice
            snow.CONST.c_o = []; % volumetric heat capacity organic
            snow.CONST.c_m = []; % volumetric heat capacity mineral
            
            snow.CONST.k_a = []; % thermal conductivity air
            snow.CONST.k_w = []; % thermal conductivity water
            snow.CONST.k_i = []; % thermal conductivity ice 
            snow.CONST.k_o = []; % thermal conductivity organic 
            snow.CONST.k_m = []; % thermal conductivity mineral 
            
            snow.CONST.rho_w = []; % water density
            snow.CONST.rho_i = []; % ice density
            snow.CONST.day_sec = []; % number of seconds in a day
        end
        
        
        function snow = finalize_init(snow, tile) 
            snow.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
%             snow.PARA.airT_height = tile.FORCING.PARA.airT_height;
            snow.PARA.latitude = tile.PARA.latitude;
            
            snow = initialize_zero_snow_BASE(snow); 
            snow.STATVAR.excessWater = 0;
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
        end
        
        
        %---time integration------
        %separate functions for CHILD pphase of snow cover
        
        function snow = get_boundary_condition_u(snow, tile) 
            % Calculates upper boundary conditions of the snow class.
            % This method is called when the SNOW class is part of the
            % regular stratigraphy (not in the CHILD phase).
            forcing = tile.FORCING;
            snow = get_ub_temperature(snow, forcing); % inherited from UB_TEMPERATURE_FORCING
            snow = get_boundary_condition_SNOW_u(snow, forcing); % inherited from SNOW
            
            % calculate snowmelt from air temperature forcing
            snow = get_boundary_condition_SNOW_meltFactor(snow, tile); % inherited from SNOW
        end
        
        
        function snow = get_boundary_condition_u_CHILD(snow, tile)
            % Calculates upper boundary conditions of the snow class.
            % This method is called when the SNOW class is in the CHILD 
            % phase.
            forcing = tile.FORCING;
            snow = get_ub_temperature(snow, forcing); % inherited from UB_TEMPERATURE_FORCING
            snow = get_boundary_condition_allSNOW_rain_u(snow, forcing); %add full snow, but rain only for snow-covered part
            
            % calculate snowmelt from air temperature forcing
            snow = get_boundary_condition_SNOW_meltFactor(snow, tile); % inherited from SNOW
        end
        
        
        function snow = get_boundary_condition_u_create_CHILD(snow, tile)
            % Calculates upper boundary conditions of the snow class.
            % This method is called right after the SNOW class is first
            % instantiated in the CHILD phase (but not on subsequent time
            % steps)
            
            forcing = tile.FORCING;
            snow = get_boundary_condition_allSNOW_u(snow, forcing); %add full snow, no rain
            snow.TEMP.d_energy = 0;
            snow.TEMP.F_ub = 0;
            snow.TEMP.F_lb = 0;
            snow.TEMP.rain_energy = 0;
            snow.TEMP.rainfall = 0;
            snow.TEMP.melt_energy = 0; % Child only instantiated when air temperature is < 0, thus no chance of melt.
            snow.TEMP.melt_factor = 0;
            snow.STATVAR.area = 0;
            snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ (snow.PARA.density ./1000); %[m] constant layerThick during CHILD phase
            snow.STATVAR.ice = 0;
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos;
        end
        
        
        function snow = get_boundary_condition_l(snow, tile)
            % Calculates lower boundary forcing. Only called for the
            % lowermost class in the stratigraphy (so only for a snow class
            % in specialized runs consisting only of a snow pack...)
            forcing = tile.FORCING;
            snow.TEMP.F_lb = forcing.PARA.heatFlux_lb .* snow.STATVAR.area(end);
            snow.TEMP.d_energy(end) = snow.TEMP.d_energy(end) + snow.TEMP.F_lb;
        end
        
        
        function snow = get_derivatives_prognostic(snow, tile)
            % calculates derivatives (fluxes) and changes in energy
            % content.
            if size(snow.STATVAR.layerThick,1) > 1
                snow = get_derivative_energy(snow);
            else
                %handled in advance_prognostic()
            end
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.rain_energy;  %add rain energy here, since it can melt snow and does not change layerThick - must be taken into account for timestep calculation
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.melt_energy;  %add melt energy here, since it can melt snow - must be taken into account for timestep calculation
        end
        
        
        function snow = get_derivatives_prognostic_CHILD(snow, tile)
            % In the child phase only handle changes in energy
            snow.TEMP.d_energy = snow.TEMP.d_energy + snow.TEMP.rain_energy;
            snow.TEMP.d_energy = snow.TEMP.d_energy + snow.TEMP.melt_energy;
        end
        
        
        function timestep = get_timestep(snow, tile) 
            % calculate the maximum timestep to be used for the snow pack
            timestep = get_timestep_SNOW(snow);
        end
        
        
        function timestep = get_timestep_CHILD(snow, tile)  
            % calculate the maximum timestep to be used for the snow pack
            % in the child phase
            timestep = get_timestep_SNOW_CHILD(snow);
        end
        
        
        function snow = advance_prognostic(snow, tile)
            % Update prognostic variables for next time step
            
            timestep = tile.timestep;
            
            % energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* snow.TEMP.d_energy;
            
            % mass
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* snow.TEMP.snow_energy;  % rainfall energy already added in get_derivatives_prognostic
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + timestep .* (snow.TEMP.snowfall + snow.TEMP.rainfall);
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* snow.TEMP.snowfall ./snow.STATVAR.area(1,1) ./ (snow.PARA.density ./ 1000);
            % snow.TEMP.snowfall and rainfall are the rates rate in m/s snow water equivalents (units converted in get_boundary_condition_u)
            % snow.PARA.density is the snow density in kg/m3
            % 1000 kg/m3 is the density of water used in conversion from swe to snow thickness
            
            % Store snow density, new snow fall is included for first cell.
            % target_density is used to ensure that e.g. snow melt does not
            % change the snow density, just reduce the layer height.
            snow.STATVAR.target_density = snow.STATVAR.ice ./ (snow.STATVAR.layerThick .* snow.STATVAR.area);   % target_density = total volume of ice / total cell volume * water-density
            % snow.STATVAR.ice is the snowpack ice content in m3 water equivalents
            % snow.STATVAR.target_density is a relative snowpack density (normalized by the density of water)
            % to get the actual dry snow pack density, multiply by 1000 kg/m3 (density of water)
            % snow.STATVAR.target_density considers only the dry snow pack,
            % any liquid water content is not taken into account, and is
            % assumed contained in the porespace of the dry snow pack.
                        
            % in first cell we include also the new snowfall
            snow.STATVAR.target_density(1) = (snow.STATVAR.ice(1) + timestep .* snow.TEMP.snowfall) ./ (snow.STATVAR.layerThick(1) .* snow.STATVAR.area(1,1)); 
            % snow.TEMP.snowfall is the snowfall rate in m/s snow water equivalents (units converted in get_boundary_condition_u)
            
        end
        
        
        function snow = advance_prognostic_CHILD(snow, tile)
            % Update prognostic variables for next time step for child phase
            
            timestep = tile.timestep;
            
            %energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy + snow.TEMP.snow_energy);
            
            %mass
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* (snow.TEMP.snowfall + snow.TEMP.rainfall);
            snow.STATVAR.area = snow.STATVAR.area + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./ 1000) ./  snow.STATVAR.layerThick ; 
            % snow.TEMP.snowfall and rainfall are the rates rate in m/s snow water equivalents (units converted in get_boundary_condition_u)
            % snow.PARA.density is the snow density in kg/m3
            % 1000 kg/m3 is the density of water
            % In the child phase, the snow layer thickness is fixed, while the area coverage of the snow changes.
            
            % Store snow density, new snow fall is included.
            % target_density is used to ensure that e.g. snow melt does not
            % change the snow density, just reduce the layer height.
            snow.STATVAR.target_density = min(1,(snow.STATVAR.ice + timestep .* snow.TEMP.snowfall) ./ snow.STATVAR.layerThick ./ snow.STATVAR.area);
            % snow.STATVAR.ice is the snowpack ice content in m3 water equivalents
            % snow.STATVAR.target_density is a relative snowpack density (normalized by the density of water)
            % to get the actual dry snow pack density, multiply by 1000 kg/m3 (density of water)
            % snow.STATVAR.target_density considers only the dry snow pack,
            % any liquid water content is not taken into account, and is
            % assumed contained in the porespace of the dry snow pack.            
        end
        
        
        function snow = compute_diagnostic_first_cell(snow, tile)
            % Computes diagnostic parameters specifically for the first
            % cell. Only called on the uppermost class in the stratigraphy.
            
            % do nothing
        end
       
        
        function snow = compute_diagnostic(snow, tile)
            % Computes diagnostic parameters for the entire grid in the class
            
            snow = get_T_water_freeW(snow); % derive temperature and water/ice contents from the energy content
            snow = subtract_water2(snow); % remove excess water from the snow pack
            
            % regrid the snowpack if necessary
            [snow, regridded_yesNo] = regrid_snow(snow, {'waterIce'; 'energy'; 'layerThick'; 'mineral'; 'organic'}, {'area'; 'target_density'}, 'ice');
            if regridded_yesNo
                snow = get_T_water_freeW(snow); % recalculate temperature and water/ice distributions after regridding
            end
            
            snow = conductivity(snow); % Calculate the conductivity of the snow
            
            % recalculate the top of the surface level of the snowpack
            snow.STATVAR.upperPos = snow.NEXT.STATVAR.upperPos + sum(snow.STATVAR.layerThick); 
            
            % reset temporary values in preparation for next iteration
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
        end
        
        
        function snow = compute_diagnostic_CHILD(snow, tile)
            %Computes diagnostic parameters for snow CHILD
            
            snow = get_T_water_freeW(snow); % derive temperature and water/ice contents from the energy content
            snow = subtract_water_CHILD(snow); % remove excess water from the snow pack
            
            snow = conductivity(snow); % Calculate the conductivity of the snow
            
            % recalculate the top of the surface level of the parent subsurface class
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos + (snow.STATVAR.layerThick .* snow.STATVAR.area ./ snow.PARENT.STATVAR.area(1,1));
            
            % reset temporary values in preparation for next iteration
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            
        end
        
        
        function snow = check_trigger(snow, tile)
            % Checks if conditions are met to reduce snowpack to child
            % phase, and if that is the case, make the conversion.
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
        
        %-----LAT_REMOVE_SURFACE_WATER-----
        function snow = lateral_push_remove_surfaceWater(snow, lateral)
            snow = lateral_push_remove_surfaceWater_simple(snow, lateral);
        end
        %-----------------------------

        
        
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE(ground);
            
            ground.PARA.class_category = 'SNOW';
            
            ground.PARA.STATVAR = {''};
            
            ground.PARA.default_value.density = {350};
            ground.PARA.comment.density = {'(initial) snow density [kg/m3]'};
            
            ground.PARA.default_value.swe_per_cell = {0.02};
            ground.PARA.comment.swe_per_cell = {'target SWE per grid cell [m]'};
            
            ground.PARA.default_value.melt_threshold = {0.5};
            ground.PARA.comment.melt_threshold = {'threshold air temperature for snow melt to occur [degC]'};
            
            ground.PARA.default_value.dt_max = {3600};
            ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
            
            ground.PARA.default_value.dE_max = {50000};
            ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end
        
    end
    
end
