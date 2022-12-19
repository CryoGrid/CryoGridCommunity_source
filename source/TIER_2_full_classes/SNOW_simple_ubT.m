%========================================================================
% CryoGrid SNOW class SNOW_simple_ubT
% simple snow scheme forced by snow depth and surface temperature
% constant density, no water flow, melt, or refreezing of melt/rain water 
% S. Westermann, October 2021
%========================================================================

classdef SNOW_simple_ubT <  HEAT_CONDUCTION & SNOW & WATER_FLUXES_LATERAL & REGRID 
    properties
        PARENT
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

        
        function snow = provide_PARA(snow)

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
            
            
            snow.CONST.rho_w = [];  % water density
            snow.CONST.rho_i = [];  %ice density
        end
        
        
        function snow = finalize_init(snow, tile)
            snow.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            
            snow = initialize_zero_snow_BASE(snow); 
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
        end
        
        
        %---time integration------
        %separate functions for CHILD pphase of snow cover
        
        function snow = get_boundary_condition_u(snow, tile) 
            forcing = tile.FORCING;
            
            snow.TEMP.F_ub = (tile.FORCING.TEMP.T_ub - snow.STATVAR.T(1,1)) .* snow.STATVAR.thermCond(1,1) ./snow.STATVAR.layerThick(1,1) .* snow.STATVAR.area(1,1);
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.F_ub;
            
            %snow = get_boundary_condition_SNOW_u(snow, forcing);
        end
        
        function snow = get_boundary_condition_u_CHILD(snow, tile)
            forcing = tile.FORCING;
            snow.TEMP.F_ub = (tile.FORCING.TEMP.T_ub - snow.STATVAR.T(1,1)) .* snow.STATVAR.thermCond(1,1) ./snow.STATVAR.layerThick(1,1) .* snow.STATVAR.area(1,1);
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.F_ub;
            
            %snow = get_boundary_condition_allSNOW_rain_u(snow, forcing); %add full snow, but rain only for snow-covered part
        end
        
        function snow = get_boundary_condition_u_create_CHILD(snow, tile)
            forcing = tile.FORCING;
            %snow = get_boundary_condition_allSNOW_u(snow, forcing); %add full snow, no rain
            snow.TEMP.d_energy = 0;
            snow.TEMP.F_ub = 0;
            snow.TEMP.F_lb = 0;
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
            end
        end
        
        function snow = get_derivatives_prognostic_CHILD(snow, tile)

        end
        
        function timestep = get_timestep(snow, tile) 
            timestep = get_timestep_SNOW(snow);
        end
        
        function timestep = get_timestep_CHILD(snow, tile)  
            timestep = get_timestep_SNOW_CHILD(snow);
        end
        
        function snow = advance_prognostic(snow, tile)
            timestep = tile.timestep;
            %energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* snow.TEMP.d_energy;

        end
        
        function snow = advance_prognostic_CHILD(snow, tile)
            timestep = tile.timestep;
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* snow.TEMP.d_energy;

        end
        
        function snow = compute_diagnostic_first_cell(snow, tile)

        end
       
        function snow = compute_diagnostic(snow, tile)
            
            snow_depth = sum(snow.STATVAR.layerThick,1);
            snow_depth_target = tile.FORCING.TEMP.snow_depth;
            if snow_depth < snow_depth_target
                difference = snow_depth_target - snow_depth; 
                snow.STATVAR.layerThick(1,1) = snow.STATVAR.layerThick(1,1) + difference;
                snow.STATVAR.energy(1,1) = snow.STATVAR.energy(1,1) + difference .* snow.PARA.density ./1000 .* (tile.FORCING.TEMP.T_ub  .* snow.CONST.c_i - snow.CONST.L_f) .* snow.STATVAR.area(1,1);
                snow.STATVAR.waterIce(1,1) = snow.STATVAR.waterIce(1,1) + difference .* snow.PARA.density ./1000 .* snow.STATVAR.area(1,1);
                
            end
            while snow_depth > snow_depth_target + 1e-6
                if  sum(snow.STATVAR.layerThick(2:end,1)) < snow_depth_target
                    remaining_fraction = (snow_depth_target - sum(snow.STATVAR.layerThick(2:end,1),1)) ./ snow.STATVAR.layerThick(1,1);
                    snow.STATVAR.layerThick(1,1) = remaining_fraction .* snow.STATVAR.layerThick(1,1);
                    snow.STATVAR.energy(1,1) = remaining_fraction .* snow.STATVAR.energy(1,1);
                    snow.STATVAR.waterIce(1,1) = remaining_fraction .* snow.STATVAR.waterIce(1,1);
                else
                    snow.STATVAR.layerThick = snow.STATVAR.layerThick(2:end,1);
                    snow.STATVAR.energy = snow.STATVAR.energy(2:end,1);
                    snow.STATVAR.waterIce = snow.STATVAR.waterIce(2:end,1);
                end
                snow_depth = sum(snow.STATVAR.layerThick,1);
            end

            [snow, ~] = regrid_snow(snow, {'waterIce'; 'energy'; 'layerThick'}, {'area'}, 'waterIce');

            snow = get_T_water_freeW(snow);

            snow = conductivity(snow);
            snow.STATVAR.upperPos = snow.NEXT.STATVAR.upperPos + sum(snow.STATVAR.layerThick);
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
        end
        
        function snow = compute_diagnostic_CHILD(snow, tile)
            
            snow_depth_target = tile.FORCING.TEMP.snow_depth;
            volume_difference = snow_depth_target .* snow.PARENT.STATVAR.area(1,1) - snow.STATVAR.layerThick .* snow.STATVAR.area;
            
            snow.STATVAR.area = snow_depth_target .* snow.PARENT.STATVAR.area(1,1) ./ snow.STATVAR.layerThick;
            snow.STATVAR.waterIce = snow_depth_target .* snow.PARA.density ./1000;

            T = double(volume_difference < 0) .* snow.STATVAR.T + double(volume_difference > 0) .* tile.FORCING.TEMP.T_ub; 
            snow.STATVAR.energy = snow.STATVAR.energy + volume_difference .* snow.PARA.density ./1000 .* (T  .* snow.CONST.c_i - snow.CONST.L_f);
            
            snow = get_T_water_freeW(snow);

            snow = conductivity(snow);
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos + (snow.STATVAR.layerThick .* snow.STATVAR.area ./ snow.PARENT.STATVAR.area(1,1));
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            
        end
        
        function snow = check_trigger(snow, tile)
            snow = make_SNOW_CHILD_ubT(snow); 
        end
        
        
        %-----non-mandatory functions-------

        
        function snow = conductivity(snow)
            snow = conductivity_snow_Yen(snow);
        end
        
        function yesNo = is_ground_surface(snow)
            yesNo = 0;
        end
                
        %-----LATERAL-------------------
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE(ground);
            
            ground.PARA.class_category = 'SNOW';
            
            ground.PARA.STATVAR = {''};
            
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
