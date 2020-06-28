

classdef SNOW_simple_bucketW_seb < SEB & HEAT_CONDUCTION & WATER_FLUXES & WATER_FLUXES_LATERAL & SNOW & INITIALIZE & REGRID

    properties
        PARENT
    end
    
    
    methods
        
        function self = SNOW_simple_bucketW_seb(index, pprovider, cprovider, forcing)  
            self@INITIALIZE(index, pprovider, cprovider, forcing);
        end
        
        function snow = provide_PARA(snow)
            
            snow.PARA.albedo = [];
            snow.PARA.max_albedo = [];
            snow.PARA.min_albedo = [];
            snow.PARA.tau_1 = [];
            snow.PARA.tau_a = [];
            snow.PARA.tau_f = [];
            
            snow.PARA.epsilon = [];
            snow.PARA.airT_height = []; %measurement height [m]
            snow.PARA.z0 = []; %roughness length [m]
            
            snow.PARA.density = []; 
            snow.PARA.field_capacity = [];
            snow.PARA.hydraulicConductivity = [];
            snow.PARA.swe_per_cell = [];
            
            snow.PARA.area = []; %initial area of the realization [m2]

            snow.PARA.heatFlux_lb = [];
            
            snow.PARA.dt_max = [] ; %[sec]
            snow.PARA.dE_max = []; %[J/m3]
        end
        
        function snow = provide_STATVAR(snow)
            
            snow.STATVAR.upperPos = [];
            snow.STATVAR.lowerPos = [];
            snow.STATVAR.layerThick = []; % [m]
            
            snow.STATVAR.waterIce = []; % [m]
            snow.STATVAR.mineral = []; % [m]
            snow.STATVAR.organic = []; % [m]
            snow.STATVAR.energy = [];  % [J/m2]
            
            snow.STATVAR.T = [];  % [degree C]
            snow.STATVAR.water = [];  % [m]
            snow.STATVAR.ice = [];
            snow.STATVAR.air = [];  % [m]
            snow.STATVAR.thermCond = [];
            snow.STATVAR.albedo = [];
            
            snow.STATVAR.Lstar = [];
            snow.STATVAR.Qh = [];
            snow.STATVAR.Qe = [];
        end
    
        function snow = provide_CONST(snow)
            
            snow.CONST.L_f = [];
            
            snow.CONST.c_w = [];
            snow.CONST.c_i = [];
            snow.CONST.c_o = [];
            snow.CONST.c_m = [];
            
            snow.CONST.k_a = [];       %air [Hillel(1982)]
            snow.CONST.k_w = [];        %water [Hillel(1982)]
            snow.CONST.k_i = [];         %ice [Hillel(1982)]
            snow.CONST.k_o = [];        %organic [Hillel(1982)]
            snow.CONST.k_m = [];
            
            snow.CONST.sigma = []; %Stefan-Boltzmann constant
            snow.CONST.kappa = [];
            snow.CONST.L_s = []; %latent heat of vaporization
            
            snow.CONST.cp = [];
            snow.CONST.g = [];
            
            snow.CONST.rho_w = [];
            snow.CONST.rho_i = [];
        end
        
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function snow = finalize_init(snow, forcing) %assign all variables, that must be calculated or assigned otherwise for initialization
            snow.PARA.heatFlux_lb = forcing.PARA.heatFlux_lb;
            snow.PARA.airT_height = forcing.PARA.airT_height;
            
            snow = initialize_zero_snow_BASE(snow);  %initialize all values to be zero
            snow.TEMP.d_energy = snow.STATVAR.energy .*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
            snow.TEMP.newSnow = 0;
            snow.STATVAR.albedo = snow.PARA.max_albedo;
        end
        
        %---time integration------
        
        function snow = get_boundary_condition_u(snow, forcing) 
            snow = surface_energy_balance(snow, forcing);
            snow = get_boundary_condition_SNOW_u(snow, forcing);
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            snow = get_sublimation(snow, forcing);
        end
        
        function snow = get_boundary_condition_u_CHILD(snow, forcing)
            snow = surface_energy_balance(snow, forcing);
            snow = get_boundary_condition_allSNOW_rain_u(snow, forcing); %add full snow, but rain only for snow-covered part
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            snow = get_sublimation(snow, forcing);
        end
        
        function snow = get_boundary_condition_u_create_CHILD(snow, forcing)
            snow = get_boundary_condition_allSNOW_u(snow, forcing); %add all snow, no rain
            
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
            snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ (snow.PARA.density ./1000); %[m] constant layerThick
            snow.STATVAR.ice = 0;
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos;
        end
        
        function snow = get_boundary_condition_l(snow, forcing)
             snow.TEMP.F_lb = forcing.PARA.heatFlux_lb .* snow.STATVAR.area(end);
             snow.TEMP.d_energy(end) = snow.TEMP.d_energy(end) + snow.TEMP.F_lb;
             snow.TEMP.F_lb_water = 0;
             snow.TEMP.F_lb_water_energy = 0;
        end
        
        function snow = get_derivatives_prognostic(snow)
            if size(snow.STATVAR.layerThick,1) > 1
                snow = get_derivative_energy(snow);
                snow = get_derivative_water_SNOW(snow);
            else
                %snow.TEMP.d_energy = snow.TEMP.d_energy + snow.TEMP.F_ub + snow.TEMP.F_lb;
                %snow.TEMP.d_water = snow.TEMP.d_water + snow.TEMP.F_ub_water + snow.TEMP.F_lb_water;
                %snow.TEMP.d_water_energy = snow.TEMP.d_water_energy + snow.TEMP.F_ub_water_energy + snow.TEMP.F_lb_water_energy;
            end
        end
        
        function snow = get_derivatives_prognostic_CHILD(snow)
            %snow.TEMP.d_energy = snow.TEMP.d_energy + snow.TEMP.F_ub + snow.TEMP.F_lb;
            %snow.TEMP.d_water = snow.TEMP.d_water + snow.TEMP.F_ub_water + snow.TEMP.F_lb_water;
            %snow.TEMP.d_water_energy = snow.TEMP.d_water_energy + snow.TEMP.F_ub_water_energy + snow.TEMP.F_lb_water_energy;
        end
        
        function timestep = get_timestep(snow) 
            timestep1 = get_timestep_heat_coduction(snow);
            timestep2 = get_timestep_SNOW_mass_balance(snow);
            timestep3 = get_timestep_water_SNOW(snow);

            timestep = min(timestep1, timestep2);
            timestep = min(timestep, timestep3);
        end
        
        function timestep = get_timestep_CHILD(snow)  
            timestep = get_timestep_SNOW_CHILD(snow);
        end
        
        function snow = advance_prognostic(snow, timestep) 
            %energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy + snow.TEMP.d_water_energy);
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* (snow.TEMP.snow_energy + snow.TEMP.sublimation_energy);  %rainfall energy already added
            %mass
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* snow.TEMP.d_water;
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + timestep .* (snow.TEMP.snowfall + snow.STATVAR.sublimation);
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* snow.STATVAR.sublimation ./snow.STATVAR.area(1,1) ./ (snow.STATVAR.ice(1) ./ snow.STATVAR.layerThick(1));
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* snow.TEMP.snowfall ./snow.STATVAR.area(1,1) ./ (snow.PARA.density ./1000);
            %store "old" density
            snow.STATVAR.target_density = snow.STATVAR.ice ./ snow.STATVAR.layerThick ./ snow.STATVAR.area;
            snow.STATVAR.target_density(1) = (snow.STATVAR.ice(1) + timestep .* (snow.TEMP.snowfall + snow.STATVAR.sublimation)) ./ snow.STATVAR.layerThick(1) ./ snow.STATVAR.area(1,1);
        
            snow = calculate_albedo_simple(snow, timestep);
        end
        
%         function snow = advance_prognostic_create_CHILD(snow, timestep) %discontinued
%             snow.STATVAR.energy = timestep .* snow.TEMP.snow_energy;
%             snow.STATVAR.waterIce = timestep .* snow.TEMP.snowfall;
%             snow.STATVAR.ice = snow.STATVAR.waterIce; %dry initially
%             snow.STATVAR.water = 0;
%             snow.STATVAR.area = snow.STATVAR.area + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./1000) ./  snow.STATVAR.layerThick ; %[m2]
% 
%             %snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ (snow.PARA.density ./1000); %[m] constant layerThick
%             %snow.STATVAR.area = timestep .* snow.TEMP.snowfall ./ (0.5.*snow.PARA.swe_per_cell) ; %[m2]
%             
%             snow.STATVAR.target_density = snow.PARA.density ./1000;
%         end
        
        function snow = advance_prognostic_CHILD(snow, timestep)

            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy + snow.TEMP.snow_energy + snow.TEMP.d_water_energy + snow.TEMP.sublimation_energy);
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* (snow.TEMP.snowfall + snow.TEMP.d_water + snow.STATVAR.sublimation);
            snow.STATVAR.area = snow.STATVAR.area + timestep .* snow.STATVAR.sublimation ./snow.STATVAR.layerThick ./ max(50, snow.STATVAR.ice ./ snow.STATVAR.layerThick);
            snow.STATVAR.area = snow.STATVAR.area + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./1000) ./  snow.STATVAR.layerThick ; %[m2]
            snow.STATVAR.target_density = min(1,(snow.STATVAR.ice + timestep .* (snow.TEMP.snowfall + snow.STATVAR.sublimation)) ./ snow.STATVAR.layerThick ./ snow.STATVAR.area);
            
            snow = calculate_albedo_simple(snow, timestep);
        end
        
        function snow = compute_diagnostic_first_cell(snow, forcing)
            snow = L_star(snow, forcing);
        end
       
        function snow = compute_diagnostic(snow, forcing)
            snow = get_T_water_freeW(snow);
            snow = subtract_water(snow);
            
            [snow, regridded_yesNo] = regrid_snow(snow, {'waterIce'; 'energy'; 'layerThick'; 'mineral'; 'organic'}, {'area'; 'target_density'}, 'waterIce');
            if regridded_yesNo
                snow = get_T_water_freeW(snow);
            end
            
            snow = conductivity(snow);
            snow.STATVAR.upperPos = snow.NEXT.STATVAR.upperPos + sum(snow.STATVAR.layerThick);
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
        end
        
        function snow = compute_diagnostic_CHILD(snow, forcing)

            snow = get_T_water_freeW(snow);
            snow = subtract_water_CHILD(snow);

            snow = conductivity(snow);
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos + (snow.STATVAR.layerThick .* snow.STATVAR.area ./ snow.PARENT.STATVAR.area(1,1));
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
            
        end
        
        function snow = check_trigger(snow, forcing)
            snow = make_SNOW_CHILD(snow); 
        end
        
        %-----non-mandatory functions-------
        function snow = surface_energy_balance(snow, forcing)
            snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
            snow.STATVAR.Sout = snow.STATVAR.albedo .*  forcing.TEMP.Sin;
            snow.STATVAR.Qh = Q_h(snow, forcing);
            snow.STATVAR.Qe = Q_eq_potET(snow, forcing);
            
            snow.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe) .* snow.STATVAR.area(1);
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.F_ub;
        end
        
        function snow = conductivity(snow)
            snow = conductivity_snow_Yen(snow);
        end
        
        
        %lateral fluxes---------------------------
        function ground = lateral_push_remove_surfaceWater(ground, lateral)
            ground = lateral_push_remove_surfaceWater_simple(ground, lateral);
        end

        
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
    end
    
end
