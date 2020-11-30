%========================================================================
% CryoGrid GROUND class GROUND_freeW_bucketW_seb
% heat conduction, free water freeze curve, surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freeW_seb < SEB & HEAT_CONDUCTION & HEAT_FLUXES_LATERAL & INITIALIZE

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function self = GROUND_freeW_seb(index, pprovider, cprovider, forcing)  
            self@INITIALIZE(index, pprovider, cprovider, forcing);
        end
        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = []; %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; %roughness length [m]
            
            ground.PARA.rs = []; %surface resistance against evapotranspiration [secm] 
            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = [];  % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = [];  % total volume of minerals [m3]
            ground.STATVAR.organic = [];  %total volume of ice [m3]
            ground.STATVAR.energy = [];  %total molar salt volume within a grid cell [mol]
            
            ground.STATVAR.T = [];   % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.ice = [];    %total volume of ice [m3]
            ground.STATVAR.air = [];   % total volume of air [m3] - NOT USED
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            
            ground.STATVAR.Lstar = [];  %Obukhov length [m]
            ground.STATVAR.Qh = []; %sensible heat flux [W/m2]
            ground.STATVAR.Qe = [];  % latent heat flux [W/m2]
        end
    
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = [];  % volumetric latent heat of fusion, freezing
            ground.CONST.c_w = [];  % volumetric heat capacity water
            ground.CONST.c_i = [];  % volumetric heat capacity ice
            ground.CONST.c_o = [];  % volumetric heat capacity organic
            ground.CONST.c_m = [];  % volumetric heat capacity mineral
            
            ground.CONST.k_a = [];   % thermal conductivity air
            ground.CONST.k_w = [];   % thermal conductivity water
            ground.CONST.k_i = [];   % thermal conductivity ice 
            ground.CONST.k_o = [];   % thermal conductivity organic 
            ground.CONST.k_m = [];   % thermal conductivity mineral 
            
            ground.CONST.sigma = []; %Stefan-Boltzmann constant
            ground.CONST.kappa = []; % von Karman constant
            ground.CONST.L_s = [];  %latent heat of sublimation, latent heat of evaporation handled in a dedicated function
            
            ground.CONST.cp = []; % specific heat capacity at constant pressure of air
            ground.CONST.g = []; % gravitational acceleration Earth surface
            
            ground.CONST.rho_w = [];   % water density
            ground.CONST.rho_i = [];   %ice density
        end
        
        function ground = finalize_init(ground, forcing)
            ground.PARA.heatFlux_lb = forcing.PARA.heatFlux_lb;
            ground.PARA.airT_height = forcing.PARA.airT_height;
            ground.STATVAR.area = forcing.PARA.area + ground.STATVAR.T .* 0;
            
            ground = get_E_freeW(ground);

            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            ground = surface_energy_balance(ground, forcing);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration 
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivative_energy(ground);

        end
        
        function timestep = get_timestep(ground, tile) 
           timestep = get_timestep_heat_coduction(ground);
        end
        
        function ground = advance_prognostic(ground, tile)           
            timestep = tile.timestep;
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)
            ground = get_T_water_freeW(ground);
            ground = conductivity(ground);
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
        end
        
        function ground = check_trigger(ground, tile)
           %do nothing 
        end
        

        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq(ground, forcing);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
        end
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        
        %-----LATERAL-------------------
        
        %-------LAT3D_HEAT-------------
        function ground = lateral3D_pull_heat(ground, lateral)
            ground = lateral3D_pull_heat_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_heat(ground, lateral)
            ground = lateral3D_push_heat_simple(ground, lateral);
        end

    end
    
end
