%========================================================================
% CryoGrid GROUND class GROUND_freezeC_seb
% heat conduction, constant water  +ice, freeze curve based on
% freezing=drying assumption, surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freezeC_ubT < SEB & HEAT_CONDUCTION & FREEZE_CURVE_KarraPainter & HEAT_FLUXES_LATERAL %& INITIALIZE & FREEZE_CURVE_DallAmico 
    
    methods
        
%         function ground = GROUND_freezeC_seb(index, pprovider, cprovider, forcing)  
%             ground@INITIALIZE(index, pprovider, cprovider, forcing);
%         end
        
        function ground = provide_PARA(ground)

            ground.PARA.LUT_size_waterIce = []; %size of lookup table for the waterIce variable [-]
            ground.PARA.LUT_size_T = [];   %size of lookup table for the (temperature) T variable [-]
            ground.PARA.min_T = []; %minimum temperature for which the LUT is calculated (modeled temperatures must be above this value) [degree C]
            ground.PARA.min_waterIce = [];  %minimum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be above this value) [-]
            ground.PARA.max_waterIce = []; %maximum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be below this value) [-]
            ground.PARA.min_mineral_organic = [];  %maximum mineral plus organic content in volumetric fraction for which the LUT is calculated (mineral plus organic content must be below this value) [-]
            
            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = [];  % upper surface elevation [m]
            ground.STATVAR.lowerPos = [];  % lower surface elevation [m]
            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = [];  % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = [];   % total volume of minerals [m3]
            ground.STATVAR.organic = [];   % total volume of organics [m3]
            ground.STATVAR.energy = [];    % total internal energy [J]
            ground.STATVAR.soil_type = []; % integer code for soil_type; 1: sand; 2: silt: 3: clay: 4: peat; 5: water (i.e. approximation of free water, very large-pore ground material).
                        
            ground.STATVAR.T = [];  % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.waterPotential = []; %soil water potential [Pa]
            ground.STATVAR.ice = [];  %total volume of ice [m3]
            ground.STATVAR.air = [];  % total volume of air [m3] - NOT USED
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            
        end
    
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = []; % volumetric latent heat of fusion, freezing
            ground.CONST.Tmfw = []; % freezing temperature of free water [K]
            
            ground.CONST.c_w = []; % volumetric heat capacity water
            ground.CONST.c_i = []; % volumetric heat capacity ice
            ground.CONST.c_o = []; % volumetric heat capacity organic
            ground.CONST.c_m = []; % volumetric heat capacity mineral
            
            ground.CONST.k_a = [];   % thermal conductivity air
            ground.CONST.k_w = [];   % thermal conductivity water
            ground.CONST.k_i = [];   % thermal conductivity ice 
            ground.CONST.k_o = [];   % thermal conductivity organic 
            ground.CONST.k_m = [];   % thermal conductivity mineral 
            
            ground.CONST.sigma = []; %Stefan-Boltzmann constant
            ground.CONST.kappa = []; % von Karman constant
            ground.CONST.L_s = [];  %latent heat of sublimation, latent heat of evaporation handled in a dedicated function
            
            ground.CONST.cp = [];  %specific heat capacity at constant pressure of air
            ground.CONST.g = [];   % gravitational acceleration Earth surface
            
            ground.CONST.rho_w = []; % water density
            ground.CONST.rho_i = []; %ice density
            
            %Mualem Van Genuchten model
            ground.CONST.alpha_water = [];  %alpha parameter for different soil types [m^-1]
            ground.CONST.alpha_sand = [];
            ground.CONST.alpha_silt = [];
            ground.CONST.alpha_clay = [];
            ground.CONST.alpha_peat = [];
            
            ground.CONST.n_water = [];  %n parameter for different soil types [-]
            ground.CONST.n_sand = [];
            ground.CONST.n_silt = [];
            ground.CONST.n_clay = [];
            ground.CONST.n_peat = [];
            
            ground.CONST.residual_wc_water = [];  %residual water content for different soil types [-]
            ground.CONST.residual_wc_sand = [];   %NOTE: this parameter is generally set to 0
            ground.CONST.residual_wc_silt = [];
            ground.CONST.residual_wc_clay = [];
            ground.CONST.residual_wc_peat = [];

        end
        
        function ground = finalize_init(ground, tile)
            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.T .* 0;
            
            ground.CONST.vanGen_alpha = [ ground.CONST.alpha_sand ground.CONST.alpha_silt ground.CONST.alpha_clay ground.CONST.alpha_peat ground.CONST.alpha_water];
            ground.CONST.vanGen_n = [ ground.CONST.n_sand ground.CONST.n_silt ground.CONST.n_clay ground.CONST.n_peat ground.CONST.n_water];
            ground.CONST.vanGen_residual_wc = [ ground.CONST.residual_wc_sand ground.CONST.residual_wc_silt ground.CONST.residual_wc_clay ground.CONST.residual_wc_peat ground.CONST.residual_wc_water];
            
            
            ground = get_E_freezeC(ground);
            ground = conductivity(ground);
            
            ground = create_LUT_freezeC(ground);

            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.F_ub = 0;
            
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile) 
            forcing = tile.FORCING;
            %ground.TEMP.F_ub = (tile.FORCING.TEMP.Tair - ground.STATVAR.T(1,1)) .* ground.STATVAR.thermCond(1,1) ./ground.STATVAR.layerThick(1,1)
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) + (tile.FORCING.TEMP.T_ub - ground.STATVAR.T(1,1)) .* ground.STATVAR.thermCond(1,1) ./ground.STATVAR.layerThick(1,1) .* ground.STATVAR.area(1,1);
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

        end
       
        function ground = compute_diagnostic(ground, tile)
            ground = get_T_water_freezeC(ground);
            ground = conductivity(ground);
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.F_ub = 0;
        end
        
        function ground = check_trigger(ground, tile)
           %do nothing 
        end
        
        
        
        %-----non-mandatory functions-------

        
        
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
        
        
        %----inherited Tier 1 functions ------------
        
        function ground = get_derivative_energy(ground)
           ground = get_derivative_energy@HEAT_CONDUCTION(ground); 
        end
        
        function ground = conductivity_mixing_squares(ground)
            ground = conductivity_mixing_squares@HEAT_CONDUCTION(ground);
        end
        
        function flux = Q_h(ground, forcing)
           flux = Q_h@SEB(ground, forcing);
        end
    
        function flux = Q_eq(ground, forcing)
            flux = Q_eq@SEB(ground, forcing);
        end
        
        function timestep = get_timestep_heat_coduction(ground)
            timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(ground);
        end
        
        function ground = L_star(ground, forcing)
           ground = L_star@SEB(ground, forcing); 
        end
        
        function [ground, S_up] = penetrate_SW_no_transmission(ground, S_down)
            [ground, S_up] = penetrate_SW_no_transmission@SEB(ground, S_down);
        end
        
        function ground = get_T_water_freeW(ground)
            ground = get_T_water_freeW@HEAT_CONDUCTION(ground);
        end
    end
    
end
