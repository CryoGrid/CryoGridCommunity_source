

classdef GROUND_freezeC_bucketW_seb < SEB & HEAT_CONDUCTION & FREEZE_CURVE & WATER_FLUXES & HEAT_FLUXES_LATERAL & WATER_FLUXES_LATERAL & INITIALIZE

    
    methods
        
        function self = GROUND_freezeC_bucketW_seb(index, pprovider, cprovider, forcing)  
            self@INITIALIZE(index, pprovider, cprovider, forcing);
        end
        
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = [];
            ground.CONST.Tmfw = [];
            
            ground.CONST.c_w = [];
            ground.CONST.c_i = [];
            ground.CONST.c_o = [];
            ground.CONST.c_m = [];
            
            ground.CONST.k_a = []; %air [Hillel(1982)]
            ground.CONST.k_w = []; %water [Hillel(1982)]
            ground.CONST.k_i = []; %ice [Hillel(1982)]
            ground.CONST.k_o = []; %organic [Hillel(1982)]
            ground.CONST.k_m = [];
            
            ground.CONST.sigma = []; %Stefan-Boltzmann constant
            ground.CONST.kappa = [];
            ground.CONST.L_s = []; %latent heat of vaporization
            
            ground.CONST.cp = [];
            ground.CONST.g = [];
            
            ground.CONST.rho_w = [];
            ground.CONST.rho_i = [];
            
            %Mualem Van Genuchten model
            ground.CONST.alpha_water = [];
            ground.CONST.alpha_sand = [];
            ground.CONST.alpha_silt = [];
            ground.CONST.alpha_clay = [];
            ground.CONST.alpha_peat = [];
            
            ground.CONST.n_water = [];
            ground.CONST.n_sand = [];
            ground.CONST.n_silt = [];
            ground.CONST.n_clay = [];
            ground.CONST.n_peat = [];
            
            ground.CONST.residual_wc_water = [];
            ground.CONST.residual_wc_sand = [];
            ground.CONST.residual_wc_silt = [];
            ground.CONST.residual_wc_clay = [];
            ground.CONST.residual_wc_peat = [];

        end
        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];
            ground.PARA.epsilon = [];
            % ground.PARA.airT_height = []; % measurement height [m] Assigned from forcing data in "finialize"
            ground.PARA.z0 = []; % roughness length [m]
            
            ground.PARA.area =[]; % initial area of the realization [m2]
            
            ground.PARA.rootDepth = [];
            ground.PARA.evaporationDepth = [];
            ground.PARA.ratioET = [];
            
            ground.PARA.hydraulicConductivity = []; %change to a prametrization later?
            % ground.PARA.externalWaterFlux = [] ; % NB: not obtained from para file, where is this assigned?
            
            % ground.PARA.heatFlux_lb = []; % Assigned from forcing data in "finialize"
            
            ground.PARA.dt_max = [] ; %[sec]
            ground.PARA.dE_max = []; %[J/m3]

            ground.PARA.LUT_size_waterIce = [];
            ground.PARA.LUT_size_T = [];
            ground.PARA.min_T = []; %minimum an maximum values for which the LUT is calculated (only roughly)
            ground.PARA.min_waterIce = [];
            ground.PARA.max_waterIce = [];
            ground.PARA.min_mineral_organic = [];
                      
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = [];
            ground.STATVAR.lowerPos = [];
            ground.STATVAR.layerThick = []; % [m]
            
            ground.STATVAR.waterIce = []; % [m]
            ground.STATVAR.mineral = []; % [m]
            ground.STATVAR.organic = []; % [m]
            ground.STATVAR.energy = [];  % [J/m2]
            ground.STATVAR.soil_type = [];
                        
            ground.STATVAR.T = [];  % [degree C]
            ground.STATVAR.water = [];  % [m]
            ground.STATVAR.waterPotential = []; % [m]
            ground.STATVAR.ice = [];
            ground.STATVAR.air = [];  % [m]
            ground.STATVAR.thermCond = [];
            ground.STATVAR.hydraulicConductivity = [];
            
            ground.STATVAR.Lstar = [];
            ground.STATVAR.Qh = [];
            ground.STATVAR.Qe = [];
            
            ground.STATVAR.field_capacity = [];
            ground.STATVAR.excessWater = 0;
            
        end
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = finalize_init(ground, forcing) %assign all variables, that must be calculated or assigned otherwise for initialization
            ground.PARA.heatFlux_lb = forcing.PARA.heatFlux_lb;
            ground.PARA.airT_height = forcing.PARA.airT_height;
            ground.STATVAR.area = ground.PARA.area + ground.STATVAR.T .* 0;
            
            ground.CONST.vanGen_alpha = [ ground.CONST.alpha_sand ground.CONST.alpha_silt ground.CONST.alpha_clay ground.CONST.alpha_peat ground.CONST.alpha_water];
            ground.CONST.vanGen_n = [ ground.CONST.n_sand ground.CONST.n_silt ground.CONST.n_clay ground.CONST.n_peat ground.CONST.n_water];
            ground.CONST.vanGen_residual_wc = [ ground.CONST.residual_wc_sand ground.CONST.residual_wc_silt ground.CONST.residual_wc_clay ground.CONST.residual_wc_peat ground.CONST.residual_wc_water];
            
            ground = get_E_freezeC(ground);
            ground = conductivity(ground);
            ground = calculate_hydraulicConductivity(ground);
            
            ground = create_LUT_freezeC(ground);

            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.STATVAR.runoff = 0;
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.TEMP.surface_runoff = 0;
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
      
            ground = surface_energy_balance(ground, forcing);
            ground = get_boundary_condition_u_water2(ground, forcing); %checked that this flux can be taken up!!
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
             ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
             ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
             ground = get_boundary_condition_l_water2(ground);  %if flux not zero, check that the water flowing out is available! Not implemented here.
        end
        
        function ground = get_derivatives_prognostic(ground)
            ground = get_derivative_energy(ground);
            ground = get_derivative_water2(ground);
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
           timestep = get_timestep_heat_coduction(ground);
           timestep = min(timestep, get_timestep_water(ground)); 
           %add check for water fluxes timestep
        end
        
        function ground = advance_prognostic(ground, timestep) 
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_water_energy; %add energy from water advection
            %water
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + timestep .* ground.TEMP.d_water; %subtract water from ET
            %ground.STATVAR.waterIce = min(ground.STATVAR.waterIce, ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic); %prevent small rounding errors 
            ground.STATVAR.excessWater = ground.STATVAR.excessWater + timestep .* ground.TEMP.surface_runoff;
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, forcing)
            
            ground = get_T_water_freezeC(ground);
            ground = conductivity(ground);
            ground = calculate_hydraulicConductivity(ground);
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.TEMP.surface_runoff = 0;
        end
        
        function ground = check_trigger(ground, forcing)
            %do nothing, but could create LAKE if there is too much surface
            %water
        end
        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe_pot = Q_eq_potET(ground, forcing);

            ground = calculateET(ground);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
        end
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        %lateral fluxes---------------------------
        function ground = lateral_push_remove_surfaceWater(ground, lateral)
            ground = lateral_push_remove_surfaceWater_simple(ground, lateral);
        end
        
        function ground = lateral_push_remove_subsurfaceWater(ground, lateral)
            ground = lateral_push_remove_subsurfaceWater_simple(ground, lateral);
        end
        
        function ground = lateral_push_remove_water_seepage(ground, lateral)
            ground = lateral_push_remove_water_seepage_simple(ground, lateral);
        end
        
        function ground = lateral_push_water_reservoir(ground, lateral)
            ground = lateral_push_water_reservoir_simple(ground, lateral);
        end

        function ground = lateral3D_pull_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_pull_water_unconfined_aquifer_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_push_water_unconfined_aquifer_simple(ground, lateral);
        end
        
        function ground = lateral3D_pull_water_general_aquifer(ground, lateral)
            ground = lateral3D_pull_water_general_aquifer_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_water_general_aquifer(ground, lateral)
            ground = lateral3D_push_water_general_aquifer_simple(ground, lateral);
        end
        
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground, lateral)
            [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_simple(ground, lateral);
        end
        
        
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
    
        function flux = Q_eq_potET(ground, forcing)
            flux = Q_eq_potET@SEB(ground, forcing);
        end
        
        function ground = calculateET(ground)
            ground = calculateET@SEB(ground);
        end
        
        function ground = get_boundary_condition_u_water2(ground, forcing)
           ground = get_boundary_condition_u_water2@WATER_FLUXES(ground, forcing);
        end
        function ground = get_derivative_water2(ground)
            ground = get_derivative_water2@WATER_FLUXES(ground);
        end
        
        function timestep = get_timestep_heat_coduction(ground)
            timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(ground);
        end
        
        function timestep = get_timestep_water(ground)
            timestep = get_timestep_water@WATER_FLUXES(ground);
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
