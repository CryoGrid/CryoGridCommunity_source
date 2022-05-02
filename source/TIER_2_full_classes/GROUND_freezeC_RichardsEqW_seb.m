%========================================================================
% CryoGrid GROUND class GROUND_freezeC_RichardsEqW_seb
% heat conduction, Richards equation water balance, Painter and Karra (2014) 
% freeze curve, surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freezeC_RichardsEqW_seb < SEB & HEAT_CONDUCTION & FREEZE_CURVE_KarraPainter & WATER_FLUXES & HEAT_FLUXES_LATERAL & WATER_FLUXES_LATERAL 

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
 
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];  %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; % roughness length [m] 

            ground.PARA.permeability = [];  %permeability for fluids/gases [m2]

            ground.PARA.conductivity_function = [];
            
            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
            ground.PARA.dWater_max = []; %%maximum possible volumteric water content change per timestep [-] 

            ground.PARA.LUT_size_waterIce = []; %size of lookup table for the waterIce variable [-]
            ground.PARA.LUT_size_T = [];   %size of lookup table for the (temperature) T variable [-]
            ground.PARA.min_T = []; %minimum temperature for which the LUT is calculated (modeled temperatures must be above this value) [degree C]
            ground.PARA.min_waterIce = [];  %minimum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be above this value) [-]
            ground.PARA.max_waterIce = []; %maximum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be below this value) [-]
            ground.PARA.min_mineral_organic = [];  %maximum mineral plus organic content in volumetric fraction for which the LUT is calculated (mineral plus organic content must be below this value) [-]
                      
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
            ground.STATVAR.hydraulicConductivity = []; % hydraulic conductivity [m/sec]
            ground.STATVAR.permeability = [];  %permeability for fluids/gases [m2]
            
            ground.STATVAR.Lstar = [];  %Obukhov length [m]
            ground.STATVAR.Qh = [];     %sensible heat flux [W/m2]
            ground.STATVAR.Qe = [];     % latent heat flux [W/m2]
            
            ground.STATVAR.field_capacity = [];  %field capacity in fraction of the total volume [-]
            ground.STATVAR.excessWater = 0;  %water volume overtopping first grid cell (i.e. surface water) [m3]
            
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
            
            ground.CONST.R = [];
            ground.CONST.molar_mass_w = [];
            
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
        
        function ground = convert_units(ground, tile)
            unit_converter = str2func(tile.PARA.unit_conversion_class);
            unit_converter = unit_converter();
            ground = convert_normal(unit_converter, ground, tile);
        end
        
        function ground = finalize_init(ground, tile)
%             ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
%             ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
%             ground.STATVAR.area = tile.PARA.area + ground.STATVAR.T .* 0;
            
            if isempty(ground.PARA.conductivity_function) || sum(isnan(ground.PARA.conductivity_function))>0
                ground.PARA.conductivity_function = 'thermalConductivity_CLM4_5';
            end
            
            ground.CONST.vanGen_alpha = [ ground.CONST.alpha_sand ground.CONST.alpha_silt ground.CONST.alpha_clay ground.CONST.alpha_peat ground.CONST.alpha_water];
            ground.CONST.vanGen_n = [ ground.CONST.n_sand ground.CONST.n_silt ground.CONST.n_clay ground.CONST.n_peat ground.CONST.n_water];
            ground.CONST.vanGen_residual_wc = [ ground.CONST.residual_wc_sand ground.CONST.residual_wc_silt ground.CONST.residual_wc_clay ground.CONST.residual_wc_peat ground.CONST.residual_wc_water];
            
            ground = get_E_freezeC(ground);
            ground = conductivity(ground);
            ground = calculate_hydraulicConductivity_RichardsEq(ground);
            
            ground = create_LUT_freezeC(ground);

            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.STATVAR.runoff = 0;
            
            ground.STATVAR.year_old = []; %str2num(datestr(tile.t, 'yyyy'));
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.TEMP.surface_runoff = 0;
        end
        
        function ground = finalize_init2(ground, tile)

            ground = get_E_freezeC(ground);
            ground = conductivity(ground);
            ground = calculate_hydraulicConductivity_RichardsEq(ground);

        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            ground = surface_energy_balance(ground, forcing);
            ground = get_boundary_condition_u_RichardsEq(ground, forcing); %checked that this flux can be taken up!!
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
            ground = get_boundary_condition_l_water2(ground);  %if flux not zero, check that the water flowing out is available! Not implemented here.
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivative_energy(ground);
            ground = get_derivative_water_RichardsEq(ground);
        end
        
        function timestep = get_timestep(ground, tile)  
           timestep = get_timestep_heat_coduction(ground);
           timestep = min(timestep, get_timestep_water_RichardsEq(ground)); 
        end
        
        function ground = advance_prognostic(ground, tile) 
            timestep = tile.timestep;
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_water_energy; %add energy from water advection
            %water
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + timestep .* ground.TEMP.d_water; %subtract water from ET
            ground.STATVAR.waterIce = max(0, ground.STATVAR.waterIce);
            
            %ground.STATVAR.waterIce = min(ground.STATVAR.waterIce, ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic); %prevent small rounding errors 
            ground.STATVAR.excessWater = ground.STATVAR.excessWater + timestep .* ground.TEMP.surface_runoff;
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)
            
            ground = get_T_water_freezeC(ground);
            ground = conductivity(ground);
            ground = calculate_hydraulicConductivity_RichardsEq(ground);
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.TEMP.surface_runoff = 0;
        end
        
        function ground = check_trigger(ground, tile)
            %do nothing at this stage - add create LAKE if there is too
            %much surface water, similar to Xice 
            % Every year, the heat reservoir is 10% further away
%             if isempty(ground.STATVAR.year_old) || str2num(datestr(tile.t, 'yyyy'))>ground.STATVAR.year_old
%                 ground.STATVAR.year_old = str2num(datestr(tile.t, 'yyyy'));
%                 tile.LATERAL.IA_CLASSES{2,1}.PARA.distance_heatReservoir = tile.LATERAL.IA_CLASSES{2,1}.PARA.distance_heatReservoir .*1.1;
%             end
        end
        
        
        %-----non-mandatory functions-------
%         function ground = surface_energy_balance(ground, forcing)
%             ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
%             ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
%             ground.STATVAR.Qh = Q_h(ground, forcing);
%             ground.STATVAR.Qe_pot = Q_eq_potET(ground, forcing);
% 
%             ground = calculateET(ground);
%             
%             ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
%             ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
%         end
        
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground = Q_evap_CLM4_5(ground, forcing);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
            
            %water -> evaporation
            ground.TEMP.d_water_ET(1,1) = ground.TEMP.d_water_ET(1,1) -  ground.STATVAR.evap.* ground.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
            ground.TEMP.d_water_ET_energy(1,1) = ground.TEMP.d_water_ET_energy(1,1) -  ground.STATVAR.evap_energy.* ground.STATVAR.area(1,1);

        end
        
        function ground = conductivity(ground)
            conductivity_function = str2func(ground.PARA.conductivity_function);
            ground = conductivity_function(ground);
            %ground = thermalConductivity_CLM4_5(ground);
            %ground = conductivity_mixing_squares(ground);
        end
        
        %-----LATERAL-------------------
        
        %-----LAT_REMOVE_SURFACE_WATER-----
        function ground = lateral_push_remove_surfaceWater(ground, lateral)
            ground = lateral_push_remove_surfaceWater_simple(ground, lateral);
        end
        
        %-----LAT_REMOVE_SUBSURFACE_WATER-----        
        function ground = lateral_push_remove_subsurfaceWater(ground, lateral)
            ground = lateral_push_remove_subsurfaceWater_simple(ground, lateral);
        end
        
        %----LAT_SEEPAGE_FACE----------            
        function ground = lateral_push_remove_water_seepage(ground, lateral)
            ground = lateral_push_remove_water_seepage_simple(ground, lateral);
        end
        
        %----LAT_WATER_RESERVOIR------------          
        function ground = lateral_push_water_reservoir(ground, lateral)
            ground = lateral_push_water_reservoir_simple(ground, lateral);
        end
        
        %----LAT_HEAT------------
        function ground = lateral_push_heat(ground, lateral)
            ground = lateral_push_heat_simple(ground, lateral);
        end

        %----LAT3D_WATER_UNCONFINED_AQUIFER------------         
        function ground = lateral3D_pull_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_pull_water_unconfined_aquifer_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_push_water_unconfined_aquifer_simple(ground, lateral);
        end
        
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground, lateral)
            [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_simple(ground, lateral);
        end
        
        %LAT3D_WATER_RESERVOIR and LAT3D_WATER_SEEPAGE_FACE do not require specific functions
         
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
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE(ground);
            
            ground.PARA.class_category = 'GROUND';
            
            %ground.PARA.options = [];
            ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'soil_type' 'field_capacity' 'permeability' 'T'};
            
            ground.PARA.default_value.albedo = {0.2};
            ground.PARA.comment.albedo = {'surface albedo [-]'};
            
            ground.PARA.default_value.epsilon = {0.99};
            ground.PARA.comment.epsilon = {'surface emissivity [-]'};
            
            ground.PARA.default_value.z0 = {0.01};
            ground.PARA.comment.z0 = {'roughness length [m]'};
            
            ground.PARA.default_value.rootDepth = {0.1};
            ground.PARA.comment.rootDepth = {'e-folding depth of transpiration reduction with depth [m]'};
            
            ground.PARA.default_value.evaporationDepth = {0.1};
            ground.PARA.comment.evaporationDepth = {'e-folding constant of evaporation reduction reduction with depth [m]'};
            
            ground.PARA.default_value.ratioET = {0.5};
            ground.PARA.comment.ratioET = {'fraction of transpiration of total evapotranspiration [-]'};
            
            ground.PARA.default_value.conductivity_function = {''};
            ground.PARA.comment.conductivity_function = {'function employed to calculate thermal conductivity, leave empty for default'};
            
            ground.PARA.default_value.dt_max = {3600};
            ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
            
            ground.PARA.default_value.dE_max = {50000};
            ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
            
            ground.PARA.default_value.dWater_max = {0.005};
            ground.PARA.comment.dWater_max = {'maximum possible volumteric water content change per timestep [-]'};
            
            ground.PARA.default_value.LUT_size_waterIce = {1000};
            ground.PARA.comment.LUT_size_waterIce = {'size of lookup table for the waterIce variable [-]'};
            
            ground.PARA.default_value.LUT_size_T = {1000};
            ground.PARA.comment.LUT_size_T = {'size of lookup table for the (temperature) T variable [-]'};
            
            ground.PARA.default_value.min_T = {-50};
            ground.PARA.comment.min_T = {'minimum temperature for which the LUT is calculated (modeled temperatures must be above this value) [degree C]'};
            
            ground.PARA.default_value.min_waterIce = {0.05};
            ground.PARA.comment.min_waterIce = {'minimum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be above this value) [-]'};
            
            ground.PARA.default_value.max_waterIce = {0.97};
            ground.PARA.comment.max_waterIce = {'maximum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be below this value) [-]'};
            
            ground.PARA.default_value.min_mineral_organic = {0.03};
            ground.PARA.comment.min_mineral_organic = {'maximum mineral plus organic content in volumetric fraction for which the LUT is calculated (mineral plus organic content must be below this value) [-]'};
            
        end
    end
    
end
