%========================================================================
% CryoGrid GROUND class GROUND_freezeC_bucketW_seb
% heat conduction, bucket water scheme, freeze curve based on
% freezing=drying assumption, surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freezeC_bucketW_seb_vegetation < GROUND_freezeC_bucketW_seb 
    
    
    properties
        VEGETATION
        IA_VEGETATION_GROUND   %exchange of water
        IA_VEGETATION_SURFACE  %surface energy balance, repositioned in case of a SNOW cover
    end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        

        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];  %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; % roughness length [m] 

            ground.PARA.rootDepth = []; %e-folding constant of transpiration reduction with depth [1/m]
            ground.PARA.evaporationDepth = []; %e-folding constant of evaporation reduction reduction with depth [1/m]
            ground.PARA.ratioET = []; %fraction of transpiration of total evapotranspiration [-]
            ground.PARA.hydraulicConductivity = [];  %saturated hydraulic conductivity [m/sec]

            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
            ground.PARA.dWater_max = [];

            ground.PARA.LUT_size_waterIce = []; %size of lookup table for the waterIce variable [-]
            ground.PARA.LUT_size_T = [];   %size of lookup table for the (temperature) T variable [-]
            ground.PARA.min_T = []; %minimum temperature for which the LUT is calculated (modeled temperatures must be above this value) [degree C]
            ground.PARA.min_waterIce = [];  %minimum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be above this value) [-]
            ground.PARA.max_waterIce = []; %maximum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be below this value) [-]
            ground.PARA.min_mineral_organic = [];  %maximum mineral plus organic content in volumetric fraction for which the LUT is calculated (mineral plus organic content must be below this value) [-]
            
            ground.PARA.VEGETATION_class = [];
            ground.PARA.VEGETATION_class_index = [];
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
            
            ground = finalize_init@GROUND_freezeC_bucketW_seb(ground, tile);
            
%             class_handle = str2func(ground.PARA.VEGETATION_class);
%         
%             ground.VEGETATION = class_handle(); 
            
            ground.VEGETATION = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(ground.PARA.VEGETATION_class){ground.PARA.VEGETATION_class_index,1});
            ground.VEGETATION.PARENT_GROUND = ground;
            ground.VEGETATION.PARENT_SURFACE = ground;  %ATTENTION: reset this pointer when making the transition to a full SNOW class

%             ground.BGC = provide_PARA(ground.BGC);
%             ground.BGC = provide_STATVAR(ground.BGC);
%             ground.BGC = provide_CONST(ground.BGC);
            ground.VEGETATION = finalize_init(ground.VEGETATION, tile);
            
            %make this a get_IA_class for vegetation in the end
%             ground.IA_VEGETATION_GROUND = IA_ml_canopy_GROUND();
%             ground.IA_VEGETATION_GROUND.GROUND = ground;
%             ground.IA_VEGETATION_GROUND.VEGETATION = ground.VEGETATION;
%             ground.VEGETATION.IA_VEGETATION_GROUND = ground.IA_VEGETATION_GROUND;
%             
%             ground.IA_VEGETATION_SURFACE = IA_ml_canopy_SURFACE();
%             ground.IA_VEGETATION_SURFACE.SURFACE = ground;
%             ground.IA_VEGETATION_SURFACE.VEGETATION = ground.VEGETATION;
%             ground.VEGETATION.IA_VEGETATION_SURFACE = ground.IA_VEGETATION_SURFACE;
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile) 
            %calculate SEB for canopy and adjusts the forcing data
            ground.VEGETATION = get_boundary_condition_u(ground.VEGETATION, tile);
            
            forcing = ground.VEGETATION.ForcingV;

            
            %radiation balance, Sout and Lout
            ground = calculate_radiation(ground, forcing);
            
            ground.VEGETATION = map_variables_no_snow(ground.VEGETATION);

%             %evaporation should eventually be handled in GROUND class - not possible for current implementation multi-layer canopy
%             ground = get_evaporation(ground, forcing);
            ground.VEGETATION = get_evaporation(ground.VEGETATION);

            ground.VEGETATION  = get_transpiration(ground.VEGETATION);
            
             %upper BC water (rainfall)
            ground = get_boundary_condition_u_water2(ground, forcing);            
        end

        
        function ground = get_boundary_condition_l(ground, tile)
            ground.VEGETATION = get_boundary_condition_l(ground.VEGETATION, tile);
            ground = get_boundary_condition_l@GROUND_freezeC_bucketW_seb(ground, tile);
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground.VEGETATION = get_derivatives_prognostic(ground.VEGETATION, tile);
            ground = get_derivatives_prognostic@GROUND_freezeC_bucketW_seb(ground, tile);
        end
        
        function timestep = get_timestep(ground, tile) 
            timestep_ground = get_timestep@GROUND_freezeC_bucketW_seb(ground, tile);
            timestep_vegetation = get_timestep(ground.VEGETATION, tile);
            
            timestep = min(timestep_ground, timestep_vegetation);
        end
        
        function ground = advance_prognostic(ground, tile) 
            ground.VEGETATION = advance_prognostic(ground.VEGETATION, tile);
            ground = advance_prognostic@GROUND_freezeC_bucketW_seb(ground, tile);
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            ground.VEGETATION = compute_diagnostic_first_cell(ground.VEGETATION, tile);
            ground = compute_diagnostic_first_cell@GROUND_freezeC_bucketW_seb(ground, tile);
        end
       
        function ground = compute_diagnostic(ground, tile)
            ground.VEGETATION = compute_diagnostic(ground.VEGETATION, tile);
            ground = compute_diagnostic@GROUND_freezeC_bucketW_seb(ground, tile);
        end
        
        function ground = check_trigger(ground, tile)
            ground.VEGETATION = check_trigger(ground.VEGETATION, tile);
            ground = check_trigger@GROUND_freezeC_bucketW_seb(ground, tile);
        end
        
%         
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
        
        function ground = calculate_radiation(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
        end
        
%         function ground = get_evaporation(ground, forcing)  %needed eventually for real vegetation scheme that is coupled to GROUND scheme
%             
%         end
        
%         function ground = conductivity(ground)
%             ground = conductivity_mixing_squares(ground);
%         end
%         
%         
%         %-----LATERAL-------------------
%         
%         %-----LAT_REMOVE_SURFACE_WATER-----
%         function ground = lateral_push_remove_surfaceWater(ground, lateral)
%             ground = lateral_push_remove_surfaceWater_simple(ground, lateral);
%         end
%         
%         %-----LAT_REMOVE_SUBSURFACE_WATER-----
%         function ground = lateral_push_remove_subsurfaceWater(ground, lateral)
%             ground = lateral_push_remove_subsurfaceWater_simple(ground, lateral);
%         end
%         
%         %----LAT_SEEPAGE_FACE----------              
%         function ground = lateral_push_remove_water_seepage(ground, lateral)
%             ground = lateral_push_remove_water_seepage_simple(ground, lateral);
%         end
%         
%         %----LAT_WATER_RESERVOIR------------  
%         function ground = lateral_push_water_reservoir(ground, lateral)
%             ground = lateral_push_water_reservoir_simple(ground, lateral);
%         end
% 
%         %----LAT3D_WATER_UNCONFINED_AQUIFER------------ 
%         function ground = lateral3D_pull_water_unconfined_aquifer(ground, lateral)
%             ground = lateral3D_pull_water_unconfined_aquifer_simple(ground, lateral);
%         end
%         
%         function ground = lateral3D_push_water_unconfined_aquifer(ground, lateral)
%             ground = lateral3D_push_water_unconfined_aquifer_simple(ground, lateral);
%         end
%         
%          %----LAT3D_WATER------------ 
%         function ground = lateral3D_pull_water_general_aquifer(ground, lateral)
%             ground = lateral3D_pull_water_general_aquifer_simple(ground, lateral);
%         end
%         
%         function ground = lateral3D_push_water_general_aquifer(ground, lateral)
%             ground = lateral3D_push_water_general_aquifer_simple(ground, lateral);
%         end
%         
%         function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground, lateral)
%             [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_simple(ground, lateral);
%         end
%         
%          %LAT3D_WATER_RESERVOIR and LAT3D_WATER_SEEPAGE_FACE do not require specific functions
%         
%         %-------LAT3D_HEAT-------------
%         function ground = lateral3D_pull_heat(ground, lateral)
%             ground = lateral3D_pull_heat_simple(ground, lateral);
%         end
%         
%         function ground = lateral3D_push_heat(ground, lateral)
%             ground = lateral3D_push_heat_simple(ground, lateral);
%         end
%         
%         
%         %----inherited Tier 1 functions ------------
%         
%         function ground = get_derivative_energy(ground)
%            ground = get_derivative_energy@HEAT_CONDUCTION(ground); 
%         end
%         
%         function ground = conductivity_mixing_squares(ground)
%             ground = conductivity_mixing_squares@HEAT_CONDUCTION(ground);
%         end
%         
%         function flux = Q_h(ground, forcing)
%            flux = Q_h@SEB(ground, forcing);
%         end
%     
%         function flux = Q_eq_potET(ground, forcing)
%             flux = Q_eq_potET@SEB(ground, forcing);
%         end
%         
%         function ground = calculateET(ground)
%             ground = calculateET@SEB(ground);
%         end
%         
%         function ground = get_boundary_condition_u_water2(ground, forcing)
%            ground = get_boundary_condition_u_water2@WATER_FLUXES(ground, forcing);
%         end
%         function ground = get_derivative_water2(ground)
%             ground = get_derivative_water2@WATER_FLUXES(ground);
%         end
%         
%         function timestep = get_timestep_heat_coduction(ground)
%             timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(ground);
%         end
%         
%         function timestep = get_timestep_water(ground)
%             timestep = get_timestep_water@WATER_FLUXES(ground);
%         end
%         
%         function ground = L_star(ground, forcing)
%            ground = L_star@SEB(ground, forcing); 
%         end
%         
%         function [ground, S_up] = penetrate_SW_no_transmission(ground, S_down)
%             [ground, S_up] = penetrate_SW_no_transmission@SEB(ground, S_down);
%         end
%         
%         function ground = get_T_water_freeW(ground)
%             ground = get_T_water_freeW@HEAT_CONDUCTION(ground);
%         end
    end
    
end
