%========================================================================
% CryoGrid GROUND class GROUND_all_constant
% initialize as for GROUND_freeW_seb, and then do nothing
% S. Westermann, October 2020
%========================================================================

classdef GROUND_all_constant < SEB & HEAT_CONDUCTION & WATER_FLUXES & FREEZE_CURVE_KarraPainter

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];  %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; %roughness length [m]
            
            ground.PARA.rootDepth = []; %e-folding constant of transpiration reduction with depth [m]
            ground.PARA.evaporationDepth = []; %e-folding constant of evaporation reduction reduction with depth [m]
            ground.PARA.ratioET = []; %fraction of transpiration of total evapotranspiration [-]
            ground.PARA.hydraulicConductivity = [];  %saturated hydraulic conductivity [m/sec]
            
            ground.PARA.dt_max = [];  %maximum possible timestep [sec]
            ground.PARA.dE_max = [];  %maximum possible energy change per timestep [J/m3]
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.XwaterIce = [];
            ground.STATVAR.Xice = [];
            ground.STATVAR.mineral = []; % total volume of minerals [m3]
            ground.STATVAR.organic = []; % total volume of organics [m3]
            ground.STATVAR.energy = [];  % total internal energy [J]
            ground.STATVAR.soil_type = [];  

            ground.STATVAR.T = [];  % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.ice = []; %total volume of ice [m3]
            ground.STATVAR.air = [];  % total volume of air [m3] - NOT USED
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            ground.STATVAR.hydraulicConductivity = []; % hydraulic conductivity [m/sec]
            
            ground.STATVAR.Lstar = []; %Obukhov length [m]
            ground.STATVAR.Qh = []; %sensible heat flux [W/m2]
            ground.STATVAR.Qe = []; % latent heat flux [W/m2]
            
            ground.STATVAR.field_capacity = []; %field capacity in fraction of the total volume [-]
            ground.STATVAR.excessWater = 0;  %water volume overtopping first grid cell (i.e. surface water [m3]
            
        end
        
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = [];  % volumetric latent heat of fusion, freezing
            ground.CONST.Tmfw = [];  % freezing temperature of free water [K]
             
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
            ground.CONST.kappa = [];  % von Karman constant
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
            
            ground.CONST.n_water = [];   %n parameter for different soil types [-]
            ground.CONST.n_sand = [];
            ground.CONST.n_silt = [];
            ground.CONST.n_clay = [];
            ground.CONST.n_peat = [];
            
            ground.CONST.residual_wc_water = [];    %residual water content for different soil types [-]
            ground.CONST.residual_wc_sand = [];    %NOTE: this parameter is generally set to 0
            ground.CONST.residual_wc_silt = [];
            ground.CONST.residual_wc_clay = [];
            ground.CONST.residual_wc_peat = [];
        end
        
        function ground = finalize_init(ground, tile) 
%             ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
%             ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.T .* 0;
            
            ground.CONST.vanGen_alpha = [ ground.CONST.alpha_sand ground.CONST.alpha_silt ground.CONST.alpha_clay ground.CONST.alpha_peat ground.CONST.alpha_water];
            ground.CONST.vanGen_n = [ ground.CONST.n_sand ground.CONST.n_silt ground.CONST.n_clay ground.CONST.n_peat ground.CONST.n_water];
            ground.CONST.vanGen_residual_wc = [ ground.CONST.residual_wc_sand ground.CONST.residual_wc_silt ground.CONST.residual_wc_clay ground.CONST.residual_wc_peat ground.CONST.residual_wc_water];
            
            
            %ground = get_E_freeW(ground);
            ground = get_E_freezeC_Xice(ground);
            
            ground = calculate_hydraulicConductivity(ground);

            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.STATVAR.runoff = 0;
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;

        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
           %do nothing 
        end
        
        function ground = get_boundary_condition_l(ground, tile)
           %do nothing 
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
           %do nothing 
        end
        
        function timestep = get_timestep(ground, tile) 
           timestep = 1e9;
        end
        
        function ground = advance_prognostic(ground, tile)           
           %do nothing 
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
           %do nothing 
        end
       
        function ground = compute_diagnostic(ground, tile)
           %do nothing 
        end
        
        function ground = check_trigger(ground, tile)
           %do nothing 
        end

        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end

    end
    
end
