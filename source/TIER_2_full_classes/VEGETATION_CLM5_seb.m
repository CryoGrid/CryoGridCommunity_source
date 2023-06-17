%========================================================================
% CryoGrid VEGETATION class VEGETATION_CLM5_seb
%
% Surface energy and water balance of vegetated surfaces as in CLM5
%
% R. B. Zweigel, September 2021
%========================================================================

classdef VEGETATION_CLM5_seb < SEB & SEB_VEGETATION & WATER_FLUXES & VEGETATION
    properties
        GROUND
        IA_GROUND
    end
    
    methods
        
        function canopy = provide_PARA(canopy)
            canopy.PARA.t_leafsprout = [];
            canopy.PARA.t_leaffall = [];
            canopy.PARA.kv = [];
            canopy.PARA.D_bh = [];
            canopy.PARA.N_tree = [];
            canopy.PARA.rho_wood = [];
            canopy.PARA.SLA = [];
            canopy.PARA.f_carbon = [];
            canopy.PARA.f_water = [];
            canopy.PARA.nir_fraction = [];
            canopy.PARA.LAI = [];
            canopy.PARA.SAI = [];
            canopy.PARA.Khi_L = [];
            canopy.PARA.alpha_leaf_vis = [];
            canopy.PARA.alpha_stem_vis = [];
            canopy.PARA.tau_leaf_vis = [];
            canopy.PARA.tau_stem_vis = [];
            canopy.PARA.alpha_leaf_nir = [];
            canopy.PARA.alpha_stem_nir = [];
            canopy.PARA.tau_leaf_nir = [];
            canopy.PARA.tau_stem_nir = [];
            canopy.PARA.dT_max = [];
            canopy.PARA.dt_max = [];
            canopy.PARA.Cv = [];
            canopy.PARA.d_leaf = [];
            canopy.PARA.Cs_dense = [];
            canopy.PARA.R_z0 = [];
            canopy.PARA.R_d = [];
            canopy.PARA.Dmax = [];
            canopy.PARA.Wmax = []; % Assume one water holding capacity for snow and water for now
            canopy.PARA.beta_root = [];
            canopy.PARA.C_leaf_max = [];
            canopy.PARA.k_shelter = [];
            canopy.PARA.psi_wilt = [];
            canopy.PARA.zeta_m = [];
            canopy.PARA.zeta_h = [];
        end
        
        function canopy = provide_STATVAR(canopy)
            canopy.STATVAR.upperPos = [];  % upper surface elevation [m]
            canopy.STATVAR.lowerPos = [];  % lower surface elevation [m]
            
            canopy.STATVAR.layerThick = []; % thickness of grid cells [m]
            canopy.STATVAR.area = [];     % grid cell area [m2]
            
            canopy.STATVAR.waterIce = [];  % total volume of water plus ice in a grid cell [m3]
            canopy.STATVAR.mineral = [];   % total volume of minerals [m3]
            canopy.STATVAR.organic = [];   % total volume of organics [m3]
            canopy.STATVAR.energy = [];    % total internal energy [J]
            
            canopy.STATVAR.T = [];  % temperature [degree C]
            canopy.STATVAR.water = [];  % total volume of water [m3]
            canopy.STATVAR.ice = [];  %total volume of ice [m3]
            canopy.STATVAR.q_s = []; % canopy specific humidity [kg/kg]
            
            canopy.STATVAR.Lstar = [];  %Obukhov length [m]
            canopy.STATVAR.u_star = [];
            canopy.STATVAR.Qh = [];     %sensible heat flux [W/m2]
            canopy.STATVAR.Qe = [];     % latent heat flux [W/m2]
            
            canopy.STATVAR.z0 = []; % roughness length [m]
            canopy.STATVAR.d = []; % displacement height
            canopy.STATVAR.f_wet = []; % wetted fraction of canopy
            canopy.STATVAR.f_dry = []; % dry and transpiring fraction of canopy
%             canopy.STATVAR.f_snow = []; % snow-covered fraction of canopy
            
        end
        
        function canopy = provide_CONST(canopy)
            canopy.CONST.Tmfw = []; % freezing temperature of free water [K]
            canopy.CONST.sigma = []; % Stefan-Boltzmann constant
            canopy.CONST.R_spec = []; % gas constant for dry air [J/K kg]
            canopy.CONST.Rw = []; %gas constant for water vapor [J/K kg]
            canopy.CONST.ypsilon = []; % Kinematic viscosity of air [kg/(m sec)]
            canopy.CONST.cp = []; % Specific heat capacity of air
            canopy.CONST.kappa = []; % von Karman constant
            canopy.CONST.Gamma_dry = []; %negative of dry adiabatic lapse rate
            canopy.CONST.g = []; % gravitational constant
            canopy.CONST.L_s = []; % Latent heat of sublimation [J/kg]
            canopy.CONST.L_f = []; % volumetric heat of subliamtion [J/m3]
            canopy.CONST.tau = []; % tortuosity of van Gneuchten model [-]
            canopy.CONST.c_w = []; % vol. heat capacity of water [J/m3/K]
            canopy.CONST.c_i = []; % vol. heat capacity of ice [J/m3/K]
            canopy.CONST.rho_w = []; % density of water [kg/m3]
        end
        
        function canopy = convert_units(canopy, tile)

        end
        
        function canopy = finalize_init(canopy, tile)
            canopy.STATVAR.area = tile.PARA.area + canopy.STATVAR.layerThick .* 0;
            canopy.STATVAR.SAI = canopy.PARA.SAI;
            
            % Construct canopy -  add/remove LAI, calc. energy, heat capacity, z0, etc. accordingly
            doy_start = tile.FORCING.PARA.start_time - datenum(year(tile.FORCING.PARA.start_time),1,1); % DayOfYear
            if doy_start >= canopy.PARA.t_leafsprout && doy_start < canopy.PARA.t_leaffall
                canopy = add_canopy(canopy);
            else
                canopy = remove_canopy(canopy);
            end

            canopy.PARA.airT_height = tile.FORCING.PARA.airT_height; % why tranfer this parameter to ground/canopy, and not use forcing directly?
            
            % assign ground for canopy-soil interactions
            canopy.GROUND = canopy.NEXT;
            canopy.IA_GROUND = IA_VEGETATION_CLM5();
            canopy.IA_GROUND.PREVIOUS = canopy;
            canopy.IA_GROUND.NEXT = canopy.NEXT;
            canopy.IA_GROUND = distribute_roots_CLM5(canopy.IA_GROUND);

%             canopy.TEMP.beta_t = 0;
            canopy.STATVAR.Lstar = -100;
            canopy.STATVAR.u_star = 0.1;
            canopy.STATVAR.q_s = .622*satPresWater(canopy,canopy.STATVAR.T+canopy.CONST.Tmfw)/101300; 
            canopy.STATVAR.Ts = canopy.STATVAR.T;
            canopy.STATVAR.Qh = 0;
            canopy.STATVAR.Qe = 0;
            canopy.STATVAR.f_wet = 0;
            canopy.STATVAR.f_dry = 0;
            
            canopy.TEMP.excess_water = 0;
            canopy.TEMP.excess_water_energy = 0;
            canopy.TEMP.d_energy = canopy.STATVAR.area.*0;
            canopy.TEMP.d_water = canopy.STATVAR.area.*0;
            canopy.TEMP.d_water_energy = canopy.STATVAR.area.*0;
            canopy.TEMP.d_water_ET = canopy.STATVAR.area.*0;
            canopy.TEMP.d_water_ET_energy = canopy.STATVAR.area.*0;
        end
        
        %--------------- time integration ---------------------------------
        %==================================================================
        
        function canopy = get_boundary_condition_u(canopy, tile)
            canopy = canopy_energy_balance(canopy,tile); % change to tile
            canopy = canopy_water_balance(canopy,tile);
        end
        
        function [canopy, L_up] = penetrate_LW(canopy, L_down)  %mandatory function when used with class that features LW penetration
            [canopy, L_up] = penetrate_LW_CLM5(canopy, L_down);
        end
        
        function [canopy, S_up] = penetrate_SW(canopy, S_down)  %mandatory function when used with class that features SW penetration
            [canopy, S_up] = penetrate_SW_CLM5(canopy, S_down);
        end
        
        function [canopy, S_up] = penetrate_SW_PARENT(canopy, S_down)  %mandatory function when used with class that features SW penetration
            [canopy, S_up] = penetrate_SW_CLM5(canopy, S_down);
        end
        
        function canopy = get_boundary_condition_l(canopy, tile)
            % MUST be used with soil class below - No lower boundary condition
        end
        
        function canopy = get_derivatives_prognostic(canopy, tile)
            canopy = get_derivate_water_canopy(canopy);
        end
        
        function timestep = get_timestep(canopy, tile)
            timestep = get_timestep_canopy_T(canopy);
            timestep = min(timestep, get_timestep_water_vegetation(canopy));
        end
        
        function canopy = advance_prognostic(canopy, tile)
            timestep = tile.timestep;
            canopy.STATVAR.energy = canopy.STATVAR.energy + timestep .* canopy.TEMP.d_energy + timestep.*canopy.TEMP.d_water_energy;
            canopy.STATVAR.waterIce = canopy.STATVAR.waterIce + timestep .* canopy.TEMP.d_water;
            canopy.STATVAR.waterIce = max(0,canopy.STATVAR.waterIce); % Avoid rounding errors
        end
        
        function canopy = compute_diagnostic_first_cell(canopy, tile)
            canopy = L_star_canopy(canopy, tile);
        end
        
        function canopy = compute_diagnostic(canopy, tile)
            canopy = get_T_water_vegetation(canopy);
            canopy = get_z0_d_vegetation(canopy);
            canopy.IA_NEXT = canopy_drip(canopy.IA_NEXT, tile);
            
            canopy.TEMP.d_energy = canopy.STATVAR.energy.*0;
            canopy.TEMP.d_water = canopy.STATVAR.energy.*0;
            canopy.TEMP.d_water_ET = canopy.STATVAR.energy.*0;
            canopy.TEMP.d_water_energy = canopy.STATVAR.energy.*0;
            canopy.TEMP.d_water_ET_energy = canopy.STATVAR.energy.*0;
        end
        
        function canopy = check_trigger(canopy, tile)
            
            doy = tile.t - datenum(year(tile.t),1,1);
            if doy >= canopy.PARA.t_leafsprout && doy < canopy.PARA.t_leaffall && canopy.STATVAR.LAI == 0
                add_canopy(canopy);
            elseif doy >= canopy.PARA.t_leaffall && canopy.STATVAR.LAI ~= 0
                remove_canopy(canopy);
            end

        end
        
        %---------------non-mandatory functions----------------------------
        %==================================================================
        
        function canopy = canopy_energy_balance(canopy,tile)
            forcing = tile.FORCING;

            % 1. Longwave penetration
            [canopy, L_up] = penetrate_LW(canopy, forcing.TEMP.Lin .* canopy.STATVAR.area(1));
%             canopy.STATVAR.Lout = L_up./canopy.STATVAR.area(1);
            
            % 2. Shortwave penetration
            canopy.TEMP.sun_angle = forcing.TEMP.sunElevation * double(forcing.TEMP.sunElevation > 0);
            canopy.TEMP.Sin_dir_fraction = forcing.TEMP.Sin_dir ./ forcing.TEMP.Sin;
            canopy.TEMP.Sin_dir_fraction(isnan(canopy.TEMP.Sin_dir_fraction)) = 0;
            [canopy, S_up] = penetrate_SW(canopy, forcing.TEMP.Sin .* canopy.STATVAR.area(1));
%             canopy.STATVAR.Sout = sum(S_up) ./ canopy.STATVAR.area(1);
%             
%             canopy.STATVAR.Sin = forcing.TEMP.Sin;
%             canopy.STATVAR.Lin = forcing.TEMP.Lin;
            
            % 3. Sensible and latent heat
            canopy = canopy_resistances_CLM5_Stewart(canopy, tile);
            canopy = Q_h_CLM5(canopy, tile);
            canopy = Q_e_CLM5_Stewart(canopy, tile);
            
            canopy.TEMP.d_energy = canopy.TEMP.d_energy - (canopy.STATVAR.Qh + canopy.STATVAR.Qe).*canopy.STATVAR.area;
            canopy.TEMP.d_water_ET = canopy.TEMP.d_water_ET - (canopy.STATVAR.evap + canopy.STATVAR.sublim).*canopy.STATVAR.area; % transpired water is removed from soil, not canopy
            canopy.TEMP.d_water_ET_energy = canopy.TEMP.d_water_ET_energy - (canopy.STATVAR.evap_energy + canopy.STATVAR.sublim_energy).*canopy.STATVAR.area;
        end
        
        function canopy = canopy_water_balance(canopy,tile)
            canopy = get_boundary_condition_u_water_canopy(canopy, tile);
            canopy.IA_GROUND = get_water_transpiration(canopy.IA_GROUND);
        end
    
        
        
        %-------------param file generation-----
        function canopy = param_file_info(canopy)
            canopy = provide_PARA(canopy);

            canopy.PARA.STATVAR = [];
            canopy.PARA.class_category = 'VEGETATION';
            canopy.PARA.options = [];

            canopy.PARA.comment.t_leafsprout = {'DayOfYear when leaves emerge (135 = 15. May)'};
            canopy.PARA.default_value.t_leafsprout = {'152'};

            canopy.PARA.comment.t_leaffall = {'DayOfYear when leaves fall (274 = 1. Oct.)'};
            canopy.PARA.default_value.t_leaffall = {'274'};

            canopy.PARA.comment.kv = {'parameter to adjust for non-cylinderness of trees'};
            canopy.PARA.default_value.kv = {0.5};

            canopy.PARA.comment.D_bh = {'mean breast-height diameter of trees'};
            canopy.PARA.default_value.D_bh = {0.28};

            canopy.PARA.comment.N_tree = {'stand density (number of trees per area)'};
            canopy.PARA.default_value.N_tree = {0.15};

            canopy.PARA.comment.rho_wood = {'Density of dry wood'};
            canopy.PARA.default_value.rho_wood = {500};

            canopy.PARA.comment.SLA = {'Specific leaf area (leaf area per unit mass of carbon)'};
            canopy.PARA.default_value.SLA = {23};
            
            canopy.PARA.comment.f_carbon = {'mass of carbon per leaf dry mass'};
            canopy.PARA.default_value.f_carbon = {0.5};

            canopy.PARA.comment.f_water = {'mass of water per leaf total mass'};
            canopy.PARA.default_value.f_water = {0.45};

            canopy.PARA.comment.nir_fraction = {'fraction of Sin that is in the NIR spectrum'};
            canopy.PARA.default_value.nir_fraction = {0.3};

            canopy.PARA.comment.LAI = {'Leaf area index'};
            canopy.PARA.default_value.LAI = {1.5};

            canopy.PARA.comment.SAI = {'Stem area index'};
            canopy.PARA.default_value.SAI = {1};

            canopy.PARA.comment.Khi_L = {'departure of leaf angles from a random distribution'};
            canopy.PARA.default_value.Khi_L = {0.01};
           
            canopy.PARA.comment.alpha_leaf_vis = {'leaf reflectence in the VIS'};
            canopy.PARA.default_value.alpha_leaf_vis = {0.07};

            canopy.PARA.comment.alpha_stem_vis = {'stem reflectence in the VIS'};
            canopy.PARA.default_value.alpha_stem_vis = {0.16};

            canopy.PARA.comment.tau_leaf_vis = {'leaf transmittance in the VIS'};
            canopy.PARA.default_value.tau_leaf_vis = {0.05};

            canopy.PARA.comment.tau_stem_vis = {'stem transmittance in the VIS'};
            canopy.PARA.default_value.tau_stem_vis = {0.001};

            canopy.PARA.comment.alpha_leaf_nir = {'leaf reflectence in the NIR'};
            canopy.PARA.default_value.alpha_leaf_nir = {0.35};

            canopy.PARA.comment.alpha_stem_nir = {'stem reflectence in the NIR'};
            canopy.PARA.default_value.alpha_stem_nir = {0.39};

            canopy.PARA.comment.tau_leaf_nir = {'leaf transmittance in the NIR'};
            canopy.PARA.default_value.tau_leaf_nir = {0.1};

            canopy.PARA.comment.tau_stem_nir = {'stem transmittance in the NIR'};
            canopy.PARA.default_value.tau_stem_nir = {0.001};

            canopy.PARA.comment.dT_max = {'max temperature change per timestep [K]'};
            canopy.PARA.default_value.dT_max = {1};

            canopy.PARA.comment.dt_max = {'maximum timestep [sec]'};
            canopy.PARA.default_value.dt_max = {3600};

            canopy.PARA.comment.Cv = {'turbulent transfer coefficient between canopy surface and canopy air'};
            canopy.PARA.default_value.Cv = {0.01};

            canopy.PARA.comment.d_leaf = {'characteristic dimension of leaves in direction of wind flow'};
            canopy.PARA.default_value.d_leaf = {0.04};
            
            canopy.PARA.comment.Cs_dense = {'dense canopy turbulent transfer coefficient (Dickinson et al. 1993)'};
            canopy.PARA.default_value.Cs_dense = {0.004};

            canopy.PARA.comment.R_z0 = {'ratio of momentum roughness length'};
            canopy.PARA.default_value.R_z0 = {0.055};

            canopy.PARA.comment.R_d = {'ratio of displacement height to canopy height'};
            canopy.PARA.default_value.R_d = {0.67};
            
            canopy.PARA.comment.Dmax = {'Maximum dry soil layer thickness, 15 mm'};
            canopy.PARA.default_value.Dmax = {0.015};

            canopy.PARA.comment.Wmax = {'water holding capacity of canopy per unit area (L+S)'}; 
            canopy.PARA.default_value.Wmax = {0.0001}; 
            
            canopy.PARA.comment.beta_root = {'root distribution parameter (CLM5)'};
            canopy.PARA.default_value.beta_root = {0.943};

            canopy.PARA.comment.C_leaf_max = {'maximum leaf conductance'};
            canopy.PARA.default_value.C_leaf_max = {0.01};

            canopy.PARA.comment.k_shelter = {'shelter factor, between 0.5 - 1 (Carlson, 1991)'};
            canopy.PARA.default_value.k_shelter = {0.5};

            canopy.PARA.comment.psi_wilt = {'water potential at permanent wilting point (-15 bar)'};
            canopy.PARA.default_value.psi_wilt = {-150};

            canopy.PARA.comment.zeta_m = {'threshold for very unstable atm. Conditions wrt. Heat/vapor fluxes'};
            canopy.PARA.default_value.zeta_m = {-0.465};

            canopy.PARA.comment.zeta_h = {'threshold for very unstable atm. Conditions wrt. mass fluxes'};
            canopy.PARA.default_value.zeta_h = {-1.574};

        end
    end
end