%========================================================================
% CryoGrid VEGETATION class VEGETATION_CLM5_seb
%
% Surface energy and water balance of vegetated surfaces as in CLM5
%
% Radiation added in September 2021
% Sensible heat exchange in October 2021
% R. B. Zweigel, September 2021
%========================================================================

classdef VEGETATION_CLM5_seb < SEB & WATER_FLUXES & VEGETATION
    
    methods
        
        function canopy = provide_PARA(canopy)
            canopy.PARA.canopy_emissivity = [];
            canopy.PARA.leaf_cp_areal = [];
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
            canopy.PARA.Cv = [];
            canopy.PARA.z_top = [];
            canopy.PARA.d_leaf = [];
            canopy.PARA.Cs_dense = [];
            canopy.PARA.R_z0 = [];
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
            
            canopy.STATVAR.Lstar = [];  %Obukhov length [m]
            canopy.STATVAR.u_star = [];
            canopy.STATVAR.Qh = [];     %sensible heat flux [W/m2]
            canopy.STATVAR.Qe = [];     % latent heat flux [W/m2]
        end
        
        function canopy = provide_CONST(canopy)
            canopy.CONST.Tmfw = []; % freezing temperature of free water [K]
            canopy.CONST.sigma = []; % Stefan-Boltzmann constant
            canopy.CONST.R_spec = []; % gas constant for dry air [J/K kg]
            canopy.CONST.ypsilon = []; % [kg/(m sec)]
            canopy.CONST.cp = []; % Specific heat capacity of air
            canopy.CONST.kappa = []; % von Karman constant
            canopy.CONST.Gamma_dry = []; %negative of dry adiabatic lapse rate
            canopy.CONST.g = []; % gravitational constant
            canopy.CONST.L_s = []; % Latent heat of sublimation
        end
        
        function canopy = finalize_init(canopy, tile)
            canopy.STATVAR.area = tile.PARA.area + canopy.STATVAR.layerThick .* 0;
            canopy.TEMP.d_energy = canopy.STATVAR.area.*0;
            canopy = get_E_simpleVegetation(canopy);
            
            canopy.PARA.airT_height = tile.FORCING.PARA.airT_height; % why tranfer this parameter to ground/canopy, and not use forcing directly?
            canopy.PARA.z0 = 0.1955; % z0_vegetation(canopy); WHY CANT I CALL THIS FUNCTION??!!
            
            canopy.STATVAR.Lstar = -100;
            canopy.STATVAR.u_star = 0.1;
            canopy.STATVAR.Qh = 0;
            canopy.STATVAR.Qe = 0;
        end
        
        %--------------- time integration ---------------------------------
        %==================================================================
        
        function canopy = get_boundary_condition_u(canopy, tile)
            forcing = tile.FORCING;
            canopy = canopy_energy_balance(canopy,forcing);
            
            canopy = canopy_water_balance(canopy,forcing);
        end
        
        function [canopy, L_up] = penetrate_LW(canopy, L_down)  %mandatory function when used with class that features LW penetration
            [canopy, L_up] = penetrate_LW_CLM5(canopy, L_down);
        end
        
        function [canopy, S_up] = penetrate_SW(canopy, S_down)  %mandatory function when used with class that features SW penetration
            [canopy, S_up] = penetrate_SW_CLM5(canopy, S_down);
        end
        
        function canopy = get_boundary_condition_l(canopy, tile)
            
        end
        
        function canopy = get_derivatives_prognostic(canopy, tile)
            
        end
        
        function timestep = get_timestep(canopy, tile)
            timestep = get_timestep_canopy_T(canopy);
        end
        
        function canopy = advance_prognostic(canopy, tile)
            timestep = tile.timestep;
            canopy.STATVAR.energy = canopy.STATVAR.energy + timestep .* canopy.TEMP.d_energy;
        end
        
        function canopy = compute_diagnostic_first_cell(canopy, tile)
            forcing = tile.FORCING;
            canopy = L_star(canopy, forcing);
        end
        
        function canopy = compute_diagnostic(canopy, tile)
            canopy = get_T_simpleVegetatation(canopy);
            canopy.TEMP.d_energy = canopy.STATVAR.energy.*0;
        end
        
        function canopy = check_trigger(canopy, tile)
            
        end
        
        %---------------non-mandatory functions----------------------------
        %==================================================================
        
        function z0v = z0_vegetation(canopy) % roughness length of vegetated surface (CLM5)
            % if L changes due to leaf/needle fall, z0v should also be
            % recalculated!!
            L = canopy.PARA.LAI; % Leaf area index
            S = canopy.PARA.SAI; % Stem area index
            z0g = canopy.NEXT.PARA.z0; % Roughness lenght of ground
            R_z0 = canopy.PARA.R_z0; % Ratio of momentum roughness length to canopy height
            z_top = sum(canopy.STATVAR.layerThick); % canopy height
            
            V = ( 1-exp(-1.*min(L+S,2)) ./ (1-exp(-2)) ); % Eq. 5.127
            z0v = exp( V.*log(z_top.*R_z0) + (1-V).*log(z0g) ); % Eq. 5.125
        end
        
        function canopy = canopy_energy_balance(canopy,forcing)
            % 1. Longwave penetration
            [canopy, L_up] = penetrate_LW(canopy, forcing.TEMP.Lin .* canopy.STATVAR.area(1));
            canopy.STATVAR.Lout = L_up./canopy.STATVAR.area(1);
            
            % 2. Shortwave penetration
            canopy.TEMP.sun_angle = forcing.TEMP.sunElevation * double(forcing.TEMP.sunElevation > 0);
            canopy.TEMP.Sin_dir_fraction = forcing.TEMP.Sin_dir ./ forcing.TEMP.Sin;
            canopy.TEMP.Sin_dir_fraction(isnan(canopy.TEMP.Sin_dir_fraction)) = 0;
            [canopy, S_up] = penetrate_SW(canopy, forcing.TEMP.Sin .* canopy.STATVAR.area(1));
            canopy.STATVAR.Sout = sum(S_up) ./ canopy.STATVAR.area(1);
            
            canopy.STATVAR.Sin = forcing.TEMP.Sin;
            canopy.STATVAR.Lin = forcing.TEMP.Lin;
            
            % 3. Sensible heat
            canopy.STATVAR.Qh = Q_h_CLM5(canopy,forcing);
            
            canopy.TEMP.d_energy = canopy.TEMP.d_energy - canopy.STATVAR.Qh;
        end
        
        function canopy = canopy_water_balance(canopy,forcing)
            % route all rainwater to next class
            canopy.NEXT = get_boundary_condition_u_water(canopy.NEXT, forcing);
        end
        %-----------------inherited Tier 1 functions ----------------------
        %==================================================================
        
        function [canopy, L_up] = penetrate_LW_CLM5(canopy,L_down)
            [canopy, L_up] = penetrate_LW_CLM5@SEB(canopy,L_down);
        end
        
        function [canopy, S_up] = penetrate_SW_CLM5(canopy,S_down)
            [canopy, S_up] = penetrate_SW_CLM5@SEB(canopy,S_down);
        end
        
        function flux = Q_h_CLM5(canopy,forcing)
            flux = Q_h_CLM5@SEB(canopy,forcing);
        end
        
        function timestep = get_timestep_canopy_T(canopy)
            timestep = get_timestep_canopy_T@VEGETATION(canopy);
        end
        
    end
end