%========================================================================
% CryoGrid VEGETATION class VEGETATION_simpleShading_seb
%
% First attempt of a VEGETATION class, providing simple shading by the
% canopy, see Link & Marks (1999), Garen & Marks (2005), and Bair et al. (2016)
%
% These parameterizations are the same as originally found in
% FORCING_slope_forest_seb_readNc
%
% R. B. Zweigel, July 2021
%========================================================================

classdef VEGETATION_simpleShading_seb < SEB & WATER_FLUXES  & VEGETATION
    
    methods
        
        function canopy = provide_PARA(canopy)
            canopy.PARA.canopy_albedo = [];
            canopy.PARA.canopy_temperature = [];
            canopy.PARA.LAI = [];
            canopy.PARA.leaf_cp_areal = [];
            canopy.PARA.canopy_transmissivity = [];
            canopy.PARA.canopy_height = [];
            canopy.PARA.canopy_extinction_coefficient = [];
            canopy.PARA.canopy_emissivity = [];
            canopy.PARA.fractional_canopy_cover = [];
            canopy.PARA.dT_max = [];
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
        end
        
        function canopy = provide_CONST(canopy)
            canopy.CONST.Tmfw = []; % freezing temperature of free water [K]
            canopy.CONST.sigma = []; %Stefan-Boltzmann constant
        end
        
        function canopy = finalize_init(canopy, tile)
            canopy.STATVAR.area = tile.PARA.area + canopy.STATVAR.layerThick .* 0;
            canopy.TEMP.d_energy = canopy.STATVAR.area.*0;
            
            if canopy.PARA.canopy_temperature == 1
                canopy = get_E_simpleVegetation(canopy);
            end
        end
        
        %--------------- time integration ---------------------------------
        %==================================================================
        
        function canopy = get_boundary_condition_u(canopy, tile)
            forcing = tile.FORCING;
            canopy = canopy_energy_balance(canopy,forcing);
            
            canopy = canopy_water_balance(canopy,forcing);
        end
        
        function [canopy, L_up] = penetrate_LW(canopy, L_down)  %mandatory function when used with class that features LW penetration
            [canopy, L_up] = penetrate_LW_simpleShading(canopy, L_down);
        end
        
        function [canopy, S_up] = penetrate_SW(canopy, S_down)  %mandatory function when used with class that features SW penetration
            [canopy, S_up] = penetrate_SW_simpleShading(canopy, S_down);
        end
        
        function canopy = get_boundary_condition_l(canopy, tile)
            
        end
        
        function canopy = get_derivatives_prognostic(canopy, tile)
            
        end
        
        function timestep = get_timestep(canopy, tile)
            if canopy.PARA.canopy_temperature == 1
                timestep = get_timestep_canopy_T(canopy);
            else
                timestep = 1e8;
            end
        end
        
        function canopy = advance_prognostic(canopy, tile)
            if canopy.PARA.canopy_temperature == 1
                timestep = tile.timestep;
                canopy.STATVAR.energy = canopy.STATVAR.energy + timestep .* canopy.TEMP.d_energy;
            end
        end
        
        function canopy = compute_diagnostic_first_cell(canopy, tile)
            forcing = tile.FORCING;
            canopy.NEXT = L_star(canopy.NEXT, forcing);
        end
        
        function canopy = compute_diagnostic(canopy, tile)
            if canopy.PARA.canopy_temperature == 1
                canopy = get_T_simpleVegetatation(canopy);
                canopy.TEMP.d_energy = canopy.STATVAR.energy.*0;
            end
        end
        
        function canopy = check_trigger(canopy, tile)
            
        end
        
        %---------------non-mandatory functions----------------------------
        %==================================================================
        
        function canopy = canopy_energy_balance(canopy,forcing)
            % 1. Longwave penetration
            if canopy.PARA.canopy_temperature == 0 % Assume canopy temperature = Tair
                canopy.STATVAR.T = forcing.TEMP.Tair; % RBZ: original code, changed when canopy temperature was introduced
            end
            
            [canopy, L_up] = penetrate_LW(canopy, forcing.TEMP.Lin .* canopy.STATVAR.area(1));
            canopy.STATVAR.Lout = sum(L_up) ./ canopy.STATVAR.area(1);
            
            % 2. Shortwave penetration
            canopy.TEMP.sun_angle = forcing.TEMP.sunElevation * double(forcing.TEMP.sunElevation > 0);
            canopy.TEMP.Sin_dir_fraction = forcing.TEMP.Sin_dir ./ forcing.TEMP.Sin;
            canopy.TEMP.Sin_dir_fraction(isnan(canopy.TEMP.Sin_dir_fraction)) = 0;
            [canopy, S_up] = penetrate_SW(canopy, forcing.TEMP.Sin .* canopy.STATVAR.area(1));
            canopy.STATVAR.Sout = sum(S_up) ./ canopy.STATVAR.area(1);
            
            canopy.STATVAR.Sin = forcing.TEMP.Sin;
            canopy.STATVAR.Lin = forcing.TEMP.Lin;
        end
        
        function canopy = canopy_water_balance(canopy,forcing)
            % route all rainwater to next class
            canopy.NEXT = get_boundary_condition_u_water(canopy.NEXT, forcing);
        end
        %-----------------inherited Tier 1 functions ----------------------
        %==================================================================
        
        function [canopy, L_up] = penetrate_LW_simpleShading(canopy,L_down)
            [canopy, L_up] = penetrate_LW_simpleShading@SEB(canopy,L_down);
        end
        
        function [canopy, S_up] = penetrate_SW_simpleShading(canopy,S_down)
            [canopy, S_up] = penetrate_SW_simpleShading@SEB(canopy,S_down);
        end
        
        function timestep = get_timestep_canopy_T(canopy)
            timestep = get_timestep_canopy_T@VEGETATION(canopy);
        end
        
    end
end