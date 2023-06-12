%========================================================================
% CryoGrid GROUND class LAKE_simple_unfrozen_bucketW_seb
% water body with heat conduction, bucket water scheme with free water evaporation, free water freeze curve, surface
% energy balance, penetration of SW radiation (single band)
% representation of unfrozen water body, works in concert with 
% LAKE_simple_bucketW_seb for frozen water body
% represents water body as single cell, assuming perfect mixing at all times
% this class can only be created by a trigger in LAKE_simple_bucketW_seb
% S. Westermann, October 2020
%========================================================================

classdef LAKE_simple_unfrozen_bucketW_seb < SEB & HEAT_CONDUCTION & WATER_FLUXES & WATER_FLUXES_LATERAL &  HEAT_FLUXES_LATERAL 

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        %initializes class when switching from frozen to unfrozen conditions
        function ground = initialize_from_LAKE_previous_season(ground, LAKE_simple_frozen)
            ground.PARA = LAKE_simple_frozen.PARA;
            ground.PARA.rs = 0;
            ground.CONST = LAKE_simple_frozen.CONST;
            ground.STATVAR = LAKE_simple_frozen.STATVAR;
            %change the STATVAR
            ground = merge_STATVAR(ground);
            ground.PARA.next_season_lake_class = class(LAKE_simple_frozen);
        end
        
        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];  %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; %roughness length [m]
            ground.PARA.SW_extinction = []; %e-folding constant of SW extinction [1/m] 
            %spectrally resolved albedo and SW_extinction are possible,
            %must be row array - NOT TESTED!
            
            ground.PARA.dt_max = [];  %maximum possible timestep [sec]
            ground.PARA.dE_max = [];  %maximum possible energy change per timestep [J/m3]
            
            ground.PARA.next_season_lake_class = [];  %LAKE class that is called by check_trigger, in this case unfozen LAKE class
            ground.PARA.threshold_water = []; %lake depth below which a trigger is called. LAKE is generally removed, depending on GROUND class below 
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = []; % total volume of minerals [m3]
            ground.STATVAR.organic = []; % total volume of organics [m3]
            ground.STATVAR.energy = [];  % total internal energy [J]
            
            ground.STATVAR.T = [];  % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.ice = []; %total volume of ice [m3]
            ground.STATVAR.air = [];  % total volume of air [m3] - NOT USED
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            
            ground.STATVAR.Lstar = []; %Obukhov length [m]
            ground.STATVAR.Qh = []; %sensible heat flux [W/m2]
            ground.STATVAR.Qe = []; % latent heat flux [W/m2]
        end
    
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = []; % volumetric latent heat of fusion, freezing
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
            
            ground.CONST.cp = []; %specific heat capacity at constant pressure of air
            ground.CONST.g = []; % gravitational acceleration Earth surface
            
            ground.CONST.rho_w = []; % water density
            ground.CONST.rho_i = []; %ice density
        end
        
        
        function ground = merge_STATVAR(ground)  %merges all cell in one 
            ground.STATVAR.layerThick = sum(ground.STATVAR.layerThick,1);
            
            ground.STATVAR.waterIce = sum(ground.STATVAR.waterIce,1);
            ground.STATVAR.mineral = sum(ground.STATVAR.mineral,1);
            ground.STATVAR.organic = sum(ground.STATVAR.organic,1);
            ground.STATVAR.energy = sum(ground.STATVAR.energy,1);
            ground.STATVAR.air = sum(ground.STATVAR.air,1);
            ground.STATVAR.area = mean(ground.STATVAR.area,1);
            
            ground = get_T_water_freeW(ground);
            ground.STATVAR.thermCond = ground.CONST.k_w;  %thermal conductivity of water
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
        end
        
        
        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_PARA(ground); 
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        

        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            ground = surface_energy_balance(ground, forcing);
            ground = get_boundary_condition_u_water_LAKE(ground, forcing);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_transmission_bulk(ground, S_down);
        end
        
        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_transmission_bulk(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
            
            ground.TEMP.F_lb_water = 0;
            ground.TEMP.F_lb_water_energy =0;
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            %do nothing, single cell only
        end
        
        function timestep = get_timestep(ground, tile)  
           timestep = get_timestep_heat_coduction(ground);
        end
        
        function ground = advance_prognostic(ground, tile) 
            timestep = tile.timestep;
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* (ground.TEMP.d_energy + ground.TEMP.d_water_energy);
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + timestep .* ground.TEMP.d_water;
            ground.STATVAR.layerThick = ground.STATVAR.layerThick + timestep .* ground.TEMP.d_water ./ ground.STATVAR.area;
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)   
            ground = get_T_water_freeW(ground);
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
        end
        
        function ground = check_trigger(ground, tile)
            trigger_yes_no = 0;
            if ground.STATVAR.energy < 0  %freezing has started
                trigger_yes_no = 1;
                ia_create_next_season_lake = get_IA_class(class(ground), ground.PARA.next_season_lake_class); %delivers IA-class that creates and initializes the next season LAKE class
                lake_next_season = create_annihilate(ia_create_next_season_lake, ground);
                
                %get the interaction classes from above and below
                lake_next_season.NEXT = ground.NEXT;
                lake_next_season.PREVIOUS = ground.PREVIOUS;
                lake_next_season.PREVIOUS.NEXT = lake_next_season;
                lake_next_season.NEXT.PREVIOUS = lake_next_season;
                %ground (CURRENT) still points to NEXT, so CURRENT.NEXT
                %will advance to the next - then it should be
                %automatically handled by the garbage collection, since
                %there is no more pointer from any active variable to it
                
                %assemble new INTERACTIONS
                if ~strcmp(class(lake_next_season.PREVIOUS), 'Top')
                    ia_class = get_IA_class(class(lake_next_season.PREVIOUS), class(lake_next_season));
                    lake_next_season.IA_PREVIOUS = ia_class;
                    lake_next_season.IA_PREVIOUS.NEXT = lake_next_season;
                    lake_next_season.IA_PREVIOUS.PREVIOUS = lake_next_season.PREVIOUS;
                    lake_next_season.PREVIOUS.IA_NEXT = ia_class;
                end
                if ~strcmp(class(lake_next_season.NEXT), 'Bottom')
                    ia_class = get_IA_class(class(lake_next_season), class(lake_next_season.NEXT));
                    lake_next_season.IA_NEXT = ia_class;
                    lake_next_season.IA_NEXT.PREVIOUS = lake_next_season;
                    lake_next_season.IA_NEXT.NEXT = lake_next_season.NEXT;
                    lake_next_season.NEXT.IA_PREVIOUS = ia_class;
                end
            end
            
            if ~trigger_yes_no & (ground.STATVAR.waterIce./ground.STATVAR.area) < ground.PARA.threshold_water
                trigger_yes_no = 1;
                trigger_remove_LAKE(ground.IA_NEXT, tile);
            end
        end
        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            [ground, S_up] = penetrate_SW(ground, forcing.TEMP.Sin .* ground.STATVAR.area(1)); %distribute SW radiation
            
            ground.STATVAR.Sout = S_up ./ ground.STATVAR.area(1);
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq(ground, forcing);
            ground = calculateET_LAKE(ground); %set F_ub_water, F_ub_water_energy
            
            ground.TEMP.F_ub = (forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
        end

        %-----LATERAL------------------
        
         %----LAT_SEEPAGE_FACE----------              
        function ground = lateral_push_remove_water_seepage(ground, lateral)
            ground = lateral_push_remove_water_seepage_lake_unfrozen(ground, lateral);
        end
        
        %----LAT_WATER_RESERVOIR-----------        
        function ground = lateral_push_water_reservoir(ground, lateral)
           ground = lateral_push_water_reservoir_lake_unfrozen(ground, lateral);
        end
        
        %---LAT_OVERLAND_FLOW----------
        function ground = lateral_push_remove_water_overland_flow(ground, lateral)
            ground = lateral_push_water_overland_flow_LAKE(ground, lateral);
        end
        
        %----LAT3D_WATER_UNCONFINED_AQUIFER------------     
        function ground = lateral3D_pull_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_pull_water_unconfined_aquifer_lake_unfrozen(ground, lateral);
        end
        
        function ground = lateral3D_push_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_push_water_unconfined_aquifer_lake_unfrozen(ground, lateral);
        end
        
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground, lateral)
            [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_lake_unfrozen(ground, lateral);
        end
        
        function ground = lateral3D_pull_water_overland_flow(ground, lateral)
            ground = lateral3D_pull_water_overland_flow_LAKE(ground, lateral);
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
        
        function ground = get_boundary_condition_u_water_LAKE(ground, forcing)
             ground = get_boundary_condition_u_water_LAKE@WATER_FLUXES(ground, forcing);
        end
        
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
        
        function ground = calculateET_LAKE(ground)
            ground = calculateET_LAKE@SEB(ground);
        end
        
        function [ground, S_up] = penetrate_SW_transmission_bulk(ground, S_down)
            [ground, S_up] = penetrate_SW_transmission_bulk@SEB(ground, S_down);
        end
            
        function timestep = get_timestep_heat_coduction(ground)
            timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(ground);
        end
        
        function ground = L_star(ground, forcing)
           ground = L_star@SEB(ground, forcing); 
        end
        
        function ground = get_T_water_freeW(ground)
            ground = get_T_water_freeW@HEAT_CONDUCTION(ground);
        end
        
        function ground = regrid_full(ground, extensive_variables)
            ground = regrid_full@REGRID(ground, extensive_variables);
        end
    end
    
end
