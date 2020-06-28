classdef LAKE_simple_bucketW_seb < SEB & HEAT_CONDUCTION & LAKE & WATER_FLUXES & REGRID & HEAT_FLUXES_LATERAL & INITIALIZE

    
    methods

        
       function self = LAKE_simple_bucketW_seb(index, pprovider, cprovider, forcing)  
            self@INITIALIZE(index, pprovider, cprovider, forcing);
        end

%          function ground = initialize_from_LAKE_unfrozen(ground, LAKE_simple_unfrozen)
%             ground.PARA = LAKE_simple_unfrozen.PARA;
%             ground.CONST = LAKE_simple_unfrozen.CONST;
%             ground.STATVAR = LAKE_simple_unfrozen.STATVAR;
%             %change the STATVAR
%             ground = create_stratigraphy_from_STATVAR(ground);
%          end
         
         function ground = initialize_from_LAKE_previous_season(ground, LAKE_simple_unfrozen)
            ground.PARA = LAKE_simple_unfrozen.PARA;
            ground.CONST = LAKE_simple_unfrozen.CONST;
            ground.STATVAR = LAKE_simple_unfrozen.STATVAR;
            %change the STATVAR
            ground = create_stratigraphy_from_STATVAR(ground);
            ground.PARA.next_season_lake_class = class(LAKE_simple_unfrozen);
         end
        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];  %spectrally resolved albedo and SW_extinction are possible, must be row array
            ground.PARA.SW_extinction = [];
            ground.PARA.epsilon = [];
            ground.PARA.airT_height = []; %measurement height [m]
            ground.PARA.z0 = []; %roughness length [m]
            ground.PARA.rs = [];
            
            ground.PARA.area =[]; %initial area of the realization [m2]
            
            ground.PARA.heatFlux_lb = [];
            
            ground.PARA.dt_max = [] ; %[sec]
            ground.PARA.dE_max = []; %[J/m3]
            
            ground.PARA.next_season_lake_class = [];
            ground.PARA.threshold_water = 0.05;
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = [];
            ground.STATVAR.lowerPos = [];
            ground.STATVAR.layerThick = []; % [m]
            
            ground.STATVAR.waterIce = []; % [m]
            ground.STATVAR.mineral = []; % [m]
            ground.STATVAR.organic = []; % [m]
            ground.STATVAR.energy = [];  % [J/m2]
            
            ground.STATVAR.T = [];  % [degree C]
            ground.STATVAR.water = [];  % [m]
            ground.STATVAR.ice = [];
            ground.STATVAR.air = [];  % [m]
            ground.STATVAR.thermCond = [];
            
            ground.STATVAR.Lstar = [];
            ground.STATVAR.Qh = [];
            ground.STATVAR.Qe = [];
        end
    
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = [];
            ground.CONST.c_m =[];
            ground.CONST.c_w = [];
            ground.CONST.c_o = [];
            ground.CONST.c_i = [];
            
            ground.CONST.k_a = [];       %air [Hillel(1982)]
            ground.CONST.k_w = [];        %water [Hillel(1982)]
            ground.CONST.k_i = [];         %ice [Hillel(1982)]
            ground.CONST.k_o = [];        %organic [Hillel(1982)]
            ground.CONST.k_m = [];
            
            ground.CONST.sigma = []; %Stefan-Boltzmann constant
            ground.CONST.kappa = [];
            ground.CONST.L_s = []; %latent heat of vaporization
            
            ground.CONST.cp = [];
            ground.CONST.g = [];
            
            ground.CONST.rho_w = [];
            ground.CONST.rho_i = [];
        end
        
        
        %----mandatory functions---------------
        %----initialization--------------------
        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_PARA(ground); 
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        
        function ground = finalize_init(ground, forcing) %assign all variables, that must be calculated or assigned otherwise for initialization
            ground.PARA.heatFlux_lb = forcing.PARA.heatFlux_lb;
            ground.PARA.airT_height = forcing.PARA.airT_height;
            
            ground.STATVAR.area = ground.PARA.area + ground.STATVAR.T .* 0;
            
            ground = get_E_freeW(ground);
            %if energy is positive, combine all grid cells in one and go to
            %the summer state
            %otherwise abort, this is nonsense to start in winter
            
            
            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
        end
        
        
        function ground = create_stratigraphy_from_STATVAR(ground)  %merges all cell in one

            ground = regrid_full(ground, {'layerThick'; 'waterIce'; 'mineral'; 'organic'});
            
            ground.STATVAR.area = ground.STATVAR.area + ground.STATVAR.layerThick .*0;
            
            energy_total = ground.STATVAR.energy; 
            ground.STATVAR.energy = ground.STATVAR.layerThick .*0;
            %assign energy top-down to separate ice from water
            i=1;
            while energy_total < 0
                ground.STATVAR.energy(i,1) = max(energy_total, -ground.STATVAR.waterIce(i,1).*ground.CONST.L_f);
                energy_total = energy_total - ground.STATVAR.energy(i,1);
                i=i+1;
            end
        
            ground = get_T_water_freeW(ground);
            ground = conductivity(ground);  %thermal conductivity of water
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
      
            ground = surface_energy_balance(ground, forcing);
            ground = get_boundary_condition_u_water_LAKE_frozen(ground, forcing);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_transmission_bulk(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
             ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
             ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;

             ground.TEMP.F_lb_water = 0;
             ground.TEMP.F_lb_water_energy =0;
        end
        
        function ground = get_derivatives_prognostic(ground)
            if size(ground.STATVAR.energy,1)>1
                ground = get_derivative_energy(ground);
            end
%             ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
%             ground.TEMP.d_water(end) = ground.TEMP.d_water(end) + ground.TEMP.F_lb_water;
%             ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
%             ground.TEMP.d_water_energy(end) = ground.TEMP.d_water_energy(end) + ground.TEMP.F_lb_water_energy;
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
           timestep = get_timestep_heat_coduction(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) 
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* (ground.TEMP.d_energy + ground.TEMP.d_water_energy);
            %water - no water added at top yet!!
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + timestep .*  ground.TEMP.d_water;
            ground.STATVAR.layerThick = ground.STATVAR.layerThick + timestep .* ground.TEMP.d_water ./ ground.STATVAR.area;
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, forcing)
            
            ground = move_ice_up(ground);
            ground = stratify(ground);
            ground = regrid_split(ground, {'layerThick'; 'waterIce'; 'mineral'; 'organic'});
            ground = get_T_water_freeW(ground);
            ground = conductivity(ground);
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
        end
        
        
        function ground = check_trigger(ground, forcing)
            trigger_yes_no = 0;
            if sum(double(ground.STATVAR.energy<0),1)==0  %all cells unfrozen, switch to LAKE_unfrozen class
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
            
            if ~trigger_yes_no & sum(ground.STATVAR.waterIce,1) < ground.PARA.threshold_water 
                trigger_yes_no = 1;
                trigger_remove_LAKE(ground.IA_NEXT, forcing);
            end

        end
        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            [ground, S_up] = penetrate_SW(ground, forcing.TEMP.Sin.* ground.STATVAR.area(1)); %distribute SW radiation
            
            ground.STATVAR.Sout = S_up ./ ground.STATVAR.area(1);
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq(ground, forcing);
            %sublimation not yet accounted for, normally handled by the
            %snow
            
            ground.TEMP.F_ub = (forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);  
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
        end
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        %--------------------
        %lateral
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
        
        function ground = stratify(ground)
            ground = stratify@LAKE(ground);
        end
        
        function ground = regrid_full(ground, variable_list)
            ground = regrid_full@REGRID(ground, variable_list);
        end
        
        function ground = regrid_split(ground, variable_list)
            ground = regrid_split@REGRID(ground, variable_list);
        end
        
    end
    
end
