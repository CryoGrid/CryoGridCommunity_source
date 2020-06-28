classdef LAKE_simple_seb < SEB & HEAT_CONDUCTION & LAKE & INITIALIZE

    
    methods

        
       function self = LAKE_simple_seb(index, pprovider, cprovider, forcing)  
            self@INITIALIZE(index, pprovider, cprovider, forcing);
        end

         function ground = initialize_from_LAKE_previous_season(ground, LAKE_simple_unfrozen)
            ground.PARA = LAKE_simple_unfrozen.PARA;
            ground.CONST = LAKE_simple_unfrozen.CONST;
            ground.STATVAR = LAKE_simple_unfrozen.STATVAR;
            %change the STATVAR
            ground = create_stratigraphy_from_STATVAR(ground);
            ground.PARA.next_season_lake_class = class(LAKE_simple_unfrozen);
         end
        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];
            ground.PARA.epsilon = [];
            ground.PARA.airT_height = []; %measurement height [m]
            ground.PARA.z0 = []; %roughness length [m]
            ground.PARA.rs = [];
            
            ground.PARA.area =[]; %initial area of the realization [m2]
            
            ground.PARA.heatFlux_lb = [];
            
            ground.PARA.dt_max = [] ; %[sec]
            ground.PARA.dE_max = []; %[J/m3]
            ground.PARA.next_season_lake_class = []; %class called by creation/annihilation IA class
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
            
        end
        
        
        function ground = create_stratigraphy_from_STATVAR(ground)  %merges all cell in one
            layerThickSum = ground.STATVAR.layerThick; 
            ground.STATVAR.layerThick = ground.STATVAR.layerThick_store; %reassign the old grid, MISSING: change grid if necessary 
            
            ground.STATVAR.waterIce = ground.STATVAR.waterIce ./ layerThickSum .* ground.STATVAR.layerThick;
            ground.STATVAR.mineral = ground.STATVAR.mineral ./ layerThickSum .* ground.STATVAR.layerThick;
            ground.STATVAR.organic = ground.STATVAR.organic ./ layerThickSum .* ground.STATVAR.layerThick;
            
            ground.STATVAR.air = ground.STATVAR.air ./ layerThickSum .* ground.STATVAR.layerThick;
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
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
      
            ground = surface_energy_balance(ground, forcing);
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
             ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
             ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_derivatives_prognostic(ground)
            ground = get_derivative_energy(ground);

        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
           timestep = get_timestep_heat_coduction(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) 
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, forcing)
            
            ground = move_ice_up(ground);
            ground = get_T_water_freeW(ground);
            ground = conductivity(ground);
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
        end
        
        %shifts to unfrozen LAKE class
        function ground = check_trigger(ground, forcing) 
            if sum(double(ground.STATVAR.energy<0),1)==0  %all cells unfrozen
                
                ia_create_next_season_lake = get_IA_class(class(ground), ground.PARA.next_season_lake_class); %delivers IA-class that creates and initializes the next season LAKE class
                lake_next_season = create_annihilate(ia_create_next_season_lake, ground);
                
                
                %get the interaction classes from above and below
                lake_next_season.NEXT = ground.NEXT;
                lake_next_season.PREVIOUS = ground.PREVIOUS;
                lake_next_season.PREVIOUS.NEXT = lake_next_season;
                lake_next_season.NEXT.PREVIOUS = lake_next_season;
                
                %ground (CURRENT) still points to NEXT, so CURRENT.NEXT
                %will advance to the next - then it is
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
        end
        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq(ground, forcing);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
        end
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
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
        
        function [ground, S_up] = penetrate_SW_no_transmission(ground, S_down)
            [ground, S_up] = penetrate_SW_no_transmission@SEB(ground, S_down);
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
    end
    
end
