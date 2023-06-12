%========================================================================
% CryoGrid GROUND class GROUND_freeW_ubft_building
% heat conduction, upper boundary temperature forcing
% J. Aga, S. Westermann, January 2023
%========================================================================


classdef GROUND_house_on_poles < HEAT_CONDUCTION &  UB_TEMPERATURE_FORCING

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = provide_PARA(ground)
            ground.PARA.thickness = []; %thickness of the space between building and ground [m]
            ground.PARA.number_of_cells = []; %number of gridcells [-]
            ground.PARA.height_above_ground = []; %height above ground [m]
            ground.PARA.epsilon = []; %surface emissivity [-]
            ground.PARA.thermCond = []; %thermal conductivity [W/mK]
            ground.PARA.heatCapacity = []; %heat capacity [J/m3K]
            ground.PARA.Tr = []; %room temperature in building [C]
            ground.PARA.exchangeAir = [];
            
            %ground.PARA.conductivity_function = [];
            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = [];  % thickness of grid cells [m]
            
            %ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            %ground.STATVAR.mineral = [];  % total volume of minerals [m3]
            %ground.STATVAR.organic = [];  % total volume of ice [m3]
            ground.STATVAR.energy = [];   % total internal energy[J]
            
            ground.STATVAR.T = [];     % temperature [degree C]
            %ground.STATVAR.water = []; % total volume of water [m3]
            %ground.STATVAR.ice = [];   % total volume of ice [m3]
            ground.STATVAR.air = [];   % total volume of air [m3] - NOT USED
            
            ground.STATVAR.thermCond = []; % thermal conductivity [W/mK]
            
        end
    
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = [];  % volumetric latent heat of fusion, freezing
            ground.CONST.c_w = [];  % volumetric heat capacity water
            ground.CONST.c_i = [];  % volumetric heat capacity ice
            ground.CONST.c_o = [];  % volumetric heat capacity organic
            ground.CONST.c_m = [];  % volumetric heat capacity mineral
            
            ground.CONST.k_a = [];  % thermal conductivity air
            ground.CONST.k_w = [];  % thermal conductivity water
            ground.CONST.k_i = [];  % thermal conductivity ice 
            ground.CONST.k_o = [];  % thermal conductivity organic 
            ground.CONST.k_m = [];  % thermal conductivity mineral 
            
            ground.CONST.sigma = []; % Stefan-Boltzmann constant
            ground.CONST.Tmfw = [];
            ground.CONST.kappa = []; % von Karman constant
            ground.CONST.L_s = [];   % latent heat of sublimation, latent heat of evaporation handled in a dedicated function
            
            ground.CONST.cp = []; % specific heat capacity at constant pressure of air
            ground.CONST.g = [];  % gravitational acceleration Earth surface
            
            ground.CONST.rho_w = [];   % water density
            ground.CONST.rho_i = [];   %ice density
        end

        
        function ground = convert_units(ground, tile)

        end


        function ground = finalize_init(ground, tile)
            
            ground.STATVAR.layerThick = zeros(ground.PARA.number_of_cells,1) + ground.PARA.thickness ./ ground.PARA.number_of_cells;

            %if isempty(ground.PARA.conductivity_function) || sum(isnan(ground.PARA.conductivity_function))>0
            %    ground.PARA.conductivity_function = 'conductivity_mixing_squares';
            %end
            ground.STATVAR.thermCond = ground.STATVAR.layerThick .* 0 + ground.PARA.thermCond;
            ground.STATVAR.heatCapacity = ground.STATVAR.layerThick .* 0 + ground.PARA.heatCapacity;

            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.layerThick .* 0;
            
            %ground = get_E_freeW(ground);
            ground.STATVAR.T = ground.STATVAR.layerThick .* 0 + ground.PARA.Tr; %initial T of 0
            energy = ground.STATVAR.T .* ground.PARA.heatCapacity; %[J/m3]
            ground.STATVAR.energy = energy .* ground.STATVAR.layerThick .* ground.STATVAR.area; %[J]

            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            
        end
                
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            %ground = get_boundary_condition_u@UB_TEMPERATURE_FORCING(ground, forcing);
            ground.TEMP.F_ub = -(ground.PARA.Tr - ground.STATVAR.T(1,1)).*ground.STATVAR.thermCond(1,1)./(ground.STATVAR.layerThick(1,1)/2); 
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) - ground.TEMP.F_ub.*ground.STATVAR.area(1,1);
        end
        
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivative_energy(ground);

        end
        
        function timestep = get_timestep(ground, tile) 
           timestep = get_timestep_heat_coduction(ground);
        end
        
        function ground = advance_prognostic(ground, tile)           
            timestep = tile.timestep;
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            % do nothing, used with surface energy ballance classes
        end
       
        function ground = compute_diagnostic(ground, tile)
            %ground = get_T_water_freeW(ground);
            ground.STATVAR.T = ground.STATVAR.energy ./ (ground.STATVAR.layerThick .* ground.STATVAR.area) ./ ground.STATVAR.heatCapacity;
            %ground = conductivity(ground);
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
        end
        
        function ground = check_trigger(ground, tile)
           %do nothing 
        end
        
        
        %-----non-mandatory functions-------
        
        
        
        %-----LATERAL-------------------
        
        %-------LAT3D_HEAT-------------
        
        
        
        %-------------param file generation-----


    end
    
end
