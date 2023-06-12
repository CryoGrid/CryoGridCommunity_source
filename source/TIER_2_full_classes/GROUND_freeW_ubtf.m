%========================================================================
% CryoGrid GROUND class GROUND_freeW_ubft
% heat conduction, free water freeze curve, upper boundary temperature
% forcing
% T.Ingeman-Nielsen, S. Westermann, December 2021
%========================================================================


classdef GROUND_freeW_ubtf < HEAT_CONDUCTION &  UB_TEMPERATURE_FORCING

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = provide_PARA(ground)
            ground.PARA.conductivity_function = [];
            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = [];  % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = [];  % total volume of minerals [m3]
            ground.STATVAR.organic = [];  % total volume of ice [m3]
            ground.STATVAR.energy = [];   % total internal energy[J]
            
            ground.STATVAR.T = [];     % temperature [degree C]
            ground.STATVAR.water = []; % total volume of water [m3]
            ground.STATVAR.ice = [];   % total volume of ice [m3]
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
            ground.CONST.kappa = []; % von Karman constant
            ground.CONST.L_s = [];   % latent heat of sublimation, latent heat of evaporation handled in a dedicated function
            
            ground.CONST.cp = []; % specific heat capacity at constant pressure of air
            ground.CONST.g = [];  % gravitational acceleration Earth surface
            
            ground.CONST.rho_w = [];   % water density
            ground.CONST.rho_i = [];   %ice density
        end

        
        function ground = convert_units(ground, tile)
                unit_converter = str2func(tile.PARA.unit_conversion_class);
                unit_converter = unit_converter();
                ground = convert_normal_ubT(unit_converter, ground, tile);
        end


        function ground = finalize_init(ground, tile)

            if isempty(ground.PARA.conductivity_function) || sum(isnan(ground.PARA.conductivity_function))>0
                ground.PARA.conductivity_function = 'conductivity_mixing_squares';
            end            

            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
%             ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.T .* 0;
            
            ground = get_E_freeW(ground);

            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            
        end
        
        function ground = finalize_init2(ground, tile)
            
            ground = get_E_freeW(ground);
            
        end
                
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            ground = get_boundary_condition_u@UB_TEMPERATURE_FORCING(ground, forcing);
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
            ground = get_T_water_freeW(ground);
            ground = conductivity(ground);
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
        end
        
        function ground = check_trigger(ground, tile)
           %do nothing 
        end
        
        
        %-----non-mandatory functions-------
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        
        %-----LATERAL-------------------
        
        %-------LAT3D_HEAT-------------
        function ground = lateral3D_pull_heat(ground, lateral)
            ground = lateral3D_pull_heat_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_heat(ground, lateral)
            ground = lateral3D_push_heat_simple(ground, lateral);
        end
        
        
        
        %-------------param file generation-----
         function ground = param_file_info(ground)
             ground = param_file_info@BASE(ground);
             
             ground.PARA.class_category = 'GROUND';
             
             ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'T'};
             
             ground.PARA.default_value.dt_max = {3600};
             ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
             
             ground.PARA.default_value.dE_max = {50000};
             ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end

    end
    
end
