%========================================================================
% CryoGrid GROUND class GROUND_freezeC_RichardsEqW_seb_pressure
% heat conduction, Richards equation water scheme, freeze curve based on
% freezing=drying assumption, surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freezeC_RichardsEqW_seb_pressure < SEB & HEAT_CONDUCTION & FREEZE_CURVE_KarraPainter & WATER_FLUXES & HEAT_FLUXES_LATERAL & WATER_FLUXES_LATERAL & SOIL_MECHANICS %& INITIALIZE & FREEZE_CURVE_Painter

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
%         function ground = GROUND_freezeC_RichardsEqW_seb(index, pprovider, cprovider, forcing)  
%             ground@INITIALIZE(index, pprovider, cprovider, forcing);
%         end

        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];  %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; % roughness length [m] 

            ground.PARA.permeability = [];  %permeability for fluids/gases [m2]
            
            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
            ground.PARA.dWater_max = []; %%maximum possible volumteric water content change per timestep [-] 

            ground.PARA.LUT_size_waterIce = []; %size of lookup table for the waterIce variable [-]
            ground.PARA.LUT_size_T = [];   %size of lookup table for the (temperature) T variable [-]
            ground.PARA.min_T = []; %minimum temperature for which the LUT is calculated (modeled temperatures must be above this value) [degree C]
            ground.PARA.min_waterIce = [];  %minimum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be above this value) [-]
            ground.PARA.max_waterIce = []; %maximum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be below this value) [-]
            ground.PARA.min_mineral_organic = [];  %maximum mineral plus organic content in volumetric fraction for which the LUT is calculated (mineral plus organic content must be below this value) [-]
            
% % %             ground.PARA.hardBottom_cutoff = [];
% % %             %ground.PARA.reference_pressure = 100; %[Pa] correpsonds to ~10cm water column
            
            ground.PARA.hardBottom_cutoff = 0.03;
%            ground.PARA.reference_pressure = 100; %[Pa] correpsonds to ~10cm water column
            ground.PARA.smoothing_factor = [];
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = [];  % upper surface elevation [m]
            ground.STATVAR.lowerPos = [];  % lower surface elevation [m]
            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = [];  % total volume of water plus ice in a grid cell [m3]
            %ground.STATVAR.saturation = [];
            ground.STATVAR.mineral = [];   % total volume of minerals [m3]
            ground.STATVAR.organic = [];   % total volume of organics [m3]
            ground.STATVAR.energy = [];    % total internal energy [J]
            ground.STATVAR.soil_type = []; % integer code for soil_type; 1: sand; 2: silt: 3: clay: 4: peat; 5: water (i.e. approximation of free water, very large-pore ground material).
            ground.STATVAR.gamma_dry = []; % dry unit weight of soil [kN/m3]
            ground.STATVAR.grain_size = []; % grain size of soil [m]
            ground.STATVAR.Cc = [];        % Compression index of soil []
            ground.STATVAR.Cs = [];        % Swelling index of soil []
            ground.STATVAR.e_ini = [];        % Initial void ratio of soil []
            
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
            
            ground.STATVAR.permeability = [];  
            
            ground.STATVAR.initial_voidRatio = []; %=por/(1-por)  
            ground.STATVAR.compression_index = [];
            ground.STATVAR.reference_pressure = []; %100; %[Pa] correpsonds to ~10cm water column
            ground.STATVAR.external_pressure = 0; %External load from above [Pa]
     
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
            
            ground.CONST.R = []; %universal gas constant
            ground.CONST.molar_mass_w = []; %molar mass of water
            
            ground.CONST.rho_w = []; % water density
            ground.CONST.rho_i = []; %ice density
            
            ground.CONST.air_pressure = 101.3; % Atmospheric pressure [kPa]
            
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
            
            %Changed Sebastian
            ground.CONST.density_water = 1000; % [kg/m3]
            ground.CONST.density_ice = 1000;
            ground.CONST.density_mineral = 2650;
            ground.CONST.density_organic = 1300;
            ground.CONST.air_pressure = 1e5; % Atmospheric pressure [Pa] -> use SI units whenever possible!
        end
        
        
        function ground = finalize_init(ground, tile)
            
            %REMOVE THIS!!
            %ground.STATVAR.initial_voidRatio = ground.STATVAR.initial_voidRatio + ground.STATVAR.T .* 0;
            %ground.STATVAR.compression_index = ground.STATVAR.compression_index + ground.STATVAR.T .* 0;
            %end remove!
            
            ground.STATVAR.initial_layerThick = ground.STATVAR.layerThick;
            
            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.T .* 0;
            
            ground.CONST.vanGen_alpha = [ ground.CONST.alpha_sand ground.CONST.alpha_silt ground.CONST.alpha_clay ground.CONST.alpha_peat ground.CONST.alpha_water];
            ground.CONST.vanGen_n = [ ground.CONST.n_sand ground.CONST.n_silt ground.CONST.n_clay ground.CONST.n_peat ground.CONST.n_water];
            ground.CONST.vanGen_residual_wc = [ ground.CONST.residual_wc_sand ground.CONST.residual_wc_silt ground.CONST.residual_wc_clay ground.CONST.residual_wc_peat ground.CONST.residual_wc_water];
            
            %Initialization for compressed soil column
            
            %get variables
            %initial_saturation = ground.STATVAR.saturation;
            initial_water = ground.STATVAR.waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area; %For initialization, soil column is completely unfrozen
            mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            organic = ground.STATVAR.organic .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            initial_porosity = ground.STATVAR.initial_voidRatio ./ (1 + ground.STATVAR.initial_voidRatio);
            %initial_water = initial_saturation .* initial_porosity .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            initial_saturation = round((initial_water ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ initial_porosity,2); %round because 100% saturation is difficult to set in excel with waterIce
            ground.STATVAR.initial_layerThick = ground.STATVAR.layerThick;
            
            %Overburden pressure
            threshold = 0.5; above_threshold = initial_saturation >= threshold;
            overburden_pressure_per_cell_normal = (mineral .* ground.CONST.density_mineral + organic .* ground.CONST.density_organic + ...
                initial_water .* ground.CONST.density_water) ./ ground.STATVAR.area .*ground.CONST.g; %[Pa]
            overburden_pressure_per_cell_normal_boyant = (mineral .* (ground.CONST.density_mineral - ground.CONST.density_water) + organic .* (ground.CONST.density_organic - ground.CONST.density_water) + ...
                + (ground.STATVAR.layerThick - mineral - organic - initial_water) .* (-ground.CONST.density_water)) ./ ground.STATVAR.area .*ground.CONST.g; %[Pa]
            overburden_pressure_per_cell = double(~above_threshold) .* overburden_pressure_per_cell_normal + ...
                double(above_threshold) .* (overburden_pressure_per_cell_normal + (initial_saturation-threshold)./ (1-threshold) .* (overburden_pressure_per_cell_normal_boyant - overburden_pressure_per_cell_normal));
            
            ground.STATVAR.overburden_pressure = cumsum(overburden_pressure_per_cell)-overburden_pressure_per_cell./2; %[Pa]
            ground.STATVAR.overburden_pressure = ground.STATVAR.overburden_pressure + ground.STATVAR.external_pressure;
            
            %Calculate porosity for given overburden pressure
            ground.STATVAR.porosity = (ground.STATVAR.initial_voidRatio - ground.STATVAR.compression_index .* log10(ground.STATVAR.overburden_pressure./ground.STATVAR.reference_pressure)) ./ ...
                (1 + ground.STATVAR.initial_voidRatio - ground.STATVAR.compression_index .* log10(ground.STATVAR.overburden_pressure./ground.STATVAR.reference_pressure));
            ground.STATVAR.porosity = min(ground.STATVAR.porosity, ground.STATVAR.initial_voidRatio ./(1 + ground.STATVAR.initial_voidRatio)); %do not allow higher porosities than initial void ration
            ground.STATVAR.void_ratio = ground.STATVAR.porosity ./ (1 - ground.STATVAR.porosity); %Void ratio
            ground.STATVAR.layerThick = ((mineral + organic) ./ (1 - ground.STATVAR.porosity)) ./ ground.STATVAR.area;
            
            ground.STATVAR.waterIce = initial_saturation .* ground.STATVAR.porosity;
            ground.STATVAR.saturation = ground.STATVAR.waterIce ./ ground.STATVAR.porosity;
            
            %Bearing capacity
            ground.STATVAR.bearing_capacity = (10.^((ground.STATVAR.initial_voidRatio - ground.STATVAR.porosity - ground.STATVAR.initial_voidRatio .* ground.STATVAR.porosity) ./ (ground.STATVAR.compression_index - ground.STATVAR.porosity .* ground.STATVAR.compression_index))) .* ground.STATVAR.reference_pressure;  %[Pa]
            ground.STATVAR.bearing_capacity = max(ground.STATVAR.reference_pressure, ground.STATVAR.bearing_capacity); %generate "Xice" (void ratio becomes smaller than initial void ratio)    

            %FreezeC
            ground = get_E_freezeC_pressure(ground);
            ground = conductivity(ground);
            ground = calculate_hydraulicConductivity_RichardsEq(ground);
            
            ground = create_LUT_freezeC(ground);
                        
            ground.STATVAR.field_capacity = (1+((ground.STATVAR.alpha./ground.CONST.g./ground.CONST.rho_w).*(ground.CONST.g.*ground.CONST.rho_w)).^ground.STATVAR.n).^(-(1-1./ground.STATVAR.n)).*ground.STATVAR.porosity;
            
%             %Update saturation for water and ice
             ground.STATVAR.saturation = (ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.porosity;
             ground.STATVAR.water_saturation = (ground.STATVAR.water ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.porosity;
%             ground.STATVAR.water_saturation = (ground.STATVAR.water ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.porosity;
%             ground.STATVAR.ice_saturation = (ground.STATVAR.ice ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.porosity;
                       
            %Find all saturated cells
            saturated = ground.STATVAR.saturation > 1-1e-6;
            %Calculate gravitational potential, assuming that all gridcells are unsaturated
            for i = size(saturated,1) : -1 : 1
                gravitationalPotential_unsaturated(i,1) = sum(ground.STATVAR.layerThick(i:size(saturated,1),1));
            end
            %Calculate gravitational potential
            gravitationalPotential = [];
            for i = 1 : size(saturated,1)
                if saturated(i,1) == 0
                    gravitationalPotential(i,1) = gravitationalPotential_unsaturated(i,1);
                    if i ~= size(saturated,1) && saturated(i+1,1) == 1 && ground.STATVAR.T(i+1,1) > 0 %&& abs(ground.STATVAR.waterPotential(i,1)-ground.STATVAR.waterPotential(i+1,1)) < gravitationalPotential_unsaturated(i,1)-gravitationalPotential_unsaturated(i+1,1)
                        %If gridcell below is saturated and waterPotential is smaller than gravitational potential
                        %--> Water would be pressed into saturated zone --> should not be the case
                        %--> Add waterPotential to gravitationalPotential so that flux into saturated zone is prevented(d_head = 0)
                        gravitationalPotential(i+1,1) = gravitationalPotential_unsaturated(i,1) + ground.STATVAR.waterPotential(i,1);
                    end
                elseif saturated(i,1) == 1
                    if size(gravitationalPotential,1) == i
                        gravitationalPotential(i+1,1) = gravitationalPotential(i,1);
                    else
                        gravitationalPotential(i,1) = gravitationalPotential_unsaturated(i,1);
                        gravitationalPotential(i+1,1) = gravitationalPotential(i,1);
                    end
                end
            end
            if size(gravitationalPotential,1) > size(saturated,1)
                gravitationalPotential(end) = [];
            end
            ground.STATVAR.gravitationalPotential = gravitationalPotential;
            
            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.STATVAR.runoff = 0;
            ground.STATVAR.Xwater = ground.STATVAR.waterIce .* 0;
            ground.STATVAR.Xice = ground.STATVAR.waterIce .* 0;
            
            ground.STATVAR.year_old = []; %str2num(datestr(tile.t, 'yyyy'));
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.TEMP.surface_runoff = 0;
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            ground = surface_energy_balance(ground, forcing);
            ground = get_boundary_condition_u_RichardsEq_pressure(ground, forcing); %checked that this flux can be taken up!!
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
            ground = get_boundary_condition_l_water2(ground);  %if flux not zero, check that the water flowing out is available! Not implemented here.
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivative_energy(ground);
            ground = get_derivative_water_RichardsEq_pressure(ground);
        end
        
        function timestep = get_timestep(ground, tile)  
           timestep = get_timestep_heat_coduction(ground);
           timestep = min(timestep, get_timestep_water_RichardsEq_pressure(ground));
        end
        
        function ground = advance_prognostic(ground, tile) 
            timestep = tile.timestep;

            prognostic_layerThick = ground.STATVAR.saturation > 1-1e-6 & (ground.TEMP.d_water > 0 | (ground.TEMP.d_water < 0 & ground.TEMP.no_air) | (ground.TEMP.d_water < 0 & ground.STATVAR.void_ratio > ground.STATVAR.initial_voidRatio));
            ground.TEMP.prognostic_layerThick = prognostic_layerThick;

            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_water_energy; %add energy from water advection
            
            %water
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + timestep .* ground.TEMP.d_water; %subtract water from ET
            ground.STATVAR.waterIce = max(0, ground.STATVAR.waterIce);
            
            %layerThick is a prognostic variable for saturated gridcells that still carry excess pressure
            ground.STATVAR.layerThick(prognostic_layerThick) = (ground.STATVAR.waterIce(prognostic_layerThick) + ground.STATVAR.mineral(prognostic_layerThick) + ground.STATVAR.organic(prognostic_layerThick)) ./ ground.STATVAR.area(prognostic_layerThick);
            
            ground.STATVAR.excessWater = ground.STATVAR.excessWater + timestep .* ground.TEMP.surface_runoff;
            
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)
            
            ground = soil_mechanics(ground, tile);
            
            ground = get_T_water_freezeC_Xice_pressure(ground);
            
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + ground.STATVAR.XwaterIce;
            ground.STATVAR.water = ground.STATVAR.water + ground.STATVAR.Xwater;
            ground.STATVAR.ice = ground.STATVAR.ice + ground.STATVAR.Xice;

            ground = gravitational_potential(ground);
            ground = conductivity(ground);
            ground = calculate_hydraulicConductivity_RichardsEq(ground);

            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.TEMP.surface_runoff = 0;
            
        end
        
         function ground = check_trigger(ground, tile)
%             % Every year, the heat reservoir is 10% further away
%             if isempty(ground.STATVAR.year_old) || str2num(datestr(tile.t, 'yyyy'))>ground.STATVAR.year_old
%                 ground.STATVAR.year_old = str2num(datestr(tile.t, 'yyyy'));
%                 tile.LATERAL.IA_CLASSES{2,1}.PARA.distance_heatReservoir = tile.LATERAL.IA_CLASSES{2,1}.PARA.distance_heatReservoir .*1.1;
%             end
         end
        
        
        %-----non-mandatory functions-------
        
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground = Q_evap_CLM4_5(ground, forcing);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
            
            %water -> evaporation
            ground.TEMP.d_water_ET(1,1) = ground.TEMP.d_water_ET(1,1) +  ground.STATVAR.evaporation.* ground.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
            ground.TEMP.d_water_ET_energy(1,1) = ground.TEMP.d_water_ET_energy(1,1) +  ground.TEMP.evaporation_energy.* ground.STATVAR.area(1,1);

        end
        
        function ground = conductivity(ground)
            ground = thermalConductivity_CLM4_5(ground);
            %ground = conductivity_mixing_squares(ground);
        end
        
        %-----LATERAL-------------------
        
        %-----LAT_REMOVE_SURFACE_WATER-----
        function ground = lateral_push_remove_surfaceWater(ground, lateral)
            ground = lateral_push_remove_surfaceWater_simple(ground, lateral);
        end
        
        %-----LAT_REMOVE_SUBSURFACE_WATER-----        
        function ground = lateral_push_remove_subsurfaceWater(ground, lateral)
            ground = lateral_push_remove_subsurfaceWater_simple(ground, lateral);
        end
        
        %----LAT_SEEPAGE_FACE----------            
        function ground = lateral_push_remove_water_seepage(ground, lateral)
            ground = lateral_push_remove_water_seepage_simple(ground, lateral);
        end
        
        %----LAT_WATER_RESERVOIR------------          
        function ground = lateral_push_water_reservoir(ground, lateral)
            ground = lateral_push_water_reservoir_RichardsEq_pressure(ground, lateral);
        end
        
        %----LAT_HEAT------------
        function ground = lateral_push_heat(ground, lateral)
            ground = lateral_push_heat_simple(ground, lateral);
        end

        %----LAT3D_WATER_UNCONFINED_AQUIFER------------         
        function ground = lateral3D_pull_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_pull_water_unconfined_aquifer_simple(ground, lateral);
        end
        
        function ground = lateral3D_push_water_unconfined_aquifer(ground, lateral)
            ground = lateral3D_push_water_unconfined_aquifer_simple(ground, lateral);
        end
        
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground, lateral)
            [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_simple(ground, lateral);
        end
        
        %LAT3D_WATER_RESERVOIR and LAT3D_WATER_SEEPAGE_FACE do not require specific functions
         
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
            
            ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'soil_type' 'field_capacity' 'permeability' 'T'}; %add further variables here
            
            ground.PARA.default_value.albedo = {0.2};
            ground.PARA.comment.albedo = {'surface albedo [-]'};
            
            ground.PARA.default_value.epsilon = {0.99};
            ground.PARA.comment.epsilon = {'surface emissivity [-]'};
            
            ground.PARA.default_value.z0 = {0.01};
            ground.PARA.comment.z0 = {'roughness length [m]'};
            
            ground.PARA.default_value.rootDepth = {0.1};
            ground.PARA.comment.rootDepth = {'e-folding depth of transpiration reduction with depth [m]'};
            
            ground.PARA.default_value.evaporationDepth = {0.1};
            ground.PARA.comment.evaporationDepth = {'e-folding constant of evaporation reduction reduction with depth [m]'};
            
            ground.PARA.default_value.ratioET = {0.5};
            ground.PARA.comment.ratioET = {'fraction of transpiration of total evapotranspiration [-]'};
            
            ground.PARA.default_value.conductivity_function = {''};
            ground.PARA.comment.conductivity_function = {'function employed to calculate thermal conductivity, leave empty for default'};
            
            ground.PARA.default_value.dt_max = {3600};
            ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
            
            ground.PARA.default_value.dE_max = {50000};
            ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
            
            ground.PARA.default_value.dWater_max = {0.005};
            ground.PARA.comment.dWater_max = {'maximum possible volumteric water content change per timestep [-]'};
            
            ground.PARA.default_value.LUT_size_waterIce = {1000};
            ground.PARA.comment.LUT_size_waterIce = {'size of lookup table for the waterIce variable [-]'};
            
            ground.PARA.default_value.LUT_size_T = {1000};
            ground.PARA.comment.LUT_size_T = {'size of lookup table for the (temperature) T variable [-]'};
            
            ground.PARA.default_value.min_T = {-50};
            ground.PARA.comment.min_T = {'minimum temperature for which the LUT is calculated (modeled temperatures must be above this value) [degree C]'};
            
            ground.PARA.default_value.min_waterIce = {0.05};
            ground.PARA.comment.min_waterIce = {'minimum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be above this value) [-]'};
            
            ground.PARA.default_value.max_waterIce = {0.97};
            ground.PARA.comment.max_waterIce = {'maximum waterIce value in volumetric fraction for which the LUT is calculated (modeled waterIce must be below this value) [-]'};
            
            ground.PARA.default_value.min_mineral_organic = {0.03};
            ground.PARA.comment.min_mineral_organic = {'maximum mineral plus organic content in volumetric fraction for which the LUT is calculated (mineral plus organic content must be below this value) [-]'};
        
            ground.PARA.default_value.hardBottom_cutoff = {0.03};
            ground.PARA.comment.hardBottom_cutoff = {'ground considered impermeable for water when water content is less than cutoff'};
            
            ground.PARA.default_value.smoothing_factor = {100};
            ground.PARA.comment.smoothing_factor = {'adjustment timescale of bearing capacity'};

        end
        
        %----inherited Tier 1 functions ------------
        
%         function ground = get_derivative_energy(ground)
%             ground = get_derivative_energy@HEAT_CONDUCTION(ground);
%         end
%         
%         function ground = conductivity_mixing_squares(ground)
%             ground = conductivity_mixing_squares@HEAT_CONDUCTION(ground);
%         end
%         
%         function flux = Q_h(ground, forcing)
%             flux = Q_h@SEB(ground, forcing);
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
%             ground = get_boundary_condition_u_water2@WATER_FLUXES(ground, forcing);
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
%             ground = L_star@SEB(ground, forcing);
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
