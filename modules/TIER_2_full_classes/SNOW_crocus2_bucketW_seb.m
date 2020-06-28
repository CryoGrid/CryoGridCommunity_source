%designed to work with Xice class, Xwater takes up excess water when the
%snow is a CHILD
%if the snow is a full class water pools up if not removed and a LAKE is
%triggered - in this case the CHILD stage is skipped and the snow energy
%and water is mixed with the LAKE

classdef SNOW_crocus2_bucketW_seb < SEB & HEAT_CONDUCTION & WATER_FLUXES & HEAT_FLUXES_LATERAL & WATER_FLUXES_LATERAL & SNOW & SNOW_FLUXES_LATERAL & INITIALIZE & REGRID

    properties
        PARENT
    end
    
    
    methods
        
        function self = SNOW_crocus2_bucketW_seb(index, pprovider, cprovider, forcing)  
            self@INITIALIZE(index, pprovider, cprovider, forcing);
        end
        
        function snow = provide_PARA(snow)

            snow.PARA.epsilon = [];
            snow.PARA.airT_height = []; %measurement height [m]
            snow.PARA.z0 = []; %roughness length [m]
            
            snow.PARA.SW_spectral_range1 = [];
            snow.PARA.SW_spectral_range2 = [];
            
            snow.PARA.field_capacity = [];
            snow.PARA.hydraulicConductivity = [];
            snow.PARA.swe_per_cell = [];
            
            snow.PARA.slope = [];
            snow.PARA.timescale_winddrift = [];
            
            snow.PARA.heatFlux_lb = [];
            
            snow.PARA.dt_max = [] ; %[sec]
            snow.PARA.dE_max = []; %[J/m3]
        end
        
        function snow = provide_STATVAR(snow)
            
            snow.STATVAR.upperPos = [];
            snow.STATVAR.lowerPos = [];
            snow.STATVAR.layerThick = []; 
            snow.STATVAR.layerThickSnowFirstCell = [];  % [m]total grid cell thicness minus the excess water phase allowed for the first cell
            snow.STATVAR.area = [];
            
            snow.STATVAR.waterIce = []; % [m]
            snow.STATVAR.mineral = []; % [m]
            snow.STATVAR.organic = []; % [m]
            snow.STATVAR.energy = [];  % [J/m2]
            
            snow.STATVAR.T = [];  % [degree C]
            snow.STATVAR.water = [];  % [m]
            snow.STATVAR.ice = [];
            snow.STATVAR.air = [];  % [m]
            snow.STATVAR.thermCond = [];
            snow.STATVAR.albedo = [];
            
            snow.STATVAR.d = [];
            snow.STATVAR.s = [];
            snow.STATVAR.gs = [];
            snow.STATVAR.time_snowfall = [];
            snow.STATVAR.target_density = [];
            snow.STATVAR.excessWater = [];
            
            snow.STATVAR.Lstar = [];
            snow.STATVAR.Qh = [];
            snow.STATVAR.Qe = [];
        end
    
        function snow = provide_CONST(snow)
            
            snow.CONST.L_f = [];
            
            snow.CONST.c_w = [];
            snow.CONST.c_i = [];
            snow.CONST.c_o = [];
            snow.CONST.c_m = [];
            
            snow.CONST.k_a = [];       %air [Hillel(1982)]
            snow.CONST.k_w = [];        %water [Hillel(1982)]
            snow.CONST.k_i = [];         %ice [Hillel(1982)]
            snow.CONST.k_o = [];        %organic [Hillel(1982)]
            snow.CONST.k_m = [];
            
            snow.CONST.sigma = []; %Stefan-Boltzmann constant
            snow.CONST.kappa = [];
            snow.CONST.L_s = []; %latent heat of vaporization
            
            snow.CONST.cp = [];
            snow.CONST.g = [];
            
            snow.CONST.rho_w = [];
            snow.CONST.rho_i = [];
        end
        
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function snow = finalize_init(snow, forcing) %assign all variables, that must be calculated or assigned otherwise for initialization
            snow.PARA.heatFlux_lb = forcing.PARA.heatFlux_lb;
            snow.PARA.airT_height = forcing.PARA.airT_height;
            
            snow = initialize_zero_snow_BASE(snow);  %initialize all values to be zero
            snow.PARA.spectral_ranges = [snow.PARA.SW_spectral_range1 snow.PARA.SW_spectral_range2 1 - snow.PARA.SW_spectral_range1 - snow.PARA.SW_spectral_range2];
            
            snow.TEMP.d_energy = snow.STATVAR.energy .*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
            
        end
        
        %---time integration------
        
        function snow = get_boundary_condition_u(snow, forcing) 
            
            snow = get_boundary_condition_SNOW_u(snow, forcing);
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            
            snow = get_snow_properties_crocus(snow,forcing); %makes a TEMP variable newSnow that contains all information on the fresh snow - which is merged in the diagnostic step
            
            snow = surface_energy_balance(snow, forcing);
            snow = get_sublimation(snow, forcing);
            
            snow.TEMP.wind = forcing.TEMP.wind;
        end
        
        function snow = get_boundary_condition_u_CHILD(snow, forcing)
            snow = get_boundary_condition_allSNOW_rain_u(snow, forcing); %add full snow, but rain only for snow-covered part
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            
            snow = get_snow_properties_crocus(snow,forcing); %makes a TEMP variable newSnow that contains all information on the fresh snow - which is merged in the diagnostic step
            
            snow = surface_energy_balance(snow, forcing); %this works inclduing penetration of SW radiation through the CHILD snow
            snow = get_sublimation(snow, forcing);
            
            snow.TEMP.wind = forcing.TEMP.wind;
        end
        
        function snow = get_boundary_condition_u_create_CHILD(snow, forcing)
            snow = get_boundary_condition_allSNOW_u(snow, forcing); %add all snow, no rain
            
            snow = get_snow_properties_crocus(snow,forcing);
            
            snow.TEMP.F_ub = 0;
            snow.TEMP.F_lb = 0;
            snow.TEMP.F_ub_water = 0;
            snow.TEMP.F_lb_water = 0;
            snow.TEMP.F_ub_water_energy = 0;
            snow.TEMP.F_lb_water_energy = 0;
            snow.STATVAR.sublimation = 0;
            snow.TEMP.sublimation_energy = 0;
            snow.TEMP.rain_energy = 0;
            snow.TEMP.rainfall = 0;
            
            snow.TEMP.d_energy = 0;
            snow.TEMP.d_water = 0;
            snow.TEMP.d_water_energy = 0;
            
            snow.STATVAR.d = 0;
            snow.STATVAR.s = 0;
            snow.STATVAR.gs = 0;
            snow.STATVAR.time_snowfall = 0;
            snow.TEMP.metam_d_d = 0;
            snow.TEMP.wind_d_d = 0;
            snow.TEMP.metam_d_s = 0;
            snow.TEMP.wind_d_s = 0;
            snow.TEMP.metam_d_gs = 0;
            snow.TEMP.compact_d_D = 0;
            snow.TEMP.wind_d_D = 0;           
            snow.TEMP.wind = forcing.TEMP.wind;
            %start with some non-zero values for area and layerThick
            snow.STATVAR.area = 1;
            snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ (snow.TEMP.newSnow.STATVAR.density ./1000); %[m] layerThick adjusted so that always 0.5 .* snow.PARA.swe_per_cell is contained
            snow.STATVAR.energy = 0;
            snow.STATVAR.waterIce = 0;
            snow.STATVAR.ice = 0;
            snow.STATVAR.excessWater = 0;
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos;
        end
        
       function [snow, S_up] = penetrate_SW(snow, S_down)  %mandatory function when used with class that features SW penetration
            [snow, S_up] = penetrate_SW_transmission_spectral(snow, S_down);
        end
        
        function snow = get_boundary_condition_l(snow, forcing)
             snow.TEMP.F_lb = forcing.PARA.heatFlux_lb .* snow.STATVAR.area(end);
             snow.TEMP.d_energy(end) = snow.TEMP.d_energy(end) + snow.TEMP.F_lb;
             
             snow.TEMP.F_lb_water = 0;
             snow.TEMP.F_lb_water_energy = 0;
        end
        
        function snow = get_derivatives_prognostic(snow)
            if size(snow.STATVAR.layerThick,1) > 1
                snow = get_derivative_energy(snow);
                snow = get_derivative_water_SNOW(snow);
                snow = get_T_gradient_snow(snow);
            else
                snow = get_T_gradient_snow_single_cell(snow);
            end
            store = snow.STATVAR.layerThick(1);
            difference = snow.STATVAR.layerThick(1) - snow.STATVAR.layerThickSnowFirstCell;
            snow.STATVAR.layerThick(1) =  snow.STATVAR.layerThickSnowFirstCell;
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) - difference .* snow.STATVAR.area(1,1);
            
            snow = prog_metamorphism(snow);
            snow = prog_wind_drift(snow); % Add blowing snow sublimation according to Gordon etAl 2006
            snow = compaction(snow);
            
            snow.STATVAR.layerThick(1) = store;
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + difference .* snow.STATVAR.area(1,1);
        end
        
        function snow = get_derivatives_prognostic_CHILD(snow)
            %snow.TEMP.d_energy = snow.TEMP.d_energy + snow.TEMP.F_ub + snow.TEMP.F_lb;
            %snow.TEMP.d_water = snow.TEMP.d_water + snow.TEMP.F_ub_water + snow.TEMP.F_lb_water;
            %snow.TEMP.d_water_energy = snow.TEMP.d_water_energy + snow.TEMP.F_ub_water_energy + snow.TEMP.F_lb_water_energy;
            
            snow = get_T_gradient_snow_single_cell(snow);
            snow = prog_metamorphism(snow);
            snow = prog_wind_drift(snow); % Add blowing snow sublimation according to Gordon etAl 2006
            snow = compaction(snow);
        end
        
        function timestep = get_timestep(snow) 
            timestep1 = get_timestep_heat_coduction(snow);
            timestep2 = get_timestep_SNOW_mass_balance(snow);
            timestep3 = get_timestep_water_SNOW(snow);

            timestep = min(timestep1, timestep2);
            timestep = min(timestep, timestep3);
        end
        
        function timestep = get_timestep_CHILD(snow)  
            timestep = get_timestep_SNOW_CHILD(snow);
        end
        
        function snow = advance_prognostic(snow, timestep)
            %energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy + snow.TEMP.d_water_energy);
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* snow.TEMP.sublimation_energy;  %snowfall energy added below, when new snow layer is merged
            %mass
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* snow.TEMP.d_water;
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + timestep .* snow.STATVAR.sublimation;
            density_ice_phase = snow.STATVAR.ice(1) ./ snow.STATVAR.layerThickSnowFirstCell ./snow.STATVAR.area(1);
            
            snow.STATVAR.ice(1) = snow.STATVAR.ice(1) + timestep .* snow.STATVAR.sublimation;

%             snow.STATVAR.layerThickSnowFirstCell = snow.STATVAR.layerThickSnowFirstCell + timestep .* snow.STATVAR.sublimation ./snow.STATVAR.area(1,1) ./ ...
%                 (snow.STATVAR.ice(1) ./ snow.STATVAR.layerThickSnowFirstCell ./snow.STATVAR.area(1));
            snow.STATVAR.layerThickSnowFirstCell = snow.STATVAR.layerThickSnowFirstCell + timestep .* snow.STATVAR.sublimation ./snow.STATVAR.area(1,1) ./ density_ice_phase;
            
            %microphysics
            snow.STATVAR.d = max(snow.STATVAR.d.*0, snow.STATVAR.d + timestep .*(snow.TEMP.metam_d_d + snow.TEMP.wind_d_d));
            snow.STATVAR.s = max(snow.STATVAR.s.*0, min(snow.STATVAR.s.*0+1, snow.STATVAR.s + timestep .*(snow.TEMP.metam_d_s + snow.TEMP.wind_d_s)));
            snow.STATVAR.gs = max(snow.STATVAR.gs, snow.STATVAR.gs + timestep .*(snow.TEMP.metam_d_gs + snow.TEMP.wind_d_gs));
            snow.STATVAR.layerThick = min(snow.STATVAR.layerThick, max(snow.STATVAR.ice ./ snow.STATVAR.area, snow.STATVAR.layerThick + timestep .*(snow.TEMP.compact_d_D + snow.TEMP.wind_d_D)));
            snow.STATVAR.layerThickSnowFirstCell = min(snow.STATVAR.layerThickSnowFirstCell, max(snow.STATVAR.ice(1) ./ snow.STATVAR.area(1), snow.STATVAR.layerThickSnowFirstCell + timestep .*(snow.TEMP.compact_d_D(1) + snow.TEMP.wind_d_D(1))));
            
            %snow.STATVAR.layerThick(1) = max(snow.STATVAR.layerThickSnowFirstCell, snow.STATVAR.waterIce(1));
            %Xwater = snow.STATVAR.layerThick(1) - snow.STATVAR.layerThickSnowFirstCell;
            %snow.STATVAR.layerThick(1) =  snow.STATVAR.layerThickSnowFirstCell;
            snow.STATVAR.layerThick(1) =  snow.STATVAR.layerThickSnowFirstCell;
            %new snow
            if snow.TEMP.snowfall >0
                snow = advance_prognostic_new_snow_crocus(snow, timestep);
                %merge with uppermost layer
                
                snow = merge_cells_intensive2(snow, 1, snow.TEMP.newSnow, 1, {'d'; 's'; 'gs'; 'time_snowfall'}, 'ice'); %account for the case that ice = 0?? Maybe not necessary
                snow = merge_cells_extensive2(snow, 1, snow.TEMP.newSnow, 1, {'waterIce'; 'energy'; 'layerThick'; 'ice'}); %CHANGE BOTH
                
                snow.STATVAR.layerThickSnowFirstCell = snow.STATVAR.layerThick(1);
            end
            
            %snow.STATVAR.layerThick(1) = max(snow.STATVAR.layerThick(1), snow.STATVAR.waterIce(1));
            
            %store "old" density - ice is updated for new snowfall and sublimation losses
            snow.STATVAR.target_density = snow.STATVAR.ice ./ snow.STATVAR.layerThick ./ snow.STATVAR.area;
            %snow.STATVAR.target_density(1) = snow.STATVAR.ice(1) ./ snow.STATVAR.layerThickSnowFirstCell ./ snow.STATVAR.area(1);
            snow.STATVAR.target_density = min(1, snow.STATVAR.target_density);  % avoids rounding errors, if >1, this might trigger a follow-up problem
            
            snow.STATVAR.layerThick(1) = max(snow.STATVAR.layerThick(1), snow.STATVAR.waterIce(1) ./ snow.STATVAR.area(1,1));
        end
        
        %--------- to here
        
        
        function snow = advance_prognostic_CHILD(snow, timestep)
            
            %energy
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* (snow.TEMP.d_energy  + snow.TEMP.d_water_energy + snow.TEMP.sublimation_energy);
            %mass
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* (snow.TEMP.d_water + snow.STATVAR.sublimation);
            snow.STATVAR.ice = snow.STATVAR.ice + timestep .* snow.STATVAR.sublimation;
            snow.STATVAR.volume = double(snow.STATVAR.ice>0) .* snow.STATVAR.layerThick .* snow.STATVAR.area;
            snow.STATVAR.volume = snow.STATVAR.volume + timestep .* snow.STATVAR.sublimation ./ max(50, snow.STATVAR.ice ./ snow.STATVAR.layerThick); 
            %microphysics
            snow.STATVAR.d = max(snow.STATVAR.d.*0, snow.STATVAR.d + timestep .*(snow.TEMP.metam_d_d + snow.TEMP.wind_d_d));
            snow.STATVAR.s = max(snow.STATVAR.s.*0, min(snow.STATVAR.s.*0+1, snow.STATVAR.s + timestep .*(snow.TEMP.metam_d_s + snow.TEMP.wind_d_s)));
            snow.STATVAR.gs = max(snow.STATVAR.gs, snow.STATVAR.gs + timestep .*(snow.TEMP.metam_d_gs + snow.TEMP.wind_d_gs));
            snow.STATVAR.volume = min(snow.STATVAR.volume, max(snow.STATVAR.ice, snow.STATVAR.volume + timestep .* snow.STATVAR.area .* (snow.TEMP.compact_d_D + snow.TEMP.wind_d_D))); %mass is conserved, reduce layerthick
            
            %new snow and merge
            if snow.TEMP.snowfall >0
                snow = advance_prognostic_new_snow_CHILD_crocus(snow, timestep);  %add new snow with the new layerThick
                %merge with uppermost layer
                snow = merge_cells_intensive2(snow, 1, snow.TEMP.newSnow, 1, {'d'; 's'; 'gs'; 'time_snowfall'}, 'ice');
                snow = merge_cells_extensive2(snow, 1, snow.TEMP.newSnow, 1, {'waterIce'; 'energy'; 'volume'; 'ice'});
            end
            
            %store "old" density - ice is updated for new snowfall and sublimation losses
            snow.STATVAR.target_density = min(1, snow.STATVAR.ice ./ snow.STATVAR.volume);
            
            %adjust layerThick, so that exactly 0.5 .* snow.PARA.swe_per_cell is contained
            snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ snow.STATVAR.target_density;
            snow.STATVAR.area = snow.STATVAR.volume ./ snow.STATVAR.layerThick;
        end
        
        function snow = compute_diagnostic_first_cell(snow, forcing)
            snow = L_star(snow, forcing);
        end
       
        function snow = compute_diagnostic(snow, forcing)

            snow = get_T_water_freeW(snow);
            

            
            %snow.STATVAR.layerThick(1) =
            %snow.STATVAR.layerThickSnowFirstCell;  %taken out SEBASTIAN 
            
            snow = subtract_water2(snow);

            
            [snow, regridded_yesNo] = regrid_snow(snow, {'waterIce'; 'energy'; 'layerThick'}, {'area'; 'target_density'; 'd'; 's'; 'gs'; 'time_snowfall'}, 'ice');
            
            if size(snow.STATVAR.ice,1) ~= size(snow.STATVAR.energy,1) || size(snow.STATVAR.T,1) ~= size(snow.STATVAR.energy,1) || size(snow.STATVAR.layerThick,1) ~= size(snow.STATVAR.energy,1)
                disp('Hallo4')
                size(snow.STATVAR.ice,1)
                size(snow.STATVAR.energy,1)
                size(snow.STATVAR.T,1)
            end
            
            if regridded_yesNo
                snow = get_T_water_freeW(snow);
            end
            
            if size(snow.STATVAR.ice,1) ~= size(snow.STATVAR.energy,1) ||size(snow.STATVAR.T,1) ~= size(snow.STATVAR.energy,1) || size(snow.STATVAR.layerThick,1) ~= size(snow.STATVAR.energy,1)
                disp('Hallo5')
                size(snow.STATVAR.ice,1)
                size(snow.STATVAR.energy,1)
                size(snow.STATVAR.T,1)
            end
            
            
            snow.STATVAR.layerThickSnowFirstCell = snow.STATVAR.layerThick(1);
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + snow.STATVAR.excessWater;
            %snow.STATVAR.water(1) = snow.STATVAR.water(1) + snow.STATVAR.excessWater;
            snow.STATVAR.water(1) = snow.STATVAR.waterIce(1) - snow.STATVAR.ice(1);
            snow.STATVAR.layerThick(1) = max(snow.STATVAR.layerThick(1), snow.STATVAR.waterIce(1) ./ snow.STATVAR.area(1,1));
            
            snow.STATVAR.excessWater = 0;
 
            snow = conductivity(snow);
            snow.STATVAR.upperPos = snow.NEXT.STATVAR.upperPos + sum(snow.STATVAR.layerThick);
            
            snow = calculate_albedo_crocus(snow, forcing); %albedo calculation is a diagnostic operation
                        
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
            

        end
        
        function snow = compute_diagnostic_CHILD(snow, forcing)

            snow = get_T_water_freeW(snow);
            snow = subtract_water_CHILD2(snow);

            snow = conductivity(snow);
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos + (snow.STATVAR.layerThick .* snow.STATVAR.area ./ snow.PARENT.STATVAR.area(1,1));
            
            snow = calculate_albedo_crocus(snow, forcing);
            
            snow.TEMP.d_energy = snow.STATVAR.energy.*0;
            snow.TEMP.d_water = snow.STATVAR.energy .*0;
            snow.TEMP.d_water_energy = snow.STATVAR.energy .*0;
            
            snow.STATVAR.layerThickSnowFirstCell = snow.STATVAR.layerThick(1); %set this to be ready for the the non-CHILD phase
            
            remove_excessWater_CHILD(snow.PARENT.IA_CHILD);  % depends on next class, coded in IA class; 

        end
        
        
        function snow = check_trigger(snow, forcing)
            trigger_yes_no = 0;
            
            if ~trigger_yes_no
                snow = make_SNOW_CHILD(snow);
            end
            
        end
        
        %-----non-mandatory functions-------
        function snow = surface_energy_balance(snow, forcing)
            snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
            
            [snow, S_up] = penetrate_SW(snow, snow.PARA.spectral_ranges .* forcing.TEMP.Sin .* snow.STATVAR.area(1)); %distribute SW radiation
            snow.STATVAR.Sout = sum(S_up) ./ snow.STATVAR.area(1);
            
            snow.STATVAR.Qh = Q_h(snow, forcing);
            snow.STATVAR.Qe = Q_eq_potET(snow, forcing);
            
            snow.TEMP.F_ub = (forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Qh - snow.STATVAR.Qe) .* snow.STATVAR.area(1);
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.F_ub;
        end
        
        
        
        function snow = conductivity(snow)
            snow = conductivity_snow_Yen(snow);
        end
        
        
        %lateral fluxes---------------------------

        
        function snow = lateral_push_remove_water_seepage(snow, lateral)
            snow = lateral_push_remove_water_seepage_snow(snow, lateral);
        end
        
        function snow = lateral3D_pull_water_unconfined_aquifer(snow, lateral)
            snow = lateral3D_pull_water_unconfined_aquifer_snow(snow, lateral);
        end
        
        function snow = lateral3D_push_water_unconfined_aquifer(snow, lateral)
            snow = lateral3D_push_water_unconfined_aquifer_snow(snow, lateral);
        end
        
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(snow, lateral)
            [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_snow(snow, lateral);
        end
        
        function snow = lateral3D_pull_heat(snow, lateral)
            snow = lateral3D_pull_heat_simple(snow, lateral);
        end
        
        function snow = lateral3D_push_heat(snow, lateral)
            snow = lateral3D_push_heat_simple(snow, lateral);
        end
        
        function snow = lateral3D_pull_snow(snow, lateral)
            snow = lateral3D_pull_snow_crocus(snow, lateral);
        end
        
        function snow = lateral3D_push_snow(snow, lateral)
            snow = lateral3D_push_snow_crocus2(snow, lateral);
        end

        
        %----inherited Tier 1 functions ------------
        
        function snow = get_derivative_energy(snow)
           snow = get_derivative_energy@HEAT_CONDUCTION(snow); 
        end
        
        function snow = conductivity_snow_Yen(snow)
            snow = conductivity_snow_Yen@HEAT_CONDUCTION(snow);
        end
        
        function flux = Q_h(snow, forcing)
           flux = Q_h@SEB(snow, forcing);
        end
    
        function flux = Q_eq_potET(snow, forcing)
            flux = Q_eq_potET@SEB(snow, forcing);
        end
        
        function timestep = get_timestep_heat_coduction(snow)
            timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(snow);
        end
        
        function timestep = get_timestep_SNOW_mass_balance(snow)
            timestep = get_timestep_SNOW_mass_balance@SNOW(snow);
        end
        
        function timestep = get_timestep_SNOW_CHILD(snow)
            timestep = get_timestep_SNOW_CHILD@SNOW(snow);
        end
        
        function snow = L_star(snow, forcing)
           snow = L_star@SEB(snow, forcing); 
        end
        
         function [snow, S_up] = penetrate_SW_transmission_spectral(snow, S_down)
             [snow, S_up] = penetrate_SW_transmission_spectral@SEB(snow, S_down);
         end
        
        function snow = get_T_water_freeW(snow)
            snow = get_T_water_freeW@HEAT_CONDUCTION(snow);
        end
        
        function snow = subtract_water(snow)
            snow = subtract_water@SNOW(snow);
        end
        
        function snow = subtract_water_CHILD(snow)
            snow = subtract_water_CHILD@SNOW(snow);
        end
        
        function snow = get_boundary_condition_SNOW_u(snow, forcing)
            snow = get_boundary_condition_SNOW_u@SNOW(snow, forcing);
        end
        
        function snow = get_boundary_condition_allSNOW_rain_u(snow, forcing)
            snow = get_boundary_condition_allSNOW_rain_u@SNOW(snow, forcing);
        end
        
        function snow = get_boundary_condition_allSNOW_u(snow, forcing) 
            snow = get_boundary_condition_allSNOW_u@SNOW(snow, forcing);
        end
        
        function snow = initialize_zero_snow_BASE(snow)
            snow = initialize_zero_snow_BASE@SNOW(snow);
        end
        
        function snow = make_SNOW_CHILD(snow)
            snow = make_SNOW_CHILD@SNOW(snow);
        end
        
        function [snow, regridded_yesNo] = regrid_snow(snow, extensive_variables, intensive_variables, intensive_scaling_variable)
            [snow, regridded_yesNo] = regrid_snow@REGRID(snow, extensive_variables, intensive_variables, intensive_scaling_variable);
        end
    end
    
end
