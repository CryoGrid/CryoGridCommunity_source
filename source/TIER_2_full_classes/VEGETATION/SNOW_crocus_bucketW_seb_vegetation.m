%========================================================================
% CryoGrid GROUND class SNOW_crocus_bucketW_seb
% CROCUS snow model Vionnet et al., 2012, but with simpler layer splitting and regridding scheme compared to CROCUS 
% temperature and windspeed-dependent initial snow density, snow microstructure (dendricity, sphericity, grain size), 
% compaction, sublimation, water flow, refreezing, variable albedo.
% Meltwater exceeding the available pore space within the snow cover is automatically removed.
% R. Zweigel, S. Westermann, October 2020
%========================================================================

classdef SNOW_crocus_bucketW_seb_vegetation < SNOW_crocus_bucketW_seb
    
    properties
        VEGETATION
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

        
        function snow = provide_PARA(snow)

            snow.PARA.epsilon = []; %surface emissivity [-]
            snow.PARA.z0 = [];      %roughness length [m]
            
            snow.PARA.SW_spectral_range1 = []; %fraction of incoming short-wave radiation in first spectral band [-], see Vionnet et al.,2012
            snow.PARA.SW_spectral_range2 = []; %fraction of incoming short-wave radiation in second spectral band [-], fraction of third spectral band calculated automatically
            
            snow.PARA.field_capacity = []; %snow field capacity in fraction of available pore space [-] NOTE: the definition is different for GROUND_XX classes
            snow.PARA.hydraulicConductivity = []; %hydraulic conductivity of snow [m/sec]
            snow.PARA.swe_per_cell = [];  %target SWE per grid cell [m]
            
            snow.PARA.slope = [];  %slope angle [-]
            snow.PARA.timescale_winddrift = []; %timescale of snow compaction for wind drift [hours!!]
            
            snow.PARA.dt_max = [];  %maximum possible timestep [sec]
            snow.PARA.dE_max = [];  %maximum possible energy change per timestep [J/m3]
        end
        
        function snow = provide_STATVAR(snow)
            
            snow.STATVAR.upperPos = [];  % upper surface elevation [m]
            snow.STATVAR.lowerPos = [];  % lower surface elevation [m]
            snow.STATVAR.layerThick = [];  % thickness of grid cells [m]
            snow.STATVAR.area = [];     % grid cell area [m2]
            
            snow.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            snow.STATVAR.mineral = [];  % total volume of minerals [m3]
            snow.STATVAR.organic = [];  % total volume of organics [m3]
            snow.STATVAR.energy = [];   % total internal energy[J]
            
            snow.STATVAR.T = [];      % temperature [degree C]
            snow.STATVAR.water = [];  % total volume of water [m3]
            snow.STATVAR.ice = [];    %total volume of ice [m3]
            snow.STATVAR.air = [];    % total volume of air [m3] - NOT USED
            snow.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            snow.STATVAR.hydraulicConductivity = []; %hydraulic conductivity of snow [m/sec]
            snow.STATVAR.albedo = []; %snow albedo [-]
            
            snow.STATVAR.d = []; %dendricity [-], range from 1 (dendric, i.e. crystal-shaped snow particles) to 0 (non-dendric, i.e. round snow particles) 
            snow.STATVAR.s = []; %sphericity [-], range from 1 (round snow particles) to 0 (elongated snow particles) 
            snow.STATVAR.gs = [];  % snow grain size [m]
            snow.STATVAR.time_snowfall = []; % average time of snowfall of a layers (i.e. grid cell) [Matlab time, days]
            snow.STATVAR.target_density = []; %ice fraction prior to  melt in diagnostic step [-]
            snow.STATVAR.excessWater = []; %water volume exceeding snow pore space [m3]
            
            snow.STATVAR.Lstar = []; %Obukhov length [m]
            snow.STATVAR.Qh = [];    %sensible heat flux [W/m2]
            snow.STATVAR.Qe = [];    %latent heat flux [W/m2]
        end
    
        function snow = provide_CONST(snow)
            
            snow.CONST.L_f = []; % volumetric latent heat of fusion, freezing
            
            snow.CONST.c_w = []; % volumetric heat capacity water
            snow.CONST.c_i = []; % volumetric heat capacity ice
            snow.CONST.c_o = []; % volumetric heat capacity organic
            snow.CONST.c_m = []; % volumetric heat capacity mineral
            
            snow.CONST.k_a = [];  % thermal conductivity air
            snow.CONST.k_w = [];  % thermal conductivity water
            snow.CONST.k_i = [];  % thermal conductivity ice 
            snow.CONST.k_o = [];  % thermal conductivity organic 
            snow.CONST.k_m = [];  % thermal conductivity mineral 
            
            snow.CONST.sigma = []; %Stefan-Boltzmann constant
            snow.CONST.kappa = []; % von Karman constant
            snow.CONST.L_s = [];  %latent heat of sublimation, evaporation handled in a dedicated function
            
            snow.CONST.cp = [];  % specific heat capacity at constant pressure of air
            snow.CONST.g = [];   % gravitational acceleration Earth surface
            
            snow.CONST.rho_w = []; % water density
            snow.CONST.rho_i = []; % ice density
        end
        
        
%         function snow = finalize_init(snow, tile) %assign all variables, that must be calculated or assigned otherwise for initialization
%             snow = finalize_init@SNOW_crocus_bucketW_seb(ground, tile);
%             
%         end
        
        %---time integration------
        %separate functions for CHILD pphase of snow cover
        
        function snow = get_boundary_condition_u(snow, tile) 
            
            snow.STATVAR.albedo4vegetation = snow.STATVAR.albedo;
            snow.STATVAR.emissivity4vegetation = snow.PARA.epsilon; %required for reflection of Lin
            snow.STATVAR.Lout4vegetation = snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+273.15).^4; %emitted part of Lout, reflected part in the vegetation
            
            %calculate SEB for canopy and adjust forcing data
            snow.VEGETATION = get_boundary_condition_u(snow.VEGETATION, tile);
                        
            forcing = snow.VEGETATION.ForcingV;
            snow.PARA.airT_height = snow.VEGETATION.STATVAR.mlcanopyinst.zs(1,2);

            %forcing = tile.FORCING;
            snow = get_boundary_condition_SNOW_u(snow, forcing);
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            
            snow = get_snow_properties_crocus(snow,forcing); %makes a TEMP variable newSnow that contains all information on the fresh snow - which is merged in the diagnostic step
            
            snow = surface_energy_balance(snow, forcing);
            snow = get_sublimation(snow, forcing);
            
            snow.TEMP.wind = tile.FORCING.TEMP.wind;
            snow.TEMP.wind_surface = forcing.TEMP.wind;
        end
        
        function snow = get_boundary_condition_u_CHILD(snow, tile)

            forcing = snow.PARENT.VEGETATION.ForcingV;
            snow = get_boundary_condition_allSNOW_rain_u(snow, forcing); %add full snow, but rain only for snow-covered part
            snow = get_boundary_condition_u_water_SNOW(snow, forcing);
            
            snow = get_snow_properties_crocus(snow,forcing); %makes a TEMP variable newSnow that contains all information on the fresh snow - which is merged in the diagnostic step
            
            snow = surface_energy_balance(snow, forcing); %this works including penetration of SW radiation through the CHILD snow
            snow = get_sublimation(snow, forcing);
            
            snow.TEMP.wind = tile.FORCING.TEMP.wind;
            snow.TEMP.wind_surface = forcing.TEMP.wind;
        end
        
        function snow = get_boundary_condition_u_create_CHILD(snow, tile)
            
            forcing = snow.PARENT.VEGETATION.ForcingV;
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
            snow.TEMP.wind = tile.FORCING.TEMP.wind;
            snow.TEMP.wind_surface = forcing.TEMP.wind;
            
            %start with  non-zero values for area and layerThick
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
        
        function snow = get_boundary_condition_l(snow, tile)
            snow.VEGETATION = get_boundary_condition_l(snow.VEGETATION, tile);
            snow = get_boundary_condition_l@SNOW_crocus_bucketW_seb(snow, tile);
        end
        
        function snow = get_derivatives_prognostic(snow, tile)
            snow.VEGETATION = get_derivatives_prognostic(snow.VEGETATION, tile);
            snow = get_derivatives_prognostic@SNOW_crocus_bucketW_seb(snow, tile);            
        end
        
        function snow = get_derivatives_prognostic_CHILD(snow, tile)
            snow = get_derivatives_prognostic_CHILD@SNOW_crocus_bucketW_seb(snow, tile);            
        end
        
        function timestep = get_timestep(snow, tile) 
            timestep_snow= get_timestep@SNOW_crocus_bucketW_seb(snow, tile);
            timestep_vegetation = get_timestep(snow.VEGETATION, tile);
            
            timestep = min(timestep_snow, timestep_vegetation);
        end
        
        function timestep = get_timestep_CHILD(snow, tile)  
            timestep = get_timestep_CHILD@SNOW_crocus_bucketW_seb(snow, tile);
        end
        
        function snow = advance_prognostic(snow, tile)
            snow.VEGETATION = advance_prognostic(snow.VEGETATION, tile);
            snow = advance_prognostic@SNOW_crocus_bucketW_seb(snow, tile);
        end
        
        
        function snow = advance_prognostic_CHILD(snow, tile)
            snow = advance_prognostic_CHILD@SNOW_crocus_bucketW_seb(snow, tile);
        end
        
        function snow = compute_diagnostic_first_cell(snow, tile)
            snow.VEGETATION = compute_diagnostic_first_cell(snow.VEGETATION, tile);
            snow = compute_diagnostic_first_cell@SNOW_crocus_bucketW_seb(snow, tile);
        end
       
        function snow = compute_diagnostic(snow, tile)
            snow.VEGETATION = compute_diagnostic(snow.VEGETATION, tile);
            snow = compute_diagnostic@SNOW_crocus_bucketW_seb(snow, tile);
        end
        
        function snow = compute_diagnostic_CHILD(snow, tile)
            snow = compute_diagnostic_CHILD@SNOW_crocus_bucketW_seb(snow, tile);
        end
        
        function snow = check_trigger(snow, tile)
            snow.VEGETATION = check_trigger(snow.VEGETATION, tile);
        
            %make SNOW a CHILD again
            if size(snow.STATVAR.layerThick,1) == 1 && snow.STATVAR.ice(1,1) ./ snow.STATVAR.area(1,1) < 0.5 .* snow.PARA.swe_per_cell
                
                ground = snow.NEXT;
                ground.VEGETATION.PARENT_SURFACE = ground;
                
                ground.PREVIOUS = snow.PREVIOUS; %reassign ground
                ground.PREVIOUS.NEXT = ground;
                ground.CHILD = snow;
                snow.PARENT = ground;
                ground.IA_CHILD = snow.IA_NEXT;
                ground.IA_CHILD.NEXT = ground;
                ground.IA_CHILD.PREVIOUS = snow;
                
                ground.IA_PREVIOUS=[]; %change to get_ia_class, if there is a possibility for another class on top of the snow cover
                
                %snow.NEXT =[];  %cut all dependencies, except for snow.NEXT which keeps being pointed to snow.PARENT, so that SW radiation can be transmitted
                snow.PREVIOUS =[];
                snow.IA_NEXT =[];
                snow.IA_PREVIOUS =[];
                
                %change to constant layerThick, variable area
                volume = snow.STATVAR.layerThick .* snow.STATVAR.area;
                snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ snow.STATVAR.target_density; %[m] constant layerThick
                snow.STATVAR.area = volume ./ snow.STATVAR.layerThick;
                
                snow = ground; %assign snow pointer to ground to return to regular stratigraphy
            end        
        
        end
        
        
        
        %-----non-mandatory functions-------
        
%         function snow = surface_energy_balance(snow, forcing)
%             snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
%             
%             [snow, S_up] = penetrate_SW(snow, snow.PARA.spectral_ranges .* forcing.TEMP.Sin .* snow.STATVAR.area(1)); %distribute SW radiation
%             snow.STATVAR.Sout = sum(S_up) ./ snow.STATVAR.area(1);
%             
%             snow.STATVAR.Qh = Q_h(snow, forcing);
%             snow.STATVAR.Qe = Q_eq_potET(snow, forcing);
%             
%             snow.TEMP.F_ub = (forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Qh - snow.STATVAR.Qe) .* snow.STATVAR.area(1);
%             snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.F_ub;
%         end
%         
%         
%         
%         function snow = conductivity(snow)
%             snow = conductivity_snow_Yen(snow);
%         end
%         
%         
%         %-----LATERAL-------------------
%         
%         %-----LAT_REMOVE_SURFACE_WATER-----
%         function snow = lateral_push_remove_surfaceWater(snow, lateral)
%             snow = lateral_push_remove_surfaceWater_simple(snow, lateral);
%         end
%         
%         %----LAT_SEEPAGE_FACE----------
%         function snow = lateral_push_remove_water_seepage(snow, lateral)
%             snow = lateral_push_remove_water_seepage_snow(snow, lateral);
%         end
%         
%         %----LAT_WATER_RESERVOIR------------
%         function snow = lateral_push_water_reservoir(snow, lateral)
%             snow = lateral_push_water_reservoir_snow(snow, lateral);
%         end
%         
%         %----LAT3D_WATER_UNCONFINED_AQUIFER------------
%         function snow = lateral3D_pull_water_unconfined_aquifer(snow, lateral)
%             snow = lateral3D_pull_water_unconfined_aquifer_snow(snow, lateral);
%         end
%         
%         function snow = lateral3D_push_water_unconfined_aquifer(snow, lateral)
%             snow = lateral3D_push_water_unconfined_aquifer_snow(snow, lateral);
%         end
%         
%         function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(snow, lateral)
%             [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_snow(snow, lateral);
%         end
%         
%         %LAT3D_WATER_RESERVOIR and LAT3D_WATER_SEEPAGE_FACE do not require specific functions
%         
%         %----LAT3D_HEAT------------
%         function snow = lateral3D_pull_heat(snow, lateral)
%             snow = lateral3D_pull_heat_simple(snow, lateral);
%         end
%         
%         function snow = lateral3D_push_heat(snow, lateral)
%             snow = lateral3D_push_heat_simple(snow, lateral);
%         end
%         
%         %----LAT3D_SNOW_CROCUS------------
%         function snow = lateral3D_pull_snow(snow, lateral)
%             snow = lateral3D_pull_snow_crocus(snow, lateral);
%         end
%         
%         function snow = lateral3D_push_snow(snow, lateral)
%             snow = lateral3D_push_snow_crocus(snow, lateral);
%         end
%         
%         
%         %----inherited Tier 1 functions ------------
%         
%         function snow = get_derivative_energy(snow)
%             snow = get_derivative_energy@HEAT_CONDUCTION(snow);
%         end
%         
%         function snow = conductivity_snow_Yen(snow)
%             snow = conductivity_snow_Yen@HEAT_CONDUCTION(snow);
%         end
%         
%         function flux = Q_h(snow, forcing)
%             flux = Q_h@SEB(snow, forcing);
%         end
%         
%         function flux = Q_eq_potET(snow, forcing)
%             flux = Q_eq_potET@SEB(snow, forcing);
%         end
%         
%         function timestep = get_timestep_heat_coduction(snow)
%             timestep = get_timestep_heat_coduction@HEAT_CONDUCTION(snow);
%         end
%         
%         function timestep = get_timestep_SNOW_mass_balance(snow)
%             timestep = get_timestep_SNOW_mass_balance@SNOW(snow);
%         end
%         
%         function timestep = get_timestep_SNOW_CHILD(snow)
%             timestep = get_timestep_SNOW_CHILD@SNOW(snow);
%         end
%         
%         function snow = L_star(snow, forcing)
%             snow = L_star@SEB(snow, forcing);
%         end
%         
%         function [snow, S_up] = penetrate_SW_transmission_spectral(snow, S_down)
%             [snow, S_up] = penetrate_SW_transmission_spectral@SEB(snow, S_down);
%         end
%         
%         function snow = get_T_water_freeW(snow)
%             snow = get_T_water_freeW@HEAT_CONDUCTION(snow);
%         end
%         
%         function snow = subtract_water(snow)
%             snow = subtract_water@SNOW(snow);
%         end
%         
%         function snow = subtract_water_CHILD(snow)
%             snow = subtract_water_CHILD@SNOW(snow);
%         end
%         
%         function snow = get_boundary_condition_SNOW_u(snow, forcing)
%             snow = get_boundary_condition_SNOW_u@SNOW(snow, forcing);
%         end
%         
%         function snow = get_boundary_condition_allSNOW_rain_u(snow, forcing)
%             snow = get_boundary_condition_allSNOW_rain_u@SNOW(snow, forcing);
%         end
%         
%         function snow = get_boundary_condition_allSNOW_u(snow, forcing)
%             snow = get_boundary_condition_allSNOW_u@SNOW(snow, forcing);
%         end
%         
%         function snow = initialize_zero_snow_BASE(snow)
%             snow = initialize_zero_snow_BASE@SNOW(snow);
%         end
%         
%         function snow = make_SNOW_CHILD(snow)
%             snow = make_SNOW_CHILD@SNOW(snow);
%         end
%         
%         function [snow, regridded_yesNo] = regrid_snow(snow, extensive_variables, intensive_variables, intensive_scaling_variable)
%             [snow, regridded_yesNo] = regrid_snow@REGRID(snow, extensive_variables, intensive_variables, intensive_scaling_variable);
%         end
    end
    
end
