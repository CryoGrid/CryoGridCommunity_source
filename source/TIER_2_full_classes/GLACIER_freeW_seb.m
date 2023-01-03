%========================================================================
% CryoGrid GLACIER CLASS class GLACIER_freeW_seb
% heat conduction, free water freeze curve, surface
% energy balance
% L. Schmidt, S. Westermann, November 2021
%========================================================================

classdef GLACIER_freeW_seb < SEB & HEAT_CONDUCTION & WATER_FLUXES & HEAT_FLUXES_LATERAL & WATER_FLUXES_LATERAL 

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        

        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = [];
            ground.PARA.epsilon = [];
            ground.PARA.z0 = []; %roughness length [m]
            
            ground.PARA.dt_max = [];  %maximum possible timestep [sec]
            ground.PARA.dE_max = [];  %maximum possible energy change per timestep [J/m3]

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

            ground.STATVAR.excessWater = 0;  %water volume overtopping first grid cell (i.e. surface water [m3]
            
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
        
        function ground = convert_units(ground, tile)
                ground.STATVAR.mineral = ground.STATVAR.T .* 0;
                ground.STATVAR.organic = ground.STATVAR.T .* 0;
                
                unit_converter = str2func(tile.PARA.unit_conversion_class);
                unit_converter = unit_converter();
                ground = convert_normal(unit_converter, ground, tile);
        end

        function ground = finalize_init(ground, tile) 
%             ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
%             ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
%             ground.STATVAR.area = tile.PARA.area + ground.STATVAR.T .* 0;
            
            ground.STATVAR.mineral = ground.STATVAR.T .* 0;
            ground.STATVAR.organic = ground.STATVAR.T .* 0;
            
            ground = get_E_freeW(ground);
            
            ground.PARA.target_layerThick = ground.STATVAR.layerThick;

            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.STATVAR.runoff = 0;
            ground.STATVAR.smb = 0; 

            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.STATVAR.sublimation = 0;
            ground.TEMP.sublimation_energy = 0; 
        end
        
        function ground = finalize_init2(ground, tile)

            ground = get_E_freeW(ground);
            ground = conductivity(ground);

        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile) 
            forcing = tile.FORCING;
           
            ground = surface_energy_balance(ground, forcing);
            ground.TEMP.F_water = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600;  %possibly add water from external source here 
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            
            ground.TEMP.d_water(1,1) = ground.TEMP.d_water(1,1) + ground.TEMP.F_water .* ground.STATVAR.area(1,1);
            ground.TEMP.d_water_energy(1,1) = ground.TEMP.d_water_energy(1,1) + ground.TEMP.F_water .* ground.TEMP.T_rainWater .* ground.CONST.c_w .* ground.STATVAR.area(1,1);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            forcing = tile.FORCING;
            ground.TEMP.F_water = 0; 
            ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
            ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivative_energy(ground);
            
        end
        
        function timestep = get_timestep(ground, tile)
 %at maximum timestep maximum half the cell is melted under melting conditions, otherwise normal heat conduction timestep
            melt_fraction = 0.9; 
            timestep = min(double(ground.TEMP.d_energy>0 & ground.STATVAR.energy >= - 0.95.* ground.STATVAR.ice.*ground.CONST.L_f) .* melt_fraction .*  (-ground.STATVAR.energy ./ ground.TEMP.d_energy) + ...
                double(ground.TEMP.d_energy<0 & ground.STATVAR.energy >= - 0.95.* ground.STATVAR.ice.*ground.CONST.L_f) .*melt_fraction .*  ((-ground.STATVAR.ice.*ground.CONST.L_f -ground.STATVAR.energy) ./ ground.TEMP.d_energy) + ...
                double( ground.STATVAR.energy < - 0.95.* ground.STATVAR.ice.*ground.CONST.L_f) .* ground.CONST.c_i .* (ground.STATVAR.ice./ground.STATVAR.layerThick./ ground.STATVAR.area)  ./ (abs(ground.TEMP.d_energy) ./ ground.STATVAR.layerThick./ ground.STATVAR.area));
            timestep(isnan(timestep)) = ground.PARA.dt_max;
            timestep = min(timestep,ground.PARA.dt_max); 
        end
        
        function ground = advance_prognostic(ground, tile) 
            timestep = tile.timestep;
            
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
%             %ground.STATVAR.energy(1) = ground.STATVAR.energy(1) + timestep .* ground.TEMP.sublimation_energy; %add sublimation energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_water_energy; %add energy from water advection
            %water
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + timestep .* ground.TEMP.d_water; 
            ground.STATVAR.layerThick = ground.STATVAR.layerThick + timestep .* ground.TEMP.d_water ./ ground.STATVAR.area(1,1); 
            %ground.STATVAR.waterIce(1) = ground.STATVAR.waterIce(1) + timestep .* ground.TEMP.sublimation;
            ground.STATVAR.layerThick(1) = ground.STATVAR.layerThick(1) + timestep .* ground.STATVAR.sublimation ./ ground.STATVAR.area(1,1);
            %smb
            %ground.STATVAR.sublimation = ground.STATVAR.sublimation + timestep .* ground.TEMP.sublimation;
            ground.STATVAR.smb = ground.STATVAR.smb +  timestep.*ground.TEMP.F_water .* ground.STATVAR.area(1,1);% + timestep .* ground.TEMP.sublimation; 
            
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)     

            ground = get_T_water_freeW(ground);
            
            ground.STATVAR.runoff = ground.STATVAR.runoff + ground.STATVAR.water(1);
            ground.STATVAR.smb = ground.STATVAR.smb - ground.STATVAR.water(1);
            ground.STATVAR.waterIce(1) = ground.STATVAR.ice(1);
%             ground.STATVAR.waterIce = max(0, ground.STATVAR.waterIce);
            ground.STATVAR.water(1) = ground.STATVAR.water(1) .* 0;
            
            
            ground.STATVAR.layerThick = (ground.STATVAR.waterIce + ground.STATVAR.organic + ground.STATVAR.mineral) ./ ground.STATVAR.area;
            
            ground = modify_grid(ground);

            ground = get_T_water_freeW(ground);
                        
            ground = conductivity(ground);
            
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_energy = ground.STATVAR.energy.*0;
            ground.TEMP.d_water_ET_energy = ground.STATVAR.energy.*0;
            ground.TEMP.sublimation_energy = 0; 
            ground.STATVAR.sublimation = 0; 
        end
        
        function ground = check_trigger(ground, tile)
            %do nothing
        end
        
        
        %-----non-mandatory functions-------
        function ground = surface_energy_balance(ground, forcing) %calculates the different fluxes of the surface energy balance and adds them up to get the upper boundary energy flux
            
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq_potET(ground, forcing);
            
            ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;

            ground.STATVAR.sublimation = -ground.STATVAR.Qe ./(ground.CONST.rho_w .* ground.CONST.L_s) .* ground.STATVAR.area(1);
            ground.TEMP.sublimation_energy =  ground.STATVAR.sublimation .* (ground.STATVAR.T(1) .* ground.CONST.c_i - ground.CONST.L_f);

        end
        
        
        function ground = modify_grid(ground)
          % divide mass onto predetermined grid of constant ice fraction (+- 10 %)

          if ground.STATVAR.ice(1) < 0.9*ground.PARA.target_layerThick(1) % move up mass
                    sf = (ground.PARA.target_layerThick(1) - ground.STATVAR.waterIce(1))./ground.STATVAR.waterIce;

                    ground.STATVAR.ice = ground.STATVAR.ice + [ground.STATVAR.ice(2:end).*sf(2:end); 0] - [0; ground.STATVAR.ice(2:end-1).*sf(2:end-1); 0];
                    ground.STATVAR.water = ground.STATVAR.water + [ground.STATVAR.water(2:end).*sf(2:end); 0] - [0; ground.STATVAR.water(2:end-1).*sf(2:end-1); 0];

                    ground.STATVAR.waterIce = ground.STATVAR.waterIce + [ground.STATVAR.waterIce(2:end).*sf(2:end); 0] - [0; ground.STATVAR.waterIce(2:end-1).*sf(2:end-1); 0];
                    ground.STATVAR.energy = ground.STATVAR.energy + [ground.STATVAR.energy(2:end).*sf(2:end); 0] - [0; ground.STATVAR.energy(2:end-1).*sf(2:end-1); 0]; 
                    ground.STATVAR.layerThick = ground.STATVAR.layerThick + [ground.STATVAR.layerThick(2:end).*sf(2:end); 0] - [0; ground.STATVAR.layerThick(2:end-1).*sf(2:end-1); 0];

                    ground.STATVAR.mineral = ground.STATVAR.mineral + [ground.STATVAR.mineral(2:end).*sf(2:end); 0] - [0; ground.STATVAR.mineral(2:end-1).*sf(2:end-1);0];
                    ground.STATVAR.organic = ground.STATVAR.organic + [ground.STATVAR.organic(2:end).*sf(2:end); 0] - [0; ground.STATVAR.organic(2:end-1).*sf(2:end-1); 0]; 
                   % ground.STATVAR.air = ground.STATVAR.air + [ground.STATVAR.air(2:end).*sf(2:end); 0] - [0; ground.STATVAR.air(2:end).*sf(2:end)]; 

          end

          if ground.STATVAR.layerThick(1) > 1.3*ground.PARA.target_layerThick(1) %move mass down
              if ground.STATVAR.ice(1) > 1.3*ground.PARA.target_layerThick(1)
                  for i=1:length(ground.STATVAR.ice)-1
                      split_fraction = ground.STATVAR.ice(i) ./ (ground.PARA.target_layerThick(i) .* ground.STATVAR.area(1));
                      sf1 = (split_fraction-1)./split_fraction;
                      sf2 = 1./split_fraction;
                      
                      if sf1 < 0 || sf1 > 1
                      else
                          ground.STATVAR.waterIce(i+1)=sf1.*ground.STATVAR.waterIce(i)+ground.STATVAR.waterIce(i+1);
                          ground.STATVAR.waterIce(i) = sf2.*ground.STATVAR.waterIce(i);
                          ground.STATVAR.energy(i+1) = sf1.*ground.STATVAR.energy(i)+ground.STATVAR.energy(i+1);
                          ground.STATVAR.energy(i) = sf2.*ground.STATVAR.energy(i);
                          ground.STATVAR.layerThick(i+1) = sf1.*ground.STATVAR.layerThick(i)+ground.STATVAR.layerThick(i+1);
                          ground.STATVAR.layerThick(i) = sf2.*ground.STATVAR.layerThick(i);

                          ground.STATVAR.mineral(i+1) = sf1.*ground.STATVAR.mineral(i)+ground.STATVAR.mineral(i+1);
                          ground.STATVAR.mineral(i) = sf2.*ground.STATVAR.mineral(i);
                          ground.STATVAR.organic(i+1) = sf1.*ground.STATVAR.organic(i)+ground.STATVAR.organic(i+1);
                          ground.STATVAR.organic(i) = sf2.*ground.STATVAR.organic(i);
 %                         ground.STATVAR.air(i+1) = sf1.*ground.STATVAR.air(i)+ground.STATVAR.air(i+1);
%                           ground.STATVAR.air(i) = sf2.*ground.STATVAR.air(i);
                      end
                  end
                end
                
          end

          
        end
        
        
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE(ground);
            %ground = provide_PARA(ground);
            
            ground.PARA.class_category = 'GLACIER';
            
            %ground.PARA.options = [];
            ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'T'};
            
            ground.PARA.default_value.albedo = {0.6};
            ground.PARA.comment.albedo = {'surface albedo [-]'};
            
            ground.PARA.default_value.epsilon = {0.99};
            ground.PARA.comment.epsilon = {'surface emissivity [-]'};
            
            ground.PARA.default_value.z0 = {0.01};
            ground.PARA.comment.z0 = {'roughness length [m]'};
            
            ground.PARA.default_value.dt_max = {3600};
            ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
            
            ground.PARA.default_value.dE_max = {50000};
            ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end
        
        %-----LATERAL-------------------
        
%         %-----LAT_REMOVE_SURFACE_WATER-----
%         function ground = lateral_push_remove_surfaceWater(ground, lateral)
%             ground = lateral_push_remove_surfaceWater_simple(ground, lateral);
%         end
%         
%         %-----LAT_REMOVE_SUBSURFACE_WATER-----
%         function ground = lateral_push_remove_subsurfaceWater(ground, lateral)
%             ground = lateral_push_remove_subsurfaceWater_simple(ground, lateral);
%         end
%         
%         %----LAT_SEEPAGE_FACE----------
%         function ground = lateral_push_remove_water_seepage(ground, lateral)
%             ground = lateral_push_remove_water_seepage_simple(ground, lateral);
%         end
% 
%         %----LAT_WATER_RESERVOIR------------           
%         function ground = lateral_push_water_reservoir(ground, lateral)
%             ground = lateral_push_water_reservoir_simple(ground, lateral);
%         end
%         
%         %----LAT3D_WATER_UNCONFINED_AQUIFER------------          
%         function ground = lateral3D_pull_water_unconfined_aquifer(ground, lateral)
%             ground = lateral3D_pull_water_unconfined_aquifer_simple(ground, lateral);
%         end
%         
%         function ground = lateral3D_push_water_unconfined_aquifer(ground, lateral)
%             ground = lateral3D_push_water_unconfined_aquifer_simple(ground, lateral);
%         end
%         
%         function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground, lateral)
%             [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_simple(ground, lateral);
%         end
%         
%         %LAT3D_WATER_RESERVOIR and LAT3D_WATER_SEEPAGE_FACE do not require specific functions
%         
%         %-------LAT3D_HEAT-------------
%         function ground = lateral3D_pull_heat(ground, lateral)
%             ground = lateral3D_pull_heat_simple(ground, lateral);
%         end
%         
%         function ground = lateral3D_push_heat(ground, lateral)
%             ground = lateral3D_push_heat_simple(ground, lateral);
%         end
%         
%         
%         %----inherited Tier 1 functions ------------
%         
%         function ground = get_derivative_energy(ground)
%            ground = get_derivative_energy@HEAT_CONDUCTION(ground); 
%         end
%         
%         function ground = conductivity_mixing_squares(ground)
%             ground = conductivity_mixing_squares@HEAT_CONDUCTION(ground);
%         end
%         
%         function flux = Q_h(ground, forcing)
%            flux = Q_h@SEB(ground, forcing);
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
%            ground = get_boundary_condition_u_water2@WATER_FLUXES(ground, forcing);
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
%            ground = L_star@SEB(ground, forcing); 
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
