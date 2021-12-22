%========================================================================
% CryoGrid GROUND class GROUND_freeW_bucketW_seb
% heat conduction, free water freeze curve, surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freeW_seb_verticalWall < SEB & HEAT_CONDUCTION & HEAT_FLUXES_LATERAL %& INITIALIZE

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = []; %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; %roughness length [m]
            
            %ground.PARA.rs = []; %surface resistance against evapotranspiration [secm] 
            ground.PARA.dt_max = []; %maximum possible timestep [sec]
            ground.PARA.dE_max = []; %maximum possible energy change per timestep [J/m3]
            
            ground.PARA.water_bucket_depth = [];
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.upperPos = []; % upper surface elevation [m]
            ground.STATVAR.lowerPos = []; % lower surface elevation [m]
            ground.STATVAR.layerThick = [];  % thickness of grid cells [m]
            
            ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = [];  % total volume of minerals [m3]
            ground.STATVAR.organic = [];  %total volume of ice [m3]
            ground.STATVAR.energy = [];  %total molar salt volume within a grid cell [mol]
            
            ground.STATVAR.T = [];   % temperature [degree C]
            ground.STATVAR.water = [];  % total volume of water [m3]
            ground.STATVAR.ice = [];    %total volume of ice [m3]
            ground.STATVAR.air = [];   % total volume of air [m3] - NOT USED
            ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
            
            ground.STATVAR.Lstar = [];  %Obukhov length [m]
            ground.STATVAR.Qh = []; %sensible heat flux [W/m2]
            ground.STATVAR.Qe = [];  % latent heat flux [W/m2]
            
            ground.STATVAR.water_bucket = 0;
            
        end
    
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = [];  % volumetric latent heat of fusion, freezing
            ground.CONST.c_w = [];  % volumetric heat capacity water
            ground.CONST.c_i = [];  % volumetric heat capacity ice
            ground.CONST.c_o = [];  % volumetric heat capacity organic
            ground.CONST.c_m = [];  % volumetric heat capacity mineral
            
            ground.CONST.k_a = [];   % thermal conductivity air
            ground.CONST.k_w = [];   % thermal conductivity water
            ground.CONST.k_i = [];   % thermal conductivity ice 
            ground.CONST.k_o = [];   % thermal conductivity organic 
            ground.CONST.k_m = [];   % thermal conductivity mineral 
            
            ground.CONST.sigma = []; %Stefan-Boltzmann constant
            ground.CONST.kappa = []; % von Karman constant
            ground.CONST.L_s = [];  %latent heat of sublimation, latent heat of evaporation handled in a dedicated function
            
            ground.CONST.cp = []; % specific heat capacity at constant pressure of air
            ground.CONST.g = []; % gravitational acceleration Earth surface
            
            ground.CONST.rho_w = [];   % water density
            ground.CONST.rho_i = [];   %ice density
        end
        
        function ground = convert_units(ground, tile)
                unit_converter = str2func(tile.PARA.unit_conversion_class);
                unit_converter = unit_converter();
                ground = convert_normal(unit_converter, ground, tile);
        end
        
        function ground = finalize_init(ground, tile)
%             ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
%             ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
%             ground.STATVAR.area = tile.PARA.area + ground.STATVAR.T .* 0;
            
            ground = get_E_freeW(ground);

            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.STATVAR.saturation = 0;
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            ground = surface_energy_balance_verticalWall(ground, forcing);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration 
            [ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
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
            forcing = tile.FORCING;
            %energy
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
            
            TForcing = ground.STATVAR.T(1);
            TForcing=TForcing+273.15;
            
            if TForcing>=273.15 %Evaporation and condensation is just applied for positive air temperatures

                L_w=1e3.*(2500.8 - 2.36.*(TForcing-273.15));  %latent heat of evaporation of water 
                rho_w = 1000; %density of the water
                
                ground.TEMP.rainfall = forcing.TEMP.rainfall ./1000 ./(24.*3600); %rainfall is in mm/day
                
                ground.TEMP.EC = - ground.STATVAR.Qe / (rho_w * L_w); %Calculation of flux that can be evaporated / condensated in this timestep
                               
                ground.STATVAR.water_bucket = ground.STATVAR.water_bucket + timestep * ground.TEMP.rainfall +  timestep * ground.TEMP.EC; %new height of water bucket after rain and evaporation / condensation

                if ground.STATVAR.water_bucket > ground.PARA.water_bucket_depth %more water is added then there is space in the bucket    
                    ground.STATVAR.water_bucket = ground.PARA.water_bucket_depth; %maximal the depth of the water bucket can be filled, surplus water is removed
                elseif ground.STATVAR.water_bucket < 0 %more water is removed than there is in the bucket
                    ground.STATVAR.water_bucket = 0; %maximal the depth of the water bucket can be removed, after that the bucket is empty
                end
                ground.STATVAR.saturation = ground.STATVAR.water_bucket / ground.PARA.water_bucket_depth;
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
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
        
        function ground = surface_energy_balance_verticalWall(ground, forcing)
            
            ground.STATVAR.Lin = forcing.TEMP.Lin;
            ground.STATVAR.Sin = forcing.TEMP.Sin;
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h_verticalWall(ground, forcing);
            ground.STATVAR.Qe = Q_eq_potET_verticalWall(ground, forcing);

            ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
            ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
            
        end
        
        function Q_e = Q_eq_potET_verticalWall(ground, forcing)
                
            uz = forcing.TEMP.wind;
            p = forcing.TEMP.p;
            q = forcing.TEMP.q;
            
            Tz = forcing.TEMP.Tair;     
            z =  ground.PARA.airT_height;
            
            z0 = ground.PARA.z0;   
            TForcing = ground.STATVAR.T(1);
            Lstar = ground.STATVAR.Lstar;   
            Tz=Tz+273.15;
            TForcing=TForcing+273.15;
            
            rho = p./(287.058.*Tz); %air density [kg m^(-3)]
            cp=1005;
            L_w=1e3.*(2500.8 - 2.36.*(TForcing-273.15));  %latent heat of evaporation of water
            L_i=1e3.*2834.1; %latent heat of sublimation
            kappa=0.4;
            g=9.81;
                    
            if isempty(ground.STATVAR.water_bucket) % only for first time step to set up saturation
                ground.STATVAR.water_bucket = ground.STATVAR.water(1) * ground.PARA.water_bucket_depth;
                ground.STATVAR.saturation = ground.STATVAR.water_bucket / ground.PARA.water_bucket_depth;
            end
                
            if TForcing<273.15 % for below-zero temperatures
                Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)).*(q-satPresIce(ground, TForcing)./p)./(log(z./z0));
                saturation = 0; %Only saturation is set to 0 for negative temperatures so that bucket water remains during frozen periods without getting changed
                Q_e = saturation * Q_e; %Always 0 for negative temperatures
            else
                Q_e = -rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)).*(q-satPresWater(ground, TForcing)./p)./(log(z./z0));
                if Q_e > 0 %in case of evaporation: Qe is reduced by multiplication of saturation
                    Q_e = ground.STATVAR.saturation * Q_e;
                elseif Q_e < 0 %in case of condensation: no reduction of Qe
                    Q_e = Q_e;
                end
            end
                
       end
            
       function Q_h = Q_h_verticalWall(ground, forcing)
                
            uz = forcing.TEMP.wind;
            z =  ground.PARA.airT_height;
            z0 = ground.PARA.z0;
            Tz = forcing.TEMP.Tair;
            TForcing = ground.STATVAR.T(1);
            Lstar = ground.STATVAR.Lstar;
            p = forcing.TEMP.p;
                
                
            Tz=Tz+273.15;
            TForcing=TForcing+273.15;
                
            rho = p./(287.058.*Tz); %air density [kg m^(-3)]
            cp=1005;
            kappa=0.4;
            g=9.81;
            sigma=5.67e-8;
                
            Q_h  =-rho.*cp.*kappa.* uz.*kappa./(log(z./z0)) .* (Tz-TForcing)./(log(z./z0));
                
        end
        
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
             
             %ground.PARA.options = [];
             ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic'  'T'};
             
             ground.PARA.default_value.albedo = {0.2};
             ground.PARA.comment.albedo = {'surface albedo [-]'};
             
             ground.PARA.default_value.epsilon = {0.99};
             ground.PARA.comment.epsilon = {'surface emissivity [-]'};
             
             ground.PARA.default_value.z0 = {0.01};
             ground.PARA.comment.z0 = {'roughness length [m]'};
             
             ground.PARA.default_value.water_bucket_depth = {0.01};
             ground.PARA.comment.water_bucket_depth = {'water holding capacity [m] of rock surface'}; 
             
             ground.PARA.default_value.dt_max = {3600};
             ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
             
             ground.PARA.default_value.dE_max = {50000};
             ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end

    end
    
end
