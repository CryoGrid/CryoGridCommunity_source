

classdef VEGETATION_ml_canopy < BASE

    properties
        ForcingV %forcing variables for ground and snow upper boundary cond

        PARENT_GROUND
        PARENT_SURFACE
        IA_VEGETATION_GROUND
        IA_VEGETATION_SURFACE
    end
    
    %         PREVIOUS
%         NEXT
%         IA_PREVIOUS
%         IA_NEXT

    methods 
        
        %mandatory functions for each class
        
        function ground = provide_PARA(ground)

            % Vegetation parameters
            ground.PARA.lai = [];
            ground.PARA.lai_winter = [];
            ground.PARA.sai = [];
            ground.PARA.ztop = [];
            ground.PARA.zref = [];
            ground.PARA.zref_old = [];
            ground.PARA.hksat = [];
            
            ground.PARA.dt_max = [] ; %[sec]
            ground.PARA.dE_max = []; %[J/m3]
        end
        
        function ground = provide_CONST(ground)
            
            ground.CONST.L_f = [];
            
            ground.CONST.c_w = [];
            ground.CONST.c_i = [];
            ground.CONST.c_o = [];
            ground.CONST.c_m = [];
            
            ground.CONST.k_a = [];       %air [Hillel(1982)]
            ground.CONST.k_w = [];        %water [Hillel(1982)]
            ground.CONST.k_i = [];         %ice [Hillel(1982)]
            ground.CONST.k_o = [];        %organic [Hillel(1982)]
            ground.CONST.k_m = [];
            
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
            
            %forcing variables need for snow and ground upper boundary
            ground.ForcingV.TEMP.tair = [];
            ground.ForcingV.TEMP.wind = [];
            ground.ForcingV.TEMP.Qh = []; % [m]
            
            ground.ForcingV.TEMP.Qe = []; % [m]
            ground.ForcingV.TEMP.Sin = []; % [m]
            ground.ForcingV.TEMP.Lin = []; % [m]
            ground.ForcingV.TEMP.p = [];  % [J/m2]
            
            ground.ForcingV.TEMP.snow_reservoir = 0;
            ground.ForcingV.TEMP.snowfall = [];  % [degree C]
            ground.ForcingV.TEMP.rainfall = [];  % [m]
            ground.ForcingV.TEMP.q = [];
            ground.ForcingV.TEMP.t = [];  % [m]
            
        end
        
        function ground = finalize_init(ground, tile)
            
            %SEBAS: is this needed?
%             ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
%             ground.PARA.airT_height = tile.FORCING.PARA.airT_height;

            ground = set_up_canopy_summer(ground);
            ground = set_up_canopy_winter(ground);            
            ground = set_up_canopy(ground);            
            %ground = finalize_STATVAR(ground, forcing); 

            % Set up the canopy structure
            % [vegetation] = SetUpCanopy();
            [ground] = initialize_mlcanopyinst(ground, tile.FORCING);
            [ground] = initialize_physcon(ground);
            [ground] = initialize_params(ground);
            [ground] = initialize_pftcon(ground);
            [ground] = initialize_leaf(ground);
            [ground] = initialize_soil(ground);
            [ground] = initialize_atmos(ground);
            [ground] = initialize_canopy(ground);
            [ground] = initialize_flux(ground);
            
            %build lookup tables here, otherwise the file will be read with each
            %iteration step - consumes alot of time
            [ground.STATVAR.vegetation.PsiLookup.dtLgridM,...
                ground.STATVAR.vegetation.PsiLookup.zdtgridM,...
                ground.STATVAR.vegetation.PsiLookup.psigridM,...
                ground.STATVAR.vegetation.PsiLookup.dtLgridH,...
                ground.STATVAR.vegetation.PsiLookup.zdtgridH,...
                ground.STATVAR.vegetation.PsiLookup.psigridH] = LookupPsihatINI();
            
            % ground.STATVAR.current_t = 0.0;
            
            %SEBAS: why is this not + 1/24 like later in the code  
            ground.STATVAR.execution_t = tile.FORCING.PARA.start_time + 0.5;
            
           % SEBAS: probably not needed?
%             ground.STATVAR.Lstar = -100;
%             ground.STATVAR.Qh = 0;
%             ground.STATVAR.Qe = 0;
            
        end
        

      
        
        
        %---------------------------------------
        
        
        function [ground] = get_boundary_condition_u(ground, tile) %functions specific for individual class, allow changing from Dirichlet to SEB
            
            forcing = tile.FORCING;
            
            ground = surface_energy_forest(ground, forcing);
            
            % Write forcing struct as the input for ground class under vegetation
            ground.ForcingV.TEMP.Tair = ground.STATVAR.vegetation.mlcanopyinst.tveg(1,2)-273.15;
            ground.ForcingV.TEMP.wind = ground.STATVAR.vegetation.mlcanopyinst.wind(1,2);
            ground.ForcingV.TEMP.wind_top = ground.STATVAR.vegetation.mlcanopyinst.wind(1,12); % wind at top of canopy (used for initial snow density calculation)
            ground.ForcingV.TEMP.Sin = ground.STATVAR.vegetation.flux.swdn(1,1) + ground.STATVAR.vegetation.flux.swdn(1,2); %ground.STATVAR.vegetation.flux.swsoi(1,1) + ground.STATVAR.vegetation.flux.swsoi(1,2); % vegetation.mlcanopyinst.sw_prof(1,2,1); %Canopy layer absorbed radiation
            ground.ForcingV.TEMP.Lin = ground.STATVAR.vegetation.flux.irdn; %ground.STATVAR.vegetation.flux.irsoi(1); %vegetation.flux.ir_source(2,1); %Longwave radiation emitted by bottom leaf layer (W/m2)
            ground.ForcingV.TEMP.p = forcing.TEMP.p; % air pressure at reference height (Pa)
            % is this really needed?
            ground.ForcingV.TEMP.snow_reservoir = ground.ForcingV.TEMP.snow_reservoir;

            ground.ForcingV.TEMP.snowfall = ground.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_snow .* (24 .*3600); % .* (24*3600); qflx_prec_grnd_snow (mm h2o/s) -> .* (24 .*3600) -> mm h2o/day
            ground.ForcingV.TEMP.rainfall = ground.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_rain .* (24 .*3600); % .* (24*3600);
            ground.ForcingV.TEMP.q = forcing.TEMP.q; %specific humidity at refernce height (kg/kg)
            ground.ForcingV.TEMP.t = forcing.TEMP.t; % time_snowfall

        end
        
        function ground = get_boundary_condition_l(ground, tile)
        %             %ground = get_boundary_condition_l@GROUND_base_class(ground, forcing);
        end
        
        
        function ground = get_derivatives_prognostic(ground, tile)
            %ground = get_derivatives_prognostic@GROUND_base_class(ground);
        end
        
        function timestep = get_timestep(ground, tile)  %could involve check for several state variables
            %timestep = get_timestep@GROUND_base_class(ground);
%             timestep = 24*60*60;
            %SEBAS: this makes more sense, since it corresponds to the main
            %timestep of the vegegation scheme
            timestep = 3600;
        end
        
        function ground = advance_prognostic(ground, tile) %real timestep derived as minimum of several classes in [sec] here!

        end

        function ground = compute_diagnostic_first_cell(ground, tile)
            %ground = L_star(ground, forcing);
        end
                
        function ground = compute_diagnostic(ground, tile)

        end
        
        function ground = check_trigger(ground, tile)
            
        end
        
        
        %non-mandatory functions
        
%         function ground = conductivity(ground)
%             ground = conductivity_mixing_squares(ground);
%         end
        
        function ground = surface_energy_forest(ground, forcing)
            
            if forcing.TEMP.t >= ground.STATVAR.execution_t || forcing.TEMP.t == forcing.PARA.start_time
                datestr(forcing.TEMP.t)
                
%                 disp(ground.STATVAR.vegetation.soilvar.transp_per_layer)
                
                %SEBAS: would make sense to make this directly
                %GROUND.STATVAR, since you only have vegetation now?
                %so vegetation = ground.STATVAR instead?
                vegetation = ground.STATVAR.vegetation;
                
                [vegetation] = set_up_forcing(vegetation, forcing);

                [vegetation] = canopy_fluxes_multilayer(vegetation);
                                
                %and ground.STATVAR = vegetation here?
                ground.STATVAR.vegetation = vegetation;
                
%                 disp(ground.STATVAR.vegetation.soilvar.transp_per_layer)
                
                ground.STATVAR.execution_t = ground.STATVAR.execution_t + 1/24;
       
            end
        end
        
        %called by GROUND, could be moved to IA class, but I guess these
        %functions could well be needed for any combination of GROUND and
        %vegetation classes
        
        function vegetation = get_transpiration(vegetation)
            ground = vegetation.PARENT_GROUND;
           
            transpiration = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;

            transp = transpiration(1) .* ground.STATVAR.area(1,1);
            ground.TEMP.d_water_ET(1,1) = ground.TEMP.d_water_ET(1,1) - transp;
            ground.TEMP.d_water_ET_energy(1,1) =  ground.TEMP.d_water_ET_energy(1,1) - transp  .* (double(ground.STATVAR.T(1,1)>=0) .* ground.CONST.c_w .* ground.STATVAR.T(1,1) + ...
                double(ground.STATVAR.T(1,1)<0) .* ground.CONST.c_i .* ground.STATVAR.T(1,1));
            
            weighting = ground.STATVAR.water(2:3) ./ sum(ground.STATVAR.water(2:3));
            transp = transpiration(2).* weighting .* ground.STATVAR.area(2:3,1);
            ground.TEMP.d_water_ET(2:3,1) = ground.TEMP.d_water_ET(2:3,1) - transp;
            ground.TEMP.d_water_ET_energy(2:3,1) =  ground.TEMP.d_water_ET_energy(2:3,1) - transp  .* (double(ground.STATVAR.T(2:3,1)>=0) .* ground.CONST.c_w .* ground.STATVAR.T(2:3,1) + ...
                double(ground.STATVAR.T(2:3,1)<0) .* ground.CONST.c_i .* ground.STATVAR.T(2:3,1));
            
            weighting = ground.STATVAR.water(4:9) ./ sum(ground.STATVAR.water(4:9));
            transp = transpiration(3).* weighting .* ground.STATVAR.area(4:9,1);
            ground.TEMP.d_water_ET(4:9,1) = ground.TEMP.d_water_ET(4:9,1) - transp;
            ground.TEMP.d_water_ET_energy(4:9,1) =  ground.TEMP.d_water_ET_energy(4:9,1) - transp  .* (double(ground.STATVAR.T(4:9,1)>=0) .* ground.CONST.c_w .* ground.STATVAR.T(4:9,1) + ...
                double(ground.STATVAR.T(4:9,1)<0) .* ground.CONST.c_i .* ground.STATVAR.T(4:9,1));            

            
            
        end

        function vegetation = get_evaporation(vegetation)
            ground = vegetation.PARENT_GROUND;
            
            evaporation = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
            
            evap = evaporation .* ground.STATVAR.area(1,1);
            ground.TEMP.d_water_ET(1,1) = ground.TEMP.d_water_ET(1,1) - evap;
            ground.TEMP.d_water_ET_energy(1,1) =  ground.TEMP.d_water_ET_energy(1,1) - evap .* (double(ground.STATVAR.T(1,1)>=0) .* ground.CONST.c_w .* ground.STATVAR.T(1,1) + ...
                double(ground.STATVAR.T(1,1)<0) .* ground.CONST.c_i .* ground.STATVAR.T(1,1));
        end
        
        %----this function is only needed for this particular vegegtation
        %scheme
        function vegetation = map_variables_no_snow(vegetation)
            ground = vegetation.PARENT_SURFACE;
           
            ground.STATVAR.Qh = vegetation.STATVAR.vegetation.mlcanopyinst.shsoi;
            ground.STATVAR.Qe = vegetation.STATVAR.vegetation.mlcanopyinst.lhsoi;
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) + vegetation.STATVAR.vegetation.mlcanopyinst.gsoi .* ground.STATVAR.area(1,1);
            
            vegetation.STATVAR.vegetation.soilvar.t_top_surfacecell = ground.STATVAR.T(1) + 273.15;
            vegetation.STATVAR.vegetation.soilvar.dz_topsurfacecell = ground.STATVAR.layerThick(1);
            vegetation.STATVAR.vegetation.soilvar.thk_topsurfacecell = ground.STATVAR.thermCond(1);
            
            
            ground = vegetation.PARENT_GROUND;
%             midPoint = cumsum([0; ground.STATVAR.layerThick]);
%             midPoint = (midPoint(2:end,1) + midPoint(1:end-1,1))./2;
%             water = ground.STATVAR.water./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
%             ice = ground.STATVAR.ice./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
%             vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(midPoint,water,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
%             vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(midPoint,water,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
%             vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(midPoint,water,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
%             vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(midPoint,ice,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
%             vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(midPoint,ice,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
%             vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(midPoint,ice,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
%             
%             vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
%             vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
%             vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
%             
%             vegetation.STATVAR.vegetation.soilvar.thk(1) = interp1(midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
%             vegetation.STATVAR.vegetation.soilvar.thk(2) = interp1(midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
%             vegetation.STATVAR.vegetation.soilvar.thk(3) = interp1(midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');

            range_cell_1 = 1;
            range_cell_2 = 2:3;
            range_cell_3 = 4:9;
            vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = sum(ground.STATVAR.water(range_cell_1))./ sum(ground.STATVAR.area(range_cell_1) .* ground.STATVAR.layerThick(range_cell_1));
            vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = sum(ground.STATVAR.water(range_cell_2))./ sum(ground.STATVAR.area(range_cell_2) .* ground.STATVAR.layerThick(range_cell_2));
            vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = sum(ground.STATVAR.water(range_cell_3))./ sum(ground.STATVAR.area(range_cell_3) .* ground.STATVAR.layerThick(range_cell_3));
            vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = sum(ground.STATVAR.ice(range_cell_1))./ sum(ground.STATVAR.area(range_cell_1) .* ground.STATVAR.layerThick(range_cell_1));
            vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = sum(ground.STATVAR.ice(range_cell_2))./ sum(ground.STATVAR.area(range_cell_2) .* ground.STATVAR.layerThick(range_cell_2));
            vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = sum(ground.STATVAR.ice(range_cell_3))./ sum(ground.STATVAR.area(range_cell_3) .* ground.STATVAR.layerThick(range_cell_3));
            
            vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = mean(ground.STATVAR.T(range_cell_1))+273.15;
            vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = mean(ground.STATVAR.T(range_cell_2))+273.15;
            vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = mean(ground.STATVAR.T(range_cell_3))+273.15;
            
            vegetation.STATVAR.vegetation.soilvar.thk(1) = mean(ground.STATVAR.thermCond(range_cell_1));
            vegetation.STATVAR.vegetation.soilvar.thk(2) = mean(ground.STATVAR.thermCond(range_cell_2));
            vegetation.STATVAR.vegetation.soilvar.thk(3) = mean(ground.STATVAR.thermCond(range_cell_3));
            
        end
        
    end  
end
