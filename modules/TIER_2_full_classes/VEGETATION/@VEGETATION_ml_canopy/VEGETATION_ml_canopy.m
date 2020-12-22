

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
            
            %ground.PARA.dt_max = [] ; %[sec]
            %ground.PARA.dE_max = []; %[J/m3]
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

%             ground.STATVAR.upperPos = [];
%             ground.STATVAR.lowerPos = [];
%             ground.STATVAR.layerThick = []; % [m]
%             
%             ground.STATVAR.waterIce = []; % [m]
%             ground.STATVAR.mineral = []; % [m]
%             ground.STATVAR.organic = []; % [m]
%             ground.STATVAR.energy = [];  % [J/m2]
%             
%             ground.STATVAR.T = [];  % [degree C]
%             ground.STATVAR.water = [];  % [m]
%             ground.STATVAR.ice = [];
%             ground.STATVAR.air = [];  % [m]
%             ground.STATVAR.thermCond = [];
            
            %forcing variables need for snow and ground upper boundary
            %ground.ForcingV.TEMP.tair = [];
            ground.ForcingV.TEMP.wind = [];
            ground.ForcingV.TEMP.Qh = []; % [m]
            
            ground.ForcingV.TEMP.Qe = []; % [m]
            ground.ForcingV.TEMP.Sin = []; % [m]
            ground.ForcingV.TEMP.Lin = []; % [m]
            ground.ForcingV.TEMP.p = [];  % [J/m2]
            
            %ground.ForcingV.TEMP.snow_reservoir = 0;
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
            
            %re-order 
            ground.STATVAR = ground.STATVAR.vegetation;
            
            %SEBAS: why is this not + 1/24 like later in the code - CHANGED
            ground.STATVAR.execution_t = tile.FORCING.PARA.start_time + 1./24;
            
            ground.STATVAR.PARENT_GROUND = ground.PARENT_GROUND;
            ground.STATVAR.PARENT_SURFACE = ground.PARENT_SURFACE;

            
        end
        

      
        
        
        %---------------------------------------
        
        
        function [ground] = get_boundary_condition_u(ground, tile) %functions specific for individual class, allow changing from Dirichlet to SEB
            
            forcing = tile.FORCING;
            
            ground = surface_energy_forest(ground, forcing);
            
            %CHECK!!!
            % Write forcing struct as the input for ground class under vegetation
            %this should be tair(1,2)!!!
            ground.ForcingV.TEMP.Tair = ground.STATVAR.mlcanopyinst.tair(1,2)-273.15;
            ground.ForcingV.TEMP.wind = ground.STATVAR.mlcanopyinst.wind(1,2);
            ground.ForcingV.TEMP.wind_top = ground.STATVAR.mlcanopyinst.wind(1,12); % wind at top of canopy (used for initial snow density calculation)
            ground.ForcingV.TEMP.Sin = ground.STATVAR.flux.swdn(1,1) + ground.STATVAR.flux.swdn(1,2); %ground.STATVAR.vegetation.flux.swsoi(1,1) + ground.STATVAR.vegetation.flux.swsoi(1,2); % vegetation.mlcanopyinst.sw_prof(1,2,1); %Canopy layer absorbed radiation
            ground.ForcingV.TEMP.Lin = ground.STATVAR.flux.irdn; %ground.STATVAR.vegetation.flux.irsoi(1); %vegetation.flux.ir_source(2,1); %Longwave radiation emitted by bottom leaf layer (W/m2)
            ground.ForcingV.TEMP.p = forcing.TEMP.p; % air pressure at reference height (Pa)
            % is this really needed?
            %ground.ForcingV.TEMP.snow_reservoir = ground.ForcingV.TEMP.snow_reservoir;

            ground.ForcingV.TEMP.snowfall = ground.STATVAR.mlcanopyinst.qflx_prec_grnd_snow .* (24 .*3600); % .* (24*3600); qflx_prec_grnd_snow (mm h2o/s) -> .* (24 .*3600) -> mm h2o/day
            ground.ForcingV.TEMP.rainfall = ground.STATVAR.mlcanopyinst.qflx_prec_grnd_rain .* (24 .*3600); % .* (24*3600);
            
            %CHECK!!!
            %this should be eair(1,2)!!! Will be needed as soon as as we
            %compute evaporation in the main model
            ground.ForcingV.TEMP.q = 0.622 .* ground.STATVAR.mlcanopyinst.eair(1,2)./forcing.TEMP.p; %specific humidity  (kg/kg) - 0.622 is the conversion from mol/mol to kg/kg!!!
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
            %SIMONE: Add change between summer and winter properties here!!!
        end
        
        
        %non-mandatory functions
        
%         function ground = conductivity(ground)
%             ground = conductivity_mixing_squares(ground);
%         end
        
        function ground = surface_energy_forest(ground, forcing)
            
            if forcing.TEMP.t >= ground.STATVAR.execution_t || forcing.TEMP.t == forcing.PARA.start_time
                disp(datestr(forcing.TEMP.t))

                
                vegetation = ground.STATVAR;
                vegetation.PARENT_SURFACE = ground.PARENT_SURFACE;
                vegetation.PARENT_GROUND = ground.PARENT_GROUND;
                
                %all could eventually become accumulated /time-averaged quantity

                
                % set parameters from GROUND models
%                 vegetation.mlcanopyinst.tg = ground.PARENT_SURFACE.STATVAR.T(1,1) + 273.15; %ground suurface temperature
%                 vegetation.soilvar.t_top_surfacecell = ground.PARENT_SURFACE.STATVAR.T(1,1) + 273.15; % kind of the same

%                 range_cell_1 = 1;
%                 range_cell_2 = 2:3;
%                 range_cell_3 = 4:9;
%                 %vegetation.soilvar.rootfr = [0.05 0.5 0.45];
%                 vegetation.soilvar.hk = [mean(ground.PARENT_GROUND.STATVAR.hydraulicConductivity(range_cell_1,1)) ...
%                     mean(ground.PARENT_GROUND.STATVAR.hydraulicConductivity(range_cell_2,1)) mean(ground.PARENT_GROUND.STATVAR.hydraulicConductivity(range_cell_3,1))];
%                 vegetation.soilvar.soil_water_matric_potential = [mean(ground.PARENT_GROUND.STATVAR.waterPotential(range_cell_1,1)) ...
%                     mean(ground.PARENT_GROUND.STATVAR.waterPotential(range_cell_2,1)) mean(ground.PARENT_GROUND.STATVAR.waterPotential(range_cell_3,1))];
% 
%                 vegetation.soilvar.h2osoi_vol(1) = sum(ground.PARENT_GROUND.STATVAR.water(range_cell_1))./ sum(ground.PARENT_GROUND.STATVAR.area(range_cell_1) .* ground.PARENT_GROUND.STATVAR.layerThick(range_cell_1));
%                 vegetation.soilvar.h2osoi_vol(2) = sum(ground.PARENT_GROUND.STATVAR.water(range_cell_2))./ sum(ground.PARENT_GROUND.STATVAR.area(range_cell_2) .* ground.PARENT_GROUND.STATVAR.layerThick(range_cell_2));
%                 vegetation.soilvar.h2osoi_vol(3) = sum(ground.PARENT_GROUND.STATVAR.water(range_cell_3))./ sum(ground.PARENT_GROUND.STATVAR.area(range_cell_3) .* ground.PARENT_GROUND.STATVAR.layerThick(range_cell_3));
%                 
% %                 vegetation.soilvar.soil_water_matric_potential = vegetation.soilvar.soil_water_matric_potential .* (double(vegetation.soilvar.h2osoi_vol < 0.1) .* (vegetation.soilvar.h2osoi_vol ./ 0.1).^-6 + ...
% %                     double(vegetation.soilvar.h2osoi_vol >= 0.1));
%                 
%                 vegetation.soilvar.h2osoi_ice(1) = sum(ground.PARENT_GROUND.STATVAR.ice(range_cell_1))./ sum(ground.PARENT_GROUND.STATVAR.area(range_cell_1) .* ground.PARENT_GROUND.STATVAR.layerThick(range_cell_1));
%                 vegetation.soilvar.h2osoi_ice(2) = sum(ground.PARENT_GROUND.STATVAR.ice(range_cell_2))./ sum(ground.PARENT_GROUND.STATVAR.area(range_cell_2) .* ground.PARENT_GROUND.STATVAR.layerThick(range_cell_2));
%                 vegetation.soilvar.h2osoi_ice(3) = sum(ground.PARENT_GROUND.STATVAR.ice(range_cell_3))./ sum(ground.PARENT_GROUND.STATVAR.area(range_cell_3) .* ground.PARENT_GROUND.STATVAR.layerThick(range_cell_3));
%                 
%                 vegetation.soilvar.t_soisno(1) = mean(ground.PARENT_GROUND.STATVAR.T(range_cell_1))+273.15;
%                 vegetation.soilvar.t_soisno(2) = mean(ground.PARENT_GROUND.STATVAR.T(range_cell_2))+273.15;
%                 vegetation.soilvar.t_soisno(3) = mean(ground.PARENT_GROUND.STATVAR.T(range_cell_3))+273.15;
%                 
%                 vegetation.soilvar.thk(1) = mean(ground.PARENT_GROUND.STATVAR.thermCond(range_cell_1));
%                 vegetation.soilvar.thk(2) = mean(ground.PARENT_GROUND.STATVAR.thermCond(range_cell_2));
%                 vegetation.soilvar.thk(3) = mean(ground.PARENT_GROUND.STATVAR.thermCond(range_cell_3));
                
                
                vegetation = set_up_forcing(vegetation, forcing);
                
                vegetation = assign_GROUND_fluxes(ground, vegetation);

                [vegetation] = canopy_fluxes_multilayer(vegetation);
                                
                vegetation = get_transpiration(ground, vegetation);
                
                %and ground.STATVAR = vegetation here?
                ground.STATVAR = vegetation;
                
%                 disp(ground.STATVAR.vegetation.soilvar.transp_per_layer)
                
                ground.STATVAR.execution_t = ground.STATVAR.execution_t + 1/24;
       
            end
        end
        
        function vegetation = assign_GROUND_fluxes(ground, vegetation)
            
            vegetation.flux.albsoib = [vegetation.PARENT_SURFACE.STATVAR.albedo4vegetation ground.PARENT_SURFACE.STATVAR.albedo4vegetation];
            vegetation.flux.albsoid = [vegetation.PARENT_SURFACE.STATVAR.albedo4vegetation ground.PARENT_SURFACE.STATVAR.albedo4vegetation];
            vegetation.flux.emissivity_ground = vegetation.PARENT_SURFACE.STATVAR.emissivity4vegetation;
            vegetation.flux.Qe_ground = vegetation.PARENT_SURFACE.STATVAR.Qe;
            vegetation.flux.Qh_ground = vegetation.PARENT_SURFACE.STATVAR.Qh;
            vegetation.flux.Lout_ground = vegetation.PARENT_SURFACE.STATVAR.Lout4vegetation;
            range = [1:vegetation.soilvar.nsoi]';
            vegetation.soilvar.h2osoi_ice = ground.PARENT_GROUND.STATVAR.ice(range)' ./ (ground.PARENT_GROUND.STATVAR.area(range)' .* ground.PARENT_GROUND.STATVAR.layerThick(range)');
            vegetation.soilvar.h2osoi_vol = ground.PARENT_GROUND.STATVAR.water(range)' ./ (ground.PARENT_GROUND.STATVAR.area(range)' .* ground.PARENT_GROUND.STATVAR.layerThick(range)');
        end

        
        
        %called by GROUND, could be moved to IA class, but I guess these
        %functions could well be needed for any combination of GROUND and
        %vegetation classes
        
        
        
        function vegetation = get_transpiration(ground, vegetation)
            
            
            %convert from mol/sec to me3 per sec using molar mass 0.018
            %kg/mol of water and the density
            %transpiration = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;

            frac_sun_shade = cat(3, vegetation.flux.fracsun, vegetation.flux.fracsha);
            leaf_et = vegetation.mlcanopyinst.trleaf .* frac_sun_shade;
            leaf_et = sum(leaf_et,3); %mol water per sec per m2 leaf per m2 ground 
            leaf_et = leaf_et .* vegetation.canopy.dlai; %mol water per sec per m2 ground per canopy layer
            leaf_et = 0.0181528e-3 .* sum(leaf_et); %in m3/m2 water/sec
            vegetation.mlcanopyinst.transpiration = leaf_et .* vegetation.mlcanopyinst.soil_et_loss; %weight with fraction assigned for each layer
            

%             transp = transpiration(1) .* ground.STATVAR.area(1,1);
%             trans = min(transp, ground.STATVAR.water(1,1) /(3600.*2));
%             ground.TEMP.d_water_ET(1,1) = ground.TEMP.d_water_ET(1,1) - transp;
%             ground.TEMP.d_water_ET_energy(1,1) =  ground.TEMP.d_water_ET_energy(1,1) - transp  .* (double(ground.STATVAR.T(1,1)>=0) .* ground.CONST.c_w .* ground.STATVAR.T(1,1) + ...
%                 double(ground.STATVAR.T(1,1)<0) .* ground.CONST.c_i .* ground.STATVAR.T(1,1));
%                         
%             
%             weighting = ground.STATVAR.water(2:3) ./ sum(ground.STATVAR.water(2:3));
%             transp = transpiration(2).* weighting .* ground.STATVAR.area(2:3,1);
%             transp = min(transp, sum(ground.STATVAR.water(2:3)) ./(3600.*2));
%             ground.TEMP.d_water_ET(2:3,1) = ground.TEMP.d_water_ET(2:3,1) - transp;
%             ground.TEMP.d_water_ET_energy(2:3,1) =  ground.TEMP.d_water_ET_energy(2:3,1) - transp  .* (double(ground.STATVAR.T(2:3,1)>=0) .* ground.CONST.c_w .* ground.STATVAR.T(2:3,1) + ...
%                 double(ground.STATVAR.T(2:3,1)<0) .* ground.CONST.c_i .* ground.STATVAR.T(2:3,1));
%             
%             weighting = ground.STATVAR.water(4:9) ./ sum(ground.STATVAR.water(4:9));
%             transp = transpiration(3).* weighting .* ground.STATVAR.area(4:9,1);
%             transp = min(transp, sum(ground.STATVAR.water(4:9)) ./(3600.*2));
%             ground.TEMP.d_water_ET(4:9,1) = ground.TEMP.d_water_ET(4:9,1) - transp;
%             ground.TEMP.d_water_ET_energy(4:9,1) =  ground.TEMP.d_water_ET_energy(4:9,1) - transp  .* (double(ground.STATVAR.T(4:9,1)>=0) .* ground.CONST.c_w .* ground.STATVAR.T(4:9,1) + ...
%                 double(ground.STATVAR.T(4:9,1)<0) .* ground.CONST.c_i .* ground.STATVAR.T(4:9,1));            


        end

%         function vegetation = get_evaporation(vegetation)
%             ground = vegetation.PARENT_GROUND;
%             
%             evaporation = (0.0181528 .* vegetation.STATVAR.mlcanopyinst.etsoi)./1000;
%             
%             evap = evaporation .* ground.STATVAR.area(1,1);
%             evap = min(evap, ground.STATVAR.water(1,1)./(3600.*2));
%             ground.TEMP.d_water_ET(1,1) = ground.TEMP.d_water_ET(1,1) - evap;
%             ground.TEMP.d_water_ET_energy(1,1) =  ground.TEMP.d_water_ET_energy(1,1) - evap .* (double(ground.STATVAR.T(1,1)>=0) .* ground.CONST.c_w .* ground.STATVAR.T(1,1) + ...
%                 double(ground.STATVAR.T(1,1)<0) .* ground.CONST.c_i .* ground.STATVAR.T(1,1));
%         end
        
        %----this function is only needed for this particular vegegtation
        %scheme
%         function vegetation = map_variables_no_snow(vegetation)
%             ground = vegetation.PARENT_SURFACE;
%            
%             ground.STATVAR.Qh = vegetation.STATVAR.mlcanopyinst.shsoi;
%             ground.STATVAR.Qe = vegetation.STATVAR.mlcanopyinst.lhsoi;
%             
%             d_energy_first_cell = vegetation.ForcingV.TEMP.Lin + vegetation.ForcingV.TEMP.Sin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
%             ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) + d_energy_first_cell .* ground.STATVAR.area(1,1);
%             
%             vegetation.STATVAR.soilvar.t_top_surfacecell = ground.STATVAR.T(1) + 273.15;
%             %vegetation.STATVAR.vegetation.soilvar.dz_topsurfacecell = ground.STATVAR.layerThick(1);
%             %vegetation.STATVAR.vegetation.soilvar.thk_topsurfacecell = ground.STATVAR.thermCond(1);
%             %not necessary, should not be used any further in vegetation
%             vegetation.STATVAR.mlcanopyinst.gsoi = d_energy_first_cell;
%             
% %             ground = vegetation.PARENT_GROUND;
% % %             midPoint = cumsum([0; ground.STATVAR.layerThick]);
% % %             midPoint = (midPoint(2:end,1) + midPoint(1:end-1,1))./2;
% % %             water = ground.STATVAR.water./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
% % %             ice = ground.STATVAR.ice./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
% % %             vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(midPoint,water,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %             vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(midPoint,water,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %             vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(midPoint,water,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %             vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(midPoint,ice,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %             vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(midPoint,ice,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %             vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(midPoint,ice,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %             
% % %             vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %             vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %             vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %             
% % %             vegetation.STATVAR.vegetation.soilvar.thk(1) = interp1(midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %             vegetation.STATVAR.vegetation.soilvar.thk(2) = interp1(midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %             vegetation.STATVAR.vegetation.soilvar.thk(3) = interp1(midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % 
% %             range_cell_1 = 1;
% %             range_cell_2 = 2:3;
% %             range_cell_3 = 4:9;
% %             vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = sum(ground.STATVAR.water(range_cell_1))./ sum(ground.STATVAR.area(range_cell_1) .* ground.STATVAR.layerThick(range_cell_1));
% %             vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = sum(ground.STATVAR.water(range_cell_2))./ sum(ground.STATVAR.area(range_cell_2) .* ground.STATVAR.layerThick(range_cell_2));
% %             vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = sum(ground.STATVAR.water(range_cell_3))./ sum(ground.STATVAR.area(range_cell_3) .* ground.STATVAR.layerThick(range_cell_3));
% %             vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = sum(ground.STATVAR.ice(range_cell_1))./ sum(ground.STATVAR.area(range_cell_1) .* ground.STATVAR.layerThick(range_cell_1));
% %             vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = sum(ground.STATVAR.ice(range_cell_2))./ sum(ground.STATVAR.area(range_cell_2) .* ground.STATVAR.layerThick(range_cell_2));
% %             vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = sum(ground.STATVAR.ice(range_cell_3))./ sum(ground.STATVAR.area(range_cell_3) .* ground.STATVAR.layerThick(range_cell_3));
% %             
% %             vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = mean(ground.STATVAR.T(range_cell_1))+273.15;
% %             vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = mean(ground.STATVAR.T(range_cell_2))+273.15;
% %             vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = mean(ground.STATVAR.T(range_cell_3))+273.15;
% %             
% %             vegetation.STATVAR.vegetation.soilvar.thk(1) = mean(ground.STATVAR.thermCond(range_cell_1));
% %             vegetation.STATVAR.vegetation.soilvar.thk(2) = mean(ground.STATVAR.thermCond(range_cell_2));
% %             vegetation.STATVAR.vegetation.soilvar.thk(3) = mean(ground.STATVAR.thermCond(range_cell_3));
% %             
%         end
        
    end  
end
