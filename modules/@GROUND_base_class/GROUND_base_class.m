%base class for a GROUND object with a free water freeze curve - upper
%boundary condition is not specified; %superclass that cannot be run alone

classdef GROUND_base_class < matlab.mixin.Copyable
    properties
        CONST %constants
        PARA %external service parameters, all other
        STATVAR  %energy, water content, etc.
        TEMP  %derivatives in prognostic timestep and optimal timestep
        PREVIOUS
        NEXT
        IA_PREVIOUS
        IA_NEXT
    end
    
    
    methods
        
        %mandatory functions for each class
        
        
        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_PARA(ground);
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        
        function variable = initialize_from_file(ground, variable, section)
            class_variables = fieldnames(variable);
            for j=1:size(class_variables,1)
                for k=1:size(section,1)
                    if strcmp(section{k,1}, class_variables(j,1))
                        variable.(class_variables{j}) = section{k,2};
                    end
                end
            end
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground.PARA.heatFlux_lb = forcing.PARA.heatFlux_lb;
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            variables = fieldnames(ground.STATVAR);
            range = (grid.MIDPOINTS > depths(1,1) & grid.MIDPOINTS <= depths(1,2));
            ground.STATVAR.layerThick = grid.LAYERTHICK(range,1);
            ground.STATVAR.upperPos = forcing.PARA.altitude - depths(1,1);
            ground.STATVAR.lowerPos = forcing.PARA.altitude - depths(1,2);
            ground.STATVAR.gridlength = grid.length(range,1);
            ground.STATVAR.midpoints_delta = grid.MIDPOINTS_delta(range,1);
            ground.STATVAR.grid = grid.GRID(range,1);
            %             ground.STATVAR.fieldCapacity (range) = ground.STATVAR.fieldCapacity (range) + soil_profile (i,7); % NC added Field C
            
            for j=1:size(variables,1)
                for i=1:size(grid.variable_names,2)
                    if strcmp(variables{j,1}, grid.variable_names{1,i})
                        ground.STATVAR.(variables{j,1}) = grid.variable_gridded(range,i);
                    end
                end
            end
            ground = finalize_STATVAR(ground); %assign all variables, that must be calculated or assigned otherwise, including energy, water and ice contents, thermal conductivity
            ground = initializeSoilThermalProperties1(ground, grid); %NC added
            ground = modify_infiltration(ground, grid); %NC added
            ground = initializeConductivityCapacity(ground, grid, forcing);
            %---- energy and water balance initialization -----------------------------
            ground = initializeBALANCE(ground);
        end
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            %empty here
            
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_F_lb(ground, forcing);
        end
        
        
        function ground = get_derivatives_prognostic(ground)
            ground = get_derivative_energy(ground);
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
            timestep = ground.PARA.dE_max ./ (max(abs(ground.TEMP.d_energy) ./ ground.STATVAR.layerThick));
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            %empty here
        end
        
        function ground = compute_diagnostic(ground, forcing)
            ground = getThermalPropertiesInfiltration(ground, forcing);
            ground = updateBALANCE(ground);
            %------ surface energy balance module ---------------------------------
            %set surface conditions (albedo, roughness length, etc.)
            ground = surfaceCondition(ground, forcing);
            %calculate the surface energy balance
            ground = surfaceEnergyBalanceInfiltration(ground, forcing);
            %------ soil module  --------------------------------------------------
            %calculate heat conduction
            ground = heatConduction(ground, forcing);
            ground = get_T_water(ground);
            ground = conductivity(ground);
            
        end
        
        function ground = compute_diagnostic_water(ground, forcing)
            ground = bucketSch(ground, forcing);
            
        end
        %non-mandatory functions -> required here so that they are usable
        %in subclasses
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        function ground = get_derivative_energy(ground)
            fluxes = (ground.STATVAR.T(1:end-1) - ground.STATVAR.T(2:end)) .* ground.STATVAR.thermCond(1:end-1) .* ground.STATVAR.thermCond(2:end) ./...
                (ground.STATVAR.thermCond(1:end-1).* ground.STATVAR.layerThick(2:end)./2 +  ground.STATVAR.thermCond(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            
            d_energy=ground.STATVAR.energy.*0;
            
            d_energy(1) = ground.TEMP.F_ub - fluxes(1);
            d_energy(2:end-1) = fluxes(1:end-1) - fluxes(2:end);
            d_energy(end) = ground.TEMP.F_lb + fluxes(end);
            
            ground.TEMP.d_energy = d_energy;
        end
        
        function ground = get_T_water(ground)
            
            Lf = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            
            E_frozen = -Lf.*ground.STATVAR.waterIce;
            
            ground.STATVAR.T = double(ground.STATVAR.energy < E_frozen) .* (ground.STATVAR.energy - E_frozen) ./ (c_i.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic) + ...
                double(ground.STATVAR.energy >0) .* ground.STATVAR.energy ./ (c_i.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic);
            ground.STATVAR.ice = double(ground.STATVAR.energy <= E_frozen) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > E_frozen & ground.STATVAR.energy < 0) .* ground.STATVAR.energy ./ (-Lf);
            ground.STATVAR.water = double(ground.STATVAR.energy >= 0) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > - Lf.*ground.STATVAR.waterIce & ground.STATVAR.energy < 0) .* (ground.STATVAR.energy + Lf.*ground.STATVAR.waterIce) ./ Lf;
            
        end
        
        %---- modification for infiltration
        function ground = modify_infiltration(ground, grid)
            ground.STATVAR.wc=ground.STATVAR.waterIce;
            ground.STATVAR.E_lb = find(ground.PARA.evaporationDepth==grid.GRID(:,1))-1;
            ground.STATVAR.T_lb= find(ground.PARA.rootDepth==grid.GRID(:,1));-1;
        end
        
        function [c_temp, k_eff, lwc_temp]=readThermalParameters(ground)
            
            T = ground.STATVAR.T;
            cT_frozen = ground.STATVAR.cT_frozen;
            cT_thawed = ground.STATVAR.cT_thawed;
            capacity = ground.STATVAR.capacity;
            K_frozen = ground.STATVAR.K_frozen;
            K_thawed = ground.STATVAR.K_thawed;
            conductivity = ground.STATVAR.conductivity;
            arraySizeT = ground.PARA.arraySize;
            liquidWaterContent = ground.STATVAR.liquidWaterContent;
            
            
            a=(T-cT_frozen)./(cT_thawed-cT_frozen)*(arraySizeT-2)+1; %T and c information live on same grid
            posT=round((a<=1).*(-a+1)+a+(a>arraySizeT-1).*(arraySizeT-a));
            posT(posT==0)=1;
            posT(isnan(posT))=arraySizeT;
            
            
            r=[1:size(capacity,1)]';
            c=posT;
            indices=(size(capacity,1))*(c-1)+r;
            c_temp=capacity(indices);
            lwc_temp=liquidWaterContent(indices);
            
            a=(T-K_frozen)./(K_thawed-K_frozen)*(arraySizeT-2)+1;
            posT=round((a<=1).*(-a+1)+a+(a>arraySizeT-1).*(arraySizeT-a));
            posT(posT==0)=1;
            posT(isnan(posT))=arraySizeT;
            
            
            r=[1:size(conductivity,1)]';
            c=posT;
            indices=(size(conductivity,1))*(c-1)+r;
            k_eff=conductivity(indices);
        end
        
        function c_temp = capacityUnfrozen(ground)
            
            c_w = ground.CONST.c_w; % 4.2*10^6; %[J/m�K]
            c_o = ground.CONST.c_o; % 2.5*10^6; %[J/m�K]
            c_m = ground.CONST.c_m; % 2*10^6; %[J/m�K]
            
            c_temp = (ground.STATVAR.mineral.*c_m + ground.STATVAR.organic.*c_o + ground.STATVAR.wc.*c_w);
            
            
            % adjust for free water
            freeWater_domain = ground.STATVAR.mineral+ground.STATVAR.organic<1e-6; % cells without soil matrix material
            c_temp(freeWater_domain) = c_w; % assume pure water for cells which consist partly of air and water
            
            % adjust for low soil matrix %% NC - actPor to field Capacity
            lowMinOrg_domain = ground.STATVAR.mineral+ground.STATVAR.organic>=1e-6 & ~ground.STATVAR.excessGroundIce & ( ground.STATVAR.field_capacity > ground.STATVAR.natPor );   % cells with lower soil matrix material than 1-natPor
            
            if sum(lowMinOrg_domain)>0
                
                water = ground.STATVAR.wc;
                mineral = ground.STATVAR.mineral;
                organic = ground.STATVAR.organic;
                matrix = mineral + organic;
                mineral(lowMinOrg_domain) = mineral(lowMinOrg_domain) ./ matrix(lowMinOrg_domain) .* (1 - GRID.soil.cT_natPor(lowMinOrg_domain)) ;
                organic(lowMinOrg_domain) = organic(lowMinOrg_domain) ./ matrix(lowMinOrg_domain) .* (1 - GRID.soil.cT_natPor(lowMinOrg_domain)) ;
                water(lowMinOrg_domain) = min( water(lowMinOrg_domain), GRID.soil.cT_natPor(lowMinOrg_domain) ) ;
                c_temp(lowMinOrg_domain) = water(lowMinOrg_domain) .* c_w + mineral(lowMinOrg_domain) .* c_m + organic(lowMinOrg_domain).* c_o;
            end
            
            
        end
        
        
        function k_temp = conductivityUnfrozen(ground)
            
            ka = ground.CONST.k_a; %0.025;       %air [Hillel(1982)]
            kw = ground.CONST.k_w; %0.57;        %water [Hillel(1982)]
            ko = ground.CONST.k_o; %0.25;        %organic [Hillel(1982)]
            km = ground.CONST.k_m; %soil.kh_bedrock;     %mineral
            
            k_freeWater = 5.0; % to ensure fast heat transfer
            
            air=1-ground.STATVAR.wc-ground.STATVAR.mineral-ground.STATVAR.organic;
            
            k_temp = (ground.STATVAR.wc.* kw.^0.5 + ground.STATVAR.mineral.* km.^0.5 + ground.STATVAR.organic.* ko.^0.5 + air.* ka.^0.5).^2 ;
            
            % adjust for free water
            freeWater_domain = ground.STATVAR.mineral+ground.STATVAR.organic<1e-6; % cells without soil matrix material
            k_temp(freeWater_domain) = k_freeWater; % assume pure water for cells which consist partly of air and water
            
            % adjust for low soil matrix
            lowMinOrg_domain = ground.STATVAR.mineral+ground.STATVAR.organic>=1e-6 & ~ground.STATVAR.excessGroundIce & ( ground.STATVAR.field_capacity > ground.STATVAR.natPor );   % cells with lower soil matrix material than 1-natPor
            
            if sum(lowMinOrg_domain)>0
                
                water = ground.STATVAR.wc;
                mineral = ground.STATVAR.mineral;
                organic = ground.STATVAR.organic;
                matrix = mineral + organic;
                mineral(lowMinOrg_domain) = mineral(lowMinOrg_domain) ./ matrix(lowMinOrg_domain) .* (1 - GRID.soil.cT_natPor(lowMinOrg_domain)) ;
                organic(lowMinOrg_domain) = organic(lowMinOrg_domain) ./ matrix(lowMinOrg_domain) .* (1 - GRID.soil.cT_natPor(lowMinOrg_domain)) ;
                water(lowMinOrg_domain) = min( water(lowMinOrg_domain), GRID.soil.cT_natPor(lowMinOrg_domain) ) ;
                air(lowMinOrg_domain) = 1-mineral(lowMinOrg_domain)-organic(lowMinOrg_domain)-water(lowMinOrg_domain);
                k_temp(lowMinOrg_domain) = (water(lowMinOrg_domain) .* kw.^0.5 + mineral(lowMinOrg_domain) .* km.^0.5 + organic(lowMinOrg_domain).* ko.^0.5 + air(lowMinOrg_domain) .* ka.^0.5).^2 ;
            end
            
        end
        
        function  ground = initializeConductivityCapacity(ground, grid, forcing) %T, wc, GRID, PARA
            
            %             [c_temp, k_temp, k_eff, lwc_temp];
            c_temp = zeros(size(grid.MIDPOINTS));
            k_temp = zeros(size(grid.MIDPOINTS));
            k_eff = zeros(size(grid.GRID));
            lwc_temp = zeros(size(grid.MIDPOINTS));
            
            
            %------- unused grid cells --------------------------------
            c_temp_air = grid.air_MIDPOINTS .*ground.CONST.c_a;    %set some value e.g. air
            k_temp_air = grid.air_MIDPOINTS .*ground.CONST.k_a;    %set some value e.g. air
            lwc_temp_air = grid.air_MIDPOINTS*0;
            
            %------- soil domain --------------------------------------
            [c_temp,k_temp,lwc_temp] = readThermalParameters(ground);
            
            %--adjust for the unfrozen part of the domain
            ground.STATVAR.c_temp= double(ground.STATVAR.T<=0).*c_temp + double(ground.STATVAR.T>0).* capacityUnfrozen(ground);
            ground.STATVAR.k_temp = double(ground.STATVAR.T<=0).*k_temp + double(ground.STATVAR.T>0).* conductivityUnfrozen(ground);
            ground.STATVAR.lwc_temp= double(ground.STATVAR.T<=0).*lwc_temp + double(ground.STATVAR.T>0).* ground.STATVAR.wc;
            
            
            %------- snow domain --------------------------------------
            %             c_temp(GRID.snow.cT_domain) = cap_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
            %                 GRID.snow.Snow_w(GRID.snow.cT_domain),...
            %                 GRID.snow.Snow_a(GRID.snow.cT_domain),...
            %                 PARA);
            %
            %             k_temp(GRID.snow.cT_domain) = cond_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
            %                 GRID.snow.Snow_w(GRID.snow.cT_domain),...
            %                 GRID.snow.Snow_a(GRID.snow.cT_domain));
            %
            %             lwc_temp(GRID.snow.cT_domain) = GRID.snow.Snow_w(GRID.snow.cT_domain)./GRID.general.K_delta(GRID.snow.cT_domain);
            %
        end
        
        function ground = initializeBALANCE(ground)
            
            % ENERGY balance
            % energy content soil domain in [J/m^2] distinguished by sensible and latent part
            ground.PARA.energy.E_soil_sens = nansum(  ( ground.CONST.c_w .* ground.STATVAR.liquidWaterContent + ground.CONST.c_i .* ....
                (ground.STATVAR.wc- ground.STATVAR.liquidWaterContent) + ground.CONST.c_m .* ground.STATVAR.mineral + ground.CONST.c_o .* ground.STATVAR.organic ) .*ground.STATVAR.T ...
                .* ground.STATVAR.layerThick );%changed K-delta to layerThick
            ground.PARA.energy.E_soil_lat = nansum( ground.CONST.rho_w .* ground.CONST.L_sl .* ground.STATVAR.liquidWaterContent .* ground.STATVAR.layerThick  );
            ground.PARA.energy.E_soil = ground.PARA.energy.E_soil_sens + ground.PARA.energy.E_soil_lat;
            % energy content snow domain in [J/m^2] distinguished by sensible and latent part
            ground.PARA.energy.E_snow_sens = 0; %nansum( c_cTgrid(GRID.snow.cT_domain) .* T(GRID.snow.cT_domain) .* GRID.general.K_delta(GRID.snow.cT_domain) );
            ground.PARA.energy.E_snow_lat = 0; %nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc_cTgrid(GRID.snow.cT_domain).* GRID.general.K_delta(GRID.snow.cT_domain) );
            ground.PARA.energy.E_snow = ground.PARA.energy.E_snow_sens + ground.PARA.energy.E_snow_lat;
            
            % accumulated changes per output timestep
            ground.PARA.energy.dE_soil_sens = 0;
            ground.PARA.energy.dE_soil_lat = 0;
            ground.PARA.energy.dE_soil = 0;
            ground.PARA.energy.dE_snow_sens = 0;
            ground.PARA.energy.dE_snow_lat = 0;
            ground.PARA.energy.dE_snow = 0;
            
            ground.STATVAR.Q_lateral = zeros( length(ground.STATVAR.midpoints) , 1 );
            
            
            % WATER ground
            % water content soil domain in [m]
            ground.PARA.water.W_soil = nansum( ground.STATVAR.wc .* ground.STATVAR.layerThick );
            % water content snow domain in [m]
            ground.PARA.water.W_snow = 0;%nansum( GRID.snow.Snow_i + GRID.snow.Snow_w ) + GRID.snow.SWEinitial;
            % accumulated changes per output timestep
            % storage
            ground.PARA.water.dW_soil = 0;
            ground.PARA.water.dW_snow = 0;
            % precipitation
            ground.PARA.water.dp_rain=0;
            ground.PARA.water.dp_snow=0; % SWE
            % evapotranspiration and sublimation
            ground.PARA.water.de=0;
            ground.PARA.water.ds=0;
            % runoff
            ground.PARA.water.dr_surface=0;
            ground.PARA.water.dr_external=0;
            ground.PARA.water.dr_snowmelt=0;
            ground.PARA.water.dr_excessSnow=0;
            ground.PARA.water.dr_lateralSnow=0;
            ground.PARA.water.dr_rain=0;  % this is only rain on frozen ground
            ground.PARA.water.dr_lateralWater=0;
            ground.PARA.water.dr_DarcyReservoir=0; % When worker is connected to a Darcy_reservoir as a boundary condition
            ground.PARA.water.dr_lateralExcess=0; % excess water when applying lateral fluxes
            % mismatch
            ground.PARA.water.dm_lacking=0;
            
            ground.PARA.water.dr_water_fluxes_out=zeros(numlabs,numlabs);
        end
        
        function ground = getThermalPropertiesInfiltration(ground, forcing)
            %[c_temp, k_temp, k_eff, lwc_temp]
            air.cT_domain = ones(80,1); %%change later ******
            %------- unused grid cells --------------------------------------------
            c_temp(air.cT_domain) = ground.CONST.c_a;
            k_temp(air.cT_domain) = ground.CONST.k_a;
            lwc_temp(air.cT_domain) = 0;
            
            %             %------- soil domain --------------------------------------------------
            [c_temp,k_temp,lwc_temp] = readThermalParameters(ground);
            
            %--adjust for the unfrozen part of the domain
            ground.STATVAR.c_temp = double(ground.STATVAR.T<=0).*c_temp + double(ground.STATVAR.T>0).* capacityUnfrozen(ground);
            ground.STATVAR.k_temp = double(ground.STATVAR.T<=0).*k_temp + double(ground.STATVAR.T>0).* conductivityUnfrozen(ground);
            ground.STATVAR.lwc_temp= double(ground.STATVAR.T<=0).*lwc_temp + double(ground.STATVAR.T>0).* ground.STATVAR.wc;
            
            
            %-------- set higher conductivity for free water ----------------------
            % now done in conductivityUnfrozen
            
            
            %             %------- snow domain --------------------------------------------------
            %             c_temp(GRID.snow.cT_domain) = cap_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
            %                 GRID.snow.Snow_w(GRID.snow.cT_domain),...
            %                 GRID.snow.Snow_a(GRID.snow.cT_domain),...
            %                 PARA);
            %
            %             k_temp(GRID.snow.cT_domain) = cond_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
            %                 GRID.snow.Snow_w(GRID.snow.cT_domain),...
            %                 GRID.snow.Snow_a(GRID.snow.cT_domain));
            %
            %             lwc_temp(GRID.snow.cT_domain) = GRID.snow.Snow_w(GRID.snow.cT_domain)./GRID.general.K_delta(GRID.snow.cT_domain);
            %
            %------- interpolate conductivity to K-grid ---------------------------
            % ground.STATVAR.k_eff(2:end-1)
            k_eff = ground.STATVAR.k_eff;
            
            k_eff = ground.STATVAR.layerThick(1:end-1)./(2.*ground.STATVAR.midpoints_delta(1:end-1)) .* (1./ground.STATVAR.k_temp(1:end-1)).^2 ...
                + ground.STATVAR.layerThick(2:end)  ./(2.*ground.STATVAR.midpoints_delta(1:end-1)) .* (1./ground.STATVAR.k_temp(2:end)).^2;
            
            k_eff(2:end-1) = k_eff(2:end-1).^(-0.5);
            %
            k_eff(1)     = ground.STATVAR.k_temp(1);
            ground.STATVAR.k_eff = k_eff;
            %                         k_eff(end)   = k_temp(end);
            %
            %                         %------ correct upper most value below air-domain ---------------------
            %                         k_eff(GRID.air.K_domain_lb+1) = k_temp(GRID.air.K_domain_lb+1);
        end
        
        
        
        function ground = updateBALANCE(ground)
            
            
            % at this point the thermal and hydrological state of the soil and snow is calculated
            % energy content soil domain in [J/m^2]
            % distinguished by sensible and latent part
            E_soil_sens_old = ground.PARA.energy.E_soil_sens;
            E_soil_lat_old =  ground.PARA.energy.E_soil_lat;
            E_soil_old =  ground.PARA.energy.E_soil;
            
            
            
            % ENERGY balance
            % energy content soil domain in [J/m^2] distinguished by sensible and latent part
            ground.PARA.energy.E_soil_sens = nansum(  ( ground.CONST.c_w .* ground.STATVAR.liquidWaterContent + ground.CONST.c_i .* ....
                (ground.STATVAR.wc- ground.STATVAR.liquidWaterContent) + ground.CONST.c_m .* ground.STATVAR.mineral + ground.CONST.c_o .* ground.STATVAR.organic ) .*ground.STATVAR.T ...
                .* ground.STATVAR.layerThick  );%changed K-delta to layerThick
            ground.PARA.energy.E_soil_lat = nansum( ground.CONST.rho_w .* ground.CONST.L_sl .* ground.STATVAR.liquidWaterContent .* ground.STATVAR.layerThick  );
            ground.PARA.energy.E_soil = ground.PARA.energy.E_soil_sens + ground.PARA.energy.E_soil_lat;
            
            
            ground.PARA.energy.dE_soil_sens = ground.PARA.energy.dE_soil_sens + ground.PARA.energy.E_soil_sens - E_soil_sens_old;
            ground.PARA.energy.dE_soil_lat = ground.PARA.energy.dE_soil_lat + ground.PARA.energy.E_soil_lat - E_soil_lat_old;
            ground.PARA.energy.dE_soil = ground.PARA.energy.dE_soil + ground.PARA.energy.E_soil - E_soil_old;
            
            
            
            % energy content snow domain in [J/m^2] distinguished by sensible and latent part
            ground.PARA.energy.E_snow_sens = 0; %nansum( c_cTgrid(GRID.snow.cT_domain) .* T(GRID.snow.cT_domain) .* GRID.general.K_delta(GRID.snow.cT_domain) );
            ground.PARA.energy.E_snow_lat = 0; %nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc_cTgrid(GRID.snow.cT_domain).* GRID.general.K_delta(GRID.snow.cT_domain) );
            ground.PARA.energy.E_snow = ground.PARA.energy.E_snow_sens + ground.PARA.energy.E_snow_lat;
            
            %             ground.PARA.energy.dE_snow_sens = ground.PARA.energy.dE_snow_sens + ground.PARA.energy.E_snow_sens - E_snow_sens_old;
            %             ground.PARA.energy.dE_snow_lat = ground.PARA.energy.dE_snow_lat + ground.PARA.energy.E_snow_lat - E_snow_lat_old;
            %             ground.PARA.energy.dE_snow = ground.PARA.energy.dE_snow + ground.PARA.energy.E_snow - E_snow_old;
            
            % WATER ground
            
            % water content soil domain in [m]
            W_soil_old = ground.PARA.water.W_soil;
            
            ground.PARA.water.W_soil = nansum( ground.STATVAR.wc .* ground.STATVAR.layerThick );
            ground.PARA.water.dW_soil = ground.PARA.water.dW_soil + (ground.PARA.water.W_soil - W_soil_old)*1000; % in [mm]
            
            % water content snow domain in [m]
            W_snow_old = ground.PARA.water.W_snow;
            ground.PARA.water.W_snow = 0;%nansum( GRID.snow.Snow_i + GRID.snow.Snow_w ) + GRID.snow.SWEinitial;
            ground.PARA.water.dW_snow =  ground.PARA.water.dW_snow + ( ground.PARA.water.W_snow - W_snow_old)*1000; % in [mm]
            
        end
        
        function ground = surfaceCondition(ground, forcing)
            
            % set surface parameters (albedo, emissivity, roughnesslength, resistance
            % to evaporation) according to the actual surface conditions
            
            ground.PARA.lake.unfrozenWaterSurface=false;
            
            %default soil surface
            ground.PARA.surf.albedo  = ground.PARA.soil.albedo;
            ground.PARA.surf.epsilon = ground.PARA.soil.epsilon;
            ground.PARA.surf.z0      = ground.PARA.soil.z0;
            ground.PARA.surf.rs      = ground.PARA.soil.rs;
            
            % check if snow cover exists  // for later when sbnow will be included - NC
            % if GRID.snow.cT_domain(GRID.air.cT_domain_lb+1)==1
            % ground.PARA.albedo  = ground.PARA.snow.albedo;
            % ground.PARA.surf.epsilon = ground.PARA.snow.epsilon;
            % ground.PARA.surf.z0      = ground.PARA.snow.z0;
            % ground.PARA.surf.rs      = ground.PARA.snow.rs;
            
            % check if water surface exists and whether it is frozen
            
            % elseif GRID.soil.cT_domain(GRID.air.cT_domain_lb+1)==1 ...
            %    && GRID.soil.cT_organic(1)+GRID.soil.cT_mineral(1)<=1e-6
            %later to modify when snow will be included
            % upper soil cell is pure water
            if ground.STATVAR.T > 0 %T(GRID.soil.cT_domain_ub)>0 % unfrozen - check later
                GRID.lake.unfrozenWaterSurface = true;
                ground.PARA.surf.albedo  = ground.PARA.water.albedo;
                ground.PARA.surf.epsilon = ground.PARA.water.epsilon;
                ground.PARA.surf.z0      = ground.PARA.water.z0;
                ground.PARA.surf.rs      = ground.PARA.water.rs;
            else %frozen
                ground.PARA.surf.albedo  = ground.PARA.ice.albedo;
                ground.PARA.surf.epsilon = ground.PARA.ice.epsilon;
                ground.PARA.surf.z0      = ground.PARA.ice.z0;
                ground.PARA.surf.rs      = ground.PARA.ice.rs;
            end
            % end
        end
        function p = satPresIce_N(T)
            
            p= 0.622.*6.112.* 100.* exp(22.46.*(T-273.15)./(272.61-273.15+T));
        end
        function Q_e = Q_eq_E(ground, forcing)
            uz = forcing.TEMP.wind;
            z =  ground.PARA.airT_height;
            z0 = ground.PARA.z0;
            Tz = forcing.TEMP.Tair;
            T_surf= ground.STATVAR.T(1);
            Lstar = ground.STATVAR.Lstar;
            p = forcing.TEMP.p;
            
            %                 Qe=real(Q_eq(forcing.DATA.wind(i), z, ground.PARA.surf.z0, forcing.DATA.q(i), forcing.DATA.Tair(i), T(1,1), Lstar, ground.PARA.surf.rs, forcing.DATA.p(i), PARA));
            Tz=Tz+273.15;
            T_surf=T_surf+273.15;
            
            rho = p./(ground.CONST.R_a*Tz); %air density [kg m^(-3)]
            kappa = ground.CONST.kappa; %0.4;
            L_w=1000.*(2500.8 - 2.36.*(T_surf-273.15));   % [J/kg] latent heat of evaporation of water
            L_i=ground.CONST.L_sg;  %[J/kg] %1e3.*2834.1; %latent heat of sublimation
            
            ground.PARA.surf.rs = 50; %change it later
            
            if T_surf<=273.15
                Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(ground, z./Lstar, z0./Lstar)).*(forcing.TEMP.q-satPresIce(ground,T_surf)./p)./(log(z./z0)- psi_H(ground, z./Lstar, z0./Lstar));
            else
                Q_e = -rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)- psi_M(ground, z./Lstar, z0./Lstar)).*(forcing.TEMP.q-satPresWater(ground,T_surf)./p)./(log(z./z0)- psi_H(ground, z./Lstar, z0./Lstar)  ...
                    + ground.PARA.surf.rs.*uz.*kappa.^2./(log(z./z0)- psi_M(ground, z./Lstar, z0./Lstar)));
            end
        end
        
        function ground =surfaceEnergyBalanceInfiltration(ground, forcing)
            
            
            Lstar=mean(ground.PARA.SEB.L_star);
            
            sigma=ground.CONST.sigma; %5.67e-8; %Stefan-Boltzmann const.
            L=ground.CONST.L_lg.*ground.CONST.rho_w;
            z=ground.PARA.technical.z;
            
            Qh=real(Q_h(ground,forcing));
            
            dwc_dt=ground.STATVAR.wc.*0;
            ground.PARA.surf.albedo = 0.2;
            ground.PARA.surf.epsilon = 0.97;
            
            %______here SW radiation is calculated_____________________________________
            dE_dt=ground.STATVAR.midpoints.*0;%GRID.general.cT_grid.*0;
            Qsolar=ground.STATVAR.midpoints.*0;%GRID.general.cT_grid.*0;
            
            %             dE_dt(GRID.air.cT_domain_lb+1)=(1-PARA.surf.albedo).*FORCING.i.Sin;
            dE_dt(1,1)=(1-ground.PARA.surf.albedo).*forcing.DATA.Sin(1,1);
            %------ snow surface (solid state green house effect)
            %--------------------- fix later - NC
            %             if ~isempty(GRID.snow.cT_domain_ub)
            %                 beta=PARA.snow.extinction;
            %                 Qsolar(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb+1) = dE_dt(GRID.snow.cT_domain_ub) .* exp(-beta.*(GRID.general.K_grid(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb+1)-GRID.general.K_grid(GRID.snow.cT_domain_ub)));
            %                 dE_dt(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb) = -Qsolar(GRID.snow.cT_domain_ub+1:GRID.snow.cT_domain_lb+1) + Qsolar(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb);
            %                 %put the rest to cell below snow
            %                 dE_dt(GRID.snow.cT_domain_lb+1) = Qsolar(GRID.snow.cT_domain_lb+1);
            %             end
            
            %__________________________________________________________________________
            %             Sout = PARA.surf.albedo*FORCING.i.Sin;
            %             Lout = PARA.surf.epsilon.*sigma.*(T(GRID.air.cT_domain_lb+1)+273.15).^4 + (1-PARA.surf.epsilon).*FORCING.i.Lin;
            %             Qnet = FORCING.i.Sin-Sout + FORCING.i.Lin - Lout ;
            
            Sout = ground.PARA.surf.albedo*forcing.DATA.Sin(1,1);
            Lout = ground.PARA.surf.epsilon.*sigma.*(ground.STATVAR.T(1,1)+273.15).^4 + (1-ground.PARA.surf.epsilon).*forcing.DATA.Lin(1,1);
            Qnet = forcing.DATA.Sin(1,1)-Sout + forcing.DATA.Lin(1,1) - Lout ;
            
            
            %calculate ET
            if ground.PARA.modules.infiltration
                
                %                 % snow cover or uppermost grid cell frozen --> no ET ; this includes the case of a frozen water body
                %                 if ~isempty(GRID.snow.cT_domain_ub) || T(GRID.soil.cT_domain_ub)<=0
                %                     Qe=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, PARA.surf.rs, FORCING.i.p, PARA));
                %                     % unfrozen water body at surface
                %                 else
                if ~ground.PARA.lake.unfrozenWaterSurface
                    Qe=real(Q_eq_E(ground, forcing));
                    dwc_dt(1)=-Qe./L; %in m water per sec, this can be evaporation or condensation
                    
                    % unfrozen soil surface
                else
                    Qe_pot=real(Q_eq(FORCING.i.wind, z, ground.PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, 0, FORCING.i.p, PARA));  %potential ET
                    if Qe_pot>0
                        % determine index of soil cell to which E and T occur
                        i_ALD = ground.PARA.location.bottomBucketSoilcTIndex; % this corresponds to the ALD or gives a reasonable maximum
                        i_Emax = GRID.soil.E_lb;
                        i_Tmax = GRID.soil.T_lb;
                        i_E = min( i_Emax, i_ALD );
                        i_T = min( i_Tmax, i_ALD );
                        r = ground.PARA.soil.ratioET;
                        % cell-wise "efficiencies"
                        fraction_E = getE_fraction( T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+i_E-1), wc(1:i_E), PARA.soil.fieldCapacity );
                        fraction_T = getT_fraction( T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+i_T-1), wc(1:i_T), PARA.soil.fieldCapacity );
                        % grid cell heights for weighting
                        K_delta_E = GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+i_E-1);
                        K_delta_T = GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+i_T-1);
                        % total efficiencies
                        efficiency_E = sum( fraction_E .* K_delta_E ) ./ sum( K_delta_E );
                        efficiency_T = sum( fraction_T .* K_delta_T ) ./ sum( K_delta_T );
                        % actual Qe
                        Qe = min( [ efficiency_E + efficiency_T * r / (1-r), 1 ] ) * Qe_pot;
                        % actual Qe_E and Qe_T partition
                        Qe_E = Qe * efficiency_E / ( efficiency_E + efficiency_T * r / (1-r) );
                        Qe_T = Qe * efficiency_T * r / (1-r) / ( efficiency_E + efficiency_T * r / (1-r) );
                        % associated changes of water amounts in [m/s]
                        dwc_dt(1:i_E) = dwc_dt(1:i_E) + -Qe_E ./ L .* fraction_E .* K_delta_E ./ sum( fraction_E .* K_delta_E );
                        dwc_dt(1:i_T) = dwc_dt(1:i_T) + -Qe_T ./ L .* fraction_T .* K_delta_T ./ sum( fraction_T .* K_delta_T );
                    else  %condensation
                        Qe=Qe_pot;
                        dwc_dt(1)=-Qe./L; %in m water per sec, put everything in uppermost grid cell
                    end
                end
            else % this is identical to case with snow cover or frozen ground
                Qe=real(Q_eq(FORCING.i.wind, z, ground.PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, PARA.surf.rs, FORCING.i.p, PARA));
            end
            %ground heat flux
            Qg   = Qnet-Qh-Qe;
            
            %surface heat flux (into upper cell, ground heat flux regards also other
            %grid cells, should be identical if no snow cover and no evapotranspiration
            %occur
            dE_dt(1,1) = dE_dt(1,1) ...
                + ground.PARA.surf.epsilon.*forcing.TEMP.Lin ...
                - ground.PARA.surf.epsilon.*sigma.*(ground.STATVAR.T(1,1)+273.15).^4 ...
                - Qh - Qe;  % Qe positive: cooling of soil => evaporation/subl. => loss of SWE
            
            % fluxes are in [ W / m^2 ]
            ground.STATVAR.Qsurf = dE_dt(1,1);
            ground.STATVAR.dE_dt_SEB = dE_dt;
            ground.STATVAR.Qnet = Qnet;
            ground.STATVAR.Qh = Qh;
            ground.STATVAR.Qe = Qe;
            ground.STATVAR.Qg = Qg;
            ground.STATVAR.Sout = Sout;
            ground.STATVAR.Lout = Lout;
        end
        
        function ground = heatConduction(ground, forcing, grid)
            
            ground.PARA.soil.Qgeo = 0.05;
            Q=ground.PARA.soil.Qgeo;
            cT_delta = ground.STATVAR.layerThick;
            cT_cellAboveSurface = 1;%GRID.air.cT_domain_lb;
            
            dE_dt=ground.STATVAR.T.*0;
            
            dE_dt(cT_cellAboveSurface)= ground.STATVAR.k_eff(cT_cellAboveSurface+1).*(ground.STATVAR.T(cT_cellAboveSurface+1)-ground.STATVAR.T(cT_cellAboveSurface))./cT_delta(cT_cellAboveSurface) ;
            
%            dE_dt(cT_cellAboveSurface+2:end-1)=  (k_eff(cT_cellAboveSurface+3:end-1).*(T(cT_cellAboveSurface+3:end)-T(cT_cellAboveSurface+2:end-1))./cT_delta(cT_cellAboveSurface+2:end) -...
     %           k_eff(cT_cellAboveSurface+2:end-2).*(T(cT_cellAboveSurface+2:end-1)-T(cT_cellAboveSurface+1:end-2))./cT_delta(cT_cellAboveSurface+1:end-1));
            
            
            % lower BC (dT_dt=geothermal heat flux)
%            dE_dt(end)= Q - k_eff(end-1).*(T(end)-T(end-1))./cT_delta(end);
            
 %           SEB.dE_dt_cond=dE_dt;
        end
        
        function ground = bucketSch(ground, forcing)
            
            ice_c = (ground.STATVAR.ice./ground.STATVAR.layerThick)*0; % ice content
            ice_c(isnan(ice_c)) = 0; % removing nan
            
            mineral_c = ground.STATVAR.mineral ./ground.STATVAR.layerThick; % mineral content
            organic_c =  ground.STATVAR.organic ./ground.STATVAR.layerThick; % organic content
            
            water_c=ground.STATVAR.water./ground.STATVAR.layerThick;% .* ground.STATVAR.layerThick; % convereted water content to total water
            waterIce_c=ground.STATVAR.waterIce ./ ground.STATVAR.layerThick;
            
            porosity = max(1 - ice_c - mineral_c - organic_c,0).* ground.STATVAR.layerThick;% porosity,
            fieldC = min(ground.STATVAR.field_capacity, porosity); % minimum feild capacity
            
            water=ground.STATVAR.water;
            water(isnan(water)) = 0;
            
            waterIce=ground.STATVAR.waterIce;
            waterIce(isnan(waterIce)) = 0;
            
            energy= ground.STATVAR.energy; %in J/m^2
            
            Troom = 25; % room temperature
            salt_solubility = (0.0692*Troom + 33.115)/ 100; %salt (NaCl)solubility at room temperature
            
            %% change in water content on the surface
            dwc_dt = 0; % include it later - NC
            
            forcing.TEMP.rainfall = 9e-1; %% my check NC
            %             forcing.TEMP.Tair = 7; %% my check NC
            
            surface_flux= forcing.TEMP.rainfall+double(dwc_dt(1,1)>0).*dwc_dt(1,1);   %add condensation to surface_flux_pool
            surface_energy_flux= forcing.TEMP.rainfall.*max(0, forcing.TEMP.Tair).*ground.CONST.c_w + double(dwc_dt(1,1)>0).*dwc_dt(1,1).*max(0,ground.STATVAR.T(1,1)).*ground.CONST.c_w;
            dwc_dt(1,1)=double(dwc_dt(1,1)<=0).*dwc_dt(1,1);
            % %%
            %% --1. take out ET------------
            dwc_dt=max(-water, dwc_dt);  %prevent from going negative
            water = water+dwc_dt;
            waterIce=waterIce+dwc_dt;
            energy= energy + dwc_dt.*ground.STATVAR.T.*ground.CONST.c_w;     %this should leave temperatures constant
            
            %% --2. find depth of water table or "hard" layer" preventing infiltration
            % i=1;
            % while water(i)<porosity(i) && porosity(i)>ground.CONST.infiltration_cutoff .* ground.STATVAR.layerThick(i) && ground.STATVAR.T(i)>0
            %     i=i+1;
            % end
            % i=i-1; %layer above saturated or hard one
            %
            % %%
            % %---3. add/subtract external fux--------------
            % while i>0 && ground.CONST.external_flux~=0 && porosity(i)>ground.CONST.infiltration_cutoff .* ground.STATVAR.layerThick(i) && ground.STATVAR.T(i)>0
            %     change=double(ground.CONST.external_flux>0).*min(max(0, porosity(i)-water(i)), ground.CONST.external_flux) + double(ground.CONST.external_flux<0).* max(min(0, fieldC(i)-water(i)), ground.CONST.external_flux );
            %     water(i) = water(i)+change;
            %     waterIce(i)=waterIce(i)+change;
            %     energy(i)=energy(i) + change.*ground.STATVAR.T(i).*ground.CONST.c_w;   %this should leave temperatures constant
            %     ground.CONST.external_flux=ground.CONST.external_flux-change;
            %     i=i-double(ground.CONST.external_flux>0)+double(ground.CONST.external_flux<0);  %move up or down
            % end
            %
            % if i>0
            % i=i-double(porosity(i)>ground.CONST.infiltration_cutoff.*ground.STATVAR.layerThick(i));
            % end
            % %i is the last unsaturated grid cell above bucket
            
            %% --2. find depth of water table or "hard" layer" preventing infiltration
            
            i=1;
            while water(i)<porosity(i) && porosity(i)> forcing.PARA.infiltration_cutoff.*ground.STATVAR.layerThick(i)%&& ground.STATVAR.T(i)>0
                i=i+1;
                if (i > length(ground.STATVAR.gridlength))
                    break;
                end
            end
            i=i-1; %layer above saturated or hard one
            
            
            %% ---4. infiltrate rain + condensation top->down
            energy_water=water.*ground.CONST.c_w.*ground.STATVAR.T;  %in J/m^2
            
            go_on=1;
            while ((i>0 && surface_flux>0) ||  (i>1 && sum(water(1:i-1)>fieldC(1:i-1))>0) )&& go_on>0
                go_on=0;
                
                j=1;
                
                change=min(porosity(j)-water(j), surface_flux);
                go_on=go_on+double(change>0);
                water(j) = water(j)+change;
                waterIce(j)=waterIce(j)+change;
                %%
                if surface_flux>0
                    energy(j)=energy(j) + change./surface_flux.*surface_energy_flux;
                    energy_water(j)=energy_water(j)+change./surface_flux.*surface_energy_flux;
                    surface_energy_flux=surface_energy_flux.*(1-change./surface_flux);
                    surface_flux=surface_flux-change;   %grid cell full
                end
                %%
                change=max(0, water(j)-fieldC(j));  %could go to next cell
                
                %                 ground.STATVAR.leaching (j) = min(ground.STATVAR.leaching_rate(j), ground.STATVAR.total_tracer(j));
                
                j=j+1;
                while j<=i
                    change=min(porosity(j)-water(j), change);  %will go to next cell
                    go_on=go_on+double(change>0);
                    water(j) = water(j)+change;
                    waterIce(j)=waterIce(j)+change;
                    energy(j)=energy(j) + min(0,change./water(j-1).*energy_water(j-1)); % NC - included min flag to remove nan problem
                    energy_water(j)=energy_water(j) + min(0,change./water(j-1).*energy_water(j-1)); % NC - included min flag to remove nan problem
                    
                    %                     ground.STATVAR.leaching (j) = max(ground.STATVAR.leaching_rate(j-1).*change.*salt_solubility, 0);
                    
                    %                     ground.STATVAR.total_tracer(j) = ground.STATVAR.total_tracer(j) + ground.STATVAR.leaching(j) .*(ground.STATVAR.total_tracer(j-1)>0);
                    
                    water(j-1)=water(j-1)-change;
                    waterIce(j-1)=waterIce(j-1)-change;
                    energy(j-1)=energy(j-1).*(1-min(0,change./water(j-1)));% NC - included min flag to remove nan problem
                    energy_water(j-1)=energy_water(j-1).*(1-min(0,change./water(j-1)));% NC - included min flag to remove nan problem
                    
                    %                     if (ground.STATVAR.leaching (j) > ground.STATVAR.total_tracer(j-1))
                    %                         ground.STATVAR.leaching(j) = ground.STATVAR.total_tracer(j-1);
                    %                     end
                    
                    %                     ground.STATVAR.total_tracer(j-1) = max(ground.STATVAR.total_tracer(j-1) - ground.STATVAR.leaching(j),0);
                    
                    change=max(0, water(j)-fieldC(j));  %could go to next cell
                    j=j+1;
                    i=i-double(water(i,1)>=porosity(i,1));  %move water table up
                end
                
                if(j>i)
                    go_on = 0; i=0;
                end
                
            end
            
            %% Consistency check
            %             ii =1;
            %             if (sum (ground.STATVAR.total_tracer )> sum (ground.STATVAR.total_tracer_fixed)+1e-3)
            %                 while(ground.STATVAR.total_tracer(ii) <= 0)
            %                     ii = ii+1;
            %                     if (ii > length(Grid_length))
            %                         break;
            %                     end
            %                 end
            %                 ground.STATVAR.total_tracer(ii) = ground.STATVAR.total_tracer(ii)-(sum (ground.STATVAR.total_tracer) - sum (ground.STATVAR.total_tracer_fixed));
            %             end
            %%
            ground.STATVAR.water=water;
            ground.STATVAR.waterIce=waterIce;
            ground.STATVAR.energy=energy;
            ground.STATVAR.porosity= porosity;%in m,
            ground.STATVAR.fieldC = min(fieldC, ground.STATVAR.porosity);
            %%
            
            set(gcf,'un','n','pos',[0.1,0.1,0.5,0.5]);
            a1= subplot(1,1, 1);
            plot(ground.STATVAR.water./ground.STATVAR.layerThick,ground.STATVAR.gridlength);
            set (a1,'Ydir','reverse')
            subplot(1,1,1);%axis([0,1,-(array_num+1),0])
            hold (a1,'on');
            title(a1,'Water content')
            %             plot(ground.STATVAR.total_tracer./sum(ground.STATVAR.total_tracer_fixed),Grid_length-length(Grid_length),'color', [0 0 0.8])%,'*')
            hold on
            pause(0.001);
            
            %%
            %             set(gcf,'un','n','pos',[0.1,0.1,0.8,0.8]);
            %             a1= subplot(1,9, 1);
            %             plot(ground.STATVAR.water./ground.STATVAR.layerThick,Grid_length-length(Grid_length))
            %             set (a1,'Ydir','reverse')
            %             subplot(1,9,1);axis([0,1,-200,0])
            %             hold (a1,'on');
            %             title(a1,'Water')
            %
            %             a2= subplot(1,9, 2);
            %             plot(ground.STATVAR.total_tracer./sum(ground.STATVAR.total_tracer_fixed),Grid_length-length(Grid_length),'color', [0 0 0.8])%,'*')
            %             %plot(ground.STATVAR.ice ./ground.STATVAR.layerThick,Grid_length-length(Grid_length))
            %             set (a2,'Ydir','reverse')
            %             hold (a2,'on');
            %             subplot(1,9,2);axis([0,1,-200,0])
            %             title(a2,'Salt')
            %
            %             a3= subplot(1,9, 3);
            %             plot(ground.STATVAR.waterIce./ground.STATVAR.layerThick,Grid_length-length(Grid_length))
            %             set (a3,'Ydir','reverse')
            %             hold (a3,'on');
            %             subplot(1,9,3);axis([0,1,-200,0])
            %             title(a3,'WaterIce')
            %
            %             a4= subplot(1,9, 4);
            %             plot(ground.STATVAR.fieldC./ground.STATVAR.layerThick,Grid_length-length(Grid_length))
            %             set (a4,'Ydir','reverse')
            %             hold (a4,'on');
            %             subplot(1,9,4);axis([0,1,-200,0])
            %             title(a4,'FC')
            %
            %             a5= subplot(1,9, 5);
            %             plot(ground.STATVAR.mineral./ground.STATVAR.layerThick,Grid_length-length(Grid_length))
            %             set (a5,'Ydir','reverse')
            %             hold (a5,'on');
            %             subplot(1,9,5);axis([0,1,-200,0])
            %             title(a5,'Mineral')
            %
            %             a6= subplot(1,9, 6);
            %             plot(ground.STATVAR.organic./ground.STATVAR.layerThick,Grid_length-length(Grid_length))
            %             set (a6,'Ydir','reverse')
            %             hold (a6,'on');
            %             subplot(1,9,6);axis([0,1,-200,0])
            %             title(a6,'Organic')
            %
            %             a7= subplot(1,9, 7);
            %             plot(ground.STATVAR.porosity./ground.STATVAR.layerThick,Grid_length-length(Grid_length))
            %             set (a7,'Ydir','reverse')
            %             hold (a7,'on');
            %             subplot(1,9,7);axis([0,1,-200,0])
            %             title(a7,'Porosity')
            %
            %             a8= subplot(1,9, 8);
            %             plot((ground.STATVAR.organic+ground.STATVAR.mineral+ground.STATVAR.ice+ground.STATVAR.water)./ground.STATVAR.layerThick,Grid_length-length(Grid_length))
            %             set (a8,'Ydir','reverse')
            %             hold (a8,'on');
            %             subplot(1,9,8);axis([0,1.0,-200,0])
            %             title(a8,'Sum')
            %
            %             a9= subplot(1,9, 9);
            %             plot(ground.STATVAR.layerThick,Grid_length-length(Grid_length))
            %             set (a9,'Ydir','reverse')
            %             % hold (a6,'on');
            %             title(a9,'Depth')
            
            pause(0.0001);
            %surface_runoff=excess_water;
            
        end
        
    end
    
    
end
