%========================================================================
% CryoGrid GROUND class GROUND_freeW_bucketW_seb
% heat conduction, free water freeze curve, constant water + ice water balance, 
% surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_LITE1 < CG_LITE
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function ground = provide_PARA(ground)
            
            ground.PARA.albedo = []; %surface albedo [-]
            ground.PARA.epsilon = []; % surface emissivity [-]
            ground.PARA.z0 = []; %roughness length [m]
            ground.PARA.rs = [];
            ground.PARA.rs_frozen = [];
            
            ground.PARA.evaporationDepth = []; %e-folding constant of evaporation reduction reduction with depth [m]
            
            ground.PARA.fieldCapacityFactor = [];
            ground.PARA.residualWaterFactor = [];

            ground.PARA.timestep = []; %timestep [day]
            ground.PARA.waterDepth = 0;
            ground.PARA.SnowDensity = [];
            ground.PARA.SnowCoverMax = [];
            ground.PARA.snowGridCellSize = [];
            
            ground.PARA.max_albedo = [];
            ground.PARA.min_albedo = [];
            ground.PARA.epsilon_snow = [];
            ground.PARA.z0_snow = [];
            ground.PARA.tau_1 = 86400;
            ground.PARA.tau_a = 0.008;
            ground.PARA.tau_f = 0.24;
   
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
            
            ground.CONST.gamma = 0.622;
            ground.CONST.R_a = 287.0580;
            ground.CONST.rho_a = 1.293;
            ground.CONST.c_a = 1057;
            ground.CONST.L_sl = 3.34e5;
            ground.CONST.L_lg = 2501000;
            ground.CONST.Pr_0 = 0.74;
            ground.CONST.beta_m = 4.7; 
            ground.CONST.beta_h = 6.3514;
            ground.CONST.gamma_m = 15;
            ground.CONST.gamma_h = 9;
            
        end

        function ground = convert_units(ground, tile)
            unit_converter = str2func(tile.PARA.unit_conversion_class);
            unit_converter = unit_converter();
            ground = convert_normal(unit_converter, ground, tile);
        end
        
        function ground = finalize_init(ground, tile)

            
            %ground = get_E_freeW(ground);
            
            %get_E_freeW
            T = ground.STATVAR.T;
            mineral= ground.STATVAR.mineral ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            organic = ground.STATVAR.organic ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            waterIce = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            layerThick = ground.STATVAR.layerThick;
            area = ground.STATVAR.area;
            
            energy = T.*(mineral .* ground.CONST.c_m + organic .* ground.CONST.c_o + double(T>=0).*(waterIce .* ground.CONST.c_w) + ...
                double(T<0).*(waterIce .* ground.CONST.c_i )) - double(T<0) .* (waterIce) .* ground.CONST.L_f;
            
            ground.STATVAR.waterIce = waterIce .* layerThick .* area; % [m3]
            ground.STATVAR.mineral = mineral .* layerThick .* area; % [m3]
            ground.STATVAR.organic = organic .* layerThick .* area; % [m3]
            ground.STATVAR.energy = energy .* layerThick .* area;  % [J]
            
            ground.STATVAR.water = double(T>=0) .* waterIce .* layerThick .* area;  % [m3]
            ground.STATVAR.ice = double(T<0) .* waterIce .* layerThick .* area; %[m3]
            ground.STATVAR.air = (1-mineral-organic-waterIce) .* layerThick .* area;  % [m3]
            %end get_E_freeW
            
            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
            ground.STATVAR.Tsurf = 0;
            ground.STATVAR.Qg = 0;
            ground.STATVAR.Qnet = 0;
            ground.TEMP.d_energy = ground.STATVAR.energy.*0;
            
            %make CGLite grids
            ground.STATVAR.layerThick = [repmat(ground.PARA.snowGridCellSize, ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.layerThick];
            ground.TEMP.Zn = [-ground.PARA.SnowCoverMax; -ground.PARA.SnowCoverMax + cumsum(ground.STATVAR.layerThick)];
            ground.TEMP.Zn = ground.TEMP.Zn(1:end-1,1);
            ground.TEMP.Zs = -ground.PARA.SnowCoverMax + cumsum(ground.STATVAR.layerThick);
            ground.TEMP.Zp = (ground.TEMP.Zs + ground.TEMP.Zn)./2;
            ground.TEMP.dxp = ground.STATVAR.layerThick;
            
            ground.STATVAR.water = [repmat(0, ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.water];
            ground.STATVAR.waterIce = [repmat(0, ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.waterIce];
            ground.STATVAR.mineral = [repmat(0, ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.mineral];
            ground.STATVAR.organic = [repmat(0, ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.organic];
            ground.STATVAR.energy = [repmat(0, ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.energy];
            ground.STATVAR.ice = [repmat(0, ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.ice];
            ground.STATVAR.air = [repmat(0, ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.air];
            ground.STATVAR.T = [repmat(0, ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.T];
            ground.STATVAR.area = [repmat(ground.STATVAR.area(1,1), ground.PARA.SnowCoverMax./ground.PARA.snowGridCellSize, 1); ground.STATVAR.area];
            
            ground.STATVAR.thermalConductivity = ThermalConductivity(ground, ground.STATVAR.waterIce./ground.STATVAR.layerThick./ground.STATVAR.area, ...
                ground.STATVAR.water./ground.STATVAR.layerThick./ground.STATVAR.area, ground.STATVAR.mineral./ground.STATVAR.layerThick./ground.STATVAR.area, ...
                ground.STATVAR.organic./ground.STATVAR.layerThick./ground.STATVAR.area);
            ground.STATVAR.heatCapacity = HeatCapacity(ground, ground.STATVAR.waterIce./ground.STATVAR.layerThick./ground.STATVAR.area, ...
                ground.STATVAR.water./ground.STATVAR.layerThick./ground.STATVAR.area, ground.STATVAR.mineral./ground.STATVAR.layerThick./ground.STATVAR.area, ...
                ground.STATVAR.organic./ground.STATVAR.layerThick./ground.STATVAR.area);
            
            [Zn, Zs, ground.TEMP.dxn, ground.TEMP.dxs, ground.TEMP.kn, ground.TEMP.ks] = MakeGrid(ground, ground.TEMP.Zp, ground.TEMP.dxp, ground.STATVAR.thermalConductivity);
            
            ground.TEMP.albedo = ground.PARA.albedo;
            ground.TEMP.emissivity = ground.PARA.epsilon;
            ground.TEMP.z0 = ground.PARA.z0;
            ground.TEMP.rs = ground.PARA.rs;
            ground.TEMP.SnowDepth = 0;
            ground.TEMP.latFlux = 0.*ground.STATVAR.layerThick;

        end
        
        function ground = finalize_init2(ground, tile)

            %ground = get_E_freeW(ground);
            
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)

        end
        
        
        function ground = get_boundary_condition_l(ground, tile)
%             forcing = tile.FORCING;
%             ground.TEMP.F_lb = forcing.PARA.heatFlux_lb .* ground.STATVAR.area(end);
%             ground.TEMP.d_energy(end) = ground.TEMP.d_energy(end) + ground.TEMP.F_lb;
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
%             ground = get_derivative_energy(ground);

        end
        
        function timestep = get_timestep(ground, tile) 
%            timestep = get_timestep_heat_coduction(ground);
           timestep = 24.*3600;
        end
        
        function ground = advance_prognostic(ground, tile)           
%             timestep = tile.timestep;
%             %energy
%             ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
              ground = CryoGridImplicit(ground, tile);
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
%             forcing = tile.FORCING;
%             ground = L_star(ground, forcing);
        end
       
        function ground = compute_diagnostic(ground, tile)
%             ground = get_T_water_freeW(ground);
%             ground = conductivity(ground);
%             ground.TEMP.d_energy = ground.STATVAR.energy.*0;
        end
        
        function ground = check_trigger(ground, tile)
           %do nothing 
        end
        

        
        
%         function ground = conductivity(ground)
%             conductivity_function = str2func(ground.PARA.conductivity_function);
%             ground = conductivity_function(ground);
%         end
        
        function ground = CryoGridImplicit(ground, tile)
            
            dt = ground.PARA.timestep .* 60*60*24;
            
            dxp = ground.TEMP.dxp;
            Zp = ground.TEMP.Zp;
            Zn = ground.TEMP.Zn;
            Zs = ground.TEMP.Zs;
            dxn = ground.TEMP.dxn;
            dxs = ground.TEMP.dxs;
            kn = ground.TEMP.kn;
            ks = ground.TEMP.ks;
            dxo = 1;
            An = ground.STATVAR.area;
            As = ground.STATVAR.area;
            Ao = 1;
            Vp = ground.STATVAR.area .* ground.STATVAR.layerThick;
            
            cp = ground.STATVAR.heatCapacity;
            kp = ground.STATVAR.thermCond;
%            lat_flux = ground.TEMP.lat_flux;
            SnowDepth = ground.TEMP.SnowDepth;
            
            
            WaterIce = ground.STATVAR.waterIce./ground.STATVAR.layerThick./ground.STATVAR.area;
            Mineral = ground.STATVAR.mineral./ground.STATVAR.layerThick./ground.STATVAR.area;
            Organic = ground.STATVAR.organic./ground.STATVAR.layerThick./ground.STATVAR.area;
            Water = ground.STATVAR.water./ground.STATVAR.layerThick./ground.STATVAR.area;
            
            WaterDepth = 0;
            WaterDensity = ground.CONST.rho_w;
            SnowDensity = ground.PARA.SnowDensity;
            Qgeo = tile.FORCING.PARA.heatFlux_lb;
            SnowDepthMax = ground.PARA.SnowCoverMax;
            
            T = ground.STATVAR.T;
            H = ground.STATVAR.energy./ground.STATVAR.layerThick./ground.STATVAR.area;
            
            N=size(T,1);
            J=size(T,2);
            
            %upper boundary temperature [°C] and snow fall rate [m/day]
            T0 = tile.FORCING.TEMP.Tair; %function SEB
            
            % snow and rain precipitation [m/day]
            snow_fall = tile.FORCING.TEMP.snowfall./1000 .* ground.PARA.timestep;
            rain_fall = tile.FORCING.TEMP.rainfall./1000 .* ground.PARA.timestep;
            
            dp = 1.0/dt;
            T_old = T;
            H_old = H;
            
            j=1;
            %Prognostic Snow Scheme
            %tsvd        [Water(:,j), WaterIce(:,j), SnowDepth(1,j), idx] = SnowCover2(Water(:,j),WaterIce(:,j),Zs,SnowDepth(1,j),SnowDepthMax(1,j),snow_fall,SnowDensity,WaterDensity);
            %[Water(:,j), WaterIce(:,j), SnowDepth(1,j), idx, ground.TEMP.albedo] = SnowCover2(Water(:,j),WaterIce(:,j),Zs,SnowDepth(1,j),SnowDepthMax(1,j),snow_fall,SnowDensity,WaterDensity,ground.TEMP.albedo,PARA);
            [Water(:,j), WaterIce(:,j), SnowDepth(1,j), idx] = SnowCover2(ground, Water(:,j),WaterIce(:,j),Zs,SnowDepth(1,j),SnowDepthMax(1,j),snow_fall,SnowDensity,WaterDensity);
            
            %update properties and state variables
            cp(:,j) = HeatCapacity(ground, WaterIce(:,j), Water(:,j), Mineral(:,j), Organic(:,j));
            kp(:,j) = ThermalConductivity(ground, WaterIce(:,j), Water(:,j), Mineral(:,j), Organic(:,j));
            [Zn, Zs, dxn, dxs, kn, ks] = MakeGrid(ground, Zp,dxp,kp);
            H(:,j) = Enthalpie(ground, T(:,j), Water(:,j), cp(:,j));
            [T(:,j), ~] = EnthalpieInv(ground, H(:,j), WaterIce(:,j), cp(:,j));
            H_old(:,j) = H(:,j);
            T_old(:,j) = T(:,j);
            
%             % Surface energy balance to get upper boundary flux
%             %temp_offset_max = 1.;
%             %iter_count = 0;
%             SEB.Tsurf(j) = T(idx,j); % update surface temperature with uppermost grid cell
%             %while abs(temp_offset_max)>1e-3 && iter_count<1000
%             %iter_count=iter_count+1;
%             %Tsurf_old = SEB.Tsurf(j);
            Porosity = 1-Mineral-Organic;
            
%             [TEMP] = SurfaceCondition(PARA, TEMP, T(:,j),  WaterIce(:,j), Mineral(:,j), Organic(:,j), Zp, j, idx, dt); %tsvd make SEB parameters surface dependent
%             SEB = SurfaceEnergyBalance( T(:,j), Water(:,j), Porosity(:,j), Zp, dxp, kp(:,j),  SEB, PARA, FORCING, TEMP, j, idx, t);
%             %temp_offset_max = SEB.Tsurf(j) - Tsurf_old;
%             %end
%             T0=SEB.Tsurf(j);
            
            %         if iter_count<1000
            %             fprintf("\t SEB for tile #%d converged after %d iterations: Tsurf=%0.2f, Qg=%0.2f\n", j, iter_count, SEB.Tsurf,SEB.Qg );
            %         else
            %             fprintf("Did not converge.\n");
            %         end
            
            %Prognostic Lake Scheme
            %get the last cell that consists of 100% liquid water
            SolidMatrix = Mineral(:,j) + Organic(:,j);
%             lic_lake = Water(idx,j)/(WaterIce(idx,j)+SolidMatrix(idx));
%             while lic_lake == 1.0
%                 %air = 1.0 - WaterIce(idx,j) - SolidMatrix(idx,j);
%                 lic_lake = Water(idx,j)/(WaterIce(idx,j)+SolidMatrix(idx));
%                 if lic_lake == 1.0
%                     idx = idx+1;
%                 end
%             end
            
            %set index of upper boundary
            ubc_idx = idx;
            
            %Main Iteration Loop
            temp_offset_max = 1;
            iter_count = 1;
            while (temp_offset_max>=1e-4 || iter_count<2) && iter_count<=1000
                if iter_count>999
                    disp('warning: did not reach convergence! Max Temperature offset:')
                    disp(temp_offset_max)
                end
                iter_count=iter_count+1;
                
                %implicit update heat capacity and thermal conductivity
                cp(:,j) = HeatCapacity(ground, WaterIce(:,j), Water(:,j), Mineral(:,j), Organic(:,j));
                kp(:,j) = ThermalConductivity(ground, WaterIce(:,j), Water(:,j), Mineral(:,j), Organic(:,j));
                
                %pre-factores acording to grid cell sizes and thermal conductivities
                kn(ubc_idx) = kp(ubc_idx);%ensures full thermal conductivity at upper boundary
                dxn(ubc_idx) = dxp(ubc_idx)/2.0; %ensures half grid cell at upper boundary
                
                anpn = An(:,j)./Vp(:,j).*kn./dxn;
                anps = As(:,j)./Vp(:,j).*ks./dxs;
                
                %Additional heat fluxes----------------------------------------
                bp = zeros(N,1);
                %flux from upper boundary
                %bp(1) = (An(:,j)./Vp(1,j).*kn(1)./dxn(1)) * T0; %[W/m³]
                %bp(1) = (An(:,j)./Vp(1,j)) * SEB.Qg(j); %[W/m³]
                
                % using surface temperature from SEB:
                bp(1:ubc_idx) = (An(1:ubc_idx,j)./Vp(1:ubc_idx,j).*kn(1:ubc_idx)./dxn(1:ubc_idx)) * ( T0 );
                
                % usign ground heat flux from SEB:
                %bp(1:ubc_idx) = (An(:,j)./Vp(1:ubc_idx,j)) * SEB.Qg(j);
                
                anps(1:ubc_idx-1) = 0.0;
                
                
                %flux from lower boundary
                bp(end) = (An(end,j)./Vp(end,j)) * Qgeo; %[W/m³]
                
                
%                 if j~=1
%                     ko = (kp(:,1) + kp(:,j))/2.0; %simple mean for testing lat flux
%                     ap = anpn + anps + Ao(:,j)./Vp(:,j).*ko./dxo(:,j);
%                     ap(end) = anpn(end) + Ao(end,j)./Vp(end,j).*ko(end)./dxo(end,j);
%                     anpn(1:ubc_idx) = 0.0;
% %                     %lateral heat flux (relative [W/m³])
% %                     bp_lat(:,j) = (Ao(:,j)./Vp(:,j).*ko./dxo(:,j)) .* T_old(:,1);
% %                     %update total lateral heat fluxes [W]
% %                     lat_flux(:,j) = (Ao(:,j).*ko./dxo(:,j)) .* (T(:,j)-T_old(:,1));
% %                     %lat_flux(:,j) = (Ao(:,j)./Vp(:,j).*ko./dxo(:,j)) .* (T(:,j)-T_old(:,1));
%                 else
                    %for tile 1 the lateral flux is assumed to be a static
                    %external flux reulting from the sum of the lateral heat
                    %fluxes from the other tiles.
                    ap = anpn + anps;
                    ap(end) = anpn(end);
                    
                    % modified for Qg
                    %ap(1:ubc_idx) = anps(ubc_idx);
                    
                    
                    anpn(1:ubc_idx) = 0.0;
%                     bp_lat(:,j) = sum(lat_flux,2)./Vp(:,j); %[W/m³]
%                 end
%                 bp = bp + bp_lat(:,j);
                
                %--------------------------------------------------------------
                [Hinv, ~] = EnthalpieInv(ground, H(:,j), WaterIce(:,j), cp(:,j)); %previous iter
                [~, dHdT] = Enthalpie(ground, Hinv, Water(:,j), cp(:,j)); %previous iter
                
                Sp = -dp*dHdT;
                Sc = dp*(H_old(:,j) - H(:,j)) - Sp.*Hinv;
                
                %TDMA solver --------------------------------------------------
                alpha = -anpn;
                beta = (ap - Sp);
                gamma = -anps;
                delta = Sc + bp;
                T(:,j) = TDMAsolver(ground, alpha,beta,gamma,delta)'; %this iter
                
                %--------------------------------------------------------------
                %update current state of H
                H(:,j) = H(:,j) + dHdT.*(T(:,j) - Hinv);%this iter
                [Hinv_check, fl] = EnthalpieInv(ground, H(:,j), WaterIce(:,j), cp(:,j));%this iter
                Water(:,j) = fl.*WaterIce(:,j); %this iter
                
                R = Hinv_check; %this iter
                B = T(:,j); %this iter
                
                temp_offset = abs(R-B);
                temp_offset_max = max(temp_offset(ubc_idx:end-1,:));
                T(:,j) = Hinv_check;
            end
            
%             OUT.T_out(:,t,j) = T(:,j);
%             OUT.Water_out(:,t,j) = Water(:,j);
%             OUT.WaterIce_out(:,t,j) = WaterIce(:,j);
%             OUT.SnowDepth_out(:,t,j) = SnowDepth(:,j);
            %OUT.Albedo_out(t,j) = PARA.albedo; % for analysis purposes
            
            ground.TEMP.dxp = dxp;
            ground.TEMP.Zp = Zp;
            ground.TEMP.Zn = Zn;
            ground.TEMP.Zs = Zs;
            ground.TEMP.dxn = dxn;
            ground.TEMP.dxs = dxs;
            ground.TEMP.kn = kn;
            ground.TEMP.ks = ks;
            %dxo = 1;
            %ground.STATVAR.area = An;
            
%             ground.STATVAR.heatCapacity = cp;
%             ground.STATVAR.thermCond = kp;
%            lat_flux = ground.TEMP.lat_flux;
            ground.TEMP.SnowDepth = SnowDepth;
            
            
            ground.STATVAR.waterIce = WaterIce .* ground.STATVAR.layerThick.*ground.STATVAR.area;
            ground.STATVAR.mineral = Mineral .* ground.STATVAR.layerThick.*ground.STATVAR.area;
            ground.STATVAR.organic = Organic .* ground.STATVAR.layerThick.*ground.STATVAR.area;
            ground.STATVAR.water = Water .* ground.STATVAR.layerThick.*ground.STATVAR.area;
            
            %WaterDepth = 0;
            %WaterDensity = ground.CONST.rho_w;
            %ground.STATVAR.SnowDensity = SnowDensity;
            %Qgeo = tile.FORCING.PARA.heatFlux_lb;
            %ground.PARA.SnowCoverMax = SnowDepthMax;
            
            ground.STATVAR.T = T;
            ground.STATVAR.energy = H .* ground.STATVAR.layerThick .* ground.STATVAR.area;

        end

        
                 
         %-------------param file generation-----
         function ground = param_file_info(ground)
             ground = param_file_info@BASE(ground);
             %ground = provide_PARA(ground);
             
             ground.PARA.class_category = 'GROUND';
             
             %ground.PARA.options = [];
             ground.PARA.STATVAR = {'waterIce' 'mineral' 'organic' 'T'};
             
             ground.PARA.default_value.albedo = {0.2};
             ground.PARA.comment.albedo = {'surface albedo [-]'};
             
             ground.PARA.default_value.epsilon = {0.99};
             ground.PARA.comment.epsilon = {'surface emissivity [-]'};
             
             ground.PARA.default_value.z0 = {0.01};
             ground.PARA.comment.z0 = {'roughness length [m]'};
             
             ground.PARA.default_value.rs = {0};
             ground.PARA.comment.rs ={'surface resistance against evapotranspiration [sec/m]'};
             
             ground.PARA.default_value.conductivity_function = {'conductivity_mixing_squares'};
             ground.PARA.comment.conductivity_function = {'function employed to calculate thermal conductivity, leave empty for default'};
             
             ground.PARA.default_value.dt_max = {3600};
             ground.PARA.comment.dt_max = {'maximum possible timestep [sec]'};
             
             ground.PARA.default_value.dE_max = {50000};
             ground.PARA.comment.dE_max = {'maximum possible energy change per timestep [J/m3]'};
        end
        

    end
    
end
