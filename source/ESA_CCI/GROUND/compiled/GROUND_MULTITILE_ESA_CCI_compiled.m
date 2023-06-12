
classdef GROUND_MULTITILE_ESA_CCI_compiled < BASE
    

    methods
        
        %-----initialize-----------------
        
        function ground = provide_PARA(ground)
            ground.PARA.virtual_gridCellSize = [];
            ground.PARA.timestep = [];
            ground.PARA.adjust_stratigraphy_date = [];
            
            ground.PARA.wind_speed_class = [];
            ground.PARA.wind_compaction_timescale = [];
            ground.PARA.water_table_depth = [];
            
            ground.PARA.speedup_size = 2100;

        end
        
        function ground = provide_CONST(ground)
            ground.CONST.L_f = []; %3.34e8;
            ground.CONST.c_w = []; %4.2e6; %[J/m³K]
            ground.CONST.c_o = []; %2.5e6; %[J/m³K]
            ground.CONST.c_m = []; %2.0e6; %[J/m³K]
            ground.CONST.c_i = []; %1.9e6;%[J/m³K]
           
            ground.CONST.k_w = []; 
            ground.CONST.k_o = []; 
            ground.CONST.k_m = [];
            ground.CONST.k_i = []; 
            ground.CONST.k_a = [];
            
            ground.CONST.day_sec = []; %24*3600;
            
            ground.CONST.g = [];
        end
        

        function ground = provide_STATVAR(ground)

            ground.STATVAR.layerThick = []; % thickness of grid cells [m]
            ground.STATVAR.layerDistance = []; % distance between midpoints of grid cells [m]
            ground.STATVAR.waterIce = [];  % total volume of water plus ice in a grid cell [m3]
            ground.STATVAR.mineral = [];   % total volume of minerals [m3]
            ground.STATVAR.organic = []; % total volume of organics [m3]
            ground.STATVAR.energy = [];   % total internal energy [J]
            ground.STATVAR.soil_type = [];  % integer code for soil_type; 1: sand; 2: silt: 3: clay: 4: peat; 5: water (i.e. approximation of free water, very large-pore ground material).
                        
            ground.STATVAR.T = [];  % temperature [degree C]
            ground.STATVAR.thermCond = [];   %thermal conductivity [W/mK]      
            
            ground.STATVAR.FT_state = [];
        end
        
        function ground = finalize_init(ground, tile)

            %multiply STATVARs with layerThick?
            ground.STATVAR.waterIce = ground.STATVAR.waterIce .* ground.STATVAR.layerThick;
            ground.STATVAR.mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick;
            ground.STATVAR.mineral(ground.STATVAR.mineral<0) = 0;
            ground.STATVAR.organic = ground.STATVAR.organic .* ground.STATVAR.layerThick;
            ground.STATVAR.organic(ground.STATVAR.organic<0) = 0;
            
            ground.STATVAR.T_onset_freezing = 0;
            ground = get_T_end_freezing(ground);
            
            ground = init_conductivity(ground);
       
            
            if size(ground.STATVAR.T,1)==1
                ground.STATVAR.T = [ground.STATVAR.T; ground.STATVAR.organic(2:end,:).*0];

                for i=1:size(ground.STATVAR.T,1)-1
                    k = double(ground.STATVAR.T(i,:) < ground.STATVAR.T_end_freezing(i,:)) .* ground.STATVAR.k_frozen(i,:) + double(ground.STATVAR.T(i,:) > ground.STATVAR.T_onset_freezing) .* ground.STATVAR.k_thawed(i,:) + ...
                        double(ground.STATVAR.T(i,:) >= ground.STATVAR.T_end_freezing(i,:) & ground.STATVAR.T(i,:) <= ground.STATVAR.T_onset_freezing) .* ground.STATVAR.k_freezing(i,:);
                    ground.STATVAR.T(i+1,:) = ground.STATVAR.T(i,:) + tile.PARA.geothermal .* ground.STATVAR.layerDistance(i,:) ./k;
                end
            end
            
            ground.PARA.snow_gridCellSize = ground.PARA.virtual_gridCellSize; %CHECK and ADAPT
            ground.STATVAR.layerThick_first_ground_cell = ground.STATVAR.layerThick(1,:);
            
            %add four cells on top, one virtual and three snow cells
            ground.STATVAR.layerThick = [repmat(ground.PARA.virtual_gridCellSize,4, size(ground.STATVAR.layerThick,2)); ground.STATVAR.layerThick];
            ground.STATVAR.layerDistance = (ground.STATVAR.layerThick(1:end-1,:) + ground.STATVAR.layerThick(2:end,:)) ./ 2;
            
            ground = calculate_E_frozen(ground);
            ground = T2E(ground);
            
            ground.STATVAR.T = [zeros(4, size(ground.STATVAR.layerThick,2)); ground.STATVAR.T];
            
            %snow
            %no snow in the beginning
            ground.TEMP.index_first_ground_cell = repmat(4, 1, size(ground.STATVAR.layerThick,2));
            ground.TEMP.snow_mat1 = [zeros(3, size(ground.STATVAR.layerThick,2)); ones(1, size(ground.STATVAR.layerThick,2))];
            ground.TEMP.snow_mat2 = ground.TEMP.snow_mat1;
            ground.TEMP.snow_base_mat = ones(4, size(ground.STATVAR.layerThick,2)); 
            ground.TEMP.snow_base_mat(2,:) = 2.* ground.TEMP.snow_base_mat(2,:); 
            ground.TEMP.snow_base_mat(3,:) = 3.* ground.TEMP.snow_base_mat(3,:); 
            ground.TEMP.snow_base_mat(4,:) = 4.* ground.TEMP.snow_base_mat(4,:); 
            
            ground.STATVAR.layerThick_snow = zeros(4, size(ground.STATVAR.layerThick,2));
            ground.STATVAR.ice_snow = ground.STATVAR.layerThick_snow;
            ground.STATVAR.upper_cell = repmat(4, 1, size(ground.STATVAR.layerThick,2));
            
            ground.STATVAR.thermCond = ground.STATVAR.layerThick.*0; % grid cell property
            ground.STATVAR.thermCond_eff = ground.STATVAR.layerThick(1:end-1,:).*0; %between grid cells property
            ground = conductivity(ground);
            
            ground.TEMP.FT_count = 0;
            ground.TEMP.count = 0;
  
            ground.STATVAR.FT_state = tile.PARA.geothermal .* NaN;
            
            ground.TEMP.d_energy = ground.STATVAR.energy .* 0;

       end
        
        
        %-----mandatory functions------------------------
        function ground = get_boundary_condition_u(ground, tile)
            
            ground.TEMP.timestep = tile.timestep;
            tile.timestep = tile.timestep .* (0.25 + 1.5*rand(1)); % max(1,tile.timestep + randn(1).*tile.timestep/2); %tile.timestep .* (0.5 + rand(1))
            
            if size(ground.STATVAR.layerThick, 2) == ground.PARA.speedup_size
                [ground.STATVAR.surf_T, ground.STATVAR.melt, ground.STATVAR.snowfall, ground.STATVAR.T, ground.TEMP.d_energy, ground.TEMP.d_new_snow_layerThick, ground.TEMP.d_new_melt_layerThick] = ...
                    get_boundary_condition_u_compiled_mex(ground.STATVAR.T, tile.FORCING.TEMP.surfT, tile.FORCING.TEMP.melt, ...
                    ground.CONST.day_sec, tile.timestep, tile.FORCING.TEMP.snowfall, ground.STATVAR.ice_snow, ...
                    ground.STATVAR.upper_cell, ground.TEMP.d_energy,ground.TEMP.snow_mat1, ground.CONST.L_f, ground.CONST.c_i, ground.STATVAR.layerThick_snow, ground.PARA.wind_speed_class);
            else
                
                
                ground.STATVAR.surf_T = tile.FORCING.TEMP.surfT;
                
                ground.STATVAR.melt = tile.FORCING.TEMP.melt ./ 1000 ./ground.CONST.day_sec .* tile.timestep;  %in [m], constant timestep
                
                ground.STATVAR.snowfall = tile.FORCING.TEMP.snowfall ./1000 ./ ground.CONST.day_sec .* tile.timestep;
                ground.STATVAR.melt = min(ground.STATVAR.melt, sum(ground.STATVAR.ice_snow,1)+ ground.STATVAR.snowfall); %limit the melt, so that it doesn't exceed the existing snow
                ground.STATVAR.surf_T = double((ground.STATVAR.snowfall - ground.STATVAR.melt + sum(ground.STATVAR.ice_snow, 1)) <=0 | (ground.STATVAR.surf_T <0)) .* ground.STATVAR.surf_T; %set to zero if there is snow and T is positive
                
                for i=1:4 %assign boundary condition T to correct cell and all cells above
                    ground.STATVAR.T(i, :) =  ground.STATVAR.T(i, :) + double(i<=ground.STATVAR.upper_cell) .* (ground.STATVAR.surf_T - ground.STATVAR.T(i, :));
                end
                
                ground.TEMP.d_energy(1:4,:) = ground.TEMP.d_energy(1:4,:) + ground.TEMP.snow_mat1 .*(repmat(ground.STATVAR.snowfall./tile.timestep .* ...
                    (-ground.CONST.L_f + ground.STATVAR.surf_T .* ground.CONST.c_i), 4, 1)  - repmat(ground.STATVAR.melt./tile.timestep, 4, 1) ...
                    .* (-ground.CONST.L_f + ground.STATVAR.T(2:5,:) .* ground.CONST.c_i));
                
                new_snow_density = get_snow_density(ground, tile);
                
                %melt > snowfall: no increase with new snow density (all new snow melts), decrease with existing snow density
                %melt < snowfall: no net melt,
                
                ground.TEMP.d_new_snow_layerThick = repmat(double(ground.STATVAR.melt < ground.STATVAR.snowfall) .* (ground.STATVAR.snowfall - ground.STATVAR.melt) .* 920 ./ new_snow_density , 4, 1) .* ground.TEMP.snow_mat1;
                ground.TEMP.d_new_melt_layerThick = repmat(double(ground.STATVAR.melt > ground.STATVAR.snowfall) .* (ground.STATVAR.melt - ground.STATVAR.snowfall), 4, 1) .* ground.TEMP.snow_mat1 .* ground.STATVAR.layerThick_snow ./ max(1e-10, ground.STATVAR.ice_snow);
                
                
                %flux is assigned in get_derivatives_prognostic
            end
        end
        
        function ground = get_boundary_condition_l(ground,  tile)
%             ground.TEMP.d_energy(end,:) = ground.TEMP.d_energy(end,:) + repmat(tile.PARA.geothermal', 1, tile.ENSEMBLE.PARA.ensemble_size); 
            ground.TEMP.d_energy(end,:) = ground.TEMP.d_energy(end,:) + tile.PARA.geothermal; 
        end
        
        %calculate spatial derivatives
        function ground = get_derivatives_prognostic(ground, tile)
            
            if size(ground.STATVAR.layerThick, 2) == ground.PARA.speedup_size
                ground.TEMP.d_energy = dE_dt_compiled_mex(ground.TEMP.d_energy, ground.STATVAR.thermCond_eff, ground.STATVAR.T, ground.STATVAR.layerDistance);
                ground.TEMP.dLayerThick_compaction = compact_windDrift_compiled_mex(ground.STATVAR.T(2:5,:), ground.STATVAR.ice_snow, ground.STATVAR.layerThick_snow, ground.TEMP.snow_mat1, ground.PARA.wind_compaction_timescale, ground.CONST.g);
            else
                %downwards flux
                ground.TEMP.d_energy = ground.TEMP.d_energy - ground.STATVAR.thermCond_eff.*(ground.STATVAR.T(2:end,:)-ground.STATVAR.T(1:end-1,:))./ground.STATVAR.layerDistance;
                %upwards flux, lower boundary already added
                ground.TEMP.d_energy(1:end-1,:) = ground.TEMP.d_energy(1:end-1,:) + ground.STATVAR.thermCond_eff(2:end,:).*(ground.STATVAR.T(3:end,:)-ground.STATVAR.T(2:end-1,:))./ground.STATVAR.layerDistance(2:end,:);
                
                ground.TEMP.dLayerThick_compaction = compact_windDrift(ground, tile);
            end
            
            ground.TEMP.d_energy(1:3,:) = ground.TEMP.d_energy(1:3,:) .*  ground.TEMP.snow_mat2(1:3,:);
            
        end
        
        %prognostic step - integrate prognostic variables in time
        function ground = advance_prognostic(ground, tile)
             
            ground.STATVAR.ice_snow = ground.STATVAR.ice_snow + ground.TEMP.snow_mat1 .*repmat(ground.STATVAR.snowfall - ground.STATVAR.melt, 4,1);
            ground.STATVAR.ice_snow(ground.STATVAR.ice_snow<0) = 0;
            ground.STATVAR.layerThick_snow =  ground.STATVAR.layerThick_snow + (ground.TEMP.d_new_snow_layerThick - ground.TEMP.d_new_melt_layerThick);
            ground.STATVAR.layerThick_snow(ground.STATVAR.layerThick_snow<0) = 0;
            ground.STATVAR.layerThick_snow(ground.STATVAR.ice_snow == 0) = 0;

            ground.STATVAR.layerThick_snow = ground.STATVAR.layerThick_snow + ground.TEMP.dLayerThick_compaction.*tile.timestep;
            ground.STATVAR.layerThick_snow(ground.STATVAR.layerThick_snow<0) = 0;
            
            if size(ground.STATVAR.layerThick, 2) == ground.PARA.speedup_size
                ground.STATVAR.energy = advance_E_compiled_mex(ground.STATVAR.energy, ground.TEMP.d_energy, tile.timestep);
            else
                ground.STATVAR.energy = ground.STATVAR.energy + ground.TEMP.d_energy .* tile.timestep;
            end
            
            for i=1:3 %set energy to 0 for unused cells
                ground.STATVAR.energy(i,:) = double(i >= ground.STATVAR.upper_cell) .* ground.STATVAR.energy(i,:);
            end
            

        end
        
        %diagnostic step - compute diagnostic variables
        function ground = compute_diagnostic(ground, tile)
            
            %diagnostic step, the 4 is the first ground cell
            ground = regrid_snow(ground, tile);
            
%             if sum(sum(ground.STATVAR.layerThick_snow == Inf))>0 || sum(sum(ground.STATVAR.layerThick_snow == -Inf))>0 || sum(sum(ground.STATVAR.layerThick_snow <0))>0
%                 'Hallo3'
%                 dhefhk
%             end

            ground.STATVAR.layerThick(2:4,:) =  ground.STATVAR.layerThick_snow(1:3,:) + double(ground.STATVAR.layerThick_snow(1:3,:) == 0) .* ground.PARA.virtual_gridCellSize; %set to snow cell size if snow, virtual_gridCellSize, otherwise
            ground.STATVAR.layerThick(5,:) = ground.STATVAR.layerThick_first_ground_cell + ground.STATVAR.layerThick_snow(4,:); %combined snow and ground cell
            ground.STATVAR.layerDistance(1:5,:) = ground.STATVAR.layerThick(2:6,:)./2 + ground.STATVAR.layerThick(1:5,:)./2; %recompute first five cells which are affected by snow; rest is constant

            ground = get_T(ground);
            
%             if sum(isnan(ground.STATVAR.energy(:)))>0 || sum(isnan(ground.STATVAR.T(:)))>0
%                 disp('T')
%             end
             
            ground = conductivity(ground);
            
            ground.TEMP.FT_count = ground.TEMP.FT_count + double(ground.STATVAR.T(5:35,:) <0); %only check until 10m for talik
            ground.TEMP.count = ground.TEMP.count + 1;
            
            ground.TEMP.d_energy = ground.STATVAR.energy .* 0;
        end
        
        %triggers
%         function ground = check_trigger(ground, tile)
%             
%             if tile.t >= ground.TEMP.adjust_stratigraphy_date
% 
%                 FT_code = ground.TEMP.FT_count ./ ground.TEMP.count;
%                 FT_code(FT_code==1)=2;  %frozen
%                 FT_code(FT_code>0 & FT_code < 1)=1;  %freeze thaw
%                 
%                 gain_loose = FT_code.*0;
%                 gain_loose(FT_code==2) = 1;  %gain when frozen
%                 %OUT.STATVAR.gain_loose(OUT.ACC.FT_isothermal==1) = 0; %no gain when isothermal
%                 gain_loose(FT_code == 1 | FT_code == 0) = -1 ;  %loose when unfrozen or FT, this depends on water table settings for the ensemble member
%                 for i=1:size(gain_loose,1)-1  %set cell above AL to gain
%                     gain_loose(i,:) = gain_loose(i,:) + double(FT_code(i,:)==1 & FT_code(i+1,:)==2) .* (1 - gain_loose(i,:));
%                 end
%                 
%                 FT_code=[ones(1,size(FT_code,2)); FT_code];
%                 FT_code = FT_code(1:end-1,:) - FT_code(2:end,:);
%                 
%                 FT_res = zeros(1, size(FT_code,2));
%                  for i = 1: size(FT_code,1)
%                      FT_res = FT_res + double (FT_res == 0 & FT_code(i,:) ==-1) .* -1 + double (FT_res == 0 & FT_code(i,:) == 1) .* 1; %initial PF yes no, -1 or 1
%                      FT_res = FT_res + double (FT_res >0 & FT_code(i,:) < 0) .* -FT_code(i,:); %switches from no PF, 1, to freeze_thaw or frozen, so 2 or 3 means talik
%                  end
%                  
%                 ground.STATVAR.FT_state = FT_res;
%                 
%                 %---------------update stratigrahy----
% 
%                 %update_stratigraphy
%                 ground = update_stratigraphy(ground, tile, gain_loose);
%                 
%                 ground.TEMP.FT_count = 0;   
%                 ground.TEMP.count = 0;
%                 ground.TEMP.adjust_stratigraphy_date = datenum([ground.PARA.adjust_stratigraphy_date num2str(str2num(datestr(tile.t, 'yyyy'))+1) ], 'dd.mm.yyyy');
%             end
%             
%         end

        %triggers
        function ground = check_trigger(ground, tile)
            
            tile.timestep = ground.TEMP.timestep;
            
            if tile.t >= ground.TEMP.adjust_stratigraphy_date

                FT_code = ground.TEMP.FT_count ./ ground.TEMP.count;
                FT_code(FT_code==1)=2;  %frozen
                FT_code(FT_code>0 & FT_code < 1)=1;  %freeze thaw
                
                gain_loose = FT_code.*0;
                gain_loose(FT_code==2) = 1;  %gain when frozen
                %OUT.STATVAR.gain_loose(OUT.ACC.FT_isothermal==1) = 0; %no gain when isothermal
                gain_loose(FT_code == 1 | FT_code == 0) = -1 ;  %loose when unfrozen or FT, this depends on water table settings for the ensemble member
%                 for i=1:size(gain_loose,1)-1  %set cell above AL to gain
%                     gain_loose(i,:) = gain_loose(i,:) + double(FT_code(i,:)==1 & FT_code(i+1,:)==2) .* (1 - gain_loose(i,:)); %& ~(OUT.ACC.FT_isothermal(i+1,:)==1)) ;
%                    % gain_loose(i,:) = gain_loose(i,:) + double(FT_code(i,:)==1 & FT_code(i+1,:)==2) .* (1 - gain_loose(i,:)); %& ~(OUT.ACC.FT_isothermal(i+1,:)==1)) ;
%                 end
                
                
                FT_code=[ones(1,size(FT_code,2)); FT_code];
                FT_code = FT_code(1:end-1,:) - FT_code(2:end,:);
                
                first_frozen_cell = 0;
                
                FT_res = zeros(1, size(FT_code,2));
                 for i = 1:size(FT_code,1)
                     first_frozen_cell = first_frozen_cell + i.* double(first_frozen_cell == 0 & FT_code(i,:) <=-1); %first frozen cell
                     FT_res = FT_res + double (FT_res == 0 & FT_code(i,:) ==-1) .* -1 + double (FT_res == 0 & FT_code(i,:) == 1) .* 1; %initial PF yes no, -1 or 1
                     FT_res = FT_res + double (FT_res >0 & FT_code(i,:) < 0) .* -FT_code(i,:); %switches from no PF, 1, to freeze_thaw or frozen, so 2 or 3 means talik
                 end
                first_frozen_cell(first_frozen_cell == 0) = i; %no frozen cell
                 
                ground.STATVAR.FT_state = FT_res;
                
                %---------------update stratigrahy----
                
                %update_stratigraphy
                ground = update_stratigraphy(ground, tile, gain_loose, first_frozen_cell);
                
                ground.TEMP.FT_count = 0;   
                ground.TEMP.count = 0;
                ground.TEMP.adjust_stratigraphy_date = datenum([ground.PARA.adjust_stratigraphy_date num2str(str2num(datestr(tile.t, 'yyyy'))+1) ], 'dd.mm.yyyy');
            end
            
        end
        
        
        
        %----non-mandatory functions
        
        
        function ground = get_T(ground) 
            
            if size(ground.STATVAR.layerThick, 2) == ground.PARA.speedup_size
                ground.STATVAR.T(6:end,:) = get_Tground_compiled_mex(ground.STATVAR.energy, ground.STATVAR.c_thawed, ground.STATVAR.c_frozen, ground.STATVAR.E_frozen, ground.STATVAR.T_end_freezing);
            else
                %1. ground without first cell
                ground.STATVAR.T(6:end,:) = double(ground.STATVAR.energy(5:end,:)>=0) .* ground.STATVAR.energy(5:end,:) ./ ground.STATVAR.c_thawed(2:end,:) + ...
                    double(ground.STATVAR.energy(5:end,:) <= ground.STATVAR.E_frozen(2:end,:)) .* ((ground.STATVAR.energy(5:end,:) - ground.STATVAR.E_frozen(2:end,:)) ./ ground.STATVAR.c_frozen(2:end,:) + ground.STATVAR.T_end_freezing(2:end,:)) + ...
                    double(ground.STATVAR.energy(5:end,:) < 0 & ground.STATVAR.energy(5:end,:) > ground.STATVAR.E_frozen(2:end,:)) .* ground.STATVAR.energy(5:end,:)./ground.STATVAR.E_frozen(2:end,:) .*(ground.STATVAR.T_end_freezing(2:end,:));
            end
            
            %2. first ground cell including initial snow - free water
            %freeze curve here, makes things much easier with the snow

            E_frozen_first_ground_cell = ground.STATVAR.E_frozen(1,:) - ground.CONST.L_f .* ground.STATVAR.ice_snow(4,:);
            %c_frozen_first_cell = PROFILE.c_frozen(1,:) +  PROFILE.const.c_i .* PROFILE.D_ice_snow(4,:);
            
            ground.STATVAR.T(5,:) = double (double(ground.STATVAR.energy(4,:)>=0)) .* ground.STATVAR.energy(4,:) ./ (ground.STATVAR.c_thawed(1,:) + ground.CONST.c_w .* ground.STATVAR.ice_snow(4,:)) + ...
                double(ground.STATVAR.energy(4,:) <= E_frozen_first_ground_cell) .* ( (ground.STATVAR.energy(4,:) - E_frozen_first_ground_cell) ./ (ground.STATVAR.c_frozen(1,:) +  ground.CONST.c_i .* ground.STATVAR.ice_snow(4,:)) ); 
            %zero degrees else

            %3. snow
            T_snow = double(ground.STATVAR.energy(1:3,:) < -ground.CONST.L_f .* ground.STATVAR.ice_snow(1:3,:)) .* (ground.STATVAR.energy(1:3,:) + ground.CONST.L_f .* ground.STATVAR.ice_snow(1:3,:)) ./ ...
                (ground.CONST.c_i .*ground.STATVAR.ice_snow(1:3,:));
            T_snow(isnan(T_snow))=0;
            T_snow(abs(T_snow)==Inf)=0; 
            ground.STATVAR.T(2:4,:) = T_snow;
            
        end
        
        
        function ground = conductivity(ground)
            
            if size(ground.STATVAR.layerThick, 2) == ground.PARA.speedup_size
                ground.STATVAR.thermCond(5:end,:) = thermCond_compiled_mex(ground.STATVAR.T, ground.STATVAR.T_end_freezing, ground.STATVAR.k_frozen, ground.STATVAR.k_freezing, ground.STATVAR.k_thawed);
                ground.STATVAR.thermCond_eff(5:end,:) = thermCond_eff_compiled_mex(ground.STATVAR.thermCond, ground.STATVAR.layerThick);
                ground.STATVAR.thermCond_snow = thermCond_snow_compiled_mex(ground.STATVAR.ice_snow, ground.STATVAR.layerThick_snow);
            else
                ground.STATVAR.thermCond(5:end,:) = double(ground.STATVAR.T(5:end,:) < ground.STATVAR.T_end_freezing) .* ground.STATVAR.k_frozen + double(ground.STATVAR.T(5:end,:) > 0) .* ground.STATVAR.k_thawed + ...
                    double(ground.STATVAR.T(5:end,:) >= ground.STATVAR.T_end_freezing & ground.STATVAR.T(5:end,:) <= 0) .* ground.STATVAR.k_freezing; % ground
                
                ground.STATVAR.thermCond_eff(5:end,:) = ground.STATVAR.thermCond(5:end-1,:).*ground.STATVAR.thermCond(6:end,:).*...
                    (ground.STATVAR.layerThick(5:end-1,:)./2 + ground.STATVAR.layerThick(6:end,:)./2) ./ (ground.STATVAR.thermCond(5:end-1,:).*ground.STATVAR.layerThick(6:end,:)./2 + ...
                    ground.STATVAR.thermCond(6:end,:).*ground.STATVAR.layerThick(5:end-1,:)./2 ); %size N
                
                %Thermal conductivity snow
                %snow_density = ground.STATVAR.ice_snow ./ max(1e-20, ground.STATVAR.layerThick_snow) .*920;
                snow_density = min(1, ground.STATVAR.ice_snow ./ max(1e-20, ground.STATVAR.layerThick_snow));
                
                %ground.STATVAR.thermCond_snow = max(5e-2, 2.3.*(snow_density./1000).^1.88);
                ground.STATVAR.thermCond_snow = max(5e-2, 2.3.*snow_density.^1.88);
                
            end
            ground.STATVAR.thermCond(2:4,:) = ground.STATVAR.thermCond_snow(1:3,:); %snow
            
            %replace conductivities above upper boundary by some high
            %value, this ensures that it is possible to divide by
            %layerDistance, and that no exception must be made for 1st
            %cell - equivalent to setting k_eff to conductivity of uppermost
            %cell and dividing by half the grid cell size.
            for i=1:4
                ground.STATVAR.thermCond(i,:) = ground.STATVAR.thermCond(i,:) + double(i <= ground.STATVAR.upper_cell) .* (100 - ground.STATVAR.thermCond(i,:));
            end

            %thermal conductivity first cell including snow, only applied
            %upwards, same as terhmal conductivity uppermost cell when
            %there is no snow
            ground.STATVAR.thermCond(5,:) = ground.STATVAR.thermCond(5,:).*ground.STATVAR.thermCond_snow(4,:).*(ground.STATVAR.layerThick_first_ground_cell./2 + ground.STATVAR.layerThick_snow(4,:)) ./ ...
                (ground.STATVAR.thermCond(5,:) .* ground.STATVAR.layerThick_snow(4,:) + ground.STATVAR.thermCond_snow(4,:) .* ground.STATVAR.layerThick_first_ground_cell./2) ; %first half cell plus snow, modified after Mamoru, error corrected
            
            ground.STATVAR.thermCond_eff(1:4,:) = ground.STATVAR.thermCond(1:4,:).*ground.STATVAR.thermCond(2:5,:).*(ground.STATVAR.layerThick(1:4,:)./2 + ground.STATVAR.layerThick(2:5,:)./2) ./ ...
                (ground.STATVAR.thermCond(1:4,:).*ground.STATVAR.layerThick(2:5,:)./2 + ground.STATVAR.thermCond(2:5,:).*ground.STATVAR.layerThick(1:4,:)./2 ); %size N
            
           %NEW
            for i=1:3 %set higher thermal conductivity for influx in first snow cell 
                ground.STATVAR.thermCond_eff(i, :) =  ground.STATVAR.thermCond_eff(i, :) + double(i<=ground.STATVAR.upper_cell) .* (2 - ground.STATVAR.thermCond_eff(i, :));
            end
            
            
%            ground.STATVAR.thermCond_eff = ground.STATVAR.thermCond_eff .* (1 + (rand(size(ground.STATVAR.thermCond_eff,1), size(ground.STATVAR.thermCond_eff,2)) - 0.5)./50);
        end
        
        
         
        function ground = get_T_end_freezing(ground)
            ground.STATVAR.T_end_freezing = double(ground.STATVAR.soil_type==1).*-0.1+ double(ground.STATVAR.soil_type==2).*-1;
            ground.STATVAR.T_end_freezing(1,:) = 0; %set first cell to zero, this makes computation of combined snow cover and ground cell easier
        end
        
        
        function ground = T2E(ground)
            E_frozen = - ground.CONST.L_f.*ground.STATVAR.waterIce + ground.STATVAR.T_end_freezing.*...
                (ground.CONST.c_w.*ground.STATVAR.waterIce./2 + ground.CONST.c_i.*ground.STATVAR.waterIce./2 + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic);
            
            ground.STATVAR.energy = double(ground.STATVAR.T>=ground.STATVAR.T_onset_freezing) .* ground.STATVAR.T .* ...
                (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic) + ...
                double(ground.STATVAR.T<=ground.STATVAR.T_end_freezing) .* ((ground.STATVAR.T-ground.STATVAR.T_end_freezing) .* (ground.CONST.c_i.*ground.STATVAR.waterIce...
                + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic) + E_frozen) + ...
                double(ground.STATVAR.T < ground.STATVAR.T_onset_freezing & ground.STATVAR.T > ground.STATVAR.T_end_freezing) .* ground.STATVAR.T./min(ground.STATVAR.T_end_freezing, -1e-12) .*E_frozen;
            %ground.E = ground.E(2:end,:);
            %ground.STATVAR.energy = ground.STATVAR.energy .* ground.STATVAR.layerThick(5:end,1); %absolute energy
            ground.STATVAR.energy=[zeros(3, size(ground.STATVAR.energy, 2)); ground.STATVAR.energy];  %add three snow cells with zero energy
            
        end
        
        function ground = init_conductivity(ground)
            waterIce = ground.STATVAR.waterIce ./ ground.STATVAR.layerThick;
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick;
            organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick;
%             
%             waterIce = ground.STATVAR.waterIce ./ ground.STATVAR.layerThick(5:end,:);
%             mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick(5:end,:);
%             organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick(5:end,:);
            air = 1 - waterIce - mineral - organic;
            ground.STATVAR.k_frozen = ( waterIce.* ground.CONST.k_i.^0.5 + mineral.* ground.CONST.k_m.^0.5 + organic.* ground.CONST.k_o.^0.5 + air.* ground.CONST.k_a.^0.5).^2;
            ground.STATVAR.k_thawed = (waterIce.* ground.CONST.k_w.^0.5 + mineral.* ground.CONST.k_m.^0.5 + organic.* ground.CONST.k_o.^0.5 + air.* ground.CONST.k_a.^0.5).^2;
            ground.STATVAR.k_freezing = (ground.STATVAR.k_frozen + ground.STATVAR.k_thawed)./2;
        end
        
        function ground = init_conductivity2(ground)
            waterIce = ground.STATVAR.waterIce ./ ground.STATVAR.layerThick(5:end,:);
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick(5:end,:);
            organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick(5:end,:);
%             
%             waterIce = ground.STATVAR.waterIce ./ ground.STATVAR.layerThick(5:end,:);
%             mineral = ground.STATVAR.mineral ./ ground.STATVAR.layerThick(5:end,:);
%             organic = ground.STATVAR.organic ./ ground.STATVAR.layerThick(5:end,:);
            air = 1 - waterIce - mineral - organic;
            ground.STATVAR.k_frozen = ( waterIce.* ground.CONST.k_i.^0.5 + mineral.* ground.CONST.k_m.^0.5 + organic.* ground.CONST.k_o.^0.5 + air.* ground.CONST.k_a.^0.5).^2;
            ground.STATVAR.k_thawed = (waterIce.* ground.CONST.k_w.^0.5 + mineral.* ground.CONST.k_m.^0.5 + organic.* ground.CONST.k_o.^0.5 + air.* ground.CONST.k_a.^0.5).^2;
            ground.STATVAR.k_freezing = (ground.STATVAR.k_frozen + ground.STATVAR.k_thawed)./2;
        end
        
        function ground = calculate_E_frozen (ground) %new function, only call this once in the beginning
            
            ground.STATVAR.T_end_freezing(1,:) = 0; %first cell freezes like free water
            ground.STATVAR.E_frozen = - ground.CONST.L_f .* ground.STATVAR.waterIce + ground.STATVAR.T_end_freezing.* ...
                (ground.CONST.c_w.*ground.STATVAR.waterIce./2 + ground.CONST.c_i.*ground.STATVAR.waterIce./2 + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic);
            
            ground.STATVAR.c_thawed = (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic); % .* ground.STATVAR.layerThick(5:end,:); %unit J/m2/K
            ground.STATVAR.c_frozen = (ground.CONST.c_i.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic); % .* ground.STATVAR.layerThick(5:end,:); %unit J/m2/K

        end
        
        %--------------

%         function ground = update_stratigraphy(ground, tile, gain_loose_in)
%         
%             layerThick = [ground.STATVAR.layerThick_first_ground_cell; ground.STATVAR.layerThick(6:end,:)];
%             depths = cumsum(layerThick); %same size as theta_w
%             gain_loose = depths.*0;
%             gain_loose(1:31,:) = gain_loose_in;
%             porosity = 1 - ground.STATVAR.mineral ./layerThick  - ground.STATVAR.organic ./layerThick;
%             field_capacity = 0.5.* porosity;   %CHANGE later
%             
%             water_table_depth = zeros(1, size(field_capacity,2));
%             for i=1:31
%                 water_table_depth = water_table_depth + double(gain_loose(i,:)==-1 & gain_loose(i+1,:)~=-1) .* (depths(i,:) - water_table_depth);
%             end
% %             water_table_depth = water_table_depth .* (1-tile.ENSEMBLE.STATVAR.water_table_depth);
%             water_table_depth = water_table_depth .* (1 - ground.PARA.water_table_depth); %this is a PARA of GROUND NOW, which si set by the ENSEMBLE class
%                         
%             for i=1:31
%                 gain_loose(i,:) = gain_loose(i,:) + double(depths(i,:) > water_table_depth) .* (1-gain_loose(i,:));
%             end
%             T_old = ground.STATVAR.T;
%             water_old = ground.STATVAR.waterIce; %in m
%             ground.STATVAR.waterIce = max(field_capacity, min(porosity, ground.STATVAR.waterIce ./layerThick + gain_loose .* 0.05)) .* layerThick;
%             water_change = ground.STATVAR.waterIce - water_old; %iin m
%             
%             ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) > 0) .* ground.STATVAR.T(5:end,:) .* ground.CONST.c_w .* water_change;
%             ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) < ground.STATVAR.T_end_freezing) .* (ground.STATVAR.T(5:end,:) .* ground.CONST.c_i - ground.CONST.L_f) .* water_change;
%             E_frozen_change = - ground.CONST.L_f .* water_change + ground.STATVAR.T_end_freezing.*(ground.CONST.c_w .* water_change./2 + ground.CONST.c_i .* water_change./2);
%             
%             ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) <=0 & ground.STATVAR.T(5:end,:)>= ground.STATVAR.T_end_freezing) .* ground.STATVAR.T(5:end,:)./min(ground.STATVAR.T_end_freezing, -1e-12) .* E_frozen_change;
%             ground.STATVAR.E_frozen = - ground.CONST.L_f.*ground.STATVAR.waterIce + ground.STATVAR.T_end_freezing.*(ground.CONST.c_w.*ground.STATVAR.waterIce./2 + ground.CONST.c_i.*ground.STATVAR.waterIce./2 + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic);
%             
%             ground.STATVAR.c_thawed = (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic) ; %unit J/m2/K
%             ground.STATVAR.c_frozen = (ground.CONST.c_i.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic); %unit J/m2/K
%             
%             ground = get_T(ground);
%             store_layerThick = ground.STATVAR.layerThick;
%             ground.STATVAR.layerThick = layerThick;
%             ground = init_conductivity(ground);
%             ground.STATVAR.layerThick = store_layerThick;
% 
%             %re-initialzie conductivities
%             ground = init_conductivity2(ground);
%             
%             ground = conductivity(ground);
%             %test= PROFILE.water_table_depth;
%             %save('save_gain_loose.mat', 'water_table_depth', 'gain_loose', 'water_change', 'T_old', 'T_new')
%             
%         end
        
         function ground = update_stratigraphy(ground, tile, gain_loose_in, first_frozen_cell)
        
            layerThick = [ground.STATVAR.layerThick_first_ground_cell; ground.STATVAR.layerThick(6:end,:)];
            depths = cumsum(layerThick); %same size as theta_w
            gain_loose = depths.*0;
            gain_loose(1:31,:) = gain_loose_in;
            porosity = 1 - ground.STATVAR.mineral ./layerThick  - ground.STATVAR.organic ./layerThick;
            field_capacity = 0.5.* porosity;   %CHANGE later
            
            water_table_depth = zeros(1, size(field_capacity,2));
%             for i=1:31
%                 water_table_depth = water_table_depth + double(gain_loose(i,:)==-1 & gain_loose(i+1,:)~=-1) .* (depths(i,:) - water_table_depth);
%             end
            for i=1:31
                water_table_depth = water_table_depth + double(i==first_frozen_cell).* depths(i,:);
            end
%             water_table_depth = water_table_depth .* (1-tile.ENSEMBLE.STATVAR.water_table_depth);
            water_table_depth = water_table_depth .* (1 - ground.PARA.water_table_depth); %this is a PARA of GROUND NOW, which si set by the ENSEMBLE class
                        
            for i=1:31
                gain_loose(i,:) = gain_loose(i,:) + double(depths(i,:) > water_table_depth) .* (1-gain_loose(i,:));
              % gain_loose(i,:) = gain_loose(i,:) + double(depths(i,:) > water_table_depth) .* (1-gain_loose(i,:));
            end
            T_old = ground.STATVAR.T;
            water_old = ground.STATVAR.waterIce; %in m
            ground.STATVAR.waterIce = max(field_capacity, min(porosity, ground.STATVAR.waterIce ./layerThick + gain_loose .* 0.025)) .* layerThick;
            water_change = ground.STATVAR.waterIce - water_old; %iin m
            
%             ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) > 0) .* ground.STATVAR.T(5:end,:) .* ground.CONST.c_w .* water_change;
%             ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) < ground.STATVAR.T_end_freezing) .* (ground.STATVAR.T(5:end,:) .* ground.CONST.c_i - ground.CONST.L_f) .* water_change;
%             E_frozen_change = - ground.CONST.L_f .* water_change + ground.STATVAR.T_end_freezing.*(ground.CONST.c_w .* water_change./2 + ground.CONST.c_i .* water_change./2);
%             
%             ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) <=0 & ground.STATVAR.T(5:end,:)>= ground.STATVAR.T_end_freezing) .* ground.STATVAR.T(5:end,:)./min(ground.STATVAR.T_end_freezing, -1e-12) .* E_frozen_change;
            ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) > 0) .* ground.STATVAR.T(5:end,:) .* ground.CONST.c_w .* water_change;
            ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) < ground.STATVAR.T_end_freezing) .* ...
                (ground.STATVAR.T_end_freezing .* (ground.CONST.c_w./2 + ground.CONST.c_i./2) + (ground.STATVAR.T(5:end,:)-ground.STATVAR.T_end_freezing) .* ground.CONST.c_i - ground.CONST.L_f) .* water_change;
            
            ice_fraction = ground.STATVAR.energy(4:end,:)./ min(-1e-12, ground.STATVAR.E_frozen);
            ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) <=0 & ground.STATVAR.T(5:end,:)>= ground.STATVAR.T_end_freezing) .*...
                (ground.STATVAR.T(5:end,:) .* water_change .* (ground.CONST.c_w ./2 + ground.CONST.c_i./2) - ground.CONST.L_f .* water_change .* ice_fraction);
            
            ground.STATVAR.E_frozen = - ground.CONST.L_f.*ground.STATVAR.waterIce + ground.STATVAR.T_end_freezing.*(ground.CONST.c_w.*ground.STATVAR.waterIce./2 + ground.CONST.c_i.*ground.STATVAR.waterIce./2 + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic);
            ground.STATVAR.c_thawed = (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic) ; %unit J/m2/K
            ground.STATVAR.c_frozen = (ground.CONST.c_i.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic); %unit J/m2/K
            
            ground = get_T(ground);
            store_layerThick = ground.STATVAR.layerThick;
            ground.STATVAR.layerThick = layerThick;
            ground = init_conductivity(ground);
            ground.STATVAR.layerThick = store_layerThick;

            %re-initialzie conductivities
            ground = init_conductivity2(ground);
            
            ground = conductivity(ground);
            
        end
        
 
        
        function rho_snow = get_snow_density(ground, tile)
            a = 109;
            b = 6;
            c = 26;
            T=min(0,ground.STATVAR.surf_T);
            rho_snow = max(50, a + b.*T + (c .* ground.PARA.wind_speed_class.^0.5));
            
%             T_air = min(0,ground.STATVAR.surf_T);
%             rho_Tair = double(T_air > -15).*(50 + 1.7.*(T_air+15).^(1.5)) ...
%                 + double(T_air <= -15).*( -3.8328.*T_air - 0.0333.*T_air.^2);
%             rho_wind = 266.861.*(0.5.*(1 + tanh(ground.PARA.wind_speed_class./5))).^8.8;
%             rho_snow = rho_Tair + rho_wind;
        end

        
      
        function dD_dt = compact_windDrift(ground, tile)

            T = min(0,ground.STATVAR.T(2:5,:));
            
            eta_0 = 7.62237e6;
            a_eta = 0.1;
            b_eta = 0.023;
            c_eta = 250;
                        
            rho_ice = 920;
            rho_max = 350;
            
            rho = ground.STATVAR.ice_snow ./ max(1e-20, ground.STATVAR.layerThick_snow) .* rho_ice;
            
            stress = ground.CONST.g .* rho_ice .* (cumsum(ground.STATVAR.ice_snow - ground.TEMP.snow_mat1 .* ground.STATVAR.ice_snow./2));             
            
            eta = eta_0 .*  rho ./ c_eta .* exp(-a_eta .* T + b_eta .* rho);
                        
            dD_dt =  - stress ./max(1e-10,eta) .* ground.STATVAR.layerThick_snow; %compaction            
            
            dD_dt = dD_dt - ground.TEMP.snow_mat1 .* rho_ice .* ground.STATVAR.ice_snow ./ max(1e-8, rho).^2.* (rho_max - min(rho, rho_max)) ./ (ground.PARA.wind_compaction_timescale .*24.*3600);
    

        end
        
        function ground = regrid_snow(ground, tile)
            
            D_ice = ground.STATVAR.ice_snow;
            D = ground.STATVAR.layerThick_snow;
            
            ice_content = D_ice./max(1e-10, D);
            
            target_SWE =0.05;
            
            D_tot=sum(D_ice,1);
            
            number_of_cells = min(3, floor(D_tot./ target_SWE));

            lower_cell = ground.TEMP.index_first_ground_cell - min(number_of_cells,1);
            ground.STATVAR.upper_cell = lower_cell - max(0, number_of_cells-1);
            for i=1:4
               ground.TEMP.snow_mat1(i,:) = double(ground.TEMP.snow_base_mat(i,:) == ground.STATVAR.upper_cell);
               ground.TEMP.snow_mat2(i,:) = double(ground.TEMP.snow_base_mat(i,:) >= ground.STATVAR.upper_cell & ground.TEMP.snow_base_mat(i,:) <= lower_cell);
            end
            
            factor = max(1, number_of_cells);
            
            test= ground.TEMP.snow_mat2.* repmat(D_tot./factor,4,1);
            
            snow_over = max(0, test(1,:) - target_SWE);
            test(1,:)= test(1,:)-snow_over;
            test(2,:)= test(2,:)+snow_over;
            
            d_D_ice = test - D_ice;
            d_D_ice_res = d_D_ice;
            %----------
            
            d_D = d_D_ice;
            
            d_D_res=zeros(4, size(ground.TEMP.index_first_ground_cell, 2));
            
            for i=4:-1:2
                
                d_D(i,:) = double(d_D_ice(i,:) > 0 & ice_content(i-1,:)>0 ) .* d_D_ice(i,:) ./ max(1e-10, ice_content(i-1,:)) + double(d_D_ice(i,:) < 0 & ice_content(i,:)>0) .* d_D_ice(i,:)./ max(1e-10, ice_content(i,:));
                d_D(i-1,:) = -d_D(i,:);
                d_D_ice(i-1,:) = d_D_ice(i-1,:) + d_D_ice(i,:);
                
                d_D_res(i,:) = d_D_res(i,:) + d_D(i,:);
                d_D_res(i-1,:) = d_D_res(i-1,:) + d_D(i-1,:);
                
            end

            d_D_res(d_D_ice_res==0)=0;
            
            ground.STATVAR.ice_snow = ground.STATVAR.ice_snow + d_D_ice_res;
            ground.STATVAR.layerThick_snow = ground.STATVAR.layerThick_snow + d_D_res;
            ground.STATVAR.energy(1:4,:) = ground.STATVAR.energy(1:4,:) + d_D_ice_res .* (-ground.CONST.L_f + ground.STATVAR.T(2:5,:) .* ground.CONST.c_i);
            
            ground.STATVAR.layerThick_snow(ground.STATVAR.layerThick_snow<0) = 0;
                        
%             if sum(sum(ground.STATVAR.layerThick_snow<0))>0
%                 'Hllo10'
%                 djkdfff
%             end
            
        end
       

    end
    
end

