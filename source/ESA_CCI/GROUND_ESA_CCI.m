
classdef GROUND_ESA_CCI < BASE
    

    methods
        
        %-----initialize-----------------
        
        function ground = provide_CONST(ground)
            ground.CONST.L_f = []; %3.34e8;
            ground.CONST.c_w = []; %4.2e6; %[J/m�K]
            ground.CONST.c_o = []; %2.5e6; %[J/m�K]
            ground.CONST.c_m = []; %2.0e6; %[J/m�K]
            ground.CONST.c_i = []; %1.9e6;%[J/m�K]
           
            ground.CONST.k_w = []; 
            ground.CONST.k_o = []; 
            ground.CONST.k_m = [];
            ground.CONST.k_i = []; 
            ground.CONST.k_a = [];
            
            ground.CONST.day_sec = []; %24*3600;
            
            ground.CONST.g = [];
        end
        
        function ground = provide_PARA(ground)
            ground.PARA.virtual_gridCellSize = [];
            ground.PARA.timestep = [];
            ground.PARA.adjust_stratigraphy_date = [];
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
            ground.TEMP.F_lb = repmat(tile.RUN_INFO.STATVAR.geothermal(tile.PARA.range,1)', 1, tile.ENSEMBLE.PARA.ensemble_size);
            
            
            ground.STATVAR.T = [ground.STATVAR.T; ground.STATVAR.organic(2:end,:).*0];
            
            for i=1:size(ground.STATVAR.T,1)-1
                k = double(ground.STATVAR.T(i,:) < ground.STATVAR.T_end_freezing(i,:)) .* ground.STATVAR.k_frozen(i,:) + double(ground.STATVAR.T(i,:) > ground.STATVAR.T_onset_freezing) .* ground.STATVAR.k_thawed(i,:) + ...
                    double(ground.STATVAR.T(i,:) >= ground.STATVAR.T_end_freezing(i,:) & ground.STATVAR.T(i,:) <= ground.STATVAR.T_onset_freezing) .* ground.STATVAR.k_freezing(i,:);
                ground.STATVAR.T(i+1,:) = ground.STATVAR.T(i,:) + ground.TEMP.F_lb .* ground.STATVAR.layerDistance(i,:) ./k;
            end

            
            ground.PARA.snow_gridCellSize = ground.PARA.virtual_gridCellSize; %CHECK and ADAPT
            ground.STATVAR.layerThick_first_ground_cell = ground.STATVAR.layerThick(1,:);
            
            %add four cells on top, one virtual and three snow cells
            ground.STATVAR.layerThick = [repmat(ground.PARA.virtual_gridCellSize,4, size(ground.STATVAR.layerThick,2)); ground.STATVAR.layerThick];
            ground.STATVAR.layerDistance = (ground.STATVAR.layerThick(1:end-1,:) + ground.STATVAR.layerThick(2:end,:)) ./ 2;
            
            %ground.STATVAR.k_mineral = repmat(ground.CONST.k_m, 1, tile.PARA.number_of_realizations);
            

            
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
  
            ground.STATVAR.FT_state = ground.TEMP.F_lb .* NaN;
            
            ground.TEMP.d_energy = ground.STATVAR.energy .* 0;
            
%             ground.TEMP.T_store = [];
%             ground.TEMP.E_store = [];
%             ground.TEMP.lt_store =[];
        end
        
        
        %-----mandatory functions------------------------
        function ground = get_boundary_condition_u(ground, tile)
            
            %has to be distributed correctly as soon as ensemble is
            %available
            ground.STATVAR.surf_T = repmat(tile.FORCING.TEMP.surfT, 1, tile.ENSEMBLE.PARA.ensemble_size);

            ground.STATVAR.melt = repmat(tile.FORCING.TEMP.melt_bare, 1, tile.ENSEMBLE.PARA.ensemble_size) .* (1 - tile.ENSEMBLE.STATVAR.melt_fraction) + ...
                repmat(tile.FORCING.TEMP.melt_forest, 1, tile.ENSEMBLE.PARA.ensemble_size) .* tile.ENSEMBLE.STATVAR.melt_fraction;
            ground.STATVAR.melt = ground.STATVAR.melt ./ 1000 ./ground.CONST.day_sec .* tile.timestep;  %in [m], constant timestep
            
            
            ground.STATVAR.snowfall = repmat(tile.FORCING.TEMP.snowfall, 1, tile.ENSEMBLE.PARA.ensemble_size) ...
                ./1000 ./ ground.CONST.day_sec .* tile.ENSEMBLE.STATVAR.snowfall_factor .* tile.timestep;
            
%             ground.STATVAR.melt = min(ground.STATVAR.melt, ground.STATVAR.snowfall + sum(ground.STATVAR.ice_snow,1)); %limit the melt, so that it doesn't exceed the existing snow
%             ground.STATVAR.surf_T = double ((ground.STATVAR.snowfall - ground.STATVAR.melt + sum(ground.STATVAR.ice_snow, 1)) <=0 | (ground.STATVAR.surf_T <0)) .* ground.STATVAR.surf_T; %set to zero if there is snow and T is positive
            
            ground.STATVAR.melt = min(ground.STATVAR.melt, sum(ground.STATVAR.ice_snow,1)+ ground.STATVAR.snowfall); %limit the melt, so that it doesn't exceed the existing snow
            ground.STATVAR.surf_T = double ((ground.STATVAR.snowfall - ground.STATVAR.melt + sum(ground.STATVAR.ice_snow, 1)) <=0 | (ground.STATVAR.surf_T <0)) .* ground.STATVAR.surf_T; %set to zero if there is snow and T is positive
            
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
            
            %ground.TEMP.dLayerThick_massBalance = ground.STATVAR.snowfall .* rho_snow ./ 917 - melt .* ground.STATVAR.ice_snow(1,:) ./ ground.STATVAR.layerThick_snow(1,:);

            %flux is assigned in get_derivatives_prognostic

%             melt = repmat(melt, 1, PROFILE.ensemble_size)./1000; %in mm
%             snowfall = repmat(snowfall, 1, PROFILE.ensemble_size)./1000 .* PROFILE.snowfall_factor;
%  
%             melt = min(melt, snowfall + sum(PROFILE.D_ice_snow,1)); %limit the melt, so that it doesn't exceed the existing snow 
%             surf_T = double (snowfall-melt + sum(PROFILE.D_ice_snow,1)<=0 | surf_T <0) .* surf_T; %set to zero if there is snow and T is positive
            
            %energy flux due to snowfall and melt - CHECK if this needs to
            %be converted to snowfall per time! - should be fine!
            %ground.TEMP.d_energy(1:4,:) = ground.TEMP.d_energy(1:4,:) + ground.TEMP.snow_mat1 .*(repmat(snowfall .* (-ground.CONST.L_f + surf_T .* ground.CONST.c_i), 4, 1)  - repmat(melt,4,1) .* (-ground.CONST.L_f + ground.STATVAR.T(2:5,:) .* ground.CONST.c_i));
            
        end
                
        function ground = get_boundary_condition_l(ground,  tile)
            %ground.TEMP.d_energy(end,:) = ground.TEMP.d_energy(end,:) + repmat(tile.RUN_INFO.STATVAR.geothermal', 1, tile.ENSEMBLE.PARA.ensemble_size); %must have been rearranged properly
            ground.TEMP.d_energy(end,:) = ground.TEMP.d_energy(end,:) + repmat(tile.PARA.geothermal', 1, tile.ENSEMBLE.PARA.ensemble_size); 
        end
        
        %calculate spatial derivatives
        function ground = get_derivatives_prognostic(ground, tile)
            
%             ground.TEMP.d_energy = ground.TEMP.d_energy + (ground.STATVAR.thermCond_eff(2:end,:).*(ground.STATVAR.T(3:end,:)-ground.STATVAR.T(2:end-1,:))./ground.STATVAR.layerDistance(2:end,:) -...
%                 ground.STATVAR.thermCond_eff(1:end-1,:).*(ground.STATVAR.T(2:end-1,:)-ground.STATVAR.T(1:end-2,:))./ground.STATVAR.layerDistance(1:end-1,:));
            
            %downwards flux
            ground.TEMP.d_energy = ground.TEMP.d_energy - ground.STATVAR.thermCond_eff.*(ground.STATVAR.T(2:end,:)-ground.STATVAR.T(1:end-1,:))./ground.STATVAR.layerDistance;
            %upwards flux, lower boundary already added
            ground.TEMP.d_energy(1:end-1,:) = ground.TEMP.d_energy(1:end-1,:) + ground.STATVAR.thermCond_eff(2:end,:).*(ground.STATVAR.T(3:end,:)-ground.STATVAR.T(2:end-1,:))./ground.STATVAR.layerDistance(2:end,:);
            
            ground.TEMP.d_energy(1:3,:) = ground.TEMP.d_energy(1:3,:) .*  ground.TEMP.snow_mat2(1:3,:);

            ground.TEMP.dLayerThick_compaction = compact_windDrift(ground, tile); 
            
        end
        
        %prognostic step - integrate prognostic variables in time
        function ground = advance_prognostic(ground, tile)
             
            ground.STATVAR.ice_snow = ground.STATVAR.ice_snow + ground.TEMP.snow_mat1 .*repmat(ground.STATVAR.snowfall - ground.STATVAR.melt, 4,1);
            ground.STATVAR.ice_snow(ground.STATVAR.ice_snow<0) = 0;
            ground.STATVAR.layerThick_snow =  ground.STATVAR.layerThick_snow + (ground.TEMP.d_new_snow_layerThick - ground.TEMP.d_new_melt_layerThick);
            ground.STATVAR.layerThick_snow(ground.STATVAR.layerThick_snow<0) = 0;
            ground.STATVAR.layerThick_snow(ground.STATVAR.ice_snow == 0) = 0;
            
%             if sum(sum(ground.STATVAR.layerThick_snow == Inf))>0 || sum(sum(ground.STATVAR.layerThick_snow == -Inf))>0 || sum(sum(ground.STATVAR.layerThick_snow <0))>0
%                 'Hallo1'
%                 dhefhk
%             end
            
            ground.STATVAR.layerThick_snow = ground.STATVAR.layerThick_snow + ground.TEMP.dLayerThick_compaction.*tile.timestep;
            ground.STATVAR.layerThick_snow(ground.STATVAR.layerThick_snow<0) = 0;

            ground.STATVAR.energy = ground.STATVAR.energy + ground.TEMP.d_energy .* tile.timestep;
            
            for i=1:3 %set energy to 0 for unused cells
                ground.STATVAR.energy(i,:) = double(i >= ground.STATVAR.upper_cell) .* ground.STATVAR.energy(i,:);
            end
            

        end
        
        %diagnostic step - compute diagnostic variables
        function ground = compute_diagnostic(ground, tile)
            
            %diagnostic step, the 4 is the first ground cell
            %[d_D_ice, d_D, PROFILE.upper_cell] = regrid_snow(PROFILE);
            ground = regrid_snow(ground, tile);
            
%             if sum(sum(ground.STATVAR.layerThick_snow == Inf))>0 || sum(sum(ground.STATVAR.layerThick_snow == -Inf))>0 || sum(sum(ground.STATVAR.layerThick_snow <0))>0
%                 'Hallo3'
%                 dhefhk
%             end
            
            %PROFILE.D_ice_snow = PROFILE.D_ice_snow + d_D_ice;
            %PROFILE.D_snow = PROFILE.D_snow + d_D;

            ground.STATVAR.layerThick(2:4,:) =  ground.STATVAR.layerThick_snow(1:3,:) + double(ground.STATVAR.layerThick_snow(1:3,:) == 0) .* ground.PARA.virtual_gridCellSize; %set to snow cell size if snow, virtual_gridCellSize, otherwise
            ground.STATVAR.layerThick(5,:) = ground.STATVAR.layerThick_first_ground_cell + ground.STATVAR.layerThick_snow(4,:); %combined snow and ground cell
            ground.STATVAR.layerDistance(1:5,:) = ground.STATVAR.layerThick(2:6,:)./2 + ground.STATVAR.layerThick(1:5,:)./2; %recompute first five cells which are affected by snow; rest is constant

            ground = get_T(ground);
            
            ground = conductivity(ground);
            
            ground.TEMP.FT_count = ground.TEMP.FT_count + double(ground.STATVAR.T(5:35,:) <0); %only check until 10m for talik
            ground.TEMP.count = ground.TEMP.count + 1;
            
            ground.TEMP.d_energy = ground.STATVAR.energy .* 0;
        end
        
        %triggers
        function ground = check_trigger(ground, tile)
            
            if tile.t >= ground.TEMP.adjust_stratigraphy_date

                FT_code = ground.TEMP.FT_count ./ ground.TEMP.count;
                FT_code(FT_code==1)=2;  %frozen
                FT_code(FT_code>0 & FT_code < 1)=1;  %freeze thaw
                
                FT_res = zeros(1, size(FT_code,2));
                 for i = 1: size(FT_code,1)
                     FT_res = FT_res + double (FT_res == 0 & FT_code(i,:) ==-1) .* -1 + double (FT_res == 0 & FT_code(i,:) == 1) .* 1; %initial PF yes no, -1 or 1
                     FT_res = FT_res + double (FT_res >0 & FT_code(i,:) < 0) .* -FT_code(i,:); %switches from no PF, 1, to freeze_thaw or frozen, so 2 or 3 means talik
                 end
                 
                ground.STATVAR.FT_state = FT_res;
                
                %---------------update stratigrahy----
                
                gain_loose = FT_code.*0;
                gain_loose(FT_code==2) = 1;  %gain when frozen
                %OUT.STATVAR.gain_loose(OUT.ACC.FT_isothermal==1) = 0; %no gain when isothermal
                gain_loose(FT_code == 1 | FT_code == 0) = -1 ;  %loose when unfrozen or FT, this depends on water table settings for the ensemble member
                for i=1:size(gain_loose,1)-1  %set cell above AL to gain
                    gain_loose(i,:) = gain_loose(i,:) + double(FT_code(i,:)==1 & FT_code(i+1,:)==2); %& ~(OUT.ACC.FT_isothermal(i+1,:)==1)) .*(1 - OUT.gain_loose(i,:));
                end
                
                %update_stratigraphy
                ground = update_stratigraphy(ground, tile, gain_loose);
                
                ground.TEMP.FT_count = 0;   
                ground.TEMP.count = 0;
                ground.TEMP.adjust_stratigraphy_date = datenum([ground.PARA.adjust_stratigraphy_date num2str(str2num(datestr(tile.t, 'yyyy'))+1) ], 'dd.mm.yyyy');
            end
            
%             ground.TEMP.T_store = [ground.TEMP.T_store ground.STATVAR.T(:,1)];
%             ground.TEMP.E_store = [ground.TEMP.E_store ground.STATVAR.energy(:,1)];
%             ground.TEMP.lt_store =[ground.TEMP.lt_store ground.STATVAR.layerThick(:,1)];
            
        end
        
        
        
        
        %----non-mandatory functions
        
        
        function ground = get_T(ground) 
            
            %1. ground without first cell
		    ground.STATVAR.T(6:end,:) = double(ground.STATVAR.energy(5:end,:)>=0) .* ground.STATVAR.energy(5:end,:) ./ ground.STATVAR.c_thawed(2:end,:) + ...
                double(ground.STATVAR.energy(5:end,:) <= ground.STATVAR.E_frozen(2:end,:)) .* ((ground.STATVAR.energy(5:end,:) - ground.STATVAR.E_frozen(2:end,:)) ./ ground.STATVAR.c_frozen(2:end,:) + ground.STATVAR.T_end_freezing(2:end,:)) + ...
                double(ground.STATVAR.energy(5:end,:) < 0 & ground.STATVAR.energy(5:end,:) > ground.STATVAR.E_frozen(2:end,:)) .* ground.STATVAR.energy(5:end,:)./ground.STATVAR.E_frozen(2:end,:) .*(ground.STATVAR.T_end_freezing(2:end,:));

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
            ground.STATVAR.thermCond(5:end,:) = double(ground.STATVAR.T(5:end,:) < ground.STATVAR.T_end_freezing) .* ground.STATVAR.k_frozen + double(ground.STATVAR.T(5:end,:) > 0) .* ground.STATVAR.k_thawed + ...
                double(ground.STATVAR.T(5:end,:) >= ground.STATVAR.T_end_freezing & ground.STATVAR.T(5:end,:) <= 0) .* ground.STATVAR.k_freezing; % ground
            
            ground.STATVAR.thermCond_eff(5:end,:) = ground.STATVAR.thermCond(5:end-1,:).*ground.STATVAR.thermCond(6:end,:).*...
                (ground.STATVAR.layerThick(5:end-1,:)./2 + ground.STATVAR.layerThick(6:end,:)./2) ./ (ground.STATVAR.thermCond(5:end-1,:).*ground.STATVAR.layerThick(6:end,:)./2 + ...
                ground.STATVAR.thermCond(6:end,:).*ground.STATVAR.layerThick(5:end-1,:)./2 ); %size N
            
            %Thermal conductivity snow
            snow_density = ground.STATVAR.ice_snow ./ max(1e-20, ground.STATVAR.layerThick_snow) .*920;
            ground.STATVAR.thermCond_snow = max(5e-3, 2.3.*(snow_density./1000).^1.88);
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
        
        function ground = calculate_E_frozen (ground) %new function, only call this once in the beginning
            
            ground.STATVAR.T_end_freezing(1,:) = 0; %first cell freezes like free water
            ground.STATVAR.E_frozen = - ground.CONST.L_f .* ground.STATVAR.waterIce + ground.STATVAR.T_end_freezing.* ...
                (ground.CONST.c_w.*ground.STATVAR.waterIce./2 + ground.CONST.c_i.*ground.STATVAR.waterIce./2 + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic);
            %ground.STATVAR.E_frozen = ground.STATVAR.E_frozen .* ground.STATVAR.layerThick(5:end,:);  %absolute energy
            
            ground.STATVAR.c_thawed = (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic); % .* ground.STATVAR.layerThick(5:end,:); %unit J/m2/K
            ground.STATVAR.c_frozen = (ground.CONST.c_i.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic); % .* ground.STATVAR.layerThick(5:end,:); %unit J/m2/K

        end
        
 %--------------

        function ground = update_stratigraphy(ground, tile, gain_loose_in)
        
            layerThick = [ground.STATVAR.layerThick_first_ground_cell; ground.STATVAR.layerThick(6:end,:)];
            depths = cumsum(layerThick); %same size as theta_w
            gain_loose = depths.*0;
            gain_loose(1:31,:) = gain_loose_in;
            porosity = 1 - ground.STATVAR.mineral ./layerThick  - ground.STATVAR.organic ./layerThick;
            field_capacity = 0.5.* porosity;   %CHANGE later
            
            water_table_depth = zeros(1, size(field_capacity,2));
            for i=1:31
                water_table_depth = water_table_depth + double(gain_loose(i,:)==-1 & gain_loose(i+1,:)~=-1) .* (depths(i,:) - water_table_depth);
            end
            water_table_depth = water_table_depth .* (1-tile.ENSEMBLE.STATVAR.water_table_depth);
            for i=1:31
                gain_loose(i,:) = gain_loose(i,:) + double(depths(i,:) > water_table_depth) .* (1-gain_loose(i,:));
            end
            T_old = ground.STATVAR.T;
            water_old = ground.STATVAR.waterIce; %in m
            ground.STATVAR.waterIce = max(field_capacity, min(porosity, ground.STATVAR.waterIce ./layerThick + gain_loose .* 0.05)) .* layerThick;
            water_change = ground.STATVAR.waterIce - water_old; %iin m
            
            ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) > 0) .* ground.STATVAR.T(5:end,:) .* ground.CONST.c_w .* water_change;
            ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) < ground.STATVAR.T_end_freezing) .* (ground.STATVAR.T(5:end,:) .* ground.CONST.c_i - ground.CONST.L_f) .* water_change;
            E_frozen_change = - ground.CONST.L_f .* water_change + ground.STATVAR.T_end_freezing.*(ground.CONST.c_w .* water_change./2 + ground.CONST.c_i .* water_change./2);
            
            ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) <=0 & ground.STATVAR.T(5:end,:)>= ground.STATVAR.T_end_freezing) .* ground.STATVAR.T(5:end,:)./min(ground.STATVAR.T_end_freezing, -1e-12) .* E_frozen_change;
            ground.STATVAR.E_frozen = - ground.CONST.L_f.*ground.STATVAR.waterIce + ground.STATVAR.T_end_freezing.*(ground.CONST.c_w.*ground.STATVAR.waterIce./2 + ground.CONST.c_i.*ground.STATVAR.waterIce./2 + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic);
            
            ground.STATVAR.c_thawed = (ground.CONST.c_w.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic) ; %unit J/m2/K
            ground.STATVAR.c_frozen = (ground.CONST.c_i.*ground.STATVAR.waterIce + ground.CONST.c_m.*ground.STATVAR.mineral + ground.CONST.c_o.*ground.STATVAR.organic); %unit J/m2/K
            
            %PROFILE = E2T_compiled(PROFILE);
            ground = get_T(ground);
            %T_new = ground.STATVAR.T;
            store_layerThick = ground.STATVAR.layerThick;
            ground.STATVAR.layerThick = layerThick;
            ground = init_conductivity(ground);
            ground.STATVAR.layerThick = store_layerThick;

            ground = conductivity(ground);
            %test= PROFILE.water_table_depth;
            %save('save_gain_loose.mat', 'water_table_depth', 'gain_loose', 'water_change', 'T_old', 'T_new')
            
        end
        

        
 
        
        function rho_snow = get_snow_density(ground, tile)
            a = 109;
            b = 6;
            c = 26;
            T=min(0,ground.STATVAR.surf_T);
            %rho_snow = max(50, a + b.*T + (c.*PROFILE.wind_speed20.^0.5)); %using 20m height wind so that snow density is not influenced by forest canopy
            rho_snow = max(50, a + b.*T + (c .* tile.ENSEMBLE.STATVAR.wind_speed_class.^0.5));
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
            
%             rho(:,2)
            

            stress = ground.CONST.g .* rho_ice .* (cumsum(ground.STATVAR.ice_snow - ground.TEMP.snow_mat1 .* ground.STATVAR.ice_snow./2)); 
            
%             stress(:,2)
            
            
            eta = eta_0 .*  rho ./ c_eta .* exp(-a_eta .* T + b_eta .* rho);
            
%             eta(:,2)
            
            dD_dt =  - stress ./max(1e-10,eta) .* ground.STATVAR.layerThick_snow; %compaction
            
%             dD_dt(:,2)
            

            %dD_dt = dD_dt - PROFILE.snow_mat1 .* rho_ice .* PROFILE.D_ice_snow ./ max(1e-8, rho).^2.* (rho_max - min(rho , rho_max)) ./ tau_0 .* S_I;
            
            dD_dt = dD_dt - ground.TEMP.snow_mat1 .* rho_ice .* ground.STATVAR.ice_snow ./ max(1e-8, rho).^2.* (rho_max - min(rho, rho_max)) ./tile.ENSEMBLE.STATVAR.wind_compaction_timescale;
    
%             dD_dt(:,2)

            
	%dD_dt = dD_dt - repmat(double(sum(PROFILE.D_snow,1)>0.05),4,1) .* PROFILE.snow_mat1 .* rho_ice .* PROFILE.D_ice_snow ./ max(1e-8, rho).^2.* (rho_max - min(rho , rho_max)) ./PROFILE.wind_compaction_timescale;

        end
        
        function ground = regrid_snow(ground, tile)
            

            
%             ground.TEMP.before = sum(ground.STATVAR.ice_snow,1);
            
            D_ice = ground.STATVAR.ice_snow;
            D = ground.STATVAR.layerThick_snow;
            
            ice_content = D_ice./max(1e-10, D);
            
            %target_SWE = 0.03;
            target_SWE =0.05;
            
            D_tot=sum(D_ice,1);
            
            number_of_cells = min(3, floor(D_tot./ target_SWE));

            lower_cell = ground.TEMP.index_first_ground_cell - min(number_of_cells,1);
            ground.STATVAR.upper_cell = lower_cell - max(0, number_of_cells-1);
            for i=1:4
               ground.TEMP.snow_mat1(i,:) = double(ground.TEMP.snow_base_mat(i,:) == ground.STATVAR.upper_cell);
               %PROFILE.snow_mat1(i,:) = double(PROFILE.snow_base_mat(i,:) >= PROFILE.upper_cell & PROFILE.snow_base_mat(i,:) <= min(PROFILE.upper_cell+1, lower_cell));
               ground.TEMP.snow_mat2(i,:) = double(ground.TEMP.snow_base_mat(i,:) >= ground.STATVAR.upper_cell & ground.TEMP.snow_base_mat(i,:) <= lower_cell);
            end
            
            factor = max(1, number_of_cells);
            
%             test=zeros(4, size(PROFILE.first_ground_cell,2));
%             for i=1:4
%                test(i,:) = double(i>=PROFILE.upper_cell & i<= lower_cell) .* D_tot./factor;
%             end
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
            
%             ground.TEMP.after = sum(ground.STATVAR.ice_snow,1);
            
%             if sum(sum(ground.STATVAR.layerThick_snow<0))>0
%                 'Hllo10'
%                 djkdfff
%             end
            
            %PROFILE.snow_mat1 .*(repmat(snowfall .* (-PROFILE.const.L + surf_T .* PROFILE.const.c_i), 4, 1)  - repmat(melt,4,1) .* (-PROFILE.const.L + PROFILE.T(2:5,:) .* PROFILE.const.c_i));

        end
       
        %-----------discontinued-----------
        
             
        function PROFILE = advance_one_timestep(PROFILE, surf_T, melt, snowfall, wind, difference_altitude)

            speedup_var = size(PROFILE.T,2) == 1400;             

            surf_T = repmat(surf_T, 1, PROFILE.ensemble_size);
            melt = repmat(melt, 1, PROFILE.ensemble_size)./1000; %in mm
            snowfall = repmat(snowfall, 1, PROFILE.ensemble_size)./1000 .* PROFILE.snowfall_factor;
 
            melt = min(melt, snowfall + sum(PROFILE.D_ice_snow,1)); %limit the melt, so that it doesn't exceed the existing snow 
            surf_T = double (snowfall-melt + sum(PROFILE.D_ice_snow,1)<=0 | surf_T <0) .* surf_T; %set to zero if there is snow and T is positive
            
            PROFILE.wind_speed = repmat(wind, 1, PROFILE.ensemble_size); % .* log(2./ PROFILE.z0) ./ log(repmat(difference_altitude, 1, PROFILE.ensemble_size) ./ PROFILE.z0) ; %downscaling of wind speeds assuming log-profile
            PROFILE.wind_speed20 = repmat(wind, 1, PROFILE.ensemble_size); % .* log(20./ PROFILE.z0) ./ log(repmat(difference_altitude, 1, PROFILE.ensemble_size) ./ PROFILE.z0) ;
           

            
            PROFILE = E2T(PROFILE); %this is actually diagnostic step form last timestep
           


            T_old=PROFILE.T(1:5,:);

            for i=1:4 %assign boundary condition T to correct cell and all cells above
                PROFILE.T(i, :) =  PROFILE.T(i, :) + double(i<=PROFILE.upper_cell) .* (surf_T - PROFILE.T(i, :));
            end
            
            
            
            %prognostic
            PROFILE.k(4:end,:) = double(PROFILE.T(5:end,:) < PROFILE.T_end_freezing) .* PROFILE.k_frozen + double(PROFILE.T(5:end,:) > 0) .* PROFILE.k_thawed + ...
                double(PROFILE.T(5:end,:) >= PROFILE.T_end_freezing & PROFILE.T(5:end,:) <= 0) .* PROFILE.k_freezing; % ground
            
	    %get k_eff for the ground, no interference with snow, valid from position 4 down, i.e. between cells 4 and 5
            %k_eff_ground = PROFILE.k(1:end-1,:).*PROFILE.k(2:end,:).*(PROFILE.K_delta(2:end-1,:)./2 + PROFILE.K_delta(3:end,:)./2) ./ (PROFILE.k(1:end-1,:).*PROFILE.K_delta(3:end,:)./2 + PROFILE.k(2:end,:).*PROFILE.K_delta(2:end-1,:)./2 ); %size N
           
if speedup_var 
 k_eff_ground = k_eff_compiled_mex(PROFILE.k, PROFILE.K_delta);
else
 k_eff_ground = PROFILE.k(4:end-1,:).*PROFILE.k(5:end,:).*(PROFILE.K_delta(5:end-1,:)./2 + PROFILE.K_delta(6:end,:)./2) ./ (PROFILE.k(4:end-1,:).*PROFILE.K_delta(6:end,:)./2 + PROFILE.k(5:end,:).*PROFILE.K_delta(5:end-1,:)./2 ); %size N
end	   



	    snow_density = PROFILE.D_ice_snow ./ max(1e-20, PROFILE.D_snow) .*920;
            PROFILE.k_snow = max(5e-3, 2.2.*(snow_density./1000).^1.88);
            PROFILE.k(1:3,:) = PROFILE.k_snow(1:3,:); %snow

            for i=1:3 %replace conductivities above upper boundary by some high value, 2 W/mK
                PROFILE.k(i,:) = PROFILE.k(i,:) + double(i < PROFILE.upper_cell) .* (2 - PROFILE.k(i,:));
            end

              
            PROFILE.k(4,:) = PROFILE.k(4,:).*PROFILE.k_snow(4,:).*(PROFILE.K_delta_first_ground_cell./2 + PROFILE.D_snow(4,:)) ./ ...
                 (PROFILE.k(4,:) .* PROFILE.D_snow(4,:) + PROFILE.k_snow(4,:) .* PROFILE.K_delta_first_ground_cell./2) ; %first half cell plus snow, modified after Mamoru, error corrected
            

         
	    k_eff_snow = PROFILE.k(1:3,:).*PROFILE.k(2:4,:).*(PROFILE.K_delta(2:4,:)./2 + PROFILE.K_delta(3:5,:)./2) ./ (PROFILE.k(1:3,:).*PROFILE.K_delta(3:5,:)./2 + PROFILE.k(2:4,:).*PROFILE.K_delta(2:4,:)./2 ); %size N
	    k_eff = [k_eff_snow; k_eff_ground];
	    
	    k_eff = [zeros(1, size(k_eff,2))+2; k_eff];



if speedup_var 
	PROFILE.dE_dt(1:end-1,:) = dE_dt_compiled_mex(k_eff, PROFILE.T, PROFILE.cT_delta);
else
            PROFILE.dE_dt(1:end-1,:) = (k_eff(2:end,:).*(PROFILE.T(3:end,:)-PROFILE.T(2:end-1,:))./PROFILE.cT_delta(2:end,:) -...
                k_eff(1:end-1,:).*(PROFILE.T(2:end-1,:)-PROFILE.T(1:end-2,:))./PROFILE.cT_delta(1:end-1,:));  %size N-1
end


            

            PROFILE.dE_dt(1:3,:) =  PROFILE.dE_dt(1:3,:) .* PROFILE.snow_mat2(1:3,:);
            
            % lower BC (dT_dt=geothermal heat flux)
            PROFILE.dE_dt(end,:) = PROFILE.F_lb - k_eff(end-1,:).*(PROFILE.T(end,:)-PROFILE.T(end-1,:))./PROFILE.cT_delta(end,:);
            

            dD_dt_compaction = compact_windDrift(PROFILE); 
            
            %add new snow
            new_density = get_snow_density(PROFILE, surf_T);
	    %new_density = get_snow_density_compiled_mex(PROFILE.wind_speed_class, surf_T);


            
            %new_snow_D = snowfall .* 920 ./ new_density ;
            %new_melt_D = melt .* diag(PROFILE.D_snow(PROFILE.upper_cell,:))'./max(1e-10, diag(PROFILE.D_ice_snow(PROFILE.upper_cell,:))');
            new_snow_D = repmat(snowfall .* 920 ./ new_density, 4, 1) .* PROFILE.snow_mat1;
            new_melt_D = repmat(melt,4,1) .* PROFILE.snow_mat1 .* PROFILE.D_snow ./ max(1e-10, PROFILE.D_ice_snow);
            
            PROFILE.D_ice_snow = PROFILE.D_ice_snow + PROFILE.snow_mat1 .*repmat(snowfall - melt, 4,1);
            PROFILE.D_snow =  PROFILE.D_snow + PROFILE.snow_mat1 .* (new_snow_D - new_melt_D);
            PROFILE.E(1:4,:) = PROFILE.E(1:4,:) + PROFILE.snow_mat1 .*(repmat(snowfall .* (-PROFILE.const.L + surf_T .* PROFILE.const.c_i), 4, 1)  - repmat(melt,4,1) .* (-PROFILE.const.L + PROFILE.T(2:5,:) .* PROFILE.const.c_i));
            

            
            PROFILE.D_snow = PROFILE.D_snow + PROFILE.timestep .* dD_dt_compaction;

if speedup_var 
            PROFILE.E = advance_E_compiled_mex(PROFILE.E , PROFILE.timestep, PROFILE.dE_dt); 
else
           PROFILE.E = PROFILE.E + PROFILE.timestep .* PROFILE.dE_dt;
end


            %book-keeping, this hopefully eliminates small rounding errors
            PROFILE.D_ice_snow(PROFILE.D_ice_snow<0) = 0;
            PROFILE.D_snow(PROFILE.D_snow<0)=0;
            PROFILE.D_snow(PROFILE.D_ice_snow == 0) = 0;
            
            for i=1:3 %set energy to 0 for unused cells
                PROFILE.E(i,:) = double(i >= PROFILE.upper_cell) .* PROFILE.E(i,:);
            end
            
            %diagnostic step, the 4 is the first ground cell
            %[d_D_ice, d_D, PROFILE.upper_cell] = regrid_snow(PROFILE);
            PROFILE = regrid_snow(PROFILE);
            
            %PROFILE.D_ice_snow = PROFILE.D_ice_snow + d_D_ice;
            %PROFILE.D_snow = PROFILE.D_snow + d_D;

            PROFILE.K_delta(2:4,:) =  PROFILE.D_snow(1:3,:) + double(PROFILE.D_snow(1:3,:) == 0) .* PROFILE.virtual_gridCellSize; %set to snow cell size if snow, virtual_gridCellSize, otherwise
            PROFILE.K_delta(5,:) = PROFILE.K_delta_first_ground_cell + PROFILE.D_snow(4,:); %combined snow and ground cell
            PROFILE.cT_delta(1:5,:) = PROFILE.K_delta(2:6,:)./2 + PROFILE.K_delta(1:5,:)./2; %recompute first five cells which are affected by snow; rest is constant
            
            PROFILE.timestamp = PROFILE.timestamp + PROFILE.timestep./PROFILE.const.day_sec;
            
        end
        
  
        
        
        function PROFILE = E2T(PROFILE) %new function
            
            %change theta_w and all other to be one entry less
            %introduce D_ground
            %1. ground without first cell
		
	speedup_var = size(PROFILE.T,2) == 1400;    

if speedup_var
           PROFILE.T(6:end,:) = get_T_compiled_mex(PROFILE.E, PROFILE.c_thawed, PROFILE.c_frozen, PROFILE.E_frozen, PROFILE.T_end_freezing);
else
            PROFILE.T(6:end,:) = double(PROFILE.E(5:end,:)>=0) .* PROFILE.E(5:end,:) ./ PROFILE.c_thawed(2:end,:) + ...
                double(PROFILE.E(5:end,:) <= PROFILE.E_frozen(2:end,:)) .* ((PROFILE.E(5:end,:) - PROFILE.E_frozen(2:end,:)) ./ PROFILE.c_frozen(2:end,:) + PROFILE.T_end_freezing(2:end,:)) + ...
                double(PROFILE.E(5:end,:) < 0 & PROFILE.E(5:end,:) > PROFILE.E_frozen(2:end,:)) .* PROFILE.E(5:end,:)./PROFILE.E_frozen(2:end,:) .*(PROFILE.T_end_freezing(2:end,:));
end

  
            %2. first ground cell including initial snow - free water
            %freeze curve here, makes things much easier with the snow

            E_frozen_first_ground_cell = PROFILE.E_frozen(1,:) - PROFILE.const.L .* PROFILE.D_ice_snow(4,:);
            %c_frozen_first_cell = PROFILE.c_frozen(1,:) +  PROFILE.const.c_i .* PROFILE.D_ice_snow(4,:);
            
            PROFILE.T(5,:) = double (double(PROFILE.E(4,:)>=0)) .* PROFILE.E(4,:) ./ (PROFILE.c_thawed(1,:) + PROFILE.const.c_w .* PROFILE.D_ice_snow(4,:)) + ...
                double(PROFILE.E(4,:) <= E_frozen_first_ground_cell) .* ( (PROFILE.E(4,:) - E_frozen_first_ground_cell) ./ (PROFILE.c_frozen(1,:) +  PROFILE.const.c_i .* PROFILE.D_ice_snow(4,:)) ); 
            %zero degrees else

            
            %3. snow
            T_snow = double(PROFILE.E(1:3,:) < -PROFILE.const.L .* PROFILE.D_ice_snow(1:3,:)) .* (PROFILE.E(1:3,:) + PROFILE.const.L .* PROFILE.D_ice_snow(1:3,:)) ./ ...
                (PROFILE.const.c_i .*PROFILE.D_ice_snow(1:3,:));
            T_snow(isnan(T_snow))=0;
            T_snow(abs(T_snow)==Inf)=0; 
            PROFILE.T(2:4,:) = T_snow;
            
        end
        

    
              
        function PROFILE = initialize_T_1(PROFILE, forcing, year_list) %generates first estimate of steady - state T profile
            PROFILE.T = PROFILE.theta_w .*0; %only ground cell
            %PROFILE.T(1,:) =  surf_T_av_estimate_initial(forcing, year_list, PROFILE.snowfall_factor, PROFILE.ensemble_size);

            %changed Sebastian, uses TTOP Jaros to initialize
            rk = mean(PROFILE.k_thawed(1:10,:),1)./ mean(PROFILE.k_frozen(1:10,:),1); %uppermost meter as proxy for AL/seasonal frost layer
            %overwrite - use the classes from Jaros
            rk(floor(PROFILE .stratigraphy) == 1 | floor(PROFILE .stratigraphy) == 7 | floor(PROFILE .stratigraphy) == 8) = 0.95;
            rk(floor(PROFILE .stratigraphy) == 2) = 0.75;
            rk(floor(PROFILE .stratigraphy) == 3) = 0.8;
            rk(floor(PROFILE .stratigraphy) == 4 | floor(PROFILE .stratigraphy) == 9) = 0.95;
            rk(floor(PROFILE .stratigraphy) == 5) = 0.9;
            rk(floor(PROFILE .stratigraphy) == 6 | floor(PROFILE .stratigraphy) == 10) = 0.55;
            
            PROFILE.T(1,:) =  surf_T_av_estimate_initial2(forcing, year_list, PROFILE.snowfall_factor, PROFILE.ensemble_size, rk');
            %----

            for i=1:size(PROFILE.T,1)-1
                k = double(PROFILE.T(i,:) < PROFILE.T_end_freezing(i,:)) .* PROFILE.k_frozen(i,:) + double(PROFILE.T(i,:) > PROFILE.T_onset_freezing) .* PROFILE.k_thawed(i,:) + ...
                    double(PROFILE.T(i,:) >= PROFILE.T_end_freezing(i,:) & PROFILE.T(i,:) <= PROFILE.T_onset_freezing) .* PROFILE.k_freezing(i,:);
                PROFILE.T(i+1,:) = PROFILE.T(i,:) + PROFILE.F_lb .* PROFILE.cT_delta(i,:) ./k;
            end
            PROFILE = T2E(PROFILE);
            PROFILE.dE_dt=PROFILE.E.*0;
            PROFILE.T=[zeros(4, size(PROFILE.T,2)); PROFILE.T];  %add three snow cells and one virtual cell with T=0
        end
        
        
        function PROFILE = initialize_T_2(PROFILE, out)
            PROFILE.T = PROFILE.theta_w .*0;
            PROFILE.T(1,:) =  surf_T_av_estimate_TTOP(PROFILE, out);
            for i=1:size(PROFILE.T,1)-1
                k = double(PROFILE.T(i,:) < PROFILE.T_end_freezing(i,:)) .* PROFILE.k_frozen(i,:) + double(PROFILE.T(i,:) > PROFILE.T_onset_freezing) .* PROFILE.k_thawed(i,:) + ...
                    double(PROFILE.T(i,:) >= PROFILE.T_end_freezing(i,:) & PROFILE.T(i,:) <= PROFILE.T_onset_freezing) .* PROFILE.k_freezing(i,:);
                PROFILE.T(i+1,:) = PROFILE.T(i,:) + PROFILE.F_lb .* PROFILE.cT_delta(i,:) ./k;
            end
            PROFILE = T2E(PROFILE);
            PROFILE.dE_dt=PROFILE.E.*0;
            PROFILE.T=[zeros(4, size(PROFILE.T,2)); PROFILE.T];  %add three snow cells and one virtual cell with T=0
            
        end
        
        
        function surf_T_av = surf_T_av_estimate_TTOP(PROFILE, out)
            rk = mean(PROFILE.k_thawed(1:10,:),1)./ mean(PROFILE.k_frozen(1:10,:),1); %uppermost meter as proxy for AL/seasonal frost layer
            PF_yes_no = (rk .* out.TDD + out.FDD)<=0;
            surf_T_av = double(PF_yes_no) .* (rk .* out.TDD + out.FDD)./out.N + double(~PF_yes_no) .* (out.TDD + 1 ./ rk .* out.FDD)./out.N;
        end

        
%         function PROFILE = set_timestamp(PROFILE, timestamp)
%             PROFILE.timestamp = timestamp;
%         end
        
        
        
    end
    
end

