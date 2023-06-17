
classdef GROUND_MULTITILE_ESA_CCI_SEB < SEB
    

    methods
        
        %-----initialize-----------------
        
        function ground = provide_PARA(ground)
            ground.PARA.albedo_ground = [];
            ground.PARA.albedo_snow_max = [];
            ground.PARA.albedo_snow_min = [];
            ground.PARA.emissivity_snow = [];
            ground.PARA.emissivity_ground = [];
            ground.PARA.turbulent_exchange_coefficient = [];
            ground.PARA.resistance_ET = [];
                        
            ground.PARA.virtual_gridCellSize = [];
            ground.PARA.timestep = [];
            ground.PARA.adjust_stratigraphy_date = [];
            
            ground.PARA.wind_speed_class = [];
            ground.PARA.wind_compaction_timescale = [];
            
            ground.PARA.field_capacity_snow = [];
%             ground.PARA.water_table_depth = [];
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
            
            ground.CONST.sigma = [];
            ground.CONST.cp = [];
            ground.CONST.Tmfw = [];
            ground.CONST.g = [];

            ground.CONST.L_f = [];
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
            ground.STATVAR.water_snow = ground.STATVAR.layerThick_snow;
            ground.STATVAR.waterIce_snow = ground.STATVAR.layerThick_snow;
            
            ground.STATVAR.upper_cell = repmat(4, 1, size(ground.STATVAR.layerThick,2));
            
            ground.STATVAR.thermCond = ground.STATVAR.layerThick.*0; % grid cell property
            ground.STATVAR.thermCond_eff = ground.STATVAR.layerThick(1:end-1,:).*0; %between grid cells property
            ground = conductivity(ground);
            
            ground.TEMP.FT_count = 0;
            ground.TEMP.count = 0;
            
            
            ground.TEMP.adjust_stratigraphy_date = datenum([ground.PARA.adjust_stratigraphy_date num2str(str2num(datestr(tile.FORCING.PARA.start_time, 'yyyy'))+1) ], 'dd.mm.yyyy');
  
            ground.STATVAR.FT_state = tile.PARA.geothermal .* NaN;
            
            ground.TEMP.d_energy = ground.STATVAR.energy .* 0;
            ground.TEMP.d_ice_snow = ground.TEMP.d_energy(1:4,:);
            
            %initialize surface variables
            ground.STATVAR.albedo = ground.STATVAR.T(1,:) .*0 + ground.PARA.albedo_ground;
            ground.STATVAR.emissivity = ground.STATVAR.T(1,:) .*0 + ground.PARA.emissivity_ground;            
            ground.STATVAR.surf_T = ground.STATVAR.T(5,:);
            ground.STATVAR.resistance_ET = ground.PARA.resistance_ET; %chnage to 0 when there is snow
            ground.TEMP.snow_yes_no = ground.STATVAR.T(1,:) .*0;
            
       end
        
        
        %-----mandatory functions------------------------
        function ground = get_boundary_condition_u(ground, tile)
            
            %has to be distributed correctly as soon as ensemble is
            %available
            
            snowfall = tile.FORCING.TEMP.snowfall ./1000  ./ ground.CONST.day_sec;
%             rainfall = tile.FORCING.TEMP.rainfall ./1000  ./ ground.CONST.day_sec;            
            
            for i=1:4 %assign boundary condition T to correct cell and all cells above, needed to get the conductive heat fluxes right
                ground.STATVAR.T(i, :) =  ground.STATVAR.T(i, :) + double(i<ground.STATVAR.upper_cell) .* (ground.STATVAR.surf_T - ground.STATVAR.T(i, :));
            end

           SEB = tile.FORCING.TEMP.Sin .* (1 - ground.STATVAR.albedo); %these variables must be initialized and then written very time in the diagnostic step
           SEB = SEB + ground.STATVAR.emissivity .* tile.FORCING.TEMP.Lin;
           SEB = SEB - ground.STATVAR.emissivity .* ground.CONST.sigma .* (ground.STATVAR.surf_T + ground.CONST.Tmfw).^4; %Lout
           SEB = SEB - ground.PARA.turbulent_exchange_coefficient .* ground.CONST.cp .* (ground.STATVAR.surf_T - tile.FORCING.TEMP.Tair); %sensible heat flux
            
           ET = 1./(1./ground.PARA.turbulent_exchange_coefficient + ground.STATVAR.resistance_ET) .* 2500e3 .* ...
               (tile.FORCING.TEMP.q - 0.622.* 6.112 .* 100 .* exp(17.62.*ground.STATVAR.surf_T./(243.12 + ground.STATVAR.surf_T))./tile.FORCING.TEMP.p);
%            ground.PARA.turbulent_exchange_coefficient = rho.*kappa.* uz.*kappa./log(z./z0).^2;
% ET = 1 .* ground.PARA.turbulent_exchange_coefficient .* 2500e3 .* ...
%     (tile.FORCING.TEMP.q - 0.622.* 6.112 .* 100 .* exp(17.62.*ground.STATVAR.surf_T./(243.12 + ground.STATVAR.surf_T))./tile.FORCING.TEMP.p);
           SEB = SEB + ET;
           
           ground.TEMP.d_energy(1:4,:) = ground.TEMP.d_energy(1:4,:) + ground.TEMP.snow_mat1 .* repmat(SEB, 4, 1);
            
%            ground.TEMP.d_water(1:4,:) = ground.TEMP.d_water(1:4,:) + repmat(rainfall, 4,1);
           ground.TEMP.d_ice_snow(1:4,:) = ground.TEMP.d_ice_snow(1:4,:) + ground.TEMP.snow_mat1 .* repmat(snowfall, 4, 1);
            
           ground.TEMP.d_energy(1:4,:) = ground.TEMP.d_energy(1:4,:) + ground.TEMP.snow_mat1 .*(repmat(snowfall .* ...
                (-ground.CONST.L_f + ground.STATVAR.surf_T .* ground.CONST.c_i), 4, 1));
            new_snow_density = get_snow_density(ground, tile);

            
            ground.TEMP.d_layerThick_snow = repmat(snowfall .* 920 ./ new_snow_density , 4, 1) .* ground.TEMP.snow_mat1;

            %flux is assigned in get_derivatives_prognostic

        end
                
        function ground = get_boundary_condition_l(ground,  tile)
            ground.TEMP.d_energy(end,:) = ground.TEMP.d_energy(end,:) + tile.PARA.geothermal; 
        end
        
        %calculate spatial derivatives
        function ground = get_derivatives_prognostic(ground, tile)
            
            %downwards flux
            ground.TEMP.d_energy = ground.TEMP.d_energy - ground.STATVAR.thermCond_eff.*(ground.STATVAR.T(2:end,:)-ground.STATVAR.T(1:end-1,:))./ground.STATVAR.layerDistance;
            %upwards flux, lower boundary already added
            ground.TEMP.d_energy(1:end-1,:) = ground.TEMP.d_energy(1:end-1,:) + ground.STATVAR.thermCond_eff(2:end,:).*(ground.STATVAR.T(3:end,:)-ground.STATVAR.T(2:end-1,:))./ground.STATVAR.layerDistance(2:end,:);
            
            ground.TEMP.d_energy(1:3,:) = ground.TEMP.d_energy(1:3,:) .*  ground.TEMP.snow_mat2(1:3,:);

            ground.TEMP.dLayerThick_compaction = compact_windDrift(ground, tile); 
            
        end
        
        %prognostic step - integrate prognostic variables in time
        function ground = advance_prognostic(ground, tile)
             
            ground.STATVAR.waterIce_snow = ground.STATVAR.waterIce_snow + ground.TEMP.d_ice_snow .* tile.timestep;
            ground.STATVAR.ice_snow = ground.STATVAR.ice_snow + ground.TEMP.d_ice_snow .* tile.timestep;

            ground.STATVAR.layerThick_snow = ground.STATVAR.layerThick_snow + ground.TEMP.d_layerThick_snow .* tile.timestep;
            ground.STATVAR.layerThick_snow = ground.STATVAR.layerThick_snow + ground.TEMP.dLayerThick_compaction.*tile.timestep;
            ground.STATVAR.layerThick_snow(ground.STATVAR.layerThick_snow<0) = 0;

            ground.TEMP.target_density = ground.STATVAR.ice_snow ./ max(ground.STATVAR.layerThick_snow, 1e-12); %ice faction
            
            ground.STATVAR.energy = ground.STATVAR.energy + ground.TEMP.d_energy .* tile.timestep;
            
            for i=1:3 %set energy to 0 for unused cells
                ground.STATVAR.energy(i,:) = double(i >= ground.STATVAR.upper_cell) .* ground.STATVAR.energy(i,:);
            end
            

        end
        
        %diagnostic step - compute diagnostic variables
        function ground = compute_diagnostic(ground, tile)
            
            %remove melted snow grid cells
            for i=1:2
                melted = double(ground.STATVAR.energy(i,:)>0);
                if sum(melted>0)
                    ground.STATVAR.energy(i+1,:) = ground.STATVAR.energy(i+1,:) + melted .* ground.STATVAR.energy(i,:);
                    ground.STATVAR.waterIce_snow(i+1,:) = ground.STATVAR.waterIce_snow(i+1,:) + melted .* ground.STATVAR.waterIce_snow(i,:);
                    
                    ground.STATVAR.upper_cell = ground.STATVAR.upper_cell + melted;
                    ground.STATVAR.waterIce_snow(i,:) = ground.STATVAR.waterIce_snow(i,:) .* double(~melted);
                    ground.STATVAR.layerThick_snow(i,:) = ground.STATVAR.layerThick_snow(i,:) .* double(~melted);
                    ground.STATVAR.energy(i,:) = ground.STATVAR.energy(i,:) .* double(~melted);
                end
            end
            
            %last grid cell melted, but snow above
            melted = ground.STATVAR.energy(3,:)>0 & ground.STATVAR.energy(2,:)<0;
            if sum(melted)>0
                ground.STATVAR.energy(3,:) = ground.STATVAR.energy(3,:) + double(melted) .* ground.STATVAR.energy(2,:);
                ground.STATVAR.waterIce_snow(3,:) = ground.STATVAR.waterIce_snow(3,:) + double(melted) .* ground.STATVAR.waterIce_snow(2,:);
                ground.STATVAR.layerThick_snow(3,:) = ground.STATVAR.layerThick_snow(3,:) + double(melted) .* (ground.STATVAR.layerThick_snow(2,:) - ground.STATVAR.layerThick_snow(3,:));
                ground.TEMP.target_density(3,:) = ground.TEMP.target_density(3,:) + double(melted) .* (ground.TEMP.target_density(2,:) - ground.TEMP.target_density(3,:));
                
                ground.STATVAR.energy(2,:) = ground.STATVAR.energy(2,:) + double(melted) .* (ground.STATVAR.energy(1,:) - ground.STATVAR.energy(2,:));
                ground.STATVAR.waterIce_snow(2,:) = ground.STATVAR.waterIce_snow(2,:) + double(melted) .* (ground.STATVAR.waterIce_snow(1,:) - ground.STATVAR.waterIce_snow(2,:));
                ground.STATVAR.layerThick_snow(2,:) = ground.STATVAR.layerThick_snow(2,:) .* double(melted) .* (ground.STATVAR.layerThick_snow(1,:) - ground.STATVAR.layerThick_snow(2,:));
                ground.TEMP.target_density(2,:) = ground.TEMP.target_density(2,:) + double(melted) .* (ground.TEMP.target_density(1,:) - ground.TEMP.target_density(2,:));
                
                ground.STATVAR.energy(1,:) = ground.STATVAR.energy(1,:) .* double(~melted);
                ground.STATVAR.waterIce_snow(1,:) = ground.STATVAR.waterIce_snow(1,:) .* double(~melted);
                ground.STATVAR.layerThick_snow(1,:) = ground.STATVAR.layerThick_snow(1,:) .* double(~melted);
                ground.TEMP.target_density(1,:) = ground.TEMP.target_density(1,:) .* double(~melted);                
                ground.STATVAR.upper_cell = ground.STATVAR.upper_cell + melted;
            end

            %last grid cell melted, remove snow completely
            melted = ground.STATVAR.energy(3,:)>0;
            if sum(melted)>0
                ground.STATVAR.energy(3,:) = ground.STATVAR.energy(3,:) .* double(~melted);
                ground.STATVAR.waterIce_snow(3,:) = ground.STATVAR.waterIce_snow(3,:) .* double(~melted);
                ground.STATVAR.layerThick_snow(3,:) = ground.STATVAR.layerThick_snow(3,:) .* double(~melted);
                ground.STATVAR.upper_cell = ground.STATVAR.upper_cell + melted;
            end
            
            %combined ground/snow grid cell melted 
            melted = ground.STATVAR.energy(4,:)>0 & ground.STATVAR.waterIce_snow(4,:)>0;
            if sum(melted)>0
                ground.STATVAR.waterIce_snow(4,:) = ground.STATVAR.waterIce_snow(4,:) .* double(~melted);
                ground.STATVAR.layerThick_snow(4,:) = ground.STATVAR.layerThick_snow(4,:) .* double(~melted);
            end
            
            %take out all water
            
            test1= sum(ground.STATVAR.waterIce_snow(1:3,:),1);

            %regrid
            %waterIce_snow = ground.STATVAR.waterIce_snow;
            %water balance
            E_frozen = -ground.CONST.L_f.*ground.STATVAR.waterIce_snow(1:3,:);
            ice_snow = double(ground.STATVAR.energy(1:3,:) <= E_frozen) .*ground.STATVAR.waterIce_snow(1:3,:) + double(ground.STATVAR.energy(1:3,:) > E_frozen) .* ground.STATVAR.energy(1:3,:) ./ (-ground.CONST.L_f);
            water_snow = double(ground.STATVAR.energy(1:3,:) > E_frozen) .* (ground.STATVAR.energy(1:3,:) + ground.CONST.L_f .* ground.STATVAR.waterIce_snow(1:3,:)) ./ ground.CONST.L_f;

            ground.STATVAR.layerThick_snow = ice_snow ./ max(1e-12, ground.TEMP.target_density(1:3,:));
            
            max_water = ground.PARA.field_capacity_snow .* (ground.STATVAR.layerThick_snow(1:3,:) - ice_snow);
            excess_water = max(0, water_snow - max_water);
            water_snow = water_snow - excess_water;
            space_left = max_water - water_snow;
            excess_water = sum(excess_water, 1);
            for i=1:2
                water_added = min(space_left(i,:), excess_water);
                water_snow(i,:) = water_snow(i,:) + water_added;
                excess_water = excess_water - water_added;
            end
            
            space_left = ground.STATVAR.layerThick_snow(1:3,:) - ice_snow;
            for i=3:-1:1
                water_added = min(space_left(i,:), excess_water);
                water_snow(i,:) = water_snow(i,:) + water_added;
                excess_water = excess_water - water_added;
            end
            
            ground.STATVAR.waterIce_snow(1:3,:) = water_snow + ice_snow;
            ground.STATVAR.ice_snow(1:3,:) = ice_snow;
            ground.STATVAR.water_snow(1:3,:) = water_snow;         


            test3= sum(ground.STATVAR.waterIce_snow(1:3,:),1) + excess_water;

            dec= (abs(test1-test3))>0.0001;

            %snow on first ground cell
            energySnow_upperGroundCell = double(ground.STATVAR.waterIce_snow(4,:)>0) .* (ground.STATVAR.energy(4,:) - ground.STATVAR.E_frozen(1,:));
            E_frozen = -ground.CONST.L_f.*ground.STATVAR.waterIce_snow(4,:);
            ice_snow = double(energySnow_upperGroundCell <= E_frozen) .*ground.STATVAR.waterIce_snow(4,:) + double(energySnow_upperGroundCell > E_frozen) .* energySnow_upperGroundCell ./ (-ground.CONST.L_f);
            ice_snow = max(0,ice_snow);
            water_snow = ground.STATVAR.waterIce_snow(4,:) - ice_snow;
            ground.STATVAR.layerThick_snow(4,:) = ice_snow ./ max(1e-12, ground.TEMP.target_density(4,:));
            space_left = ground.STATVAR.layerThick_snow(4,:) - ice_snow;
            water_snow = min(water_snow, space_left);
            ground.STATVAR.waterIce_snow(4,:) = ice_snow + water_snow;
            ground.STATVAR.ice_snow(4,:) = ice_snow;
            ground.STATVAR.water_snow(4,:) = water_snow;     


            
            %diagnostic step, the 4 is the first ground cell
            ground = regrid_snow(ground, tile);

            ground.STATVAR.layerThick(2:4,:) =  ground.STATVAR.layerThick_snow(1:3,:) + double(ground.STATVAR.layerThick_snow(1:3,:) == 0) .* ground.PARA.virtual_gridCellSize; %set to snow cell size if snow, virtual_gridCellSize, otherwise
            ground.STATVAR.layerThick(5,:) = ground.STATVAR.layerThick_first_ground_cell + ground.STATVAR.layerThick_snow(4,:); %combined snow and ground cell
            ground.STATVAR.layerDistance(1:5,:) = ground.STATVAR.layerThick(2:6,:)./2 + ground.STATVAR.layerThick(1:5,:)./2; %recompute first five cells which are affected by snow; rest is constant

            ground = get_T(ground);
            
            ground = conductivity(ground);
            
%             ground.TEMP.FT_count = ground.TEMP.FT_count + double(ground.STATVAR.T(5:35,:) <0); %only check until 10m for talik
%             ground.TEMP.count = ground.TEMP.count + 1;
            
            ground.TEMP.snow_yes_no = sum(ground.STATVAR.waterIce_snow, 1) >1e-10;
            ground.STATVAR.albedo = double(~ground.TEMP.snow_yes_no) .* ground.PARA.albedo_ground + double(ground.TEMP.snow_yes_no) .* ground.PARA.albedo_snow_max;
            ground.STATVAR.emissivity = double(~ground.TEMP.snow_yes_no) .* ground.PARA.emissivity_ground + double(ground.TEMP.snow_yes_no) .* ground.PARA.emissivity_snow;            
            for i=1:4
                ground.STATVAR.surf_T = double(ground.STATVAR.upper_cell == i) .* ground.STATVAR.T(i+1,:);
            end
            ground.STATVAR.resistance_ET = double(~ground.TEMP.snow_yes_no) .*ground.PARA.resistance_ET; %chnage to 0 when there is snow

            ground.TEMP.d_energy = ground.STATVAR.energy .* 0;
            ground.TEMP.d_ice_snow = ground.TEMP.d_ice_snow .* 0;
        
            
            ground.TEMP.FT_count = ground.TEMP.FT_count + double(ground.STATVAR.T(5:35,:) <0); %only check until 10m for talik
            ground.TEMP.count = ground.TEMP.count + 1;
        
        
        end
        
        %triggers
        function ground = check_trigger(ground, tile)
            
            if tile.t >= ground.TEMP.adjust_stratigraphy_date

                FT_code = ground.TEMP.FT_count ./ ground.TEMP.count;
                FT_code(FT_code==1)=2;  %frozen
                FT_code(FT_code>0 & FT_code < 1)=1;  %freeze thaw
                
                gain_loose = FT_code.*0;
                gain_loose(FT_code==2) = 1;  %gain when frozen
                %OUT.STATVAR.gain_loose(OUT.ACC.FT_isothermal==1) = 0; %no gain when isothermal
                gain_loose(FT_code == 1 | FT_code == 0) = -1 ;  %loose when unfrozen or FT, this depends on water table settings for the ensemble member
                for i=1:size(gain_loose,1)-1  %set cell above AL to gain
                    gain_loose(i,:) = gain_loose(i,:) + double(FT_code(i,:)==1 & FT_code(i+1,:)==2) .* (1 - gain_loose(i,:)); %& ~(OUT.ACC.FT_isothermal(i+1,:)==1)) ;
                   % gain_loose(i,:) = gain_loose(i,:) + double(FT_code(i,:)==1 & FT_code(i+1,:)==2) .* (1 - gain_loose(i,:)); %& ~(OUT.ACC.FT_isothermal(i+1,:)==1)) ;
                end
                
                
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
            
            %1. ground without first cell
		    ground.STATVAR.T(6:end,:) = double(ground.STATVAR.energy(5:end,:)>=0) .* ground.STATVAR.energy(5:end,:) ./ ground.STATVAR.c_thawed(2:end,:) + ...
                double(ground.STATVAR.energy(5:end,:) <= ground.STATVAR.E_frozen(2:end,:)) .* ((ground.STATVAR.energy(5:end,:) - ground.STATVAR.E_frozen(2:end,:)) ./ ground.STATVAR.c_frozen(2:end,:) + ground.STATVAR.T_end_freezing(2:end,:)) + ...
                double(ground.STATVAR.energy(5:end,:) < 0 & ground.STATVAR.energy(5:end,:) > ground.STATVAR.E_frozen(2:end,:)) .* ground.STATVAR.energy(5:end,:)./ground.STATVAR.E_frozen(2:end,:) .*(ground.STATVAR.T_end_freezing(2:end,:));

            %2. first ground cell including initial snow - free water
            %freeze curve here, makes things much easier with the snow

            E_frozen_first_ground_cell = ground.STATVAR.E_frozen(1,:) - ground.CONST.L_f .* ground.STATVAR.waterIce_snow(4,:);
            
            ground.STATVAR.T(5,:) = double (double(ground.STATVAR.energy(4,:)>=0)) .* ground.STATVAR.energy(4,:) ./ (ground.STATVAR.c_thawed(1,:) + ground.CONST.c_w .* ground.STATVAR.waterIce_snow(4,:)) + ...
                double(ground.STATVAR.energy(4,:) <= E_frozen_first_ground_cell) .* ( (ground.STATVAR.energy(4,:) - E_frozen_first_ground_cell) ./ (ground.STATVAR.c_frozen(1,:) +  ground.CONST.c_i .* ground.STATVAR.waterIce_snow(4,:)) ); 
            %zero degrees else

            %3. snow
            T_snow = double(ground.STATVAR.energy(1:3,:) < -ground.CONST.L_f .* ground.STATVAR.waterIce_snow(1:3,:)) .* (ground.STATVAR.energy(1:3,:) + ground.CONST.L_f .* ground.STATVAR.waterIce_snow(1:3,:)) ./ ...
                (ground.CONST.c_i .*ground.STATVAR.waterIce_snow(1:3,:));
            T_snow(isnan(T_snow))=0;
            T_snow(abs(T_snow)==Inf)=0; 
            ground.STATVAR.T(2:4,:) = T_snow;
            
            E_frozen = -ground.CONST.L_f.*ground.STATVAR.waterIce_snow(1:3,:);
            ground.STATVAR.ice_snow(1:3,:) = double(ground.STATVAR.energy(1:3,:) <= E_frozen) .*ground.STATVAR.waterIce_snow(1:3,:) + double(ground.STATVAR.energy(1:3,:) > E_frozen) .* ground.STATVAR.energy(1:3,:) ./ (-ground.CONST.L_f);
            ground.STATVAR.water_snow(1:3,:) = double(ground.STATVAR.energy(1:3,:) > E_frozen) .* (ground.STATVAR.energy(1:3,:) + ground.CONST.L_f .* ground.STATVAR.waterIce_snow(1:3,:)) ./ ground.CONST.L_f;

            energySnow_upperGroundCell = double(ground.STATVAR.waterIce_snow(4,:)>0) .* (ground.STATVAR.energy(4,:) - ground.STATVAR.E_frozen(1,:));
            E_frozen = -ground.CONST.L_f.*ground.STATVAR.waterIce_snow(4,:);
            ice_snow = double(energySnow_upperGroundCell <= E_frozen) .*ground.STATVAR.waterIce_snow(4,:) + double(energySnow_upperGroundCell > E_frozen) .* energySnow_upperGroundCell ./ (-ground.CONST.L_f);
            ice_snow = max(0,ice_snow);
            water_snow = ground.STATVAR.waterIce_snow(4,:) - ice_snow;

            ground.STATVAR.ice_snow(4,:) = ice_snow;
            ground.STATVAR.water_snow(4,:) = water_snow;
        end
        
        
        function ground = conductivity(ground)
            ground.STATVAR.thermCond(5:end,:) = double(ground.STATVAR.T(5:end,:) < ground.STATVAR.T_end_freezing) .* ground.STATVAR.k_frozen + double(ground.STATVAR.T(5:end,:) > 0) .* ground.STATVAR.k_thawed + ...
                double(ground.STATVAR.T(5:end,:) >= ground.STATVAR.T_end_freezing & ground.STATVAR.T(5:end,:) <= 0) .* ground.STATVAR.k_freezing; % ground
            
            ground.STATVAR.thermCond_eff(5:end,:) = ground.STATVAR.thermCond(5:end-1,:).*ground.STATVAR.thermCond(6:end,:).*...
                (ground.STATVAR.layerThick(5:end-1,:)./2 + ground.STATVAR.layerThick(6:end,:)./2) ./ (ground.STATVAR.thermCond(5:end-1,:).*ground.STATVAR.layerThick(6:end,:)./2 + ...
                ground.STATVAR.thermCond(6:end,:).*ground.STATVAR.layerThick(5:end-1,:)./2 ); %size N
            
            %Thermal conductivity snow
            snow_density = ground.STATVAR.waterIce_snow ./ max(1e-20, ground.STATVAR.layerThick_snow) .*920;
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
            
           %NEW
            for i=1:3 %set higher thermal conductivity for influx in first snow cell 
                ground.STATVAR.thermCond_eff(i, :) =  ground.STATVAR.thermCond_eff(i, :) + double(i<=ground.STATVAR.upper_cell) .* (2 - ground.STATVAR.thermCond_eff(i, :));
            end
            
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
            ground.STATVAR.waterIce = max(field_capacity, min(porosity, ground.STATVAR.waterIce ./layerThick + gain_loose .* 0.05)) .* layerThick;
            water_change = ground.STATVAR.waterIce - water_old; %iin m
            
%             ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) > 0) .* ground.STATVAR.T(5:end,:) .* ground.CONST.c_w .* water_change;
%             ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) < ground.STATVAR.T_end_freezing) .* (ground.STATVAR.T(5:end,:) .* ground.CONST.c_i - ground.CONST.L_f) .* water_change;
%             E_frozen_change = - ground.CONST.L_f .* water_change + ground.STATVAR.T_end_freezing.*(ground.CONST.c_w .* water_change./2 + ground.CONST.c_i .* water_change./2);
%             
%             ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) <=0 & ground.STATVAR.T(5:end,:)>= ground.STATVAR.T_end_freezing) .* ground.STATVAR.T(5:end,:)./min(ground.STATVAR.T_end_freezing, -1e-12) .* E_frozen_change;
            ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) > 0) .* ground.STATVAR.T(5:end,:) .* ground.CONST.c_w .* water_change;
            ground.STATVAR.energy(4:end,:) = ground.STATVAR.energy(4:end,:) + double(ground.STATVAR.T(5:end,:) < ground.STATVAR.T_end_freezing) .* ...
                (ground.STATVAR.T_end_freezing .* (ground.CONST.c_w./2 + ground.CONST.c_i./2) + (ground.STATVAR.T(5:end,:)-ground.STATVAR.T_end_freezing) .* ground.CONST.c_i - ground.CONST.L_f) .* water_change;
            
            ice_fraction = ground.STATVAR.energy(4:end,:)./ground.STATVAR.E_frozen;
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
%             a = 109;
%             b = 6;
%             c = 26;
%             T=min(0,ground.STATVAR.surf_T);
%             %rho_snow = max(50, a + b.*T + (c.*PROFILE.wind_speed20.^0.5)); %using 20m height wind so that snow density is not influenced by forest canopy
            %             rho_snow = max(50, a + b.*T + (c .* tile.ENSEMBLE.STATVAR.wind_speed_class.^0.5));
            T_air = min(0,ground.STATVAR.surf_T);
            rho_Tair = double(T_air > -15).*(50 + 1.7.*(T_air+15).^(1.5)) ...
                + double(T_air <= -15).*( -3.8328.*T_air - 0.0333.*T_air.^2);
            rho_wind = 266.861.*(0.5.*(1 + tanh(ground.PARA.wind_speed_class./5))).^8.8;
            rho_snow = rho_Tair + rho_wind;
        end

        
      
        function dD_dt = compact_windDrift(ground, tile)

            T = min(0,ground.STATVAR.T(2:5,:));
            
            eta_0 = 7.62237e6;
            a_eta = 0.1;
            b_eta = 0.023;
            c_eta = 250;
                        
            rho_ice = 920;
            rho_max = 350;
            
            rho = ground.STATVAR.waterIce_snow ./ max(1e-20, ground.STATVAR.layerThick_snow) .* rho_ice;
            
            stress = ground.CONST.g .* rho_ice .* (cumsum(ground.STATVAR.waterIce_snow - ground.TEMP.snow_mat1 .* ground.STATVAR.waterIce_snow./2));             
            
            eta = eta_0 .*  rho ./ c_eta .* exp(-a_eta .* T + b_eta .* rho);
                        
            dD_dt =  - stress ./max(1e-10,eta) .* ground.STATVAR.layerThick_snow; %compaction            
            
            dD_dt = dD_dt - ground.TEMP.snow_mat1 .* rho_ice .* ground.STATVAR.waterIce_snow ./ max(1e-8, rho).^2.* (rho_max - min(rho, rho_max)) ./ repmat(ground.PARA.wind_compaction_timescale .*24.*3600,4,1);
    

        end
        
        function ground = regrid_snow(ground, tile)
            
            test1= ground.STATVAR.energy(1:4,:);

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
            
            target_ice = ground.TEMP.snow_mat2.* repmat(D_tot./factor,4,1);
            
            snow_over = max(0, target_ice(1,:) - target_SWE);
            target_ice(1,:)= target_ice(1,:)-snow_over;
            target_ice(2,:)= target_ice(2,:)+snow_over;
            
            d_D_ice = target_ice - D_ice;
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
            
            
            %energy and water
            ice_down = d_D_ice_res(1:3,:) .* 0;
            ice_up = ice_down;
            void_up = ice_up;
            void_down= ice_up;
            
            temp_d_D_ice_res = d_D_ice_res;
            temp_d_void_res = d_D_res - d_D_ice_res;
            
            for i=1:3
                ice_down(i,:) = double(temp_d_D_ice_res(i,:)<0) .* -temp_d_D_ice_res(i,:);
                temp_d_D_ice_res(i,:) = temp_d_D_ice_res(i,:) + ice_down(i,:);
                temp_d_D_ice_res(i+1,:) = temp_d_D_ice_res(i+1,:) - ice_down(i,:);
                
                ice_up(i,:) = double(temp_d_D_ice_res(i,:)>0) .* temp_d_D_ice_res(i,:);
                temp_d_D_ice_res(i,:) = temp_d_D_ice_res(i,:) - ice_up(i,:);
                temp_d_D_ice_res(i+1,:) = temp_d_D_ice_res(i+1,:) + ice_up(i,:);
                
                void_down(i,:) = double(temp_d_void_res(i,:)<0) .* -temp_d_void_res(i,:);
                temp_d_void_res(i,:) = temp_d_void_res(i,:) + void_down(i,:);
                temp_d_void_res(i+1,:) = temp_d_void_res(i+1,:) - void_down(i,:);
                
                void_up(i,:) = double(temp_d_void_res(i,:)>0) .* temp_d_void_res(i,:);
                temp_d_void_res(i,:) = temp_d_void_res(i,:) - void_up(i,:);
                temp_d_void_res(i+1,:) = temp_d_void_res(i+1,:) + void_up(i,:);
            end
            
            ice_up = ice_up .* double(ground.STATVAR.ice_snow(2:4,:)>0);
            ice_down = ice_down .* double(ground.STATVAR.ice_snow(1:3,:)>0);
            ice_down(3,:) = ice_down(3,:) .* double(ground.STATVAR.ice_snow(3,:)<target_SWE);

            void_up = void_up .* double(ground.STATVAR.ice_snow(2:4,:)>0);
            void_down = void_down .* double(ground.STATVAR.ice_snow(1:3,:)>0);
            void_down(3,:) = void_down(3,:) .* double(ground.STATVAR.ice_snow(3,:)<target_SWE);


            void_space = ground.STATVAR.layerThick_snow - ground.STATVAR.ice_snow;
            
            energy_prior = ground.STATVAR.energy(1:4,:);
            energy_prior(4,:) = double(ground.STATVAR.waterIce_snow(4,:)>0) .* double(ground.STATVAR.energy(4,:) >= ground.STATVAR.E_frozen(1,:) - ground.CONST.L_f .* ground.STATVAR.waterIce_snow(4,:)) .* (ground.STATVAR.energy(4,:) - ground.STATVAR.E_frozen(1,:)); %melting snow
            energy_prior(4,:) = energy_prior(4,:) + double(ground.STATVAR.waterIce_snow(4,:)>0) .* double(ground.STATVAR.energy(4,:) < ground.STATVAR.E_frozen(1,:) - ground.CONST.L_f .* ground.STATVAR.waterIce_snow(4,:)) .* ...
                (-ground.CONST.L_f .* ground.STATVAR.waterIce_snow(4,:) + (ground.CONST.c_i .* ground.STATVAR.waterIce_snow(4,:)) ./ (ground.STATVAR.c_frozen(1,:) + ground.CONST.c_i .* ground.STATVAR.waterIce_snow(4,:)) .* (ground.STATVAR.energy(4,:) - ground.STATVAR.E_frozen(1,:) + ground.CONST.L_f .* ground.STATVAR.waterIce_snow(4,:)));

            %regrid starts here!
            ground.STATVAR.waterIce_snow = ground.STATVAR.waterIce_snow + d_D_ice_res;
            ground.STATVAR.layerThick_snow = ground.STATVAR.layerThick_snow + d_D_res;
            
            ground.STATVAR.energy(1:3,:) = ground.STATVAR.energy(1:3,:) - ice_down ./ max(1e-12, ground.STATVAR.ice_snow(1:3,:)) .* energy_prior(1:3,:) + ...
                ice_up ./ max(1e-12, ground.STATVAR.ice_snow(2:4,:)) .* energy_prior(2:4,:);
            ground.STATVAR.energy(2:4,:) = ground.STATVAR.energy(2:4,:) + ice_down ./ max(1e-12, ground.STATVAR.ice_snow(1:3,:)) .* energy_prior(1:3,:) - ...
                ice_up ./ max(1e-12, ground.STATVAR.ice_snow(2:4,:)) .* energy_prior(2:4,:);
            
            water_prior = ground.STATVAR.water_snow;
            ground.STATVAR.waterIce_snow(1:3,:) = ground.STATVAR.waterIce_snow(1:3,:) - void_down ./ max(1e-12, void_space(1:3,:)) .* water_prior(1:3,:) + ...
                void_up ./ max(1e-12, void_space(2:4,:)) .* water_prior(2:4,:);
            ground.STATVAR.waterIce_snow(2:4,:) = ground.STATVAR.waterIce_snow(2:4,:) + void_down ./ max(1e-12, void_space(1:3,:)) .* water_prior(1:3,:) - ...
                void_up ./ max(1e-12, void_space(2:4,:)) .* water_prior(2:4,:);
            
            ground.STATVAR.layerThick_snow(ground.STATVAR.layerThick_snow<0) = 0;
            ground.STATVAR.waterIce_snow(ground.STATVAR.waterIce_snow<0) = 0;

            ground.STATVAR.waterIce_snow(4,:) = ground.STATVAR.waterIce_snow(4,:) .* double(ground.STATVAR.waterIce_snow(3,:) ==0);
            ground.STATVAR.layerThick_snow(4,:) = ground.STATVAR.layerThick_snow(4,:) .* double(ground.STATVAR.layerThick_snow(3,:) ==0);


            test = (abs(sum(test1,1) - sum(ground.STATVAR.energy(1:4,:),1)))>100;

        end
       

        
        
    end
    
end

