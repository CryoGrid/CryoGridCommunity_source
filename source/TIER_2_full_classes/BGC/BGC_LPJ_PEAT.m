%========================================================================
% CryoGrid GROUND class BGC_CLASS_TEMPLATE

% S. Westermann, November 2020
%========================================================================

classdef BGC_LPJ_PEAT < PEAT_ACCUMULATION & PEAT_DECOMPOSE
    
    properties
        PARENT
        IA_BGC
    end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = provide_PARA(ground)%, tile)
            %initialize the PARA here
            
            %sebastian can you pass tile here
            %                       ground.PARA.start_year = datevec(tile.FORCING.PARA.start_time);
            %             ground.PARA.end_year = datevec(tile.FORCING.PARA.end_time);
            
            % % Initializing variable
            
%             ground.PARA.noy = 20; %end_year(1,1)-start_year(1,1) + 1; % no. of years
%             ground.PARA.nop = 1; % no. of patches
%             ground.PARA.noi = 1; % no. of iteration
%             ground.PARA.tyears = ground.PARA.noy*ground.PARA.noi; % total no. of years to run
%             ground.PARA.totaldays = ground.PARA.noy*ground.PARA.noi*365; % total no. of years to run
%             ground.PARA.total_water=  2500;%initialize total water (mm)
%             ground.PARA.snow = 0; % snow initialization
%             ground.PARA.boundry_condition = 200; % boundary condition in case to run only for limit years
%             ground.PARA.atmfrac = 0.0; % amount of NPP lost to the atm
%             ground.PARA.rec= 1;
%             ground.PARA.rept = 1; % paramater to adjust parallel runs (np. of time model can run)
%             ground.PARA.npp_frac_shrub = 0.5; % adjusting NPP based on shrub productivity
%             ground.PARA.npp_frac_graminoid  = 2.0; % adjusting NPP based on graminoid productivity
%             ground.PARA.heightOfWater = 200; % height of water above ground (cm)
%             ground.PARA.snow_days = 0;% initialising snow days
%             ground.PARA.precipitaion_days = 0;% initialising ppt. days
%             ground.PARA.adjusting_SH = (0.3)*1.0; % andjusting uneven ground
            
%             ground.PARA.max_evaporation = 2.0; %maximum evaporationoration
%             ground.PARA.c6 = 0.5;
%             ground.PARA.c7 = (200)^-1;
%             
%             ground.PARA.ids = 20e-2/365;
%             ground.PARA.idg = 30e-2/365;
%             ground.PARA.idm = 10e-2/365;
            
%             ground.PARA.water_in = 0;
%             ground.PARA.snow_in = 0;
%             ground.PARA.T_in = 15;
            
%             ground.PARA.t = 1;
%             ground.PARA.t_old = 0;
%             ground.PARA.dayofyear = 1;
%             ground.PARA.endtime = ground.PARA.totaldays+1;
            
%             ground.PARA.lastyear = 0;
%             ground.PARA.flag = 1;
%             ground.PARA.flag2 = 1;
%             ground.PARA.flag3 = 1;
%             ground.PARA.flag4 = 1;
            
%             ground.PARA.initial_topo = 200;
%             ground.PARA.adjusting_SH = 0.3;
            
%             ground.TEMP.timestep_PEAT = 0;
%             ground.PARA.target_timestep_PEAT = 3600;
%             ground.TEMP.C_derivative = 0;
%             ground.TEMP.dailytimestep = 1;
%             
%             
%             ground.PARA.minbulkDensity = 40;
%             ground.PARA.diffbulkDensity = 80;
%             
%             ground.PARA.porosity = 0.9; % porosity in the soil
%             ground.PARA.mineral_bulkDensity =160;% bulk density (in kgC/m3)
%             
%             
%             ground.PARA.initialDecomposition = 10e-2/365;
%             ground.PARA.decompo = 1.0;
%             
%             ground.PARA.accumulation_month = 1;
%             ground.PARA.accumulation_day = 1;

            ground.PARA.Q10 = 2;

            ground.PARA.initialDecomposition = 0.055/365; %10e-2/365;
            ground.PARA.decompo = 1.0;

            ground.PARA.fieldCapacity = 0.75;% soil water at field capacity
            ground.PARA.bulkDensity = 105; %is this kg/m3
            
            ground.PARA.decompose_timestep = 1.25; %[days]
            %ground.PARA.accumulate_timestep = 365; %[days]
            
            ground.PARA.number_of_growing_days = 60;
            ground.PARA.TDD_onset_growing_season = 100; %NPP allocation starts of the first cell has exceeded this number of thawing degree days
            
            ground.PARA.annual_NPP = 0.3; %kg/m2;
        end
        
        
        function ground = provide_STATVAR(ground)
            %initialize the STATVAR here
%             ground.STATVAR.total_whc = 0;%zeros(ground.PARA.tyears,ground.PARA.nop);%intialize total water holding capacity for every patch
%             ground.STATVAR.mx= 0;%zeros(1,ground.PARA.nop);%intialize total water holding capacity for every patch
%             ground.STATVAR.peatD_all= 0;%intialize peatdepth of moss
%             ground.STATVAR.twater= 0;%zeros(ground.PARA.tyears,ground.PARA.nop); %intialize total water in peat
%             ground.STATVAR.tpeat= 0;%zeros(ground.PARA.tyears,ground.PARA.nop);%zeros(1,ground.PARA.nop); %intialize total peat
%             
%             ground.STATVAR.peat_moss = 0;%zeros(ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.peat_shrub = 0;%zeros(ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.peat_graminoid = 0;%zeros(ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.new_peat = 0;
%             ground.STATVAR.total_peat = 0;%zeros(ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.max_waterHoldingCapacity= 0;%zeros(ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.peat_depth = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
            
            %             ground.STATVAR.water_influx = zeros(1,ground.PARA.totaldays);
            %             ground.STATVAR.snow_influx = zeros(1,ground.PARA.totaldays);
            %             ground.STATVAR.daily_melt = zeros(1,ground.PARA.totaldays);
            %             ground.STATVAR.twas = zeros(1,ground.PARA.totaldays);
            %             ground.STATVAR.snowa = zeros(1,ground.PARA.totaldays);
            
%             ground.STATVAR.mxx= 0;%zeros(1,ground.PARA.nop);
%             ground.STATVAR.wtde = 0;%zeros(ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.difcuu = 0;%zeros(1,ground.PARA.tyears);
%             ground.STATVAR.total_cumulative_peat = 0;%zeros(ground.PARA.tyears,ground.PARA.nop);
            %ground.STATVAR.runoff = zeros(1,ground.PARA.nop);
            %ground.STATVAR.evaporation= zeros(ground.PARA.totaldays,ground.PARA.nop);
            %ground.STATVAR.temp_modifier  = zeros (ground.PARA.totaldays,ground.PARA.nop);
            %ground.STATVAR.water_modifier = zeros(ground.PARA.totaldays,ground.PARA.nop);
            %ground.STATVAR.water_content = ones(ground.PARA.totaldays,ground.PARA.nop);
            %ground.STATVAR.dwater_depth = zeros(365,ground.PARA.nop);
%             %ground.STATVAR.water_depth = zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.catos= 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.cats = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.catom= 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.catm = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.catog= 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.catg = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
            %ground.STATVAR.decdepth= 0;%zeros(ground.PARA.tyears,ground.PARA.nop);
            %ground.STATVAR.soiltd = 0;%zeros (ground.PARA.totaldays,ground.PARA.nop);
%             %ground.STATVAR.mw = 0;%zeros (ground.PARA.totaldays,ground.PARA.nop);
%             ground.STATVAR.peataccu = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.awater_depth = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.fpcgg = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.fpcss  = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.fpcmm  = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
            
%             ground.STATVAR.fpm = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.fpg = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.fps  = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.T  = ones (ground.PARA.tyears,ground.PARA.nop).*5;
            %ground.STATVAR.water  = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
            
%             ground.STATVAR.totalpeatC = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             %ground.STATVAR.totalpeatC_originalMass = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.massRemain = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.tempModifier = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             
%             
%             ground.STATVAR.waterContent = ones(ground.PARA.tyears,ground.PARA.nop)*0.5;
%             ground.STATVAR.waterModifier = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.cato= 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             %             ground.STATVAR.bulkDensity= 0;%ones (ground.PARA.tyears,ground.PARA.nop)*105;
%             
             ground.STATVAR.layerThick = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);%[0.01; 0.2; 0.2];
%             ground.STATVAR.upperPosition_gridCell = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);%[0; 0.01; 0.21]; %upper position of each grid cell, not used now, but will get relevant for Xice
%             ground.STATVAR.area = 0;%ones (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.XwaterIce = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.organic = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.mineral = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.energy = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
%             ground.STATVAR.waterIce = 0;%zeros (ground.PARA.tyears,ground.PARA.nop);
            
            
            %initialize with some random values, is overwritten in first timestep
            %             ia_BGC.BGC.STATVAR.vol_water = ground.STATVAR.layerThick .*0 + 0.5;
            %             ia_BGC.BGC.STATVAR.vol_mineral = ground.STATVAR.layerThick .*0 + 0.3;
            %             ia_BGC.BGC.STATVAR.porosity = ground.STATVAR.layerThick .*0 + 0.5;
            %             ia_BGC.BGC.STATVAR.field_capacity = ground.STATVAR.layerThick .*0 + 0.3;
            %             ia_BGC.BGC.STATVAR.T = ground.STATVAR.layerThick .*0 + 5;
        end
        
        function ground = provide_CONST(ground)
            %initialize the CONST here
%             ground.CONST.mtomm = 1000;
%             ground.CONST.mtocm = 100;
%             ground.CONST.cmtomm = 10;
%             ground.CONST.mmtocm = 0.1;
%             ground.CONST.mmtom = 0.001;
%             ground.CONST.cmtom = 0.01;
            ground.CONST.day_sec = 24.*3600;
            
            % constant based on literature values
%             ground.CONST.SOC_fraction_new_peat = 0.5; %weight fraction of C of the total material -> convert kgC to kg new material
%             ground.CONST.buk_density_new_peat = 0.08 .* 1000; % bulk density = dry organic material mass / dry organic material volume, unit [kg/m3]
%             ground.CONST.porosity_new_peat = 0.9;
            
            ground.CONST.organicDensity = 1300; %in kg/m3, initial porosity of peat is 0.91 with bulk density of 105 kg/m3
            
        end
        
        function ground = finalize_init(ground, tile)
            
            ground = provide_PARA(ground);
            ground = provide_STATVAR(ground);
            ground = provide_CONST(ground);
            
%             ground = readingNppData(ground);  %reading NPP data
%             ground = distributing_NPP(ground);
            
            ground.STATVAR.total_peat = [];
            ground.STATVAR.totalpeatC_originalMass = [];
            

            ground.STATVAR.next_decompose_timestamp = ceil(tile.FORCING.PARA.start_time + ground.PARA.decompose_timestep);
            ground.STATVAR.next_accumulate_timestamp = ground.STATVAR.next_decompose_timestamp;
            
            ground.TEMP.year_count = 0;
            ground.TEMP.year = -100; %str2num(datestr(tile.FORCING.PARA.startt,'yyyy'));
            ground.TEMP.TDD = 0;
            ground.TEMP.growing_days_count = 0;
            ground.TEMP.accumulation_new_cell_flag = 1;
            ground.TEMP.regrid_flag = 0;
            
            ground.TEMP.d_layerThick = 0;
            ground.TEMP.d_organic = 0;

        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)

        end
        
        
        function ground = get_boundary_condition_l(ground, tile)

        end
        
        function ground = get_derivatives_prognostic(ground, tile)

            if tile.t == ground.STATVAR.next_decompose_timestamp
                get_ground_variables(ground.IA_BGC, tile); %get the physical variables over to the BGC class, so that it starts with correct state from the last timestep
                
                ground = temp_modifier(ground);
                ground = water_modifier(ground);
                ground = peat_decompose(ground);
                

                
                
                % extracting date, month and year from the tile
                Date_vector = datevec(tile.t);
                accumulation_year = Date_vector(:, 1);
                accumulation_month = Date_vector(:, 2);
                accumulation_day = Date_vector(:,3);
                
                %reset FDD on 1 Jan
                if accumulation_year ~= ground.TEMP.year
                    ground.TEMP.year = accumulation_year;
                    ground.TEMP.TDD = 0;
                    ground.TEMP.accumulation_new_cell_flag = 1; 
                    ground.TEMP.growing_days_count = 0;
                    ground.TEMP.regrid_flag = 1;
                end
                
                if ~isempty(ground.STATVAR.T)
                    ground.TEMP.TDD = ground.TEMP.TDD + double(ground.STATVAR.T(1,1)>0) .* ground.STATVAR.T(1,1).* ground.PARA.decompose_timestep;
                end
                % peat accumulates starts
                if ground.TEMP.TDD > ground.PARA.TDD_onset_growing_season && ground.TEMP.growing_days_count <= ground.PARA.number_of_growing_days && ground.STATVAR.T(1,1) > 0
                    ground.TEMP.growing_days_count = ground.TEMP.growing_days_count + ground.PARA.decompose_timestep;
                    
                
                    
                    if ground.TEMP.accumulation_new_cell_flag == 1
                        %ground = get_peatC(ground);
                        ground = peat_accumulation_newCell(ground);
                        ground.TEMP.year_count = ground.TEMP.year_count+1;
                        ground.TEMP.accumulation_new_cell_flag = 0;
                        add_grid_cell(ground.IA_BGC, tile);
                    else
                        ground = peat_accumulation(ground);
                    end
                    
                end
                
                %               ground = updatebulkD(ground);
                

            end
            
        end
        
        function timestep = get_timestep(ground, tile)
            %modify
%             if tile.t == ground.STATVAR.next_decompose_timestamp
%                 ground.STATVAR.next_decompose_timestamp = ground.STATVAR.next_decompose_timestamp + ground.PARA.decompose_timestep;
%             end
            
            timestep = (ground.STATVAR.next_decompose_timestamp - tile.t) .* ground.CONST.day_sec;
            
            
            %@NITIN: update the long peat timestep
            %ground.TEMP.timestep_PEAT = ground.TEMP.timestep_PEAT + timestep;
        end
        
        function ground = advance_prognostic(ground, tile)
            %add here
            
            
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            %add here
        end
        
        function ground = compute_diagnostic(ground, tile)
          
            if tile.t == ground.STATVAR.next_decompose_timestamp
                send_BGC_variables(ground.IA_BGC, tile); %send the BGC variables over to the GROUND class
                
                ground.STATVAR.next_decompose_timestamp = ground.STATVAR.next_decompose_timestamp + ground.PARA.decompose_timestep;
                ground.TEMP.d_layerThick = 0;
                ground.TEMP.d_organic = 0;
            end
            
        end
        
        function ground = check_trigger(ground, tile)
            %regrid the pysical stratigraphy
            if ground.TEMP.regrid_flag == 1
                regrid_stratigraphy(ground.IA_BGC, tile);
                ground.TEMP.regrid_flag = 0;
            end
        end
        
        
        function ground = reset_time(ground, tile) %used e.g. with TILE_BUILDER update_forcing_out
            ground.STATVAR.next_decompose_timestamp = tile.t + ground.PARA.decompose_timestep;
            ground.TEMP.year = str2num(datestr(tile.t, 'yyyy'));
        end

    end
    
end
