%========================================================================
% CryoGrid BGC class based on Frolking et al. 2010, https://esd.copernicus.org/articles/1/1/2010/esd-1-1-2010.pdf
% The model is combined with a dependence of NPP on solar radiation, ET and
% near-surface temperature
% S. Westermann, November 2021
%========================================================================

classdef BGC_Frolking_PEAT < PEAT_ACCUMULATION & PEAT_DECOMPOSE
    
    properties
        PARENT
        IA_BGC
    end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = provide_PARA(ground)%, tile)

            ground.PARA.Q10 = 2;
            ground.PARA.initialDecomposition = 0.055/365; %10e-2/365;

            ground.PARA.fieldCapacity = 0.75;% soil water at field capacity
            ground.PARA.bulkDensity = 105; %is this kg/m3
            
            ground.PARA.BGC_timestep = 1.25; %[days]
            
            %accumulation
            ground.PARA.start_PhotSyn_T = 0;
            ground.PARA.start_fullPhotSyn_T = 10;
            ground.PARA.end_fullPhotSyn_T = 25;
            ground.PARA.end_PhotSyn_T = 35;

        end
        
        
        function ground = provide_STATVAR(ground)
             
             ground.STATVAR.layerThick = 0;

        end
        
        function ground = provide_CONST(ground)

            ground.CONST.day_sec = 24.*3600;
            
            ground.CONST.organicDensity = 1300; %in kg/m3, initial porosity of peat is 0.91 with bulk density of 105 kg/m3
            
        end
        
        function ground = finalize_init(ground, tile)
            
            ground = provide_PARA(ground);
            ground = provide_STATVAR(ground);
            ground = provide_CONST(ground);
            
            ground = define_NPP_variables_Frolking(ground);

            
            ground.STATVAR.total_peat = [];
            ground.STATVAR.total_peat_PFT = zeros(0,size(ground.PARA.initial_decomposability,2));
            ground.STATVAR.totalpeatC_originalMass = zeros(0,size(ground.PARA.initial_decomposability,2));
            

            ground.STATVAR.next_decompose_timestamp = ceil(tile.FORCING.PARA.start_time + ground.PARA.BGC_timestep);
            ground.STATVAR.next_accumulate_timestamp = ground.STATVAR.next_decompose_timestamp;
            
            ground.TEMP.year_count = 0;
            ground.TEMP.year = -100; %str2num(datestr(tile.FORCING.PARA.startt,'yyyy'));
            ground.TEMP.TDD = 0;
            ground.TEMP.growing_days_count = 0;
            ground.TEMP.accumulation_new_cell_flag = 1;
            ground.TEMP.regrid_flag = 0;
            
            ground.TEMP.d_layerThick = 0;
            ground.TEMP.d_organic = 0;
            
            ground.TEMP.GPP_acc = 0;
            ground.TEMP.GPP_acc2 = 0;
            ground.TEMP.water_table_depth_acc = 0;
            
            ground.STATVAR.water_table_depth = 0;
            ground.STATVAR.peat_depth = 0;
            ground = get_annual_NPP_Frolking(ground);
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            get_NPP_Frolking(ground.IA_BGC, tile);
        end
        
        
        function ground = get_boundary_condition_l(ground, tile)

        end
        
        function ground = get_derivatives_prognostic(ground, tile)
                        
            if tile.t == ground.STATVAR.next_decompose_timestamp
                get_ground_variables(ground.IA_BGC, tile); %get the physical variables over to the BGC class, so that it starts with correct state from the last timestep
               
%                 ground = temp_modifier_Frolking(ground);
%                 ground = water_modifier_Frolking(ground);
%                 ground = peat_decompose_Frolking(ground);
                
                ground = temp_modifier(ground);
                ground = water_modifier(ground);
                ground = peat_decompose_Frolking(ground);
                
                % extracting date, month and year from the tile
                Date_vector = datevec(tile.t);
                accumulation_year = Date_vector(:, 1);
                
                
                %reset  on 1 Jan
                if accumulation_year ~= ground.TEMP.year
                    ground.TEMP.year = accumulation_year;

                    ground.TEMP.accumulation_new_cell_flag = 1; 

                    ground.TEMP.regrid_flag = 1;
                    
                    
                end
                
                % peat accumulates 
                if ground.TEMP.accumulation_new_cell_flag == 1
                     %compute average annual water table 
                    ground.STATVAR.water_table_depth =[ground.STATVAR.water_table_depth; ground.TEMP.water_table_depth_acc./ground.TEMP.GPP_acc2];
                    if size(ground.STATVAR.water_table_depth,1)>5
                        ground.STATVAR.water_table_depth = ground.STATVAR.water_table_depth(2:end,1);
                    end
                     
                    ground.TEMP.water_table_depth_acc = 0;                   
                    ground.TEMP.GPP_acc2 = 0;
                    
                    %compute peat depth
                    ground.STATVAR.peat_depth = [ground.STATVAR.peat_depth; get_peat_depth(ground.IA_BGC, tile)];
                    if size(ground.STATVAR.peat_depth,1)>5
                        ground.STATVAR.peat_depth = ground.STATVAR.peat_depth(2:end,1);
                    end
                    
                    %get annual NPP for the different PFTs
                    ground = get_annual_NPP_Frolking(ground);
                    
                    ground = peat_accumulation_Frolking_newCell(ground);

                    ground.TEMP.year_count = ground.TEMP.year_count+1;
                    ground.TEMP.accumulation_new_cell_flag = 0;
                    add_grid_cell(ground.IA_BGC, tile);
                    

                else
                    ground = peat_accumulation_Frolking(ground);
                    %ground = peat_accumulation(ground);
                end

                %               ground = updatebulkD(ground);
                
                ground.TEMP.GPP_acc = 0;
                ground.STATVAR.GPP_acc2 = ground.TEMP.GPP_acc2;
                
            end
            
        end
        
        function timestep = get_timestep(ground, tile)
            
            timestep = (ground.STATVAR.next_decompose_timestamp - tile.t) .* ground.CONST.day_sec;

        end
        
        function ground = advance_prognostic(ground, tile)
            %add here
            ground.TEMP.GPP_acc = ground.TEMP.GPP_acc + ground.TEMP.GPP  .* tile.timestep./ground.CONST.day_sec;
            ground.TEMP.GPP_acc2 = ground.TEMP.GPP_acc2 + ground.TEMP.GPP  .* tile.timestep./ground.CONST.day_sec;
            ground.TEMP.water_table_depth_acc = ground.TEMP.water_table_depth_acc + ground.TEMP.water_table_depth .* ground.TEMP.GPP  .* tile.timestep./ground.CONST.day_sec;
            
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            %add here
        end
        
        function ground = compute_diagnostic(ground, tile)
          
            if tile.t == ground.STATVAR.next_decompose_timestamp
                send_BGC_variables(ground.IA_BGC, tile); %send the BGC variables over to the GROUND class
                
                ground.STATVAR.next_decompose_timestamp = ground.STATVAR.next_decompose_timestamp + ground.PARA.BGC_timestep;
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
            ground.STATVAR.next_decompose_timestamp = tile.t + ground.PARA.BGC_timestep;
            [year, ~,~] = datevec(tile.t);
            ground.TEMP.year = year; %str2num(datestr(tile.t, 'yyyy'));
        end

    end
    
end
