%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction between two GROUND
% classes without water cycle
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_BGC_Xice <  IA_BGC
    
    properties
        GROUND
        BGC
    end
    
    methods
        
        function finalize_init(ia_BGC, tile)
            ia_BGC.GROUND.STATVAR.non_BGC_layerThick = ia_BGC.GROUND.STATVAR.layerThick; %layerThick of each grid celll before BGC accumulation
            ia_BGC.GROUND.STATVAR.non_BGC_mineral = ia_BGC.GROUND.STATVAR.mineral; %layerThick of each grid celll before BGC accumulation
            ia_BGC.GROUND.STATVAR.non_BGC_organic = ia_BGC.GROUND.STATVAR.organic; %layerThick of each grid celll before BGC accumulation
            
            ia_BGC.GROUND.STATVAR.first_mineral_cell = 1;
            
            %ia_BGC.BGC.STATVAR.GROUND_overlap_vector = []; %same size as part of GROUND vector that contains BGC grid cells
            ia_BGC.BGC.STATVAR.BGC_overlap_vector = []; %same size as BGC vector
        end
        

        
        function get_ground_variables(ia_BGC, tile)
            
            %vol_water_ground = ia_BGC.GROUND.STATVAR.waterIce ./ (ia_BGC.GROUND.STATVAR.layerThick .* ia_BGC.GROUND.STATVAR.area - ia_BGC.GROUND.STATVAR.XwaterIce); %must be water plus ice, not water, water_modifer would otherwise lead to strong degradation when frozen (= little water)
            vol_water_ground = ia_BGC.GROUND.STATVAR.waterIce ./ (ia_BGC.GROUND.STATVAR.layerThick .* ia_BGC.GROUND.STATVAR.area - ia_BGC.GROUND.STATVAR.XwaterIce - ia_BGC.GROUND.STATVAR.mineral - ia_BGC.GROUND.STATVAR.organic); %must be water plus ice, not water, water_modifer would otherwise lead to strong degradation when frozen (= little water)
            %vol_water is actually saturation in LPJ_peat!!!
            if isempty(ia_BGC.BGC.STATVAR.BGC_overlap_vector) %before 1st BGC grid cell is present
                ia_BGC.BGC.STATVAR.T = ia_BGC.GROUND.STATVAR.T(1,1);
                ia_BGC.BGC.STATVAR.vol_water  = vol_water_ground(1,1);
                ia_BGC.BGC.STATVAR.field_capacity = ia_BGC.GROUND.STATVAR.field_capacity(1,1);
            else %normal BGC operations
                for i=1:size(ia_BGC.BGC.STATVAR.BGC_overlap_vector,1)
                    ia_BGC.BGC.STATVAR.T(i,1) = ia_BGC.GROUND.STATVAR.T(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1), 1);
                    ia_BGC.BGC.STATVAR.vol_water(i,1) = vol_water_ground(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1), 1);
                    ia_BGC.BGC.STATVAR.field_capacity(i,1) = ia_BGC.GROUND.STATVAR.field_capacity(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1);
                end
            end
        end
           
        
        function send_BGC_variables(ia_BGC, tile)
%             disp('Hallo1')
%             disp(ia_BGC.BGC.STATVAR.BGC_overlap_vector)
%             disp('Hallo2')
%             disp(ia_BGC.BGC.TEMP.d_organic)

            d_layerThick = ia_BGC.GROUND.STATVAR.layerThick .* 0;
            d_organic = ia_BGC.GROUND.STATVAR.layerThick .* 0;
            for i=1:size(ia_BGC.BGC.STATVAR.BGC_overlap_vector,1)
                %change in layerThick, organic and energy
                index = ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1);
                
                d_organic(index, 1) = d_organic(index, 1) + ia_BGC.BGC.TEMP.d_organic(i,1) .* ia_BGC.GROUND.STATVAR.area(index, 1);
                ia_BGC.GROUND.STATVAR.energy(index, 1) = ia_BGC.GROUND.STATVAR.energy(index, 1) + ia_BGC.BGC.TEMP.d_organic(i,1) .* ia_BGC.GROUND.STATVAR.T(index, 1) .* ...
                    ia_BGC.GROUND.CONST.c_o .* ia_BGC.GROUND.STATVAR.area(index, 1);
                
                d_layerThick(index, 1) = d_layerThick(index, 1) + ia_BGC.BGC.TEMP.d_layerThick(i,1);

            end


            plus_LT = d_layerThick > 0;
            minus_LT = d_layerThick < 0;

            %plus_LT, layerThick increases, "void space" (air) of
            %d_layerThic - d_organic created, which can be filled by
            %XwaterIce if available. 1. determine d_air; 2. redistribute XwaterIce to waterIce up to d_air 3. increase layerThick if necessary 
            
            d_air_pot = d_layerThick(plus_LT) .* ia_BGC.GROUND.STATVAR.area(plus_LT) - d_organic(plus_LT);
            d_air = max(0, d_air_pot - ia_BGC.GROUND.STATVAR.XwaterIce(plus_LT));
            water_change = d_air_pot-d_air;
            ia_BGC.GROUND.STATVAR.XwaterIce(plus_LT) = ia_BGC.GROUND.STATVAR.XwaterIce(plus_LT) - water_change;
            ia_BGC.GROUND.STATVAR.waterIce(plus_LT) = ia_BGC.GROUND.STATVAR.waterIce(plus_LT) + water_change;
            ia_BGC.GROUND.STATVAR.layerThick(plus_LT) = ia_BGC.GROUND.STATVAR.layerThick(plus_LT) + (d_air + d_organic(plus_LT)) ./ ia_BGC.GROUND.STATVAR.area(plus_LT);
            
            %minus_LT, layerThick decreases, air is pressed out, oif no air
            %available, waterIce is converted to XwaterIce  
            
            d_air_pot = -(d_layerThick(minus_LT) .* ia_BGC.GROUND.STATVAR.area(minus_LT) - d_organic(minus_LT));
            d_air_max = max(0, ia_BGC.GROUND.STATVAR.layerThick(minus_LT) .* ia_BGC.GROUND.STATVAR.area(minus_LT) - ia_BGC.GROUND.STATVAR.XwaterIce(minus_LT) - ...
                 ia_BGC.GROUND.STATVAR.waterIce(minus_LT) - ia_BGC.GROUND.STATVAR.mineral(minus_LT) - ia_BGC.GROUND.STATVAR.organic(minus_LT) );
             water_change = max(0, d_air_pot - d_air_max);
             d_air = -(d_air_pot - water_change);
             
             ia_BGC.GROUND.STATVAR.XwaterIce(minus_LT) = ia_BGC.GROUND.STATVAR.XwaterIce(minus_LT) + water_change;
             ia_BGC.GROUND.STATVAR.waterIce(minus_LT) = ia_BGC.GROUND.STATVAR.waterIce(minus_LT) - water_change;
             ia_BGC.GROUND.STATVAR.layerThick(minus_LT) = ia_BGC.GROUND.STATVAR.layerThick(minus_LT) + (d_air + d_organic(minus_LT)) ./ ia_BGC.GROUND.STATVAR.area(minus_LT);

             %update organic
             ia_BGC.GROUND.STATVAR.organic = ia_BGC.GROUND.STATVAR.organic + d_organic;
            
%             ia_BGC.GROUND.STATVAR.layerThick(minus_LT) = ia_BGC.GROUND.STATVAR.layerThick(minus_LT) + d_layerThick(minus_LT);
%             %adjust waterIce and XwaterIce
%             residual_waterIce =  max(0, ia_BGC.GROUND.STATVAR.XwaterIce(minus_LT) + ia_BGC.GROUND.STATVAR.waterIce + ia_BGC.GROUND.STATVAR.mineral + ia_BGC.GROUND.STATVAR.organic - ia_BGC.GROUND.STATVAR.layerThick .* ia_BGC.GROUND.STATVAR.area);
%             ia_BGC.GROUND.STATVAR.layerThick = max(ia_BGC.GROUND.STATVAR.layerThick, ...
%                 (ia_BGC.GROUND.STATVAR.XwaterIce + ia_BGC.GROUND.STATVAR.waterIce + ia_BGC.GROUND.STATVAR.mineral + ia_BGC.GROUND.STATVAR.organic)./ia_BGC.GROUND.STATVAR.area); 
%             ia_BGC.GROUND.STATVAR.waterIce = ia_BGC.GROUND.STATVAR.waterIce - residual_waterIce;
%             ia_BGC.GROUND.STATVAR.XwaterIce = ia_BGC.GROUND.STATVAR.XwaterIce + residual_waterIce;
            
        end
        
        function add_grid_cell(ia_BGC, tile)
            ia_BGC.BGC.STATVAR.BGC_overlap_vector = [1; ia_BGC.BGC.STATVAR.BGC_overlap_vector];
            
        end
        
        function regrid_stratigraphy(ia_BGC, tile)
            %split up the extensive prognostic variables
            if ~isempty(ia_BGC.BGC.STATVAR.BGC_overlap_vector)
                BGC_organic = ia_BGC.GROUND.STATVAR.organic - ia_BGC.GROUND.STATVAR.non_BGC_organic;
                non_BGC_organic = ia_BGC.GROUND.STATVAR.non_BGC_organic;
                BGC_layerThick = ia_BGC.GROUND.STATVAR.layerThick - ia_BGC.GROUND.STATVAR.non_BGC_layerThick;
                
                saturation_waterIce = ia_BGC.GROUND.STATVAR.waterIce ./ (ia_BGC.GROUND.STATVAR.layerThick .* ia_BGC.GROUND.STATVAR.area - ia_BGC.GROUND.STATVAR.XwaterIce - ia_BGC.GROUND.STATVAR.mineral - ia_BGC.GROUND.STATVAR.organic);
                fraction_XwaterIce = ia_BGC.GROUND.STATVAR.XwaterIce ./ (ia_BGC.GROUND.STATVAR.layerThick .* ia_BGC.GROUND.STATVAR.area - ia_BGC.GROUND.STATVAR.XwaterIce);
                
                %split up energy
                energy_BGC_organic = ia_BGC.GROUND.STATVAR.T .* BGC_organic .* ia_BGC.GROUND.CONST.c_o;
                energy_non_BGC_organic = ia_BGC.GROUND.STATVAR.T .* non_BGC_organic .* ia_BGC.GROUND.CONST.c_o;
                energy_mineral = ia_BGC.GROUND.STATVAR.T .* ia_BGC.GROUND.STATVAR.mineral .* ia_BGC.GROUND.CONST.c_m;
                energy_waterIce = ia_BGC.GROUND.STATVAR.energy - energy_BGC_organic - energy_non_BGC_organic - energy_mineral;
                
                
                waterIce_BGC = ia_BGC.BGC.STATVAR.layerThick .*0;
                XwaterIce_BGC = ia_BGC.BGC.STATVAR.layerThick .*0;
                organic_BGC = ia_BGC.BGC.STATVAR.layerThick .*0;
                energy_BGC = ia_BGC.BGC.STATVAR.layerThick .*0;
                
                i=1;
                j=1;
                while i<=size(ia_BGC.BGC.STATVAR.BGC_overlap_vector,1)
                    while j<=size(ia_BGC.BGC.STATVAR.BGC_overlap_vector,1) && ia_BGC.BGC.STATVAR.BGC_overlap_vector(j,1) == ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1)
                        j=j+1;
                    end
                    index = i:j-1;
                    %assign to BGC class grid cells
                    organic_BGC(index,1) = ia_BGC.BGC.STATVAR.total_peat(index,1) ./ sum(ia_BGC.BGC.STATVAR.total_peat(index,1)) .* BGC_organic(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1);
                    waterIce_BGC(index,1) = (ia_BGC.GROUND.STATVAR.area(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) .* ia_BGC.BGC.STATVAR.layerThick(index,1) - organic_BGC(index,1)) .* saturation_waterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1);
                    XwaterIce_BGC(index,1) = ia_BGC.BGC.STATVAR.layerThick(index,1) .* ia_BGC.GROUND.STATVAR.area(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) .* fraction_XwaterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1);
                    XwaterIce_thickness_BGC(index,1) = ia_BGC.BGC.STATVAR.layerThick(index,1) .* fraction_XwaterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1);
                    area_BGC(index,1) = ia_BGC.GROUND.STATVAR.area(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1);
                    energy_BGC(index,1) = organic_BGC(index,1)./sum(organic_BGC(index,1),1) .* energy_BGC_organic(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) + ...
                        (waterIce_BGC(index,1) + XwaterIce_BGC(index,1)) ./ max(1e-40, ia_BGC.GROUND.STATVAR.waterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) + ia_BGC.GROUND.STATVAR.XwaterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1)) .* energy_waterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1);
%                     energy_BGC(index,1) = organic_BGC(index,1)./sum(organic_BGC(index,1),1) .* energy_BGC_organic(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) + ...
%                         (waterIce_BGC(index,1) + XwaterIce_BGC(index,1)) ./ max(1e-40, sum(waterIce_BGC(index,1) + XwaterIce_BGC(index,1),1)) .* energy_waterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1);
                    %subtract from GROUND class grid cells
                    ia_BGC.GROUND.STATVAR.layerThick(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) = ia_BGC.GROUND.STATVAR.layerThick(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) - sum(ia_BGC.BGC.STATVAR.layerThick(index,1),1) - sum(XwaterIce_thickness_BGC(index,1),1);
                    ia_BGC.GROUND.STATVAR.organic(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) = ia_BGC.GROUND.STATVAR.organic(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) - sum(organic_BGC(index,1),1);
                    ia_BGC.GROUND.STATVAR.waterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) = ia_BGC.GROUND.STATVAR.waterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) - sum(waterIce_BGC(index,1),1);
                    ia_BGC.GROUND.STATVAR.XwaterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) = ia_BGC.GROUND.STATVAR.XwaterIce(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) - sum(XwaterIce_BGC(index,1),1);
                    ia_BGC.GROUND.STATVAR.energy(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) = ia_BGC.GROUND.STATVAR.energy(ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1),1) - sum(energy_BGC(index,1),1);
                    
                    i=j;
                end
                
                %reassemble the stratigraphy top_down
                %- make new ia_BGC.BGC.STATVAR.BGC_overlap_vector according to
                %LayerThick_BGC
                layerThick_phys = [];
                organic_phys = [];
                waterIce_phys = [];
                XwaterIce_phys = [];
                area_phys=[];
                energy_phys = [];
                
                %layerThick_BGC_acc = cumsum(ia_BGC.BGC.STATVAR.layerThick + XwaterIce_thickness_BGC);
                layerThick_BGC_acc = cumsum(ia_BGC.BGC.STATVAR.layerThick);
                 
                i=1;
                start_pos =1;
                while i.*ia_BGC.GROUND.PARA.target_grid_cell_size<=layerThick_BGC_acc(end,1)
                    [~, pos] = min((i.*ia_BGC.GROUND.PARA.target_grid_cell_size - layerThick_BGC_acc).^2);
                    layerThick_phys = [layerThick_phys; sum(ia_BGC.BGC.STATVAR.layerThick(start_pos:pos,1) + XwaterIce_thickness_BGC(start_pos:pos,1))];
                    organic_phys = [organic_phys; sum(organic_BGC(start_pos:pos,1))];
                    waterIce_phys = [waterIce_phys; sum(waterIce_BGC(start_pos:pos,1))];
                    XwaterIce_phys = [XwaterIce_phys; sum(XwaterIce_BGC(start_pos:pos,1))];
                    area_phys=[area_phys; mean(area_BGC(start_pos:pos,1))];
                    energy_phys = [energy_phys; sum(energy_BGC(start_pos:pos,1))];
                    ia_BGC.BGC.STATVAR.BGC_overlap_vector(start_pos:pos,1) = i;
                    start_pos = pos+1;
                    i=i+1;
                end
                ia_BGC.BGC.STATVAR.BGC_overlap_vector(start_pos:end,1) = i;
                
                %reassemble stratigraphy
                first_mineral_cell_old = ia_BGC.GROUND.STATVAR.first_mineral_cell;
                ia_BGC.GROUND.STATVAR.first_mineral_cell = i;
                
                ia_BGC.GROUND.STATVAR.layerThick = [layerThick_phys; ia_BGC.GROUND.STATVAR.layerThick(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.organic = [organic_phys; ia_BGC.GROUND.STATVAR.organic(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.waterIce = [waterIce_phys; ia_BGC.GROUND.STATVAR.waterIce(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.XwaterIce = [XwaterIce_phys; ia_BGC.GROUND.STATVAR.XwaterIce(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.area = [area_phys; ia_BGC.GROUND.STATVAR.area(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.energy = [energy_phys; ia_BGC.GROUND.STATVAR.energy(first_mineral_cell_old:end,1)];
                
                ia_BGC.GROUND.STATVAR.mineral = [energy_phys.*0; ia_BGC.GROUND.STATVAR.mineral(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.soil_type = [energy_phys.*0+4; ia_BGC.GROUND.STATVAR.soil_type(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.n = [energy_phys.*0 + ia_BGC.GROUND.CONST.vanGen_n(1,4); ia_BGC.GROUND.STATVAR.n(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.alpha = [energy_phys.*0 + ia_BGC.GROUND.CONST.vanGen_alpha(1,4); ia_BGC.GROUND.STATVAR.alpha(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.field_capacity = [energy_phys.*0 + 0.5; ia_BGC.GROUND.STATVAR.field_capacity(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.satHydraulicConductivity = [energy_phys.*0 + 1e-5; ia_BGC.GROUND.STATVAR.satHydraulicConductivity(first_mineral_cell_old:end,1)];
                
                ia_BGC.GROUND.STATVAR.non_BGC_layerThick = [energy_phys.*0; ia_BGC.GROUND.STATVAR.non_BGC_layerThick(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.non_BGC_mineral = [energy_phys.*0; ia_BGC.GROUND.STATVAR.non_BGC_mineral(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.non_BGC_organic = [energy_phys.*0; ia_BGC.GROUND.STATVAR.non_BGC_organic(first_mineral_cell_old:end,1)];

                
                %add remaining BGC cells to the first mineral grid cell
                pos = size(ia_BGC.BGC.STATVAR.layerThick,1);
                ia_BGC.GROUND.STATVAR.layerThick(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.layerThick(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(ia_BGC.BGC.STATVAR.layerThick(start_pos:pos,1) + XwaterIce_thickness_BGC(start_pos:pos,1));
                ia_BGC.GROUND.STATVAR.organic(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.organic(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(organic_BGC(start_pos:pos,1));
                ia_BGC.GROUND.STATVAR.waterIce(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.waterIce(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(waterIce_BGC(start_pos:pos,1));
                ia_BGC.GROUND.STATVAR.XwaterIce(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.XwaterIce(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(XwaterIce_BGC(start_pos:pos,1));
                %ia_BGC.GROUND.STATVAR.area(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.area(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(area_BGC(start_pos:pos,1));
                ia_BGC.GROUND.STATVAR.energy(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.energy(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(energy_BGC(start_pos:pos,1));
                
                ia_BGC.GROUND = get_T_water_freezeC_Xice(ia_BGC.GROUND);
                
                ia_BGC.GROUND = conductivity(ia_BGC.GROUND);
                ia_BGC.GROUND = calculate_hydraulicConductivity_Xice(ia_BGC.GROUND);
                
                ia_BGC.GROUND = set_TEMP_2zero(ia_BGC.GROUND);
            end
            %ia_BGC.GROUND.STATVAR.first_mineral_cell = 1;
            %-check the correct length of the GROUND vector, maybe make
            %varibale for the 1st mineral grid cell position
            %-go through the new BGC_overlap_vector and reassemble GROUND
            %from BGC
            %-compute_diagnostic for the GROUND
        end
        
        function add_new_BGC_cells_from_READ(ia_BGC, tile)  %add the newly gained peat cells to first GROUND grid cell after switching from READ
           %CHNAGE THIS COMPLETELY, FOR EACH CELL COMPARE BEFORE AN AFTER
           %and add/subtract - make loop over entire overlap_vector
           % if cell has shrunk, route water to Xwater, if it increases,
           % add water an energy according to gain, so that saturation and
           % XwaterFraction are conserved
           
           i=1;
           while ~isempty(find(ia_BGC.BGC.STATVAR.BGC_overlap_vector == i))
            
               new_cells = find(ia_BGC.BGC.STATVAR.BGC_overlap_vector == i);
               BGC_layerThick = sum(ia_BGC.BGC.STATVAR.layerThick(new_cells),1);
               GROUND_layerThick = ia_BGC.GROUND.STATVAR.layerThick(i,1) - ia_BGC.GROUND.STATVAR.XwaterIce(i,1) ./ ia_BGC.GROUND.STATVAR.area(i,1) - ia_BGC.GROUND.STATVAR.non_BGC_layerThick(i,1);
               
               if BGC_layerThick > GROUND_layerThick %cell size has increased
                   
                   T_GROUND = ia_BGC.GROUND.STATVAR.T(i,1);
                   saturation_waterIce = ia_BGC.GROUND.STATVAR.waterIce(i,1) ./ (ia_BGC.GROUND.STATVAR.layerThick(i,1) .* ia_BGC.GROUND.STATVAR.area(i,1) - ia_BGC.GROUND.STATVAR.XwaterIce(i,1) - ia_BGC.GROUND.STATVAR.mineral(i,1) - ia_BGC.GROUND.STATVAR.organic(i,1));
                   fraction_XwaterIce = ia_BGC.GROUND.STATVAR.XwaterIce(i,1) ./ (ia_BGC.GROUND.STATVAR.layerThick(i,1) .* ia_BGC.GROUND.STATVAR.area(i,1) - ia_BGC.GROUND.STATVAR.XwaterIce(i,1));
                   energy_density_waterIce = (ia_BGC.GROUND.STATVAR.energy(i,1) - T_GROUND .* ia_BGC.GROUND.STATVAR.organic(i,1) .* ia_BGC.GROUND.CONST.c_o - ...
                       - T_GROUND .* ia_BGC.GROUND.STATVAR.mineral(i,1) .* ia_BGC.GROUND.CONST.c_m) ./ (ia_BGC.GROUND.STATVAR.waterIce(i,1) + ia_BGC.GROUND.STATVAR.XwaterIce(i,1));
                   
                   ia_BGC.GROUND.STATVAR.layerThick(i,1) = (sum(ia_BGC.BGC.STATVAR.layerThick(new_cells),1) + ia_BGC.GROUND.STATVAR.non_BGC_layerThick(i,1)) .* (1 + fraction_XwaterIce);
                   ia_BGC.GROUND.STATVAR.organic(i,1) = sum(ia_BGC.BGC.STATVAR.total_peat(new_cells),1) ./ ia_BGC.BGC.CONST.organicDensity .* ia_BGC.GROUND.STATVAR.area(i,1) + ia_BGC.GROUND.STATVAR.non_BGC_organic(i,1);
                   ia_BGC.GROUND.STATVAR.XwaterIce(i,1) = ia_BGC.GROUND.STATVAR.layerThick(i,1) .* fraction_XwaterIce ./ (1 + fraction_XwaterIce) .* ia_BGC.GROUND.STATVAR.area(i,1);
                   ia_BGC.GROUND.STATVAR.waterIce(i,1) =(ia_BGC.GROUND.STATVAR.layerThick(i,1) .* ia_BGC.GROUND.STATVAR.area(i,1) - ia_BGC.GROUND.STATVAR.XwaterIce(i,1)- ia_BGC.GROUND.STATVAR.organic(i,1) - ia_BGC.GROUND.STATVAR.mineral(i,1)) .* saturation_waterIce; 
                   %ia_BGC.GROUND.STATVAR.waterIce(i,1) = (sum(ia_BGC.BGC.STATVAR.layerThick(new_cells),1) - sum(ia_BGC.BGC.STATVAR.total_peat(new_cells),1) ./ ia_BGC.BGC.CONST.organicDensity) .* saturation_waterIce .* ia_BGC.GROUND.STATVAR.area(1,1);
                  
                   %ia_BGC.GROUND.STATVAR.XwaterIce(i,1) = sum(ia_BGC.BGC.STATVAR.layerThick(new_cells),1) .* fraction_XwaterIce .* ia_BGC.GROUND.STATVAR.area(i,1);
                   ia_BGC.GROUND.STATVAR.energy(i,1) =  T_GROUND .* (ia_BGC.GROUND.STATVAR.organic(i,1) .* ia_BGC.GROUND.CONST.c_o + ia_BGC.GROUND.STATVAR.mineral(i,1) .* ia_BGC.GROUND.CONST.c_m);
                   ia_BGC.GROUND.STATVAR.energy(i,1) = ia_BGC.GROUND.STATVAR.energy(i,1) + energy_density_waterIce .* (ia_BGC.GROUND.STATVAR.waterIce(i,1) + ia_BGC.GROUND.STATVAR.XwaterIce(i,1));
                   
               elseif BGC_layerThick < GROUND_layerThick %cell size has decreased
                   d_layerThick = GROUND_layerThick - BGC_layerThick;
                   d_organic = ia_BGC.GROUND.STATVAR.organic(i,1) -  ia_BGC.GROUND.STATVAR.non_BGC_organic(i,1) - sum(ia_BGC.BGC.STATVAR.total_peat(new_cells),1) ./ ia_BGC.BGC.CONST.organicDensity .* ia_BGC.GROUND.STATVAR.area(i,1);
                   
                   ia_BGC.GROUND.STATVAR.layerThick(i,1) = ia_BGC.GROUND.STATVAR.layerThick(i,1) - d_layerThick ;
                   ia_BGC.GROUND.STATVAR.organic(i,1) = ia_BGC.GROUND.STATVAR.organic(i,1) - d_organic;
                   ia_BGC.GROUND.STATVAR.energy(i,1) = ia_BGC.GROUND.STATVAR.energy(i,1) - d_organic .* T_GROUND .* ia_BGC.GROUND.CONST.c_o;
                   new_waterIce = min(ia_BGC.GROUND.STATVAR.waterIce(i,1), ia_BGC.GROUND.STATVAR.layerThick(i,1) .* ia_BGC.GROUND.STATVAR.area(i,1) - ia_BGC.GROUND.STATVAR.organic(i,1) - ia_BGC.GROUND.STATVAR.mineral(i,1));
                   difference_waterIce = ia_BGC.GROUND.STATVAR.waterIce(i,1) - new_waterIce;
                   ia_BGC.GROUND.STATVAR.waterIce(i,1) = new_waterIce;
                   ia_BGC.GROUND.STATVAR.XwaterIce(i,1) = ia_BGC.GROUND.STATVAR.XwaterIce(i,1) + difference_waterIce;
                   ia_BGC.GROUND.STATVAR.layerThick(i,1) = ia_BGC.GROUND.STATVAR.layerThick(i,1) +  ia_BGC.GROUND.STATVAR.XwaterIce(i,1)  ./  ia_BGC.GROUND.STATVAR.area(i,1); 
                   %no change in energy
                   
               end
               i=i+1;
           end
           
           regrid_stratigraphy(ia_BGC, tile);
        end
        
        
        function get_NPP_Frolking(ia_BGC, tile)
            T = ia_BGC.GROUND.STATVAR.T(1,1);
            temp_modifier = max(0,min(1,min((T - ia_BGC.BGC.PARA.start_PhotSyn_T)./(ia_BGC.BGC.PARA.start_fullPhotSyn_T - ia_BGC.BGC.PARA.start_PhotSyn_T), ...
                (T - ia_BGC.BGC.PARA.end_PhotSyn_T)./ (ia_BGC.BGC.PARA.end_fullPhotSyn_T - ia_BGC.BGC.PARA.end_PhotSyn_T))));
            %disp([T temp_modifier])
            radiation_modifier = tile.FORCING.TEMP.Sin;
            ET_modifier = max(0, ia_BGC.GROUND.STATVAR.Qe ./ ia_BGC.GROUND.STATVAR.Qe_pot);
            
            ia_BGC.BGC.TEMP.GPP = max(0, temp_modifier .* radiation_modifier .* ET_modifier);
            
            %get water table depth
            saturation = ia_BGC.GROUND.STATVAR.waterIce ./ (ia_BGC.GROUND.STATVAR.layerThick .* ia_BGC.GROUND.STATVAR.area - ia_BGC.GROUND.STATVAR.XwaterIce - ia_BGC.GROUND.STATVAR.mineral - ia_BGC.GROUND.STATVAR.organic);
            ia_BGC.BGC.TEMP.water_table_depth = sum(ia_BGC.GROUND.STATVAR.layerThick(1:find(saturation>0.95,1)-1,1),1) - ia_BGC.GROUND.STATVAR.XwaterIce(1,1)./ ia_BGC.GROUND.STATVAR.area(1,1);
        end
        
        function peat_depth = get_peat_depth(ia_BGC, tile)
        
            peat_depth = sum(ia_BGC.GROUND.STATVAR.layerThick,1) - sum(ia_BGC.GROUND.STATVAR.non_BGC_layerThick,1) - ia_BGC.GROUND.STATVAR.XwaterIce(1,1)./ ia_BGC.GROUND.STATVAR.area(1,1);
        end
        
%         function send_BGC_variables(ia_BGC, tile)
%             %do nothing at this point - for fully coupled runs, this must
%             %must include transfer of mass to and within
%             %the physics stratigraphy
%         end
%         
%         function get_ground_variables(ia_BGC, tile)
%             %map the physical variables required by the BGC module to the BGC grid
%             depths_ground =  cumsum([0; ia_BGC.GROUND.STATVAR.layerThick]);
%             depths_BGC =  cumsum([0; ia_BGC.BGC.STATVAR.layerThick]);
%             %set everything to zero
%             ia_BGC.BGC.STATVAR.vol_water = ia_BGC.BGC.STATVAR.layerThick .*0;
%             ia_BGC.BGC.STATVAR.vol_mineral = ia_BGC.BGC.STATVAR.layerThick .*0;
%             ia_BGC.BGC.STATVAR.porosity = ia_BGC.BGC.STATVAR.layerThick .*0;
%             if ~isempty(ia_BGC.BGC.STATVAR.layerThick)
%                 ia_BGC.BGC.STATVAR.T = ia_BGC.BGC.STATVAR.layerThick .*0;
%             else
%                 ia_BGC.BGC.STATVAR.T = ia_BGC.GROUND.STATVAR.T(1,1);
%             end
%             ia_BGC.BGC.STATVAR.field_capacity = ia_BGC.BGC.STATVAR.layerThick .*0;
%             vol_water_ground = ia_BGC.GROUND.STATVAR.water ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;
%             vol_mineral_ground = ia_BGC.GROUND.STATVAR.mineral ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;
%             porosity_ground = 1- (ia_BGC.GROUND.STATVAR.mineral - ia_BGC.GROUND.STATVAR.organic) ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;
%             
%             
%             overlap = get_overlap_cells(ia_BGC, depths_BGC, depths_ground);
%             ia_BGC.BGC.STATVAR.overlap = overlap;
%             
%             for i=1:size(overlap,1)%                 ia_BGC.BGC.STATVAR.T(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.T(overlap(i, 1),1) + overlap(i,3) .* ia_BGC.GROUND.STATVAR.T(overlap(i, 2),1);
%                 ia_BGC.BGC.STATVAR.vol_water(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.vol_water(overlap(i, 1),1) + overlap(i,3) .* vol_water_ground(overlap(i, 2),1);
%                 ia_BGC.BGC.STATVAR.vol_mineral(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.vol_mineral(overlap(i, 1),1) + overlap(i,3) .* vol_mineral_ground(overlap(i, 2),1);
%                 ia_BGC.BGC.STATVAR.porosity(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.porosity(overlap(i, 1),1) + overlap(i,3) .* porosity_ground(overlap(i, 2),1);
%                 ia_BGC.BGC.STATVAR.field_capacity(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.field_capacity(overlap(i, 1),1) + overlap(i,3) .* ia_BGC.GROUND.STATVAR.field_capacity(overlap(i, 2),1);
%             end
%             
%         end
    end
end