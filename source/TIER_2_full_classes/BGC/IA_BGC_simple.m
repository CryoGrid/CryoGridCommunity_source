%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction between two GROUND
% classes without water cycle
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_BGC_simple <  IA_BGC
    
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
            
            vol_water_ground = ia_BGC.GROUND.STATVAR.water ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;

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
            for i=1:size(ia_BGC.BGC.STATVAR.BGC_overlap_vector,1)
                %change in layerThick, organic and energy
                index = ia_BGC.BGC.STATVAR.BGC_overlap_vector(i,1);
                
                ia_BGC.GROUND.STATVAR.organic(index, 1) = ia_BGC.GROUND.STATVAR.organic(index, 1) + ia_BGC.BGC.TEMP.d_organic(i,1) .* ia_BGC.GROUND.STATVAR.area(index, 1);
                ia_BGC.GROUND.STATVAR.energy(index, 1) = ia_BGC.GROUND.STATVAR.energy(index, 1) + ia_BGC.BGC.TEMP.d_organic(i,1) .* ia_BGC.GROUND.STATVAR.T(index, 1) .* ...
                    ia_BGC.GROUND.CONST.c_o .* ia_BGC.GROUND.STATVAR.area(index, 1);
                
                ia_BGC.GROUND.STATVAR.layerThick(index, 1) = ia_BGC.GROUND.STATVAR.layerThick(index, 1) + ia_BGC.BGC.TEMP.d_layerThick(i,1);

            end
            
            %adjust waterIce and XwaterIce
            residual_waterIce =  max(0, ia_BGC.GROUND.STATVAR.XwaterIce + ia_BGC.GROUND.STATVAR.waterIce + ia_BGC.GROUND.STATVAR.mineral + ia_BGC.GROUND.STATVAR.organic - ia_BGC.GROUND.STATVAR.layerThick .* ia_BGC.GROUND.STATVAR.area);
            ia_BGC.GROUND.STATVAR.layerThick = max(ia_BGC.GROUND.STATVAR.layerThick, ...
                (ia_BGC.GROUND.STATVAR.XwaterIce + ia_BGC.GROUND.STATVAR.waterIce + ia_BGC.GROUND.STATVAR.mineral + ia_BGC.GROUND.STATVAR.organic)./ia_BGC.GROUND.STATVAR.area); 
            ia_BGC.GROUND.STATVAR.waterIce = ia_BGC.GROUND.STATVAR.waterIce - residual_waterIce;
            ia_BGC.GROUND.STATVAR.XwaterIce = ia_BGC.GROUND.STATVAR.XwaterIce + residual_waterIce;
            
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
                
                layerThick_BGC_acc = cumsum(ia_BGC.BGC.STATVAR.layerThick + XwaterIce_thickness_BGC);
                 
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
                
                first_mineral_cell_old = ia_BGC.GROUND.STATVAR.first_mineral_cell;
                ia_BGC.GROUND.STATVAR.first_mineral_cell = i;
                
                ia_BGC.GROUND.STATVAR.layerThick = [layerThick_phys; ia_BGC.GROUND.STATVAR.layerThick(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.organic = [organic_phys; ia_BGC.GROUND.STATVAR.organic(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.waterIce = [waterIce_phys; ia_BGC.GROUND.STATVAR.waterIce(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.XwaterIce = [XwaterIce_phys; ia_BGC.GROUND.STATVAR.XwaterIce(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.area = [area_phys; ia_BGC.GROUND.STATVAR.area(first_mineral_cell_old:end,1)];
                ia_BGC.GROUND.STATVAR.energy = [energy_phys; ia_BGC.GROUND.STATVAR.energy(first_mineral_cell_old:end,1)];
                
                %add remaining BGC cells to the first mineral grid cell
                pos = size(ia_BGC.BGC.STATVAR.layerThick,1);
                ia_BGC.GROUND.STATVAR.layerThick(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.layerThick(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(ia_BGC.BGC.STATVAR.layerThick(start_pos:pos,1) + XwaterIce_thickness_BGC(start_pos:pos,1));
                ia_BGC.GROUND.STATVAR.organic(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.organic(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(organic_BGC(start_pos:pos,1));
                ia_BGC.GROUND.STATVAR.waterIce(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.waterIce(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(waterIce_BGC(start_pos:pos,1));
                ia_BGC.GROUND.STATVAR.XwaterIce(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.XwaterIce(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(XwaterIce_BGC(start_pos:pos,1));
                %ia_BGC.GROUND.STATVAR.area(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.area(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(area_BGC(start_pos:pos,1));
                ia_BGC.GROUND.STATVAR.energy(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) = ia_BGC.GROUND.STATVAR.energy(ia_BGC.GROUND.STATVAR.first_mineral_cell,1) + sum(energy_BGC(start_pos:pos,1));
                
                
                ia_BGC.GROUND = compute_diagnostic(ia_BGC.GROUND, tile);
            end
            %ia_BGC.GROUND.STATVAR.first_mineral_cell = 1;
            %-check the correct length of the GROUND vector, maybe make
            %varibale for the 1st mineral grid cell position
            %-go through the new BGC_overlap_vector and reassemble GROUND
            %from BGC
            %-compute_diagnostic for the GROUND
            
            
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
%             for i=1:size(overlap,1)
%                 ia_BGC.BGC.STATVAR.T(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.T(overlap(i, 1),1) + overlap(i,3) .* ia_BGC.GROUND.STATVAR.T(overlap(i, 2),1);
%                 ia_BGC.BGC.STATVAR.vol_water(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.vol_water(overlap(i, 1),1) + overlap(i,3) .* vol_water_ground(overlap(i, 2),1);
%                 ia_BGC.BGC.STATVAR.vol_mineral(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.vol_mineral(overlap(i, 1),1) + overlap(i,3) .* vol_mineral_ground(overlap(i, 2),1);
%                 ia_BGC.BGC.STATVAR.porosity(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.porosity(overlap(i, 1),1) + overlap(i,3) .* porosity_ground(overlap(i, 2),1);
%                 ia_BGC.BGC.STATVAR.field_capacity(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.field_capacity(overlap(i, 1),1) + overlap(i,3) .* ia_BGC.GROUND.STATVAR.field_capacity(overlap(i, 2),1);
%             end
%             
%         end
    end
end