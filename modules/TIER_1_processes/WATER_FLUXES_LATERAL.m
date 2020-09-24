classdef WATER_FLUXES_LATERAL < BASE


    methods

        
        %remove surface water
        function ground = lateral_push_remove_surfaceWater_simple(ground, lateral)
            lateral.STATVAR.surface_run_off = lateral.STATVAR.surface_run_off + ground.STATVAR.excessWater;
            ground.STATVAR.excessWater = 0;
        end
        
        function ground = lateral_push_remove_surfaceWater_Xice(ground, lateral)
            lateral.STATVAR.surface_run_off = lateral.STATVAR.surface_run_off + ground.STATVAR.Xwater(1);
            ground.STATVAR.XwaterIce(1) = ground.STATVAR.XwaterIce(1) - ground.STATVAR.Xwater(1);
            ground.STATVAR.energy(1) = ground.STATVAR.energy(1) - ground.STATVAR.Xwater(1) .*  ground.CONST.c_w .* ground.STATVAR.T(1);
            ground.STATVAR.layerThick(1) = ground.STATVAR.layerThick(1) - ground.STATVAR.Xwater(1) ./ ground.STATVAR.area(1);
            ground.STATVAR.Xwater(1) = 0;
        end
        
        %remove subsurface water, function also works for Xice
        function ground = lateral_push_remove_subsurfaceWater_simple(ground, lateral)
            exceeding_field_capacity = ground.STATVAR.water > ground.STATVAR.field_capacity .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            excess_water = ground.STATVAR.water(exceeding_field_capacity) - ground.STATVAR.field_capacity(exceeding_field_capacity) .* ground.STATVAR.layerThick(exceeding_field_capacity) .* ground.STATVAR.area(exceeding_field_capacity);
            ground.STATVAR.waterIce(exceeding_field_capacity) = ground.STATVAR.waterIce(exceeding_field_capacity) - excess_water;
            ground.STATVAR.energy(exceeding_field_capacity) = ground.STATVAR.energy(exceeding_field_capacity) - excess_water .* ground.STATVAR.T(exceeding_field_capacity).* ...
                (ground.CONST.c_w .* double(ground.STATVAR.T(exceeding_field_capacity)>=0) + ground.CONST.c_i .* double(ground.STATVAR.T(exceeding_field_capacity)<0)) ;
            ground.STATVAR.water(exceeding_field_capacity) = ground.STATVAR.water(exceeding_field_capacity) - excess_water;
            
            lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off + sum(excess_water);
        end
        
        
        %seepage face
        function ground = lateral_push_remove_water_seepage_simple(ground, lateral)
            water_volumetric = ground.STATVAR.water ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            mobile = water_volumetric  > ground.STATVAR.field_capacity;
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area) > 0.999; %avoid rounding errors
            saturated = [saturated; 0]; %change later, so that first cell of next class is checked -> next class assumed to be hard-bottom now
            hardBottom = water_volumetric <= lateral.PARA.hardBottom_cutoff;
            hardBottom = [hardBottom; 1];
            
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            fraction_above_out = min(1,max(0,(depths(1:end-1,1) - lateral.PARA.upperElevation) ./ ground.STATVAR.layerThick)); %fraction of cell above or below the seepage face
            fraction_below_out = min(1,max(0, (lateral.PARA.lowerElevation - depths(2:end,1))  ./ ground.STATVAR.layerThick));

            for i=1:size(ground.STATVAR.layerThick,1)
                
                %conditions: water in cell must be mobile AND next cell must be either saturated or hardBottom or threshold is reached
                %head
                if (mobile(i,1) && ~hardBottom(i,1)) && (saturated(i+1,1) || hardBottom(i+1,1) || fraction_below_out(i,1) > 0)
                    saturated_height = (water_volumetric(i,1) - ground.STATVAR.field_capacity(i,1)) ./ (porosity(i,1) - ground.STATVAR.field_capacity(i,1)) .* ground.STATVAR.layerThick(i,1);
                    saturated_height = max(0, saturated_height - fraction_below_out(i,1) .* ground.STATVAR.layerThick(i,1)); %seepage face threshold below
                    lateral.TEMP.head = lateral.TEMP.head + saturated_height;                 
                else
                   lateral.TEMP.head = 0;
                end

                
                %flow
                if lateral.TEMP.head >0  %avoid unnecessary computation 
                    saturated_height = min(saturated_height , (1-fraction_above_out(i,1)) .* ground.STATVAR.layerThick(i,1)); %seepage face threshold above
                    cross_section = lateral.PARA.seepage_contact_length .* saturated_height;
                    
                    %replace by better hydraulic conductivity formulation later
                    flux = ground.STATVAR.hydraulicConductivity(i,1) .* (lateral.TEMP.head ./ lateral.PARA.distance_seepageFace .* cross_section);
                    water_removed = min(max(0,ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* ground.STATVAR.layerThick(i,1) .* ground.STATVAR.area(i,1)), flux .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec);
                    ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) - water_removed;
                    ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) - water_removed .* ground.STATVAR.T(i,1).* ...
                        (ground.CONST.c_w .* double(ground.STATVAR.T(i,1)>=0) + ground.CONST.c_i .* double(ground.STATVAR.T(i,1)<0)) ;
                    ground.STATVAR.water(i,1) = ground.STATVAR.water(i,1) - water_removed;
                    
                    lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off + water_removed;
                end
            end
        end
        
        function ground = lateral_push_remove_water_seepage_Xice(ground, lateral)
            water_volumetric = ground.STATVAR.water ./ (ground.STATVAR.layerThick.* ground.STATVAR.area  - ground.STATVAR.XwaterIce);
            porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
            mobile = water_volumetric  > ground.STATVAR.field_capacity;
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce) > 0.999; %avoid rounding errors
            saturated = [saturated; 0]; %change later, so that first cell of next class is checked -> next class assumed to be hard-bottom now
            hardBottom = (water_volumetric <= lateral.PARA.hardBottom_cutoff) | (ground.STATVAR.Xice > 0);
            hardBottom = [hardBottom; 1];
            
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            fraction_above_out = min(1,max(0,(depths(1:end-1,1) - lateral.PARA.upperElevation) ./ (ground.STATVAR.layerThick - ground.STATVAR.XwaterIce./ground.STATVAR.area) )); %fraction of cell above or below the seepage face
            fraction_below_out = min(1,max(0, (lateral.PARA.lowerElevation - depths(2:end,1))  ./ (ground.STATVAR.layerThick - ground.STATVAR.XwaterIce./ground.STATVAR.area) ));

            for i=1:size(ground.STATVAR.layerThick,1)
                
                %conditions: water in cell must be mobile AND next cell must be either saturated or hardBottom or threshold is reached
                %head
                if (mobile(i,1) && ~hardBottom(i,1)) && (saturated(i+1,1) || hardBottom(i+1,1) || fraction_below_out(i,1) > 0)
                    saturated_height = (water_volumetric(i,1) - ground.STATVAR.field_capacity(i,1)) ./ (porosity(i,1) - ground.STATVAR.field_capacity(i,1)) .* ground.STATVAR.layerThick(i,1);
                    saturated_height = max(0, saturated_height - fraction_below_out(i,1) .* ground.STATVAR.layerThick(i,1) ); %seepage face threshold below
                    lateral.TEMP.head = lateral.TEMP.head + saturated_height;                 
                end
                if hardBottom(i,1) %set head to zero if hard layer is reached
                    lateral.TEMP.head = 0;
                end

                %flow
                if lateral.TEMP.head >0  %avoid unnecessary computation 
                    saturated_height = min(saturated_height , (1-fraction_above_out(i,1)) .* ground.STATVAR.layerThick(i,1)); %seepage face threshold above
                    cross_section = lateral.PARA.seepage_contact_length .* saturated_height;
                    
                    %replace by better hydraulic conductivity formulation later
                    flux = ground.STATVAR.hydraulicConductivity(i,1) .* (lateral.TEMP.head ./ lateral.PARA.distance_seepageFace .* cross_section);
                    water_removed = min(max(0, ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* (ground.STATVAR.layerThick(i,1) .* ground.STATVAR.area(i,1) - ground.STATVAR.XwaterIce(i,1))),...
                        flux .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec);
                    ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) - water_removed;
                    ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) - water_removed .* ground.STATVAR.T(i,1).* ...
                        (ground.CONST.c_w .* double(ground.STATVAR.T(i,1)>=0) + ground.CONST.c_i .* double(ground.STATVAR.T(i,1)<0)) ;
                    ground.STATVAR.water(i,1) = ground.STATVAR.water(i,1) - water_removed;
                    
                    lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off + water_removed;
                end
            end

        end
        
        function ground = lateral_push_remove_water_seepage_snow(ground, lateral)
            water_volumetric = ground.STATVAR.water ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic + ground.STATVAR.ice)./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            mobile = water_volumetric  > ground.PARA.field_capacity .* porosity;
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area) > 0.999; %avoid rounding errors
            hardBottom = (water_volumetric <= lateral.PARA.hardBottom_cutoff) & saturated;
            
            saturated = [saturated; 0]; %change later, so that first cell of next class is checked -> next class assumed to be hard-bottom now
            hardBottom = [hardBottom; 1];
            
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            fraction_above_out = min(1,max(0,(depths(1:end-1,1) - lateral.PARA.upperElevation) ./ ground.STATVAR.layerThick)); %fraction of cell above or below the seepage face
            fraction_below_out = min(1,max(0, (lateral.PARA.lowerElevation - depths(2:end,1))  ./ ground.STATVAR.layerThick));

            for i=1:size(ground.STATVAR.layerThick,1)
                
                %conditions: water in cell must be mobile AND next cell must be either saturated or hardBottom or threshold is reached
                %head
                if (mobile(i,1) && ~hardBottom(i,1)) && (saturated(i+1,1) || hardBottom(i+1,1) || fraction_below_out(i,1) > 0)
                    saturated_height = (water_volumetric(i,1) - ground.PARA.field_capacity.*porosity(i,1)) ./ (porosity(i,1) - ground.PARA.field_capacity.*porosity(i,1)) .* ground.STATVAR.layerThick(i,1);
                    saturated_height = max(0, saturated_height - fraction_below_out(i,1) .* ground.STATVAR.layerThick(i,1)); %seepage face threshold below
                    lateral.TEMP.head = lateral.TEMP.head + saturated_height;                 
                end
                if hardBottom(i,1) %set head to zero if hard layer is reached
                    lateral.TEMP.head = 0;
                end

                %flow
                if lateral.TEMP.head >0  %avoid unnecessary computation 
                    saturated_height = min(saturated_height , (1-fraction_above_out(i,1)) .* ground.STATVAR.layerThick(i,1)); %seepage face threshold above
                    cross_section = lateral.PARA.seepage_contact_length .* saturated_height;
                    
                    
                    flux = ground.STATVAR.hydraulicConductivity(i,1) .* (lateral.TEMP.head ./ lateral.PARA.distance_seepageFace .* cross_section);
                    water_removed = min(max(0,ground.STATVAR.water(i,1) - ground.PARA.field_capacity .* porosity(i,1) .* ground.STATVAR.layerThick(i,1) .* ground.STATVAR.area(i,1)), flux .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec);
                    ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) - water_removed;
                    ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) - water_removed .* ground.STATVAR.T(i,1).* ground.CONST.c_w;
                    ground.STATVAR.water(i,1) = ground.STATVAR.water(i,1) - water_removed;
                    
                    lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off + water_removed;
                end
            end
        end
        
        %---------------------------
        %water reservoir
        %-------------------------
        %remove water from each cell
        function ground = lateral_push_water_reservoir_simple(ground, lateral)
            water_volumetric = ground.STATVAR.water ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            %mobile = water_volumetric  > ground.STATVAR.field_capacity;
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area) > 0.999; %avoid rounding errors
            saturated = [saturated; 0]; %change later, so that first cell of next class is checked -> next class assumed to be hard-bottom now
            hardBottom = water_volumetric <= lateral.PARA.hardBottom_cutoff;
            hardBottom = [hardBottom; 1];
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            
            %go down from top
            %1. identify water table height
            %2. if external water table is higher, identify the space left
            %above to fill up 
            %3. calculate the cross sectional area 
            %4. calculate flux, 
            %5. fill up or take away from the water table height.
            
            pore_space_below_reservoir = 0; %m3 - water than can be added 

            %go down until hitting the water table or a hard bottom
            i = 1;
            while i <= size(ground.STATVAR.layerThick,1)
                if ~hardBottom(i,1)
                    fraction_below = min(1, max(0, (lateral.PARA.reservoir_elevation - depths(i+1,1)) ./ ground.STATVAR.layerThick(i,1)));
                    target_water = fraction_below .* (ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.mineral(i,1) - ground.STATVAR.organic(i,1)) + ...
                        (1-fraction_below) .* (ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.mineral(i,1) - ground.STATVAR.organic(i,1)) .* ground.STATVAR.field_capacity(i,1);
                    %water if the cell was filled up right to the reservoir elevation
                    pore_space_below_reservoir = pore_space_below_reservoir + double(fraction_below>0) .* max(0, target_water - ground.STATVAR.waterIce(i,1));
                    
                    if (saturated(i+1,1) || hardBottom(i+1,1))  %water table identified
                        water_table_cell = i;
                        height_saturated_zone = max(0, (water_volumetric(i,1) - ground.STATVAR.field_capacity(i,1)) ./ ...
                            (porosity(i,1) - ground.STATVAR.field_capacity(i,1)) .* ground.STATVAR.layerThick(i,1));
                        water_table_elevation = depths(i+1,1) + height_saturated_zone;  %absolute elevation of the water table

                        delta_head = water_table_elevation - lateral.PARA.reservoir_elevation; %positive flux = outflow
                        
                        cross_section = lateral.PARA.reservoir_contact_length .* max(height_saturated_zone, 0.05);
                        flux = ground.STATVAR.hydraulicConductivity(i,1) .* (delta_head ./ lateral.PARA.distance_reservoir .* cross_section);
                        mobile_water_saturated_zone = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1)));
                        
                        i=i+1;
                        while (saturated(i,1) && ~hardBottom(i,1)) % && (saturated(i+1,1) || hardBottom(i+1,1))
                            cross_section = lateral.PARA.reservoir_contact_length .* ground.STATVAR.layerThick(i,1);
                            flux = flux + ground.STATVAR.hydraulicConductivity(i,1) .* (delta_head ./ lateral.PARA.distance_reservoir .* cross_section);
                            mobile_water_saturated_zone = mobile_water_saturated_zone + max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1)));

                            i=i+1;
                        end
                        i=i-1; %counter at the last cell of the saturated zone
                        
                        %add/subtract flux bottom-up or top-down -> for
                        %Xice make this cell-wise
                        flux = flux .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec;
                        flux = max(min(flux, mobile_water_saturated_zone), -pore_space_below_reservoir); %limit the fluxes to accounto for water or pore space limitation
                        
                        lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off + flux;
                        
                        j= water_table_cell;
                        if flux > 0
                            while flux>0
                                water_removed = min(max(0,ground.STATVAR.water(j,1) - ground.STATVAR.field_capacity(j,1) .* ground.STATVAR.layerThick(j,1) .* ground.STATVAR.area(j,1)), ...
                                flux);
                                flux = flux - water_removed;
                                                    
                                ground.STATVAR.waterIce(j,1) = ground.STATVAR.waterIce(j,1) - water_removed;
                                ground.STATVAR.energy(j,1) = ground.STATVAR.energy(j,1) - water_removed .* ground.STATVAR.T(j,1).* ...
                                    (ground.CONST.c_w .* double(ground.STATVAR.T(j,1)>=0) + ground.CONST.c_i .* double(ground.STATVAR.T(j,1)<0)) ;
                                ground.STATVAR.water(j,1) = ground.STATVAR.water(j,1) - water_removed;

                                j=j+1;
                            end
                        elseif flux <0
                            while flux<0
                                water_added = min(max(0,ground.STATVAR.layerThick(j,1) .* ground.STATVAR.area(j,1)) - ground.STATVAR.waterIce(j,1) - ground.STATVAR.mineral(j,1) - ground.STATVAR.organic(j,1), ...
                                    -flux);
                                flux = flux + water_added; %flux negative
                                ground.STATVAR.waterIce(j,1) = ground.STATVAR.waterIce(j,1) + water_added;
                                ground.STATVAR.energy(j,1) = ground.STATVAR.energy(j,1) + water_added .* ground.STATVAR.T(j,1).* ...
                                    (ground.CONST.c_w .* double(ground.STATVAR.T(j,1)>=0) + ground.CONST.c_i .* double(ground.STATVAR.T(j,1)<0)) ;
                                ground.STATVAR.water(j,1) = ground.STATVAR.water(j,1) + water_added;
                                
                                j=j-1;
                            end
                        end
                        pore_space_below_reservoir = 0;
                    end
                else    
                    pore_space_below_reservoir = 0;
                end
                i = i+1;
            end
        end
        
        
        %remove water from each cell
        function ground = lateral_push_water_reservoir_Xice(ground, lateral) %causes problems when columen is unfrozen, check
            water_volumetric = ground.STATVAR.water ./ (ground.STATVAR.layerThick.* ground.STATVAR.area  - ground.STATVAR.XwaterIce);
            porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce) > 0.999 %avoid rounding errors
            saturated = [saturated; 0]; %change later, so that first cell of next class is checked -> next class assumed to be hard-bottom now
            hardBottom = (water_volumetric <= lateral.PARA.hardBottom_cutoff) | (ground.STATVAR.Xice > 0);
            hardBottom = [hardBottom; 1];
            
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            
            pore_space_below_reservoir = 0; %m3 - water than can be added
            %go down until hitting the water table or a hard bottom
            i = 1;
            while i <= size(ground.STATVAR.layerThick,1)
                if ~hardBottom(i,1)
                    fraction_below = min(1, max(0, (lateral.PARA.reservoir_elevation - depths(i+1,1)) ./ ground.STATVAR.layerThick(i,1)));
                    target_water = fraction_below .* (ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.mineral(i,1) - ground.STATVAR.organic(i,1) - ground.STATVAR.XwaterIce(i,1)) + ...
                        (1-fraction_below) .* (ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.mineral(i,1) - ground.STATVAR.organic(i,1) - ground.STATVAR.XwaterIce(i,1)) .* ground.STATVAR.field_capacity(i,1);
                    %water if the cell was filled up right to the reservoir elevation
                    pore_space_below_reservoir = pore_space_below_reservoir + double(fraction_below>0) .* max(0, target_water - ground.STATVAR.waterIce(i,1));
                    
                    if (saturated(i+1,1) || hardBottom(i+1,1))  %water table identified
                        flux=[];
                        water_table_cell = i;
                        height_saturated_zone = max(0, (water_volumetric(i,1) - ground.STATVAR.field_capacity(i,1)) ./ ...
                            (porosity(i,1) - ground.STATVAR.field_capacity(i,1)) .* (ground.STATVAR.layerThick(i,1) - ground.STATVAR.XwaterIce(i,1)./ground.STATVAR.area(i,1)));
                        water_table_elevation = depths(i+1,1) + height_saturated_zone + ground.STATVAR.XwaterIce(i,1);  %absolute elevation of the water table
                        delta_head = water_table_elevation - lateral.PARA.reservoir_elevation; %positive flux = outflow
                        
                        cross_section = lateral.PARA.reservoir_contact_length .* max(height_saturated_zone, 0.05); %5cm minimum contact height, otherwise no inflow if completely dry
                        maximum_flux = ground.STATVAR.hydraulicConductivity(i,1) .* (delta_head ./ lateral.PARA.distance_reservoir .* cross_section) .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec;
                        mobile_water = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* (ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.XwaterIce(i,1))));
                        flux = [flux; min(maximum_flux, mobile_water)]; %full flux stored in case of inflow!
                        i=i+1;
                        %go down to bottom of saturated zone
                        while (saturated(i,1) && ~hardBottom(i,1)) % && (saturated(i+1,1) || hardBottom(i+1,1))
                            cross_section = lateral.PARA.reservoir_contact_length .* ground.STATVAR.layerThick(i,1);
                            maximum_flux = ground.STATVAR.hydraulicConductivity(i,1) .* (delta_head ./ lateral.PARA.distance_reservoir .* cross_section) .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec;
                            mobile_water = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* (ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.XwaterIce(i,1))));
                            flux = [flux; min(maximum_flux, mobile_water)]; %full flux stored in case of inflow!
                            i=i+1;
                        end
                        i=i-1; %counter at the last cell of the saturated zone
                        
                        
                        %add/subtract flux cell-wise for Xice
                        if sum(flux,1) > 0 %outflow
                            ground.STATVAR.waterIce(water_table_cell:i,1) = ground.STATVAR.waterIce(water_table_cell:i,1) - flux;
                            ground.STATVAR.energy(water_table_cell:i,1) = ground.STATVAR.energy(water_table_cell:i,1) - flux .* ground.STATVAR.T(water_table_cell:i,1).* ...
                                (ground.CONST.c_w .* double(ground.STATVAR.T(water_table_cell:i,1)>=0) + ground.CONST.c_i .* double(ground.STATVAR.T(water_table_cell:i,1)<0));
                            ground.STATVAR.water(water_table_cell:i,1) = ground.STATVAR.water(water_table_cell:i,1) - flux;
                        elseif sum(flux,1) < 0 %inflow
                            flux_correction_factor = double(~lateral.TEMP.open_system) .* min(1, pore_space_below_reservoir./(-sum(flux,1))) + double(lateral.TEMP.open_system);
                            %no flux corretcion for open system, otherwise correct fluxes so that only the available pore pace gets filled
                            flux =  flux .* flux_correction_factor;
                            pore_space = max(0, (ground.STATVAR.layerThick(water_table_cell:i,1).*ground.STATVAR.area(water_table_cell:i,1) - ground.STATVAR.mineral(water_table_cell:i,1) - ...
                                ground.STATVAR.organic(water_table_cell:i,1) - ground.STATVAR.waterIce(water_table_cell:i,1) - ground.STATVAR.XwaterIce(water_table_cell:i,1)));
                            normal_flux = min(-flux, pore_space);
                            X_flux = -flux - normal_flux;
                            ground.STATVAR.waterIce(water_table_cell:i,1) = ground.STATVAR.waterIce(water_table_cell:i,1) + normal_flux;
                            ground.STATVAR.water(water_table_cell:i,1) = ground.STATVAR.water(water_table_cell:i,1) + normal_flux;
                            ground.STATVAR.XwaterIce(water_table_cell:i,1) = ground.STATVAR.XwaterIce(water_table_cell:i,1) + X_flux;
                            ground.STATVAR.Xwater(water_table_cell:i,1) = ground.STATVAR.Xwater(water_table_cell:i,1) + X_flux;
                            ground.STATVAR.layerThick(water_table_cell:i,1) = ground.STATVAR.layerThick(water_table_cell:i,1) + X_flux ./ ground.STATVAR.area(water_table_cell:i,1);
                            %disp(ground.STATVAR.XwaterIce(1:5,1))
                            % allow for advection of heat
                            if isempty(lateral.PARA.reservoir_temperature)
                                inflow_temperature = ground.STATVAR.T(water_table_cell:i,1);
                            else
                                inflow_temperature = lateral.PARA.reservoir_temperature;
                            end
                            ground.STATVAR.energy(water_table_cell:i,1) = ground.STATVAR.energy(water_table_cell:i,1) - flux .* inflow_temperature .* ...
                                (ground.CONST.c_w .* double(inflow_temperature>=0) + ground.CONST.c_i .* double(inflow_temperature<0));
                        end
                        
                        lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off + sum(flux,1);
                        
                        lateral.TEMP.open_system = 0;
                        pore_space_below_reservoir = 0;
                    end
                else
                    lateral.TEMP.open_system = 0;
                    pore_space_below_reservoir = 0;
                end
                i = i+1;
            end
        end
        
        function ground = lateral_push_water_reservoir_RichardsEq_Xice(ground, lateral)
            depths = ground.STATVAR.upperPos - cumsum(ground.STATVAR.layerThick) + ground.STATVAR.layerThick./2; %midpoints!

            water_volumetric = ground.STATVAR.water ./ (ground.STATVAR.layerThick.* ground.STATVAR.area  - ground.STATVAR.XwaterIce);
            %porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
            pore_space = max(0,ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce - ground.STATVAR.waterIce - ground.STATVAR.mineral - ground.STATVAR.organic);
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce) > 0.999; %avoid rounding errors
             %change later, so that first cell of next class is checked -> next class assumed to be hard-bottom now
            hardBottom = (water_volumetric <= lateral.PARA.hardBottom_cutoff) | (ground.STATVAR.Xice > 0);
            hardBottom = hardBottom;
            Xwater_mobile = ~hardBottom & ground.STATVAR.Xwater > 0;
            
            head =  depths + ground.STATVAR.waterPotential;
            %hydrostatic pressure
            i=1;
            head0=depths(1);
            while i< size(ground.STATVAR.layerThick,1)
                if saturated(i,1) & ~hardBottom(i,1)
                    head(i) = head0;
                else
                   head0 = depths(i+1,1);
                end
                i=i+1;
            end
            if saturated(i,1) %last cell
                head(i) = head0;
            end
            %overburden pressure
            overburden_pressure = ((cumsum(ground.STATVAR.mineral) + ground.STATVAR.mineral./2) .* ground.CONST.rho_m ./ground.CONST.rho_w + ...
                (cumsum(ground.STATVAR.organic) + ground.STATVAR.organic./2) .* ground.CONST.rho_o ./ground.CONST.rho_w + ...
                (cumsum(ground.STATVAR.waterIce) + ground.STATVAR.waterIce./2 + cumsum(ground.STATVAR.XwaterIce) + ground.STATVAR.XwaterIce./2)) ./ ground.STATVAR.area;
            head(Xwater_mobile) = depths(Xwater_mobile) + ground.STATVAR.waterPotential(Xwater_mobile) + overburden_pressure(Xwater_mobile);
            
            cross_section = lateral.PARA.reservoir_contact_length .* ground.STATVAR.layerThick;    
            flux = ground.STATVAR.hydraulicConductivity  .* (head - lateral.PARA.reservoir_elevation) ./ lateral.PARA.distance_reservoir .* cross_section .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec;
            
            %old code no Xice formation laterally
%             flux(hardBottom) = 0;
%             outflow_normal = flux>0 & ~Xwater_mobile; 
%             outflow_Xwater = flux>0 & Xwater_mobile;
%             
%             inflow_normal = flux<0 & ~Xwater_mobile;
%             inflow_Xwater = flux<0 & Xwater_mobile; 
%             
%             flux(outflow_normal) = min(flux(outflow_normal), ground.STATVAR.water(outflow_normal)./2); % maximum half of water in a cell can drain
%             flux(outflow_Xwater) = min(flux(outflow_Xwater), ground.STATVAR.Xwater(outflow_Xwater)); %only Xwater can drain
%             flux(inflow_normal) = max(flux(inflow_normal), -pore_space(inflow_normal)); %fill only available porespace
%                       
%             flux(inflow_Xwater) = 0; %do not add further Xwater
%             
%             ground.STATVAR.waterIce(~outflow_Xwater) = ground.STATVAR.waterIce(~outflow_Xwater) - flux(~outflow_Xwater);
%             ground.STATVAR.water(~outflow_Xwater) = ground.STATVAR.water(~outflow_Xwater) - flux(~outflow_Xwater);
%             
%             ground.STATVAR.XwaterIce(outflow_Xwater) = ground.STATVAR.XwaterIce(outflow_Xwater) - flux(outflow_Xwater);
%             ground.STATVAR.Xwater(outflow_Xwater) = ground.STATVAR.Xwater(outflow_Xwater) - flux(outflow_Xwater);
%             ground.STATVAR.layerThick(outflow_Xwater) = ground.STATVAR.layerThick(outflow_Xwater) - flux(outflow_Xwater) ./ ground.STATVAR.area(outflow_Xwater);
            
            
            %new code Xice formation laterally allowed

            outflow_normal = flux>0 & ~Xwater_mobile; 
            outflow_Xwater = flux>0 & Xwater_mobile;
            
            inflow_normal = flux<0 & ~Xwater_mobile;
            X_water_generation = flux .*0;
            X_water_generation(inflow_normal) = flux(inflow_normal) > pore_space(inflow_normal);
            time2full = flux .*0;
            if sum(X_water_generation)>=1
                time2full(X_water_generation) = pore_space(X_water_generation)./ flux(X_water_generation);
            end
            flux(inflow_normal) = max(flux(inflow_normal), -pore_space(inflow_normal)); %fill available porespace first
            

            flux(outflow_normal) = min(flux(outflow_normal), ground.STATVAR.water(outflow_normal)./2); % maximum half of water in a cell can drain
            flux(outflow_Xwater) = min(flux(outflow_Xwater), ground.STATVAR.Xwater(outflow_Xwater)); %only Xwater can drain
           
            flow_normal = outflow_normal | inflow_normal;
            ground.STATVAR.waterIce(flow_normal) = ground.STATVAR.waterIce(flow_normal) - flux(flow_normal);
            ground.STATVAR.water(flow_normal) = ground.STATVAR.water(flow_normal) - flux(flow_normal);
            
            if sum(X_water_generation)>=1
                head(X_water_generation) = depths(X_water_generation) + ground.STATVAR.waterPotential(X_water_generation) + overburden_pressure(X_water_generation);
                flux(X_water_generation) = ground.STATVAR.hydraulicConductivity(X_water_generation)  .* (head(X_water_generation) - lateral.PARA.reservoir_elevation) ./ lateral.PARA.distance_reservoir .* ...
                    cross_section(X_water_generation) .* lateral.PARENT.IA_TIME_INCREMENT .* (1 - time2full(X_water_generation));
                flux(X_water_generation) = double(flux(X_water_generation)<0) .* flux(X_water_generation) ;
            end
            
            flow_Xwater = flux<0 & Xwater_mobile | X_water_generation | outflow_Xwater; 
                       
            ground.STATVAR.XwaterIce(flow_Xwater) = ground.STATVAR.XwaterIce(flow_Xwater) - flux(flow_Xwater);
            ground.STATVAR.Xwater(flow_Xwater) = ground.STATVAR.Xwater(flow_Xwater) - flux(flow_Xwater);
            ground.STATVAR.layerThick(flow_Xwater) = ground.STATVAR.layerThick(flow_Xwater) - flux(flow_Xwater) ./ ground.STATVAR.area(flow_Xwater);
            
            if sum(X_water_generation)>=1
                flux(X_water_generation) = flux(X_water_generation) -pore_space(X_water_generation);
            end
            %old code continueing
            
            if isempty(lateral.PARA.reservoir_temperature)
                inflow_temperature = ground.STATVAR.T;
            else
                inflow_temperature = lateral.PARA.reservoir_temperature;
            end
            ground.STATVAR.energy = ground.STATVAR.energy - flux .* inflow_temperature .* ...
                (ground.CONST.c_w .* double(inflow_temperature>=0) + ground.CONST.c_i .* double(inflow_temperature<0));

            lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off + sum(flux,1);
            lateral.TEMP.open_system = 0;
          
        end
        
        
        
        function ground = lateral_push_water_reservoir_lake_unfrozen(ground, lateral)
            %water removed and added instantaneously
            delta_layerThick = min((ground.STATVAR.upperPos - lateral.PARA.reservoir_elevation), ground.STATVAR.layerThick./2);
            delta_water = delta_layerThick .* ground.STATVAR.area(1,1);
            
            
            ground.STATVAR.waterIce(1,1) = ground.STATVAR.waterIce(1,1) - delta_water;
            ground.STATVAR.water(1,1) = ground.STATVAR.water(1,1) - delta_water;
            ground.STATVAR.layerThick(1,1) = ground.STATVAR.layerThick(1,1) - delta_layerThick;
            
            % allow for advection of heat
            if isempty(lateral.PARA.reservoir_temperature) || delta_layerThick>0
                inflow_temperature = ground.STATVAR.T(1,1);
            else
                inflow_temperature = lateral.PARA.reservoir_temperature;
            end
            ground.STATVAR.energy(1,1) = ground.STATVAR.energy(1,1) - delta_water .* inflow_temperature .* ...
                (ground.CONST.c_w .* double(inflow_temperature>=0) + ground.CONST.c_i .* double(inflow_temperature<0));
            
            lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off + delta_water;
            
            lateral.TEMP.open_system = 0; %closed system below, all water exchange handled through lake
        end
        
        
        %---------3D-coupled-fluxes-------
        %flow between uppermost unconfined aquifers
        function ground = lateral3D_pull_water_unconfined_aquifer_simple(ground, lateral)
            water_volumetric = ground.STATVAR.water ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            waterIce_volumetric = ground.STATVAR.waterIce ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area) > 0.999; %avoid rounding errors
            
            if ~strcmp(class(ground.NEXT), 'Bottom') %get the values for the first cell of the next class
                [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground.NEXT, lateral);
            else
                saturated_next = 0; 
                hardBottom_next = 1;
            end
            saturated = [saturated; saturated_next];
            hardBottom = water_volumetric <= lateral.PARA.hardBottom_cutoff;
            hardBottom = [hardBottom; hardBottom_next];
            %add one cell based on NEXT 
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            %lateral.PARENT.STATVAR.water_available = 0;
            
            mobile_water = depths .* 0; %m3 - water than can be added 
            hydraulicConductivity = depths .* 0;
            T_water = depths .* 0;

            %go down until hitting the water table or a hard bottom
            i = 1;
            j=0;
            while i <= size(ground.STATVAR.layerThick,1) && lateral.TEMP.open_system
                if ~hardBottom(i,1)
                    
                    mobile_water(i,1) =  -max(0, ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.mineral(i,1) - ground.STATVAR.organic(i,1) - ground.STATVAR.waterIce(i,1));
                    
                    if (saturated(i+1,1) || hardBottom(i+1,1))  %water table identified
                                      
                        water_table_top_cell = i;
%                         height_saturated_zone = max(0, (water_volumetric(i,1) - ground.STATVAR.field_capacity(i,1)) ./ ...
%                             (porosity(i,1) - ground.STATVAR.field_capacity(i,1)) .* ground.STATVAR.layerThick(i,1));
                        height_saturated_zone = max(0, (waterIce_volumetric(i,1) - ground.STATVAR.field_capacity(i,1)) ./ ...
                            (porosity(i,1) - ground.STATVAR.field_capacity(i,1)) .* ground.STATVAR.layerThick(i,1));
                        
                        water_table_elevation = depths(i+1,1) + height_saturated_zone;  %absolute elevation of the water table
                        
                        if i==1 && saturated(i,1) % && height_saturated_zone >= ground.STATVAR.layerThick(i,1) %first cell saturated
                            hydraulicConductivity(i,1) = ground.STATVAR.hydraulicConductivity(i,1);
                            T_water(i,1) = ground.STATVAR.T(i,1);
                             mobile_water(i,1) = max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1));
                            
                        end
                        if height_saturated_zone > 0 && ~saturated(i,1)% && height_saturated_zone < ground.STATVAR.layerThick(i,1) %water level in between grid cells
                             depths = [depths(1:i,1); water_table_elevation; depths(i+1:end,1)];
                             hydraulicConductivity(i+1,1) = ground.STATVAR.hydraulicConductivity(i,1);
                             T_water(i+1,1) = ground.STATVAR.T(i,1);
                             mobile_water(i+1,1) = max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1));
                             j = 1;
                             lateral.PARENT.STATVAR.water_table_top_cell =  water_table_top_cell;
                        end
                        
                        i=i+1;
                        
                        while (saturated(i,1) && ~hardBottom(i,1)) &&  i <= size(ground.STATVAR.layerThick,1)
                            mobile_water(i+j,1) = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1)));
                            hydraulicConductivity(i+j,1) = ground.STATVAR.hydraulicConductivity(i,1);
                            T_water(i+j,1) = ground.STATVAR.T(i,1);
                            i=i+1;
                        end
                        i=i-1; %counter at the last cell of the saturated zone
                        water_table_bottom_cell = i;
                        lateral.TEMP.open_system = 0;
                        if i == size(ground.STATVAR.layerThick,1) && ~hardBottom(i+1,1)
                            lateral.TEMP.open_system = 1;
                        end
                    else
                        water_table_bottom_cell = i;
                        water_table_elevation = [];
                    end
                else
                    lateral.TEMP.open_system = 0;
                    if i==1   %first cell hardBottom
                        j=0;
                        water_table_bottom_cell=0;
                        mobile_water = 0;
                        hydraulicConductivity=0;
                        T_water = 0;
                        water_table_elevation = depths(1,1);
                    end
                end
                i = i+1;
            end
            if isempty(lateral.PARENT.STATVAR.depths)
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths(1:water_table_bottom_cell+j+1,1)];
            else
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths(2:water_table_bottom_cell+j+1,1)];
            end
            %water_table_elevation
            lateral.PARENT.STATVAR.water_status = [lateral.PARENT.STATVAR.water_status; mobile_water(1:water_table_bottom_cell+j,1)];
            lateral.PARENT.STATVAR.hydraulicConductivity = [lateral.PARENT.STATVAR.hydraulicConductivity;  hydraulicConductivity(1:water_table_bottom_cell+j,1)];
            if isempty(lateral.PARENT.STATVAR.water_table_elevation)
                lateral.PARENT.STATVAR.water_table_elevation = water_table_elevation;
            else
                lateral.PARENT.STATVAR.water_table_elevation = max(water_table_elevation, lateral.PARENT.STATVAR.water_table_elevation);
            end
            lateral.PARENT.STATVAR.water_available = lateral.PARENT.STATVAR.water_available || sum(double(mobile_water>0) .* mobile_water,1) > 0;
            lateral.PARENT.STATVAR.T_water = [lateral.PARENT.STATVAR.T_water; T_water(1:water_table_bottom_cell+j,1)];
        end
        
                
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_simple(ground, lateral)
            water_volumetric = ground.STATVAR.water(1,1) ./ ground.STATVAR.layerThick(1,1) ./ ground.STATVAR.area(1,1);
            saturated_next = (ground.STATVAR.waterIce(1,1) + ground.STATVAR.mineral(1,1) + ground.STATVAR.organic(1,1)) ./ (ground.STATVAR.layerThick(1,1) .* ground.STATVAR.area(1,1)) > 0.999; %avoid rounding errors;
            hardBottom_next = water_volumetric <= lateral.PARA.hardBottom_cutoff;
        end
        
        
        function ground = lateral3D_push_water_unconfined_aquifer_simple(ground, lateral)
            %move bottom up and allocate excess water to next cell 
            %if ~isempty(lateral.PARENT.STATVAR.water_flux)
            bottom_cell = size(lateral.PARENT.STATVAR.water_flux,1);
            CURRENT = ground.PREVIOUS;
            while ~strcmp(class(CURRENT), 'Top')
                bottom_cell = bottom_cell - size(CURRENT.STATVAR.energy,1);
                CURRENT = CURRENT.PREVIOUS;
            end
            for i = bottom_cell:-1:1
                    ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) + lateral.PARENT.STATVAR.water_flux(end,1) + lateral.PARENT.STATVAR.water_up;
                    lateral.PARENT.STATVAR.water_up = max(0, ground.STATVAR.waterIce(i,1) + ground.STATVAR.mineral(i,1) + ground.STATVAR.organic(i,1) - ground.STATVAR.layerThick(i,1) .* ground.STATVAR.area(i,1));
                    ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) - lateral.PARENT.STATVAR.water_up;
                    
                    ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) + lateral.PARENT.STATVAR.water_flux_energy(end,1) + lateral.PARENT.STATVAR.water_up_energy;
                    lateral.PARENT.STATVAR.water_up_energy = lateral.PARENT.STATVAR.water_up .* ground.CONST.c_w .* ground.STATVAR.T(i,1);
                    ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) - lateral.PARENT.STATVAR.water_up_energy;
                    
                    lateral.PARENT.STATVAR.water_flux(end,:) = [];
                    lateral.PARENT.STATVAR.water_flux_energy(end,:) = [];
            end
        end
        %------------------------------
        
        
        %Xice
        function ground = lateral3D_pull_water_unconfined_aquifer_Xice(ground, lateral)
            water_volumetric = ground.STATVAR.water ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
            waterIce_volumetric = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
            porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce) > 0.999; %avoid rounding errors
            
            if ~strcmp(class(ground.NEXT), 'Bottom') %get the values for the first cell of the next class
                [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground.NEXT, lateral);
            else
                saturated_next = 0; 
                hardBottom_next = 1;
            end
            saturated = [saturated; saturated_next];
            hardBottom = (water_volumetric <= lateral.PARA.hardBottom_cutoff) | ground.STATVAR.Xice > 0;  %main difference to simple: layer also considered hard when it has Xice
            hardBottom = [hardBottom; hardBottom_next];
            %add one cell based on NEXT 
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            %lateral.PARENT.STATVAR.water_available = 0;
            
            mobile_water = depths .* 0; %m3 - water than can be added 
            hydraulicConductivity = depths .* 0;
            T_water = depths .* 0;

            %go down until hitting the water table or a hard bottom
            i = 1;
            j=0;
            while i <= size(ground.STATVAR.layerThick,1) && lateral.TEMP.open_system
                
                if ~hardBottom(i,1)
                    
                    mobile_water(i,1) =  -max(0, ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.XwaterIce(i,1) - ground.STATVAR.mineral(i,1) - ground.STATVAR.organic(i,1) - ground.STATVAR.waterIce(i,1));
                    
                    if (saturated(i+1,1) || hardBottom(i+1,1))  %water table identified
                                      
                        water_table_top_cell = i;

                        height_saturated_zone = max(0, (waterIce_volumetric(i,1) - ground.STATVAR.field_capacity(i,1)) ./ ...
                            (porosity(i,1) - ground.STATVAR.field_capacity(i,1)) .* (ground.STATVAR.layerThick(i,1) - ground.STATVAR.XwaterIce(i,1)./ground.STATVAR.area(i,1)));
                        
                        water_table_elevation = depths(i+1,1) + height_saturated_zone;  %absolute elevation of the water table
                        
                        if i==1 && saturated(i,1) % && height_saturated_zone >= ground.STATVAR.layerThick(i,1) %first cell saturated
                            hydraulicConductivity(i,1) = ground.STATVAR.hydraulicConductivity(i,1);
                            T_water(i,1) = ground.STATVAR.T(i,1);
                             mobile_water(i,1) = max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1));
                            
                        end
                        if height_saturated_zone > 0 && ~saturated(i,1)% && height_saturated_zone < ground.STATVAR.layerThick(i,1) %water level in between grid cells
                             depths = [depths(1:i,1); water_table_elevation; depths(i+1:end,1)];
                             hydraulicConductivity(i+1,1) = ground.STATVAR.hydraulicConductivity(i,1);
                             T_water(i+1,1) = ground.STATVAR.T(i,1);
                             mobile_water(i+1,1) = max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1));
                             j = 1;
                             lateral.PARENT.STATVAR.water_table_top_cell =  water_table_top_cell;
                        end
                        
                        i=i+1;
                        
                        while (saturated(i,1) && ~hardBottom(i,1)) &&  i <= size(ground.STATVAR.layerThick,1)
                            mobile_water(i+j,1) = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* (ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.XwaterIce(i,1))));
                            hydraulicConductivity(i+j,1) = ground.STATVAR.hydraulicConductivity(i,1);
                            T_water(i+j,1) = ground.STATVAR.T(i,1);
                            i=i+1;
                        end
                        i=i-1; %counter at the last cell of the saturated zone
                        water_table_bottom_cell = i;
                        lateral.TEMP.open_system = 0;
                        if i == size(ground.STATVAR.layerThick,1) && ~hardBottom(i+1,1)
                            lateral.TEMP.open_system = 1;
                        end
                    else
                        water_table_bottom_cell = i;
                        water_table_elevation = [];
                    end
                else
                    lateral.TEMP.open_system = 0;
                    if i==1   %first cell hardBottom
                        j=0;
                        water_table_bottom_cell=0;
                        mobile_water = 0;
                        hydraulicConductivity=0;
                        T_water = 0;
                        water_table_elevation = depths(1,1);
                    end
                end
                i = i+1;
            end
            
            if isempty(lateral.PARENT.STATVAR.depths)
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths(1:water_table_bottom_cell+j+1,1)];
            else
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths(2:water_table_bottom_cell+j+1,1)];
            end
            %water_table_elevation
            lateral.PARENT.STATVAR.water_status = [lateral.PARENT.STATVAR.water_status; mobile_water(1:water_table_bottom_cell+j,1)];
            lateral.PARENT.STATVAR.hydraulicConductivity = [lateral.PARENT.STATVAR.hydraulicConductivity;  hydraulicConductivity(1:water_table_bottom_cell+j,1)];
            if isempty(lateral.PARENT.STATVAR.water_table_elevation)
                lateral.PARENT.STATVAR.water_table_elevation = water_table_elevation;
            else
                lateral.PARENT.STATVAR.water_table_elevation = max(water_table_elevation, lateral.PARENT.STATVAR.water_table_elevation);
            end
            lateral.PARENT.STATVAR.water_available = lateral.PARENT.STATVAR.water_available || sum(double(mobile_water>0) .* mobile_water,1) > 0;
            lateral.PARENT.STATVAR.T_water = [lateral.PARENT.STATVAR.T_water; T_water(1:water_table_bottom_cell+j,1)];
        end
        
                
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_Xice(ground, lateral)
            water_volumetric = ground.STATVAR.water(1,1) ./ (ground.STATVAR.layerThick(1,1) .* ground.STATVAR.area(1,1) - ground.STATVAR.XwaterIce(1,1));
            saturated_next = (ground.STATVAR.waterIce(1,1) + ground.STATVAR.mineral(1,1) + ground.STATVAR.organic(1,1)) ./ (ground.STATVAR.layerThick(1,1) .* ground.STATVAR.area(1,1) - ground.STATVAR.XwaterIce(1,1)) > 0.999; %avoid rounding errors
            hardBottom_next = (water_volumetric <= lateral.PARA.hardBottom_cutoff) | ground.STATVAR.Xice(1,1) > 0;
        end
        
        
        function ground = lateral3D_push_water_unconfined_aquifer_Xice(ground, lateral)
            %move bottom up and allocate excess water to next cell 
            %if ~isempty(lateral.PARENT.STATVAR.water_flux)
            bottom_cell = size(lateral.PARENT.STATVAR.water_flux,1);
            CURRENT = ground.PREVIOUS;
            while ~strcmp(class(CURRENT), 'Top')
                bottom_cell = bottom_cell - size(CURRENT.STATVAR.energy,1);
                CURRENT = CURRENT.PREVIOUS;
            end
            for i = bottom_cell:-1:1
                ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) + lateral.PARENT.STATVAR.water_flux(end,1) + lateral.PARENT.STATVAR.water_up;
                lateral.PARENT.STATVAR.water_up = max(0, ground.STATVAR.waterIce(i,1) + ground.STATVAR.XwaterIce(i,1) + ground.STATVAR.mineral(i,1) + ground.STATVAR.organic(i,1) - ground.STATVAR.layerThick(i,1) .* ground.STATVAR.area(i,1));
                ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) - lateral.PARENT.STATVAR.water_up;
                
                ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) + lateral.PARENT.STATVAR.water_flux_energy(end,1) + lateral.PARENT.STATVAR.water_up_energy;
                lateral.PARENT.STATVAR.water_up_energy = lateral.PARENT.STATVAR.water_up .* ground.CONST.c_w .* ground.STATVAR.T(i,1);
                ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) - lateral.PARENT.STATVAR.water_up_energy;
                
                lateral.PARENT.STATVAR.water_flux(end,:) = [];
                lateral.PARENT.STATVAR.water_flux_energy(end,:) = [];
            end
            
            if strcmp(class(ground.PREVIOUS), 'Top')
               ground.STATVAR.XwaterIce(1,1) = ground.STATVAR.XwaterIce(1,1) + lateral.PARENT.STATVAR.water_up;
               %lateral.PARENT.STATVAR.water_up = 0;
               ground.STATVAR.layerThick(1,1) = ground.STATVAR.layerThick(1,1) + lateral.PARENT.STATVAR.water_up ./ ground.STATVAR.area(1,1);
               lateral.PARENT.STATVAR.water_up = 0;
               ground.STATVAR.energy(1,1) = ground.STATVAR.energy(1,1) + lateral.PARENT.STATVAR.water_up_energy;
               lateral.PARENT.STATVAR.water_up_energy = 0;
            end
        end
        
        
        
        
        %------------------------------
        %snow
         function snow = lateral3D_pull_water_unconfined_aquifer_snow(snow, lateral)
            water_volumetric = snow.STATVAR.water ./ snow.STATVAR.layerThick ./ snow.STATVAR.area;
            porosity = 1 - (snow.STATVAR.ice)./ (snow.STATVAR.layerThick .* snow.STATVAR.area);
            saturated = (snow.STATVAR.waterIce) ./ (snow.STATVAR.layerThick .* snow.STATVAR.area) > 0.999; %avoid rounding errors
            
            if ~strcmp(class(snow.NEXT), 'Bottom') %get the values for the first cell of the next class
                [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(snow.NEXT, lateral);
            else
                saturated_next = 0; 
                hardBottom_next = 1;
            end
            saturated = [saturated; saturated_next];
            hardBottom = porosity <= lateral.PARA.hardBottom_cutoff; %DIFFERENT FOR SNOW!!!
            hardBottom = [hardBottom; hardBottom_next];
            %add one cell based on NEXT 
            depths = snow.STATVAR.upperPos - cumsum([0; snow.STATVAR.layerThick]);
                        
            mobile_water = depths .* 0; %m3 - water than can be added 
            hydraulicConductivity = depths .* 0;
            T_water = depths .* 0;

            %go down until hitting the water table or a hard bottom
            i = 1;
            j=0;
            while i <= size(snow.STATVAR.layerThick,1) && lateral.TEMP.open_system
                if ~hardBottom(i,1)

                    mobile_water(i,1) =  -max(0, snow.STATVAR.layerThick(i,1).*snow.STATVAR.area(i,1) - snow.STATVAR.waterIce(i,1));  %space left in each cell
                    
                    if (saturated(i+1,1) || hardBottom(i+1,1))  %water table identified
                                                
                        water_table_top_cell = i;
%                         height_saturated_zone = max(0, (water_volumetric(i,1) - snow.STATVAR.field_capacity(i,1)) ./ ...
%                             (porosity(i,1) - snow.STATVAR.field_capacity(i,1)) .* snow.STATVAR.layerThick(i,1));
                        height_saturated_zone = max(0, (snow.STATVAR.water(i,1) ./ snow.STATVAR.area(i,1) ./ porosity(i,1) - snow.STATVAR.layerThick(i,1) .* snow.PARA.field_capacity) ./ ...
                            (1 - snow.PARA.field_capacity));

                        water_table_elevation = depths(i+1,1) + height_saturated_zone;  %absolute elevation of the water table
                        if i==1 && saturated(i,1) %first cell saturated
                            hydraulicConductivity(i,1) = snow.STATVAR.hydraulicConductivity(i,1);
                            T_water(i,1) = snow.STATVAR.T(i,1);
%                              mobile_water(i,1) = max(0, (snow.STATVAR.water(i,1) ./ snow.STATVAR.waterIce(i,1) .* porosity(i,1) ...
%                                  - snow.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* snow.STATVAR.area(i,1));
                             mobile_water(i,1) = max(0, porosity(i,1) .* (1 - snow.PARA.field_capacity).* height_saturated_zone .* snow.STATVAR.area(i,1));
                        end
                        if height_saturated_zone > 0 && ~saturated(i,1)    %height_saturated_zone < snow.STATVAR.layerThick(i,1) %water level in between grid cells
                             depths = [depths(1:i,1); water_table_elevation; depths(i+1:end,1)];
                             hydraulicConductivity(i+1,1) = snow.STATVAR.hydraulicConductivity(i,1);
                             T_water(i+1,1) = snow.STATVAR.T(i,1);
                             mobile_water(i+1,1) = max(0, porosity(i,1) .* (1 - snow.PARA.field_capacity) .* height_saturated_zone .* snow.STATVAR.area(i,1));
                             j = 1;
                             lateral.PARENT.STATVAR.water_table_top_cell =  water_table_top_cell;
                        end
                        
                        i=i+1;
                        
                        while (saturated(i,1) && ~hardBottom(i,1)) &&  i <= size(snow.STATVAR.layerThick,1)
                            %mobile_water(i+j,1) = max(0,(snow.STATVAR.water(i,1) - snow.STATVAR.field_capacity(i,1) .* snow.STATVAR.layerThick(i,1).*snow.STATVAR.area(i,1)));
                            mobile_water(i+1,1) = max(0, porosity(i,1) .* (1 - snow.PARA.field_capacity) .*  height_saturated_zone .* snow.STATVAR.area(i,1));
                            hydraulicConductivity(i+j,1) = snow.STATVAR.hydraulicConductivity(i,1);
                            T_water(i+j,1) = snow.STATVAR.T(i,1);
                            i=i+1;
                        end
                        i=i-1; %counter at the last cell of the saturated zone
                        water_table_bottom_cell = i;
                        lateral.TEMP.open_system = 0;
                        if i == size(snow.STATVAR.layerThick,1) && ~hardBottom(i+1,1)
                            lateral.TEMP.open_system = 1;
                        end
                    else
                        water_table_bottom_cell = i;
                        water_table_elevation = [];
                    end
                else
                    lateral.TEMP.open_system = 0;
                    if i==1   %first cell hardBottom
                        j=0;
                        water_table_bottom_cell=0;
                        mobile_water = 0;
                        hydraulicConductivity=0;
                        T_water = 0;
                        water_table_elevation = depths(1,1);
                    end
                end
                i = i+1;
            end
            if isempty(lateral.PARENT.STATVAR.depths)
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths(1:water_table_bottom_cell+j+1,1)];
            else
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths(2:water_table_bottom_cell+j+1,1)];
            end
            %water_table_elevation
            lateral.PARENT.STATVAR.water_status = [lateral.PARENT.STATVAR.water_status; mobile_water(1:water_table_bottom_cell+j,1)];
            lateral.PARENT.STATVAR.hydraulicConductivity = [lateral.PARENT.STATVAR.hydraulicConductivity;  hydraulicConductivity(1:water_table_bottom_cell+j,1)];
            if isempty(lateral.PARENT.STATVAR.water_table_elevation)
                lateral.PARENT.STATVAR.water_table_elevation = water_table_elevation;
            else
                lateral.PARENT.STATVAR.water_table_elevation = max(water_table_elevation, lateral.PARENT.STATVAR.water_table_elevation);
            end
            lateral.PARENT.STATVAR.water_available = lateral.PARENT.STATVAR.water_available || sum(double(mobile_water>0) .* mobile_water,1) > 0;
            lateral.PARENT.STATVAR.T_water = [lateral.PARENT.STATVAR.T_water; T_water(1:water_table_bottom_cell+j,1)];
        end
        
                
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_snow(snow, lateral)
            porosity = 1 - (snow.STATVAR.ice(1,1))./ (snow.STATVAR.layerThick(1,1) .* snow.STATVAR.area(1,1));
            saturated_next = (snow.STATVAR.waterIce(1,1)) ./ (snow.STATVAR.layerThick(1,1) .* snow.STATVAR.area(1,1)) > 0.999; %avoid rounding errors;
            hardBottom = porosity <= lateral.PARA.hardBottom_cutoff; %DIFFERENT FOR SNOW!!!
        end
        
        
        function snow = lateral3D_push_water_unconfined_aquifer_snow(snow, lateral)
            %move bottom up and allocate excess water to next cell 
            %if ~isempty(lateral.PARENT.STATVAR.water_flux)
            bottom_cell = size(lateral.PARENT.STATVAR.water_flux,1);
            CURRENT = snow.PREVIOUS;
            while ~strcmp(class(CURRENT), 'Top')
                bottom_cell = bottom_cell - size(CURRENT.STATVAR.energy,1);
                CURRENT = CURRENT.PREVIOUS;
            end
            
            for i = bottom_cell:-1:1
                snow.STATVAR.waterIce(i,1) = snow.STATVAR.waterIce(i,1) + lateral.PARENT.STATVAR.water_flux(end,1) + lateral.PARENT.STATVAR.water_up;
                lateral.PARENT.STATVAR.water_up = max(0, snow.STATVAR.waterIce(i,1) - snow.STATVAR.layerThick(i,1) .* snow.STATVAR.area(i,1));
                snow.STATVAR.waterIce(i,1) = snow.STATVAR.waterIce(i,1) - lateral.PARENT.STATVAR.water_up;
                
                snow.STATVAR.energy(i,1) = snow.STATVAR.energy(i,1) + lateral.PARENT.STATVAR.water_flux_energy(end,1) + lateral.PARENT.STATVAR.water_up_energy;
                lateral.PARENT.STATVAR.water_up_energy = lateral.PARENT.STATVAR.water_up .* snow.CONST.c_w .* snow.STATVAR.T(i,1);
                snow.STATVAR.energy(i,1) = snow.STATVAR.energy(i,1) - lateral.PARENT.STATVAR.water_up_energy;
                
                lateral.PARENT.STATVAR.water_flux(end,:) = [];
                lateral.PARENT.STATVAR.water_flux_energy(end,:) = [];
            end
        end
        
        function snow = lateral3D_push_water_unconfined_aquifer_snow2(snow, lateral)
            %move bottom up and allocate excess water to next cell 
            %if ~isempty(lateral.PARENT.STATVAR.water_flux)
            bottom_cell = size(lateral.PARENT.STATVAR.water_flux,1);
            CURRENT = snow.PREVIOUS;
            while ~strcmp(class(CURRENT), 'Top')
                bottom_cell = bottom_cell - size(CURRENT.STATVAR.energy,1);
                CURRENT = CURRENT.PREVIOUS;
            end
            
            for i = bottom_cell:-1:1
                snow.STATVAR.waterIce(i,1) = snow.STATVAR.waterIce(i,1) + lateral.PARENT.STATVAR.water_flux(end,1) + lateral.PARENT.STATVAR.water_up;
                lateral.PARENT.STATVAR.water_up = max(0, snow.STATVAR.waterIce(i,1) - snow.STATVAR.layerThick(i,1) .* snow.STATVAR.area(i,1));
                snow.STATVAR.waterIce(i,1) = snow.STATVAR.waterIce(i,1) - lateral.PARENT.STATVAR.water_up;
                
                snow.STATVAR.energy(i,1) = snow.STATVAR.energy(i,1) + lateral.PARENT.STATVAR.water_flux_energy(end,1) + lateral.PARENT.STATVAR.water_up_energy;
                lateral.PARENT.STATVAR.water_up_energy = lateral.PARENT.STATVAR.water_up .* snow.CONST.c_w .* snow.STATVAR.T(i,1);
                snow.STATVAR.energy(i,1) = snow.STATVAR.energy(i,1) - lateral.PARENT.STATVAR.water_up_energy;
                
                lateral.PARENT.STATVAR.water_flux(end,:) = [];
                lateral.PARENT.STATVAR.water_flux_energy(end,:) = [];
            end
            %add remaining water to uppermost cell
            if isempty(lateral.PARENT.STATVAR.water_flux)
              
                snow.STATVAR.waterIce(1,1) = snow.STATVAR.waterIce(1,1) + lateral.PARENT.STATVAR.water_up;
                snow.STATVAR.layerThick(1,1) =  snow.STATVAR.layerThick(1,1) +  lateral.PARENT.STATVAR.water_up ./ snow.STATVAR.area(1,1);
                snow.STATVAR.energy(1,1) = snow.STATVAR.energy(1,1) + lateral.PARENT.STATVAR.water_up_energy;
                lateral.PARENT.STATVAR.water_up = 0;
                lateral.PARENT.STATVAR.water_up_energy = 0;
            end
            
        end
        
        
        %lake unfrozen

        function ground = lateral3D_pull_water_unconfined_aquifer_lake_unfrozen(ground, lateral)
            hydraulic_conductivity_water = 5.*1e-4; %set a higher value than  ground -> maybe problematic when two lakes are communicating! 
            saturated = 1; 
            hardBottom = 0;
            if ~strcmp(class(ground.NEXT), 'Bottom') %get the values for the first cell of the next class
                [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground.NEXT, lateral);
            else
                saturated_next = 0; 
                hardBottom_next = 1;
            end
            saturated = [saturated; saturated_next];
            hardBottom = [hardBottom; hardBottom_next];

            
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);

            lateral.TEMP.open_system = 1;
            lateral.PARENT.STATVAR.water_table_top_cell =  0; %water level is not split between cells!
            
            if isempty(lateral.PARENT.STATVAR.depths)
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths];
            else
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths(2,1)];
            end
            lateral.PARENT.STATVAR.water_status = [lateral.PARENT.STATVAR.water_status; ground.STATVAR.waterIce];
            lateral.PARENT.STATVAR.hydraulicConductivity = [lateral.PARENT.STATVAR.hydraulicConductivity;  hydraulic_conductivity_water];
            if isempty(lateral.PARENT.STATVAR.water_table_elevation)
                lateral.PARENT.STATVAR.water_table_elevation = ground.STATVAR.upperPos;
            else
                lateral.PARENT.STATVAR.water_table_elevation = max(water_table_elevation, ground.STATVAR.upperPos);
            end
            lateral.PARENT.STATVAR.water_available = 1;
            lateral.PARENT.STATVAR.T_water = [lateral.PARENT.STATVAR.T_water; ground.STATVAR.T];
        end
        
                
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell_lake_unfrozen(ground, lateral)
            saturated_next = 1;
            hardBottom_next = 0;
        end
        
        
        function ground = lateral3D_push_water_unconfined_aquifer_lake_unfrozen(ground, lateral)
            %move bottom up and allocate excess water to next cell 
            %if ~isempty(lateral.PARENT.STATVAR.water_flux)
            bottom_cell = size(lateral.PARENT.STATVAR.water_flux,1);
            CURRENT = ground.PREVIOUS;
            while ~strcmp(class(CURRENT), 'Top')
                bottom_cell = bottom_cell - size(CURRENT.STATVAR.energy,1);
                CURRENT = CURRENT.PREVIOUS;
            end
            for i = bottom_cell:-1:1
                    ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) + lateral.PARENT.STATVAR.water_flux(end,1) + lateral.PARENT.STATVAR.water_up;

                    
                    ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) + lateral.PARENT.STATVAR.water_flux_energy(end,1) + lateral.PARENT.STATVAR.water_up_energy;
                    
                    
                    ground.STATVAR.layerThick(i,1) = ground.STATVAR.layerThick(i,1) + (lateral.PARENT.STATVAR.water_flux(end,1) + lateral.PARENT.STATVAR.water_up) ./ ground.STATVAR.area(i,1);
                                        
                    lateral.PARENT.STATVAR.water_up_energy = 0;
                    lateral.PARENT.STATVAR.water_up = 0; %all water taken up by lake
                    lateral.PARENT.STATVAR.water_flux(end,:) = [];
                    lateral.PARENT.STATVAR.water_flux_energy(end,:) = [];
            end
        end
        
        %---
        %genernal flow between confined and unconfined aquifers
        function ground = lateral3D_pull_water_general_aquifer_simple(ground, lateral)
            water_volumetric = ground.STATVAR.water ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            waterIce_volumetric = ground.STATVAR.waterIce ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area) > 0.999; %avoid rounding errors
            hardBottom = (water_volumetric <= lateral.PARA.hardBottom_cutoff | ground.STATVAR.T < 0) & saturated;
            
            if ~strcmp(class(ground.NEXT), 'Bottom') %get the values for the first cell of the next class
                [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground.NEXT, lateral);
            else
                saturated_next = 0; 
                hardBottom_next = 1;
            end
            saturated = [saturated; saturated_next];
            
            hardBottom = [hardBottom; hardBottom_next];
            %add one cell based on NEXT 
            depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            %lateral.PARENT.STATVAR.water_available = 0;
            aquifer_index = water_volumetric .* 0;
            
            mobile_water = water_volumetric .* 0; %m3 - water than can be added 
            hydraulicConductivity = water_volumetric .* 0;
            T_water = water_volumetric .* 0;
            head = water_volumetric .* 0;
            head_unknown = water_volumetric .* 0;

            %go down until hitting the water table or a hard bottom
            i = 1;
            j=0;
                         %same size as lateral.PARENT.depths - uppermost aquifer gets index 1, then 2, etc. - 0 for hard layers
%             lateral.STATVAR.head = []; %same size as number of aquifers - 1D head 
%             lateral.STATVAR.effective_hydraulic_conductivity =[]; %same size as number of aquifers sum(K_eff/d_eff .* A)
%             lateral.STATVAR.empty_volume_left = []; 
%             lateral.STATVAR.available_water_volume = [];
%             lateral.STATVAR.saturated = [];
            j = 0;
            
            while i <= size(ground.STATVAR.layerThick,1)
                if ~hardBottom(i,1)
                    
                    mobile_water(i+j,1) =  -max(0, ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.mineral(i,1) - ground.STATVAR.organic(i,1) - ground.STATVAR.waterIce(i,1));
                    aquifer_index(i+j,1) = lateral.TEMP.aquifer_index_count;
                    head(i+j,1) = (depths(i+j,1) + depths(i+j+1,1)) ./ 2;
                    head_unknown(i+j,1) = double(lateral.TEMP.open_system == 0) + double(lateral.TEMP.open_system == 0 && lateral.TEMP.head_space_available == 0);
                    hydraulicConductivity(i+j,1) = ground.STATVAR.hydraulicConductivity(i,1);
                    T_water(i+j,1) = ground.STATVAR.T(i,1);
                    
                    if (saturated(i+1,1) || hardBottom(i+1,1))  %water table identified
                                      
                        water_table_top_cell = i;
                        aquifer_info(1,2) = water_table_top_cell;
                        height_saturated_zone = max(0, (waterIce_volumetric(i,1) - ground.STATVAR.field_capacity(i,1)) ./ ...
                            (porosity(i,1) - ground.STATVAR.field_capacity(i,1)) .* ground.STATVAR.layerThick(i,1));
                        
                        water_table_elevation = depths(i+j+1,1) + height_saturated_zone;  %absolute elevation of the water table
                                                
                        if saturated(i,1) %first cell saturated
                            hydraulicConductivity(i+j,1) = ground.STATVAR.hydraulicConductivity(i,1);
                            T_water(i+j,1) = ground.STATVAR.T(i,1);
                             mobile_water(i+j,1) = max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1));
                             head(i+j,1) = water_table_elevation;
                             head_unknown(i+j,1) = double(lateral.TEMP.open_system == 0) + double(lateral.TEMP.open_system == 0 && lateral.TEMP.head_space_available == 0);
                            
                        end
                        if height_saturated_zone > 0 && ~saturated(i,1)% water level in between grid cells
                             depths = [depths(1:i+j,1); water_table_elevation; depths(i+j+1:end,1)];
                             aquifer_index = [aquifer_index(1:i+j,1); lateral.TEMP.aquifer_index_count; aquifer_index(i+j+1:end,1)];
                             hydraulicConductivity = [hydraulicConductivity(1:i+j,1); ground.STATVAR.hydraulicConductivity(i,1); hydraulicConductivity(i+j+1:end,1)];
                             %hydraulicConductivity(i+j+1,1) = ground.STATVAR.hydraulicConductivity(i,1);
                             %T_water(i+j+1,1) = ground.STATVAR.T(i,1);
                             T_water = [T_water(1:i+j,1); ground.STATVAR.T(i,1); T_water(i+j+1:end,1)];
                             mobile_water = [mobile_water(1:i+j,1); max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1)); mobile_water(i+j+1:end,1)];

                             
                             lateral.TEMP.head_space_available = 1;
                             head(i+j,1) = (depths(i+j,1) + depths(i+j+1,1)) ./ 2;
                             head = [head(1:i+j,1); water_table_elevation; head(i+j+1:end,1)];
                             head_unknown(i+j,1) = double(lateral.TEMP.open_system == 0) + double(lateral.TEMP.open_system == 0 && lateral.TEMP.head_space_available == 0);
                             head_unknown = [head_unknown(1:i+j,1); head_unknown(i+j,1) ; head_unknown(i+j+1:end,1)];
                             j = j+1;
                             lateral.PARENT.STATVAR.water_table_top_cell =  water_table_top_cell;

                        end
                        
                        i=i+1;
                        
                        while (saturated(i,1) && ~hardBottom(i,1)) &&  i <= size(ground.STATVAR.layerThick,1)
                            mobile_water(i+j,1) = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1)));
                            hydraulicConductivity(i+j,1) = ground.STATVAR.hydraulicConductivity(i,1);
                            T_water(i+j,1) = ground.STATVAR.T(i,1);
                            aquifer_index(i+j,1) = lateral.TEMP.aquifer_index_count;
                            
                            head(i+j,1) = water_table_elevation;
                            head_unknown(i+j,1) = double(lateral.TEMP.open_system == 0) + double(lateral.TEMP.open_system == 0 && lateral.TEMP.head_space_available == 0);
                            
                            i=i+1;
                            
                        end
                        i=i-1; %counter at the last cell of the saturated zone
                        water_table_bottom_cell = i;
                        lateral.TEMP.open_system = 0;
                        if i == size(ground.STATVAR.layerThick,1) && ~hardBottom(i+1,1)
                            lateral.TEMP.open_system = 1;
                        end
                        lateral.TEMP.aquifer_index_count = lateral.TEMP.aquifer_index_count+1;
                    else
                        water_table_bottom_cell = i;
                        water_table_elevation = [];
                        lateral.TEMP.head_space_available = 1;
                       
                        head(i+j,1) = (depths(i+j,1) + depths(i+j+1,1)) ./ 2;
                        head_unknown(i+j,1) = double(lateral.TEMP.open_system == 0) + double(lateral.TEMP.open_system == 0 && lateral.TEMP.head_space_available == 0);
                        hydraulicConductivity(i+j,1) = ground.STATVAR.hydraulicConductivity(i,1);
                    end
                else
                    lateral.TEMP.open_system = 0;
                    lateral.TEMP.head_space_available = 0;
                    head_unknown(i+j,1) = -1;
                end
                i = i+1;
            end

            
            if isempty(lateral.PARENT.STATVAR.depths)
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths];
            else
                lateral.PARENT.STATVAR.depths = [lateral.PARENT.STATVAR.depths; depths(2:end,1)];
            end
            
            lateral.PARENT.STATVAR.water_status = [lateral.PARENT.STATVAR.water_status; mobile_water];

            lateral.PARENT.STATVAR.hydraulicConductivity = [lateral.PARENT.STATVAR.hydraulicConductivity;  hydraulicConductivity];

            lateral.PARENT.STATVAR.aquifer_index = [lateral.PARENT.STATVAR.aquifer_index; aquifer_index];

            lateral.PARENT.STATVAR.water_available = lateral.PARENT.STATVAR.water_available || sum(double(mobile_water>0) .* mobile_water,1) > 0;

            lateral.PARENT.STATVAR.T_water = [lateral.PARENT.STATVAR.T_water; T_water];
            lateral.PARENT.STATVAR.head = [lateral.PARENT.STATVAR.head; head];
            
            lateral.PARENT.STATVAR.head_unknown = [lateral.PARENT.STATVAR.head_unknown; head_unknown];
        end
        
        
        function ground = lateral3D_push_water_general_aquifer_simple(ground, lateral)
            %move bottom up and allocate excess water to next cell
            %if ~isempty(lateral.PARENT.STATVAR.water_flux)
            %             bottom_cell = size(lateral.PARENT.STATVAR.water_flux,1);
            %             CURRENT = ground.PREVIOUS;
            %             while ~strcmp(class(CURRENT), 'Top')
            %                 bottom_cell = bottom_cell - size(CURRENT.STATVAR.energy,1);
            %                 CURRENT = CURRENT.PREVIOUS;
            %             end
            %             for i = bottom_cell:-1:1
            %                 ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) + lateral.PARENT.STATVAR.water_flux(end,1) + lateral.PARENT.STATVAR.water_up;
            %                 lateral.PARENT.STATVAR.water_up = max(0, ground.STATVAR.waterIce(i,1) + ground.STATVAR.mineral(i,1) + ground.STATVAR.organic(i,1) - ground.STATVAR.layerThick(i,1) .* ground.STATVAR.area(i,1));
            %                 ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) - lateral.PARENT.STATVAR.water_up;
            %
            %                 ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) + lateral.PARENT.STATVAR.water_flux_energy(end,1) + lateral.PARENT.STATVAR.water_up_energy;
            %                 lateral.PARENT.STATVAR.water_up_energy = lateral.PARENT.STATVAR.water_up  .* ground.STATVAR.T(i,1).* ...
            %                     (double(ground.STATVAR.T(i,1)>0) .* ground.CONST.c_w +  double(ground.STATVAR.T(i,1)<0) .* ground.CONST.c_i);
            %                 ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) - lateral.PARENT.STATVAR.water_up_energy;
            %
            %                 lateral.PARENT.STATVAR.water_flux(end,:) = [];
            %                 lateral.PARENT.STATVAR.water_flux_energy(end,:) = [];
            %             end
            %         end
            
            %add the water flux, allowing overflow for individual cells
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + lateral.PARENT.STATVAR.water_flux(end-size(ground.STATVAR.waterIce,1)+1:end,1);
            ground.STATVAR.waterIce(end,1) = ground.STATVAR.waterIce(end,1) + lateral.PARENT.STATVAR.water_up;
            ground.STATVAR.energy = ground.STATVAR.energy + lateral.PARENT.STATVAR.water_flux_energy(end-size(ground.STATVAR.waterIce,1)+1:end,1);
            ground.STATVAR.energy(end,1) = ground.STATVAR.energy(end,1)  + lateral.PARENT.STATVAR.water_up_energy;
            water_max = ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral + ground.STATVAR.organic;
            water_down = 0;
            
            
            pos=size(ground.STATVAR.waterIce,1);
            while pos>0
                if ground.STATVAR.waterIce(pos,1)>=water_max(pos,1)
                    water_up = max(0, ground.STATVAR.waterIce(pos,1) - water_max(pos,1));
                    water_up_energy =  water_up  .* ground.STATVAR.T(pos,1).* (double(ground.STATVAR.T(pos,1)>0) .* ground.CONST.c_w +  double(ground.STATVAR.T(pos,1)<0) .* ground.CONST.c_i);
                    ground.STATVAR.waterIce(pos,1) = ground.STATVAR.waterIce(pos,1) - water_up;
                    ground.STATVAR.energy(pos,1) = ground.STATVAR.energy(pos,1) - water_up_energy;
                    if pos>1
                        ground.STATVAR.waterIce(pos-1,1) = ground.STATVAR.waterIce(pos-1,1) + water_up;
                        ground.STATVAR.energy(pos-1,1) = ground.STATVAR.energy(pos-1,1) + water_up_energy;
                        water_up = 0;
                        water_up_energy = 0;
                    end
                    pos = pos-1;
                elseif ground.STATVAR.waterIce(pos,1) < water_max(pos,1)  %find the next cell with excess water above
                    water_found = 0;
                    pos2=pos-1;
                    while pos2>0 && ~water_found
                        if ground.STATVAR.waterIce(pos2,1) > water_max(pos2,1)
                            water_found = 1;
                        else
                            pos2=pos2-1;
                        end
                    end
                    if water_found  %this can be done better
                        water_down = min( ground.STATVAR.waterIce(pos2,1)-water_max(pos2,1), sum(water_max(pos2+1:pos,1)-ground.STATVAR.waterIce(pos2+1:pos,1),1));
                        water_down_energy = water_down .* ground.STATVAR.T(pos2,1).* (double(ground.STATVAR.T(pos2,1)>0) .* ground.CONST.c_w +  double(ground.STATVAR.T(pos2,1)<0) .* ground.CONST.c_i);
                        ground.STATVAR.waterIce(pos2,1) = ground.STATVAR.waterIce(pos2,1) - water_down;
                        ground.STATVAR.energy(pos2,1) = ground.STATVAR.energy(pos2,1) - water_down_energy;
                        j=pos2+1;
                        while j<=pos
                            ground.STATVAR.waterIce(j,1) = ground.STATVAR.waterIce(j,1) + water_down;
                            ground.STATVAR.energy(j,1) = ground.STATVAR.energy(j,1) + water_down_energy;
                            water_down = max(0, ground.STATVAR.waterIce(j,1) - water_max(j,1));
                            water_down_energy = water_down .* ground.STATVAR.T(j,1).* (double(ground.STATVAR.T(j,1)>0) .* ground.CONST.c_w +  double(ground.STATVAR.T(j,1)<0) .* ground.CONST.c_i);
                            
                            ground.STATVAR.waterIce(j,1) = ground.STATVAR.waterIce(j,1) - water_down;
                            ground.STATVAR.energy(j,1) = ground.STATVAR.energy(j,1) - water_down_energy;
                            j=j+1;
                        end
                    else
                        pos=0; %terminate
                    end
                end
            end
            lateral.PARENT.STATVAR.water_up = water_up;
            lateral.PARENT.STATVAR.water_up_energy = water_up_energy;
            
        end
    end
end

