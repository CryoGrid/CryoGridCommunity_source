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
                    flux = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1)) .* (lateral.TEMP.head ./ lateral.PARA.distance_seepageFace .* cross_section);
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
                    flux = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1)) .* (lateral.TEMP.head ./ lateral.PARA.distance_seepageFace .* cross_section);
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
                    
                    %replace by better hydraulic conductivity formulation later
                    flux = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1)) .* (lateral.TEMP.head ./ lateral.PARA.distance_seepageFace .* cross_section);
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
                        %disp(water_table_elevation - 1222);
                        delta_head = water_table_elevation - lateral.PARA.reservoir_elevation; %positive flux = outflow
                        
                        cross_section = lateral.PARA.reservoir_contact_length .* max(height_saturated_zone, 0.05);
                        flux = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1)) .* (delta_head ./ lateral.PARA.distance_reservoir .* cross_section);
                        mobile_water_saturated_zone = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1)));
                        
                        i=i+1;
                        while (saturated(i,1) && ~hardBottom(i,1)) % && (saturated(i+1,1) || hardBottom(i+1,1))
                            cross_section = lateral.PARA.reservoir_contact_length .* ground.STATVAR.layerThick(i,1);
                            flux = flux + ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1)) .* (delta_head ./ lateral.PARA.distance_reservoir .* cross_section);
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
        function ground = lateral_push_water_reservoir_Xice(ground, lateral)
            water_volumetric = ground.STATVAR.water ./ (ground.STATVAR.layerThick.* ground.STATVAR.area  - ground.STATVAR.XwaterIce);
            porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic)./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce);
            saturated = (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce) > 0.999; %avoid rounding errors
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
                        maximum_flux = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1)) .* (delta_head ./ lateral.PARA.distance_reservoir .* cross_section) .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec;
                        mobile_water = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* (ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.XwaterIce(i,1))));
                        flux = [flux; min(maximum_flux, mobile_water)]; %full flux stored in case of inflow!
                        i=i+1;
                        %go down to bottom of saturated zone
                        while (saturated(i,1) && ~hardBottom(i,1)) % && (saturated(i+1,1) || hardBottom(i+1,1))
                            cross_section = lateral.PARA.reservoir_contact_length .* ground.STATVAR.layerThick(i,1);
                            maximum_flux = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1)) .* (delta_head ./ lateral.PARA.distance_reservoir .* cross_section) .* lateral.PARENT.IA_TIME_INCREMENT .* lateral.CONST.day_sec;
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
                            hydraulicConductivity(i,1) = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1));
                            T_water(i,1) = ground.STATVAR.T(i,1);
                             mobile_water(i,1) = max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1));
                            
                        end
                        if height_saturated_zone > 0 && ~saturated(i,1)% && height_saturated_zone < ground.STATVAR.layerThick(i,1) %water level in between grid cells
                             depths = [depths(1:i,1); water_table_elevation; depths(i+1:end,1)];
                             hydraulicConductivity(i+1,1) = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1));
                             T_water(i+1,1) = ground.STATVAR.T(i,1);
                             mobile_water(i+1,1) = max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1));
                             j = 1;
                             lateral.PARENT.STATVAR.water_table_top_cell =  water_table_top_cell;
                        end
                        
                        i=i+1;
                        
                        while (saturated(i,1) && ~hardBottom(i,1)) &&  i <= size(ground.STATVAR.layerThick,1)
                            mobile_water(i+j,1) = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1)));
                            hydraulicConductivity(i+j,1) = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1));
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
            %end
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
                            (porosity(i,1) - ground.STATVAR.field_capacity(i,1)) .* (ground.STATVAR.layerThick(i,1) - ground.STATVAR.XwaterIce./ground.STATVAR.area(i,1)));
                        
                        water_table_elevation = depths(i+1,1) + height_saturated_zone;  %absolute elevation of the water table
                        
                        if i==1 && saturated(i,1) % && height_saturated_zone >= ground.STATVAR.layerThick(i,1) %first cell saturated
                            hydraulicConductivity(i,1) = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1));
                            T_water(i,1) = ground.STATVAR.T(i,1);
                             mobile_water(i,1) = max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1));
                            
                        end
                        if height_saturated_zone > 0 && ~saturated(i,1)% && height_saturated_zone < ground.STATVAR.layerThick(i,1) %water level in between grid cells
                             depths = [depths(1:i,1); water_table_elevation; depths(i+1:end,1)];
                             hydraulicConductivity(i+1,1) = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1));
                             T_water(i+1,1) = ground.STATVAR.T(i,1);
                             mobile_water(i+1,1) = max(0, (ground.STATVAR.water(i,1) ./ ground.STATVAR.waterIce(i,1) .* porosity(i,1) ...
                                 - ground.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* ground.STATVAR.area(i,1));
                             j = 1;
                             lateral.PARENT.STATVAR.water_table_top_cell =  water_table_top_cell;
                        end
                        
                        i=i+1;
                        
                        while (saturated(i,1) && ~hardBottom(i,1)) &&  i <= size(ground.STATVAR.layerThick,1)
                            mobile_water(i+j,1) = max(0,(ground.STATVAR.water(i,1) - ground.STATVAR.field_capacity(i,1) .* (ground.STATVAR.layerThick(i,1).*ground.STATVAR.area(i,1) - ground.STATVAR.XwaterIce)));
                            hydraulicConductivity(i+j,1) = ground.PARA.hydraulicConductivity .* (ground.STATVAR.water(i,1)./ground.STATVAR.waterIce(i,1));
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
                lateral.PARENT.STATVAR.water_up = max(0, ground.STATVAR.waterIce(i,1) + ground.STATVAR.mineral(i,1) + ground.STATVAR.organic(i,1) - ground.STATVAR.layerThick(i,1) .* ground.STATVAR.area(i,1));
                ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) - lateral.PARENT.STATVAR.water_up;
                
                ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) + lateral.PARENT.STATVAR.water_flux_energy(end,1) + lateral.PARENT.STATVAR.water_up_energy;
                lateral.PARENT.STATVAR.water_up_energy = lateral.PARENT.STATVAR.water_up .* ground.CONST.c_w .* ground.STATVAR.T(i,1);
                ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) - lateral.PARENT.STATVAR.water_up_energy;
                
                lateral.PARENT.STATVAR.water_flux(end,:) = [];
                lateral.PARENT.STATVAR.water_flux_energy(end,:) = [];
            end
            
            if strcmp(class(ground.PREVIOUS), 'Top')
               ground.STATVAR.XwaterIce(1,1) = ground.STATVAR.XwaterIce(1,1) + lateral.PARENT.STATVAR.water_up;
               lateral.PARENT.STATVAR.water_up = 0;
               ground.STATVAR.layerThick(1,1) = ground.STATVAR.layerThick(1,1) + lateral.PARENT.STATVAR.water_up ./ ground.STATVAR.area(1,1);
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
                            hydraulicConductivity(i,1) = snow.PARA.hydraulicConductivity .* (snow.STATVAR.water(i,1)./snow.STATVAR.waterIce(i,1));
                            T_water(i,1) = snow.STATVAR.T(i,1);
%                              mobile_water(i,1) = max(0, (snow.STATVAR.water(i,1) ./ snow.STATVAR.waterIce(i,1) .* porosity(i,1) ...
%                                  - snow.STATVAR.field_capacity(i,1) ).* height_saturated_zone .* snow.STATVAR.area(i,1));
                             mobile_water(i,1) = max(0, porosity(i,1) .* (1 - snow.PARA.field_capacity).* height_saturated_zone .* snow.STATVAR.area(i,1));
                        end
                        if height_saturated_zone > 0 && ~saturated(i,1)    %height_saturated_zone < snow.STATVAR.layerThick(i,1) %water level in between grid cells
                             depths = [depths(1:i,1); water_table_elevation; depths(i+1:end,1)];
                             hydraulicConductivity(i+1,1) = snow.PARA.hydraulicConductivity .* (snow.STATVAR.water(i,1)./snow.STATVAR.waterIce(i,1));
                             T_water(i+1,1) = snow.STATVAR.T(i,1);
                             mobile_water(i+1,1) = max(0, porosity(i,1) .* (1 - snow.PARA.field_capacity) .* height_saturated_zone .* snow.STATVAR.area(i,1));
                             j = 1;
                             lateral.PARENT.STATVAR.water_table_top_cell =  water_table_top_cell;
                        end
                        
                        i=i+1;
                        
                        while (saturated(i,1) && ~hardBottom(i,1)) &&  i <= size(snow.STATVAR.layerThick,1)
                            %mobile_water(i+j,1) = max(0,(snow.STATVAR.water(i,1) - snow.STATVAR.field_capacity(i,1) .* snow.STATVAR.layerThick(i,1).*snow.STATVAR.area(i,1)));
                            mobile_water(i+1,1) = max(0, porosity(i,1) .* (1 - snow.PARA.field_capacity) .*  height_saturated_zone .* snow.STATVAR.area(i,1));
                            hydraulicConductivity(i+j,1) = snow.PARA.hydraulicConductivity .* (snow.STATVAR.water(i,1)./snow.STATVAR.waterIce(i,1));
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
            %end
        end
        
        
        %lake unfrozen

        function ground = lateral3D_pull_water_unconfined_aquifer_lake_unfrozen(ground, lateral)
            hydraulic_conductivity_water = 1e-3; %set a large value compared to ground -> can be problematic when two lakes are communicating! 
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
                    lateral.PARENT.STATVAR.water_up = 0; %all water taken up by lake
                    
                    ground.STATVAR.energy(i,1) = ground.STATVAR.energy(i,1) + lateral.PARENT.STATVAR.water_flux_energy(end,1) + lateral.PARENT.STATVAR.water_up_energy;
                    lateral.PARENT.STATVAR.water_up_energy = 0;
                    
                    ground.STATVAR.layerThick(i,1) = ground.STATVAR.layerThick(i,1) + (lateral.PARENT.STATVAR.water_flux(end,1) + lateral.PARENT.STATVAR.water_up) ./ ground.STATVAR.area(i,1);
                    
                    lateral.PARENT.STATVAR.water_flux(end,:) = [];
                    lateral.PARENT.STATVAR.water_flux_energy(end,:) = [];
            end
        end
        
        
    end
end

