%========================================================================
% CryoGrid TIER1 library class, functions related to vertical fluxes of water
% S. Westermann, October 2020
%========================================================================

classdef WATER_FLUXES < BASE

    methods
        
        %normal upper boundary condition for water fluxes for bucketW water scheme (the "2" has no meaning!)
        function ground = get_boundary_condition_u_water2(ground, forcing)
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  %possibly add water from external source here 
            
            %partition already here in infiltration and surface runoff,
            %considering ET losses and potentially external fluxes
            saturation_first_cell = (ground.STATVAR.waterIce(1) - ground.STATVAR.field_capacity(1) .* ground.STATVAR.layerThick(1).* ground.STATVAR.area(1))./ (ground.STATVAR.layerThick(1).*ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.field_capacity(1).*ground.STATVAR.layerThick(1).*ground.STATVAR.area(1));
            saturation_first_cell = max(0,min(1,saturation_first_cell)); % 0 water at field capacity, 1: water at saturation
            saturation_first_cell(saturation_first_cell >= (1 - 1e-9)) = 1;
            
            evap = double(ground.TEMP.d_water_ET(1)<0).*ground.TEMP.d_water_ET(1);
            condensation = double(ground.TEMP.d_water_ET(1)>0).*ground.TEMP.d_water_ET(1);
            
            rainfall = rainfall + condensation; %add condensation to rainfall to avoid overflowing of grid cell
            ground.TEMP.d_water_ET(1) = evap; %evaporation (water loss) subrtacted in get_derivative
            
            ground.TEMP.F_ub_water = double(rainfall <= -evap) .* rainfall + ...
                double(rainfall > -evap) .* (-evap + (rainfall + evap) .* reduction_factor_in(saturation_first_cell, ground));
            ground.TEMP.surface_runoff = rainfall - ground.TEMP.F_ub_water;
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
        
        function ground = get_boundary_condition_u_RichardsEq(ground, forcing)
            
            max_infiltration = max(0, ground.STATVAR.hydraulicConductivity(1,1).* ((0 - ground.STATVAR.waterPotential(1,1)) ./ (ground.STATVAR.layerThick(1,1) ./ 2) + 1) .* ground.STATVAR.area(1,1));
            
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  %possibly add water from external source here 

            %partition already here in infiltration and surface runoff,
            %considering ET losses and potentially external fluxes
            saturation_first_cell = ground.STATVAR.waterIce(1)./ (ground.STATVAR.layerThick(1).*ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1));
            saturation_first_cell = max(0,min(1,saturation_first_cell)); % 0 water at field capacity, 1: water at saturation
            saturation_first_cell(saturation_first_cell >= (1 - 1e-9)) = 1;
            
            evap = double(ground.TEMP.d_water_ET(1)<0).*ground.TEMP.d_water_ET(1);
            condensation = double(ground.TEMP.d_water_ET(1)>0).*ground.TEMP.d_water_ET(1);
            
            rainfall = rainfall + condensation; %add condensation to rainfall to avoid overflowing of grid cell
            excessRain = max(0, rainfall-max_infiltration);
            rainfall = min(rainfall, max_infiltration);
            
            ground.TEMP.d_water_ET(1) = evap; %evaporation (water loss) subrtacted in get_derivative
            
            ground.TEMP.F_ub_water = double(rainfall <= -evap) .* rainfall + ...
                double(rainfall > -evap) .* (-evap + (rainfall + evap) .* reduction_factor_in(saturation_first_cell, ground));
            ground.TEMP.surface_runoff = rainfall - ground.TEMP.F_ub_water + excessRain;
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
        
        function ground = get_boundary_condition_u_RichardsEq_pressure(ground, forcing)
            
            max_infiltration = max(0, ground.STATVAR.hydraulicConductivity(1,1).* ((0 - ground.STATVAR.waterPotential(1,1)) ./ (ground.STATVAR.layerThick(1,1) ./ 2) + 1) .* ground.STATVAR.area(1,1));
            
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  
            
            %partition already here in infiltration and surface runoff,
            %considering ET losses and potentially external fluxes
            saturation_first_cell = ground.STATVAR.waterIce(1)./ (ground.STATVAR.layerThick(1).*ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1));
            saturation_first_cell = max(0,min(1,saturation_first_cell)); % 0 water at field capacity, 1: water at saturation
            
            evap = double(ground.TEMP.d_water_ET(1)<0).*ground.TEMP.d_water_ET(1);
            condensation = double(ground.TEMP.d_water_ET(1)>0).*ground.TEMP.d_water_ET(1);
            
            rainfall = rainfall + condensation; %add condensation to rainfall to avoid overflowing of grid cell
            excessRain = max(0, rainfall-max_infiltration);
            rainfall = min(rainfall, max_infiltration);
            
            ground.TEMP.d_water_ET(1) = evap; %evaporation (water loss) subtracted in get_derivative
            
            ground.TEMP.F_ub_water = rainfall;
            ground.TEMP.surface_runoff = excessRain;
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
        
        %bucketW hydrological scheme
        function ground = get_boundary_condition_u_water_Xice(ground, forcing)  %simply add the water to first grid cell, excess taken up by Xwater, no checks needed
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  %possibly add water from external source here 
            
            %partition already here in infiltration and surface runoff, considering ET losses and potentially external fluxes
            volume_matrix = ground.STATVAR.layerThick(1) .* ground.STATVAR.area(1) - ground.STATVAR.XwaterIce(1);
            saturation_first_cell = (ground.STATVAR.waterIce(1)  - ground.STATVAR.field_capacity(1) .* volume_matrix)./...
                (volume_matrix - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.field_capacity(1) .* volume_matrix);
            saturation_first_cell = max(0,min(1,saturation_first_cell)); % 0 water at field capacity, 1: water at saturation
            %NEW SW
            saturation_first_cell(saturation_first_cell >= (1 - 1e-9)) = 1;
            
            evap = double(ground.TEMP.d_water_ET(1)<0).*ground.TEMP.d_water_ET(1); %negative
            condensation = double(ground.TEMP.d_water_ET(1)>0).*ground.TEMP.d_water_ET(1);
            
            rainfall = rainfall + condensation; %add condensation to rainfall to avoid overflowing of grid cell
            ground.TEMP.d_water_ET(1) = evap; %evaporation (water loss) subrtacted in get_derivative
            
            ground.TEMP.F_ub_water = double(rainfall <= -evap) .* rainfall + ...
                double(rainfall > -evap) .* (-evap + (rainfall + evap) .* reduction_factor_in(saturation_first_cell, ground));
            ground.TEMP.F_ub_Xwater = rainfall - ground.TEMP.F_ub_water;
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            ground.TEMP.F_ub_Xwater_energy = ground.TEMP.F_ub_Xwater .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
            ground.TEMP.d_Xwater(1) = ground.TEMP.d_Xwater(1) + ground.TEMP.F_ub_Xwater;
            ground.TEMP.d_Xwater_energy(1) = ground.TEMP.d_Xwater_energy(1) + ground.TEMP.F_ub_Xwater_energy;
            
        end
        
        %DISCONTINUED
        function ground = get_boundary_condition_u_water_RichardsEq_Xice(ground, forcing)  %simply add the water to first grid cell, excess taken up by Xwater, no checks needed
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  %possibly add water from external source here 
            
            evap = double(ground.TEMP.d_water_ET(1)<0).*ground.TEMP.d_water_ET(1); %negative
            condensation = double(ground.TEMP.d_water_ET(1)>0).*ground.TEMP.d_water_ET(1);
            
            rainfall = rainfall + condensation; %add condensation to rainfall to avoid overflowing of grid cell
            ground.TEMP.d_water_ET(1) = evap; %evaporation (water loss) subrtacted in get_derivative
            
            ground.TEMP.F_ub_water = rainfall;
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
        
                
        function ground = get_boundary_condition_u_water_RichardsEq_Xice2(ground, forcing)  %simply add the water to first grid cell, excess taken up by Xwater, no checks needed
             max_infiltration = max(0, ground.STATVAR.hydraulicConductivity(1,1).* ((0 - ground.STATVAR.waterPotential(1,1)) ./ (ground.STATVAR.layerThick(1,1) ./ 2) + 1) .* ground.STATVAR.area(1,1));
            
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  
            
            %partition already here in infiltration and surface runoff, considering ET losses and potentially external fluxes
%             volume_matrix = ground.STATVAR.layerThick(1) .* ground.STATVAR.area(1) - ground.STATVAR.XwaterIce(1);
%             saturation_first_cell = (ground.STATVAR.waterIce(1)  - ground.STATVAR.field_capacity(1) .* volume_matrix)./...
%                 (volume_matrix - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.field_capacity(1) .* volume_matrix);
            
            saturation_first_cell = ground.STATVAR.waterIce(1)./ (ground.STATVAR.layerThick(1).*ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.XwaterIce(1));

            saturation_first_cell = max(0,min(1,saturation_first_cell)); % 0 water at field capacity, 1: water at saturation
            saturation_first_cell(saturation_first_cell >= (1 - 1e-6 .* rand(1))) = 1;
            
            evap = double(ground.TEMP.d_water_ET(1)<0).*ground.TEMP.d_water_ET(1);
            condensation = double(ground.TEMP.d_water_ET(1)>0).*ground.TEMP.d_water_ET(1);
            
            rainfall = rainfall + condensation; %add condensation to rainfall to avoid overflowing of grid cell
            excessRain = max(0, rainfall-max_infiltration);
            rainfall = min(rainfall, max_infiltration);
            
            ground.TEMP.d_water_ET(1) = evap; %evaporation (water loss) subtracted in get_derivative
            
            ground.TEMP.F_ub_water = double(rainfall <= -evap) .* rainfall + ...
                double(rainfall > -evap) .* (-evap + (rainfall + evap) .* reduction_factor_in(saturation_first_cell, ground));
            ground.TEMP.F_ub_Xwater = rainfall - ground.TEMP.F_ub_water + excessRain;
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            ground.TEMP.F_ub_Xwater_energy = ground.TEMP.F_ub_Xwater .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
            ground.TEMP.d_Xwater(1) = ground.TEMP.d_Xwater(1) + ground.TEMP.F_ub_Xwater;
            ground.TEMP.d_Xwater_energy(1) = ground.TEMP.d_Xwater_energy(1) + ground.TEMP.F_ub_Xwater_energy;

        end
        
        
        %bucektW for snow classes
        function ground = get_boundary_condition_u_water_SNOW(ground, forcing)
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  
            
            %partition already here in infiltration and surface runoff,
            %considering ET losses and potentially external fluxes
            remaining_pore_space = ground.STATVAR.layerThick(1).* ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.ice(1);
            saturation_first_cell = (ground.STATVAR.waterIce(1) - ground.PARA.field_capacity .* remaining_pore_space) ./ ...
                (ground.STATVAR.layerThick(1).*ground.STATVAR.area(1) - remaining_pore_space); 
            saturation_first_cell = max(0,min(1,saturation_first_cell)); % 0 water at field capacity, 1: water at saturation
            %NEW SW
            saturation_first_cell(saturation_first_cell >= (1 - 1e-9)) = 1;
            
            ground.TEMP.F_ub_water = rainfall .* reduction_factor_in(saturation_first_cell, ground);
            ground.TEMP.surface_runoff = rainfall - ground.TEMP.F_ub_water;  %route this to surface pool
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
        
        
        function ground = get_boundary_condition_u_water_SNOW2(ground, forcing) %infiltrate all water
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  

            ground.TEMP.F_ub_water = rainfall ;
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
        
        
        %LAKE unfrozen
        function ground = get_boundary_condition_u_water_LAKE(ground, forcing)
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  
            snowfall = forcing.TEMP.snowfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);
            ground.TEMP.F_ub_water = rainfall + snowfall;
            
            T_rainWater =  max(0,forcing.TEMP.Tair);
            T_snow = min(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = rainfall .* ground.CONST.c_w .* T_rainWater + snowfall .* (ground.CONST.c_i .* T_snow - ground.CONST.L_f); 
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
        
        %LAKE frozen
        function ground = get_boundary_condition_u_water_LAKE_frozen(ground, forcing)
            rainfall = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600 .* ground.STATVAR.area(1);  
            ground.TEMP.F_ub_water = rainfall;
            
            T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = rainfall .* ground.CONST.c_w .* T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
        
        %zero water flux lower boundary condition
        function ground = get_boundary_condition_l_water2(ground)
            ground.TEMP.F_lb_water = 0;
            ground.TEMP.F_lb_water_energy = 0;
            
            ground.TEMP.d_water(end) = ground.TEMP.d_water(end) + ground.TEMP.F_lb_water;
            ground.TEMP.d_water_energy(end) = ground.TEMP.d_water_energy(end) + ground.TEMP.F_lb_water_energy;
        end
        
        %water fluxes for bucketW hydrological scheme (the "2" has no meaning)
        function ground = get_derivative_water2(ground) %adapts the fluxes automatically so that no checks are necessary when advancing the prognostic variable 
            saturation = (ground.STATVAR.waterIce - ground.STATVAR.field_capacity .* ground.STATVAR.layerThick.*ground.STATVAR.area)./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.field_capacity.*ground.STATVAR.layerThick.*ground.STATVAR.area);
            saturation = max(0,min(1,saturation)); % 0 water at field capacity, 1: water at saturation
            saturation(saturation >= (1 - 1e-6)) = 1;
            saturation(saturation <= (0 + 1e-6)) = 0;
            
            guaranteed_flow = ground.TEMP.d_water_ET;  %add other external fluxes here
            guaranteed_flow_energy = ground.TEMP.d_water_ET_energy;
            
            %outflow
            %saturation_water = ground.STATVAR.water./ground.STATVAR.waterIce;
            vol_water = ground.STATVAR.water ./ (ground.STATVAR.layerThick.*ground.STATVAR.area);
            d_water_out = max(0, ground.STATVAR.hydraulicConductivity .* ground.STATVAR.area ); %.* saturation_water
            d_water_out(vol_water <= ground.STATVAR.field_capacity) = 0;

            guaranteed_inflow = guaranteed_flow.* double(guaranteed_flow > 0); 
            d_water_out = double(guaranteed_inflow >= d_water_out) .* d_water_out + double(guaranteed_inflow < d_water_out) .* ...
                 (guaranteed_inflow + (d_water_out - guaranteed_inflow) .* reduction_factor_out(saturation, ground)); %this is positive when flowing out
            d_water_out(end,1) = 0; % lower boundary handled elsewhere
             %d_water_out(end,1) = -ground.TEMP.F_lb_water; %positive
             
            %inflow
            d_water_in = d_water_out .*0;
            d_water_in(2:end) = d_water_out(1:end-1);
            guaranteed_outflow = guaranteed_flow.* double(guaranteed_flow < 0);
            d_water_in = double(-guaranteed_outflow >= d_water_in) .* d_water_in + double(-guaranteed_outflow < d_water_in) .* ...
                (-guaranteed_outflow + (d_water_in + guaranteed_outflow).* reduction_factor_in(saturation, ground));
            %d_water_in(1) = ground.TEMP.F_ub_water; %already checked in UB, that space is available
            
            %avoid rounding errors
            %saturated = ground.STATVAR.layerThick.*ground.STATVAR.area <= ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic;
            %d_water_in(saturated) = min(d_water_in(saturated), -guaranteed_outflow(saturated));
            
            %readjust outflow
            d_water_out(1:end-1) = d_water_in(2:end); %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* (double(ground.STATVAR.T>=0) .* ground.CONST.c_w + double(ground.STATVAR.T<0) .* ground.CONST.c_i) .* ground.STATVAR.T;
            d_water_in_energy = d_water_out.*0;
            d_water_in_energy(2:end,1) = d_water_out_energy(1:end-1,1); 
            %d_water_in_energy(1) = ground.TEMP.F_ub_water_energy;
            %d_water_out_energy(end) = -ground.TEMP.F_lb_water_energy;
            
            %sum up               
            ground.TEMP.d_water = ground.TEMP.d_water + guaranteed_flow - d_water_out + d_water_in; 
            ground.TEMP.d_water_energy = ground.TEMP.d_water_energy + guaranteed_flow_energy - d_water_out_energy + d_water_in_energy;
            
            ground.TEMP.d_water_in = d_water_in; % at this stage nice-to-have variables, good for troubleshooting
            ground.TEMP.d_water_out = d_water_out;
        end
        
        
        function ground = get_derivative_water_Xice(ground) %adapts the fluxes automatically so that no checks are necessary when advancing the prognostic variable
            volume_matrix = ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce;
            saturation = (ground.STATVAR.waterIce  - ground.STATVAR.field_capacity .* volume_matrix)./...
                (volume_matrix - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.field_capacity .* volume_matrix);
            saturation = max(0,min(1,saturation)); % 0 water at field capacity, 1: water at saturation
            %test to avoid problems with too small timesteps
            saturation(saturation >= (1 - 1e-9)) = 1;
            
            guaranteed_flow = ground.TEMP.d_water_ET;  %add other external fluxes here
            guaranteed_flow_energy = ground.TEMP.d_water_ET_energy;
            
            %saturation(10)-1
            
            %outflow
            d_water_out = ground.STATVAR.hydraulicConductivity .* ground.STATVAR.area; % area cancels out; make this depended on both involved cells? 
            guaranteed_inflow = guaranteed_flow.* double(guaranteed_flow > 0); 
            d_water_out = double(guaranteed_inflow >= d_water_out) .* d_water_out + double(guaranteed_inflow < d_water_out) .* ...
                 (guaranteed_inflow + (d_water_out - guaranteed_inflow) .* reduction_factor_out(saturation, ground)); %this is positive when flowing out
            d_water_out(end,1) = 0; % lower boundary handled elsewhere
             %d_water_out(end,1) = -ground.TEMP.F_lb_water; %positive
             
            %inflow
            d_water_in = d_water_out .*0;
            d_water_in(2:end) = d_water_out(1:end-1);
            guaranteed_outflow = guaranteed_flow.* double(guaranteed_flow < 0);
            %d_water_in(10)
            %test= reduction_factor_in(saturation, ground);
            %test(10)
            d_water_in = double(-guaranteed_outflow >= d_water_in) .* d_water_in + double(-guaranteed_outflow < d_water_in) .* ...
                (-guaranteed_outflow + (d_water_in + guaranteed_outflow).* reduction_factor_in(saturation, ground));
            %d_water_in(1) = ground.TEMP.F_ub_water; %already checked in UB, that space is available
            %d_water_in(10)
            %readjust outflow
            d_water_out(1:end-1) = d_water_in(2:end); %reduce outflow if inflow is impossible
            
            
            %energy advection
            d_water_out_energy = d_water_out .* (double(ground.STATVAR.T>=0) .* ground.CONST.c_w + double(ground.STATVAR.T<0) .* ground.CONST.c_i) .* ground.STATVAR.T;
            d_water_in_energy = d_water_out.*0;
            d_water_in_energy(2:end,1) = d_water_out_energy(1:end-1,1); 
            %d_water_in_energy(1) = ground.TEMP.F_ub_water_energy;
            %d_water_out_energy(end) = -ground.TEMP.F_lb_water_energy;
            
            %sum up               
            ground.TEMP.d_water = ground.TEMP.d_water + guaranteed_flow - d_water_out + d_water_in; 
            ground.TEMP.d_water_energy = ground.TEMP.d_water_energy + guaranteed_flow_energy - d_water_out_energy + d_water_in_energy;
            
            ground.TEMP.d_water_in = d_water_in; % at this stage nice-to-have variables, good for troubleshooting, later necessary to route solutes
            ground.TEMP.d_water_out = d_water_out;
            
            
%             timestep = double(ground.TEMP.d_water(1) > 0) .* (ground.STATVAR.layerThick(1) .* ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.waterIce(1) - ground.STATVAR.XwaterIce(1)) ...
%                 ./ ground.TEMP.d_water(1);
%             if timestep>0 && timestep<1e-13
% %                 disp(saturation_first_cell)
%                 disp(ground.TEMP.F_ub_water)
% %                 disp(double(rainfall <= -evap) .* rainfall)
% %                 disp(double(rainfall > -evap) .* (-evap + (rainfall + evap) .* reduction_factor_in(saturation_first_cell, ground)))
% %                 disp(reduction_factor_in(saturation_first_cell, ground))
%                 disp('next2')
%             end
        end
        
        %upward fluxes of excess water in Xice classes
        function ground = get_derivative_Xwater(ground)  %routes Xwater up when Xice has melted
            %saturation = ground.STATVAR.Xwater ./ ground.STATVAR.area ./ (ground.PARA.hydraulicConductivity .* ground.PARA.dt_max);
            %saturation = max(0,min(1,saturation)); % 0 no Xwater, 1: water routed up within maximum timestep
            d_Xwater_out = ground.STATVAR.hydraulicConductivity .* ground.STATVAR.area; %  .* reduction_factor_out(saturation, ground); 
            d_Xwater_out(1,1) = 0; %Xwater stays in uppermost cell, must be removed elesewhere
            
            d_Xwater_out(d_Xwater_out > 0) = min(d_Xwater_out(d_Xwater_out > 0), ground.STATVAR.Xwater(d_Xwater_out > 0) ./ ground.PARA.dt_max); %makes explicit timestep check unnecessary

            d_Xwater_in = d_Xwater_out .*0;
            d_Xwater_in(1:end-1) = d_Xwater_out(2:end) .* double(ground.STATVAR.T(1:end-1)>0); % water can only be taken up by unfrozen cells, important for lateral Xice melt (e.g. palsa case)
            
            d_Xwater_out(2:end) = d_Xwater_in(1:end-1); %reduce outflow if inflow is impossible
            
            d_Xwater_out_energy = d_Xwater_out .* ground.CONST.c_w .* ground.STATVAR.T;
            d_Xwater_in_energy = d_Xwater_out .*0;
            d_Xwater_in_energy(1:end-1) = d_Xwater_out_energy(2:end);
            
            ground.TEMP.d_Xwater = ground.TEMP.d_Xwater - d_Xwater_out + d_Xwater_in; 
            %ground.TEMP.d_Xwater(1) = ground.TEMP.d_Xwater(1) + ground.TEMP.F_ub_Xwater;
            ground.TEMP.d_Xwater_energy = ground.TEMP.d_Xwater_energy - d_Xwater_out_energy + d_Xwater_in_energy;
            %ground.TEMP.d_Xwater_energy(1) = ground.TEMP.d_Xwater_energy(1) + ground.TEMP.F_ub_Xwater_energy;
        end
        
        %bucketW hydrological scheme for SNOW classes
        function ground = get_derivative_water_SNOW(ground) %adapts the fluxes automatically so that no checks are necessary when advancing the prognostic variable
            remaining_pore_space = ground.STATVAR.layerThick.* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.ice;
            %saturation = (ground.STATVAR.waterIce - ground.PARA.field_capacity .* remaining_pore_space) ./ ...
            %    (ground.STATVAR.layerThick.*ground.STATVAR.area - remaining_pore_space); 
            saturation = (ground.STATVAR.water - ground.PARA.field_capacity .* remaining_pore_space) ./ ...
                (remaining_pore_space - ground.PARA.field_capacity .* remaining_pore_space); 
            
            saturation = max(0,min(1,saturation)); % 0 water at field capacity, 1: water at saturation
            saturation(saturation >= (1 - 1e-9)) = 1;
            
            %outflow
            d_water_out = ground.STATVAR.hydraulicConductivity .* ground.STATVAR.area; % area cancels out; make this depended on both involved cells?
            d_water_out = d_water_out .* reduction_factor_out(saturation, ground); %this is positive when flowing out
            d_water_out(end,1) = 0; % lower boundary handled elsewhere
            %d_water_out(end,1) = -ground.TEMP.F_lb_water; %positive
            
            %inflow
            d_water_in = d_water_out .*0;
            d_water_in(2:end) = d_water_out(1:end-1);
            d_water_in = d_water_in .* reduction_factor_in(saturation, ground);
            %d_water_in(1) = ground.TEMP.F_ub_water; %already checked in UB, that space is available
            
            %readjust outflow
            d_water_out(1:end-1) = d_water_in(2:end); %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* ground.CONST.c_w .* ground.STATVAR.T;
            d_water_in_energy = d_water_out.*0;
            d_water_in_energy(2:end,1) = d_water_out_energy(1:end-1,1);
            %d_water_in_energy(1) = ground.TEMP.F_ub_water_energy;
            %d_water_out_energy(end) = -ground.TEMP.F_lb_water_energy;
            
            %sum up
            ground.TEMP.d_water = ground.TEMP.d_water - d_water_out + d_water_in;
            ground.TEMP.d_water_energy = ground.TEMP.d_water_energy - d_water_out_energy + d_water_in_energy;
            
            ground.TEMP.d_water_in = d_water_in; % at this stage nice-to-have variables, good for troubleshooting
            ground.TEMP.d_water_out = d_water_out;
        end
        
     
        %Richards equation         
        function ground = get_derivative_water_RichardsEq(ground) %adapts the fluxes automatically so that no checks are necessary when advancing the prognostic variable

            waterIce_saturation = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            waterIce_saturation = max(0,min(1,waterIce_saturation));
            water_saturation = (ground.STATVAR.waterIce - ground.PARA.min_waterIce .* ground.STATVAR.layerThick.*ground.STATVAR.area) ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            %test Sebastian
            %water_saturation = ground.STATVAR.water ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            water_saturation = max(0,min(1,water_saturation));
            water_saturation(water_saturation >= (1 - 1e-9)) = 1;
            
            guaranteed_flow = ground.TEMP.d_water_ET;  %add other external fluxes here
            guaranteed_flow_energy = ground.TEMP.d_water_ET_energy;
            
            k_eff = ground.STATVAR.hydraulicConductivity(1:end-1,1).*ground.STATVAR.hydraulicConductivity(2:end,1) ./ ...
                (ground.STATVAR.hydraulicConductivity(1:end-1,1).*ground.STATVAR.layerThick(2:end,1)./2 + ground.STATVAR.hydraulicConductivity(2:end,1).*ground.STATVAR.layerThick(1:end-1,1)./2);
            
            fluxes = -k_eff.* (ground.STATVAR.waterPotential(1:end-1,1) - ground.STATVAR.waterPotential(2:end,1) + (ground.STATVAR.layerThick(2:end,1)./2 + ground.STATVAR.layerThick(1:end-1,1)./2)) .* 0.5.* (ground.STATVAR.area(1:end-1,1) + ground.STATVAR.area(2:end,1));
            %minus means flux downwards
            %ground.TEMP.fluxes_prev = fluxes;
            %guaranteed_flow = ground.TEMP.d_water_ET(end) + ground.TEMP.d_water(end,1);
            
            %1. reduce fluxes due to lack of water than can flow out
            guaranteed_inflow = guaranteed_flow.* double(guaranteed_flow > 0); 
            flux_out_down = ground.STATVAR.hydraulicConductivity .* 0;
            flux_out_down(1:end-1,1) = -fluxes .* double(fluxes<0);
            flux_out_up = ground.STATVAR.hydraulicConductivity .* 0;
            flux_out_up(2:end,1) = fluxes .* double(fluxes>0);
            flux_out = flux_out_down + flux_out_up;
            
            flux_out = double(guaranteed_inflow >= flux_out) .* flux_out + double(guaranteed_inflow < flux_out) .* ...
                 (guaranteed_inflow + (flux_out - guaranteed_inflow) .* reduction_factor_out(water_saturation, ground)); %this is positive when flowing out
            flux_out_down = flux_out ./ (flux_out_down + flux_out_up) .* flux_out_down;
            flux_out_down(isnan(flux_out_down)) = 0;
            flux_out_up = flux_out ./ (flux_out_down + flux_out_up) .* flux_out_up;
            flux_out_up(isnan(flux_out_up)) = 0;
            
            %2.inflow             
            guaranteed_outflow = guaranteed_flow.* double(guaranteed_flow < 0); 
            flux_in_from_above = ground.STATVAR.hydraulicConductivity .* 0;
            flux_in_from_above(2:end,1) = -fluxes .* double(fluxes<0);
            flux_in_from_below = ground.STATVAR.hydraulicConductivity .* 0;
            flux_in_from_below(1:end-1,1) = fluxes .* double(fluxes>0);
            flux_in = flux_in_from_above + flux_in_from_below;
            
            flux_in = double(-guaranteed_outflow >= flux_in) .* flux_in + double(-guaranteed_outflow < flux_in) .* ...
                (-guaranteed_outflow + (flux_in + guaranteed_outflow).* reduction_factor_in(waterIce_saturation, ground));
            flux_in_from_above = flux_in ./ (flux_in_from_above + flux_in_from_below) .* flux_in_from_above;
            flux_in_from_above(isnan(flux_in_from_above)) = 0;
            
            flux_in_from_below = flux_in ./ (flux_in_from_above + flux_in_from_below) .* flux_in_from_below;
            flux_in_from_below(isnan(flux_in_from_below)) = 0;
            
            fluxes = max(fluxes, -flux_out_down(1:end-1,1));
            fluxes = max(fluxes, -flux_in_from_above(2:end,1));
            fluxes = min(fluxes, flux_out_up(2:end,1));
            fluxes = min(fluxes, flux_in_from_below(1:end-1,1));
            
            ground.TEMP.fluxes = fluxes;
            
            %same as for bucketW
            d_water_out = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_out(1:end-1,1) = -fluxes .* double(fluxes <0);
            d_water_out(2:end,1) = d_water_out(2:end,1)  + fluxes .* double(fluxes >0);
            
            d_water_in_from_above = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_in_from_above(2:end,1) = -fluxes .* double(fluxes<0);
            d_water_in_from_below = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_in_from_below(1:end-1,1) = fluxes .* double(fluxes>0);
            
            %energy advection
            d_water_out_energy = d_water_out .* (double(ground.STATVAR.T>=0) .* ground.CONST.c_w + double(ground.STATVAR.T<0) .* ground.CONST.c_i) .* ground.STATVAR.T;
            d_water_in_energy = d_water_out_energy.*0;
            d_water_in_energy(2:end,1) = d_water_in_energy(2:end,1) + d_water_in_from_above(2:end,1) .* ...
                (double(ground.STATVAR.T(1:end-1,1)>=0) .* ground.CONST.c_w + double(ground.STATVAR.T(1:end-1,1)<0) .* ground.CONST.c_i) .* ground.STATVAR.T(1:end-1,1);
            d_water_in_energy(1:end-1,1) = d_water_in_energy(1:end-1,1) + d_water_in_from_below(1:end-1,1) .* ...
                (double(ground.STATVAR.T(2:end,1)>=0) .* ground.CONST.c_w + double(ground.STATVAR.T(2:end,1)<0) .* ground.CONST.c_i) .* ground.STATVAR.T(2:end,1);
            
            
            %sum up               
            ground.TEMP.d_water = ground.TEMP.d_water + guaranteed_flow - d_water_out + d_water_in_from_above + d_water_in_from_below; 
            ground.TEMP.d_water_energy = ground.TEMP.d_water_energy + guaranteed_flow_energy - d_water_out_energy + d_water_in_energy;
            
            ground.TEMP.d_water_in_from_above = d_water_in_from_above; % at this stage nice-to-have variables, good for troubleshooting
            ground.TEMP.d_water_in_from_below = d_water_in_from_below;
            ground.TEMP.d_water_out = d_water_out;

        end
        
        %Richards equation with pressure
        function ground = get_derivative_water_RichardsEq_pressure(ground) %adapts the fluxes automatically so that no checks are necessary when advancing the prognostic variable
                                   
            k_eff = ground.STATVAR.hydraulicConductivity(1:end-1,1).*ground.STATVAR.hydraulicConductivity(2:end,1) ./ ...
                (ground.STATVAR.hydraulicConductivity(1:end-1,1).*ground.STATVAR.layerThick(2:end,1)./2 + ground.STATVAR.hydraulicConductivity(2:end,1).*ground.STATVAR.layerThick(1:end-1,1)./2);

            ground.TEMP.d_head_waterPotential = ground.STATVAR.waterPotential(1:end-1,1) - ground.STATVAR.waterPotential(2:end,1);
            ground.TEMP.d_head_gravitationalPotential = ground.STATVAR.gravitationalPotential(1:end-1,1) - ground.STATVAR.gravitationalPotential(2:end,1);
            %ground.TEMP.d_head_gravity = (double(ground.STATVAR.saturation(1:end-1,1) <= 1-1e-6) .* ground.STATVAR.layerThick(1:end-1,1) + double(ground.STATVAR.saturation(2:end,1) <= 1-1e-6) .* ground.STATVAR.layerThick(2:end,1))./2;
            %ground.TEMP.d_head_hydrostatic_pressure = double(ground.STATVAR.saturation(1:end-1,1) >= 1-1e-6 & ground.STATVAR.saturation(2:end,1) >= 1-1e-6) .* ((ground.STATVAR.hydrostatic_pressure(1:end-1,1) - ground.STATVAR.hydrostatic_pressure(2:end,1))./ground.CONST.density_water./ground.CONST.g);
            %ground.TEMP.d_head_soilMechanics = (double(ground.STATVAR.saturation(1:end-1,1) > 1-1e-6) .* double(ground.STATVAR.saturation(2:end,1) > 1-1e-6) .* ((ground.STATVAR.overburden_pressure(1:end-1,1) - ground.STATVAR.bearing_capacity(1:end-1,1)) ...
            %    - ( ground.STATVAR.overburden_pressure(2:end,1) - ground.STATVAR.bearing_capacity(2:end,1)))) ./ ground.CONST.density_water./ground.CONST.g;           
            %ground.TEMP.d_head_soilMechanics = (double(ground.STATVAR.saturation(2:end,1) > 1-1e-6) .* ((ground.STATVAR.overburden_pressure(1:end-1,1) - ground.STATVAR.bearing_capacity(1:end-1,1)) ...
            %    - ( ground.STATVAR.overburden_pressure(2:end,1) - ground.STATVAR.bearing_capacity(2:end,1)))) ./ ground.CONST.density_water./ground.CONST.g;           
            %%%ground.TEMP.d_head_soilMechanics = ((double(ground.STATVAR.saturation(1:end-1,1) > 1-1e-6) .* (ground.STATVAR.overburden_pressure(1:end-1,1) - ground.STATVAR.bearing_capacity(1:end-1,1))) ...
            %%%    - (double(ground.STATVAR.saturation(2:end,1) > 1-1e-6) .* (ground.STATVAR.overburden_pressure(2:end,1) - ground.STATVAR.bearing_capacity(2:end,1)))) ./ ground.CONST.density_water./ground.CONST.g;
            
            ground.TEMP.d_head_soilMechanics = ((double(ground.STATVAR.saturation(1:end-1,1) > 1-1e-6) .* (ground.STATVAR.overburden_pressure(1:end-1,1) - ground.STATVAR.bearing_capacity(1:end-1,1))) ...
                - (double(ground.STATVAR.saturation(2:end,1) > 1-1e-6) .*(ground.STATVAR.overburden_pressure(2:end,1) - ground.STATVAR.bearing_capacity(2:end,1)))) ./ ground.CONST.density_water./ground.CONST.g;
            
            ground.TEMP.d_head = ground.TEMP.d_head_waterPotential + ground.TEMP.d_head_gravitationalPotential + ground.TEMP.d_head_soilMechanics;
            
            if ground.TEMP.d_head(5,1) < 0.01
                xxx = 0;
            end
            
            %fluxes from one cell into another
            fluxes = -k_eff.* ground.TEMP.d_head .* 0.5.* (ground.STATVAR.area(1:end-1,1) + ground.STATVAR.area(2:end,1));
            %fluxes to the surface --> has to be positive as water can be pressed out of soil but no water available to be drawn into the soil
            %flux_out_surface_soilMechanics = max(0,double(ground.STATVAR.saturation(1,1) >= 1-1e-6) .* ((ground.STATVAR.hydraulicConductivity(1,1) ./ (ground.STATVAR.layerThick(1,1)./2)) .* (ground.STATVAR.overburden_pressure(1,1) - ground.STATVAR.bearing_capacity(1,1) - 0) .* ground.STATVAR.area(1,1)));
           
            ground.TEMP.fluxes = fluxes;
         
            for i = 2 : size(ground.STATVAR.porosity,1)
                if ground.STATVAR.porosity(i,1) > 0.49
                    xxxx = 0;
                end
            end
            for i = 7
                if ground.STATVAR.porosity(i,1) < 0.4
                    xxxx = 0;
                end
            end
            
            %same as for bucketW
            d_water_out = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_out(1:end-1,1) = -fluxes .* double(fluxes <0);
            d_water_out(2:end,1) = d_water_out(2:end,1)  + fluxes .* double(fluxes >0);
            
            %d_water_out(1,1) = d_water_out(1,1) + flux_out_surface_soilMechanics;
            Xwater_out = 0;
            if ground.STATVAR.Xwater(1,1) / ground.STATVAR.waterIce(1,1) > 0.05
                Xwater_out = ground.STATVAR.Xwater(1,1) - 0.05 * ground.STATVAR.waterIce(1,1);
            end
            d_water_out(1,1) = d_water_out(1,1) + Xwater_out;
            
            d_water_in_from_above = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_in_from_above(2:end,1) = -fluxes .* double(fluxes<0);
            d_water_in_from_below = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_in_from_below(1:end-1,1) = fluxes .* double(fluxes>0);

            %energy advection
            d_water_out_energy = d_water_out .* (double(ground.STATVAR.T>=0) .* ground.CONST.c_w + double(ground.STATVAR.T<0) .* ground.CONST.c_i) .* ground.STATVAR.T;
            d_water_in_energy = d_water_out_energy.*0;
            d_water_in_energy(2:end,1) = d_water_in_energy(2:end,1) + d_water_in_from_above(2:end,1) .* ...
                (double(ground.STATVAR.T(1:end-1,1)>=0) .* ground.CONST.c_w + double(ground.STATVAR.T(1:end-1,1)<0) .* ground.CONST.c_i) .* ground.STATVAR.T(1:end-1,1);
            d_water_in_energy(1:end-1,1) = d_water_in_energy(1:end-1,1) + d_water_in_from_below(1:end-1,1) .* ...
                (double(ground.STATVAR.T(2:end,1)>=0) .* ground.CONST.c_w + double(ground.STATVAR.T(2:end,1)<0) .* ground.CONST.c_i) .* ground.STATVAR.T(2:end,1);
 
            %sum up
            ground.TEMP.d_water = ground.TEMP.d_water + ground.TEMP.d_water_ET - d_water_out + d_water_in_from_above + d_water_in_from_below;
            ground.TEMP.d_water_energy = ground.TEMP.d_water_energy + ground.TEMP.d_water_ET_energy - d_water_out_energy + d_water_in_energy;
            
            ground.TEMP.d_water_in_from_above = d_water_in_from_above; % at this stage nice-to-have variables, good for troubleshooting
            ground.TEMP.d_water_in_from_below = d_water_in_from_below;
            ground.TEMP.d_water_out = d_water_out;
            
            %this assumes that air can always be sucked in from "the side"
            %if the cell above or below is not suaturated.
            air_available_above = ground.STATVAR.saturation(1:end-1,1) < 1-1e-6;
            air_available_above = [1; air_available_above]; %add condition if air can be drawn from above this class (e.g. in case of lake)
            air_available_below = ground.STATVAR.saturation(2:end,1) < 1-1e-6;
            air_available_below = [air_available_below; 0]; % no air reservoir below last cell
            ground.TEMP.no_air = ground.TEMP.d_water <0 & ~air_available_above & ~air_available_below;

            dBearingCapacity_dPorosity = ground.STATVAR.reference_pressure .* log(10) .* 10.^((ground.STATVAR.initial_voidRatio - ground.STATVAR.porosity - ground.STATVAR.initial_voidRatio .* ground.STATVAR.porosity) ./ ...
                (ground.STATVAR.compression_index - ground.STATVAR.porosity .* ground.STATVAR.compression_index)) .* ...
                ((-1-ground.STATVAR.initial_voidRatio).*(ground.STATVAR.compression_index - ground.STATVAR.porosity.*ground.STATVAR.compression_index) + ...
                ground.STATVAR.compression_index.*(ground.STATVAR.initial_voidRatio - ground.STATVAR.porosity - ground.STATVAR.initial_voidRatio .* ground.STATVAR.porosity)) ./ ...
                ((ground.STATVAR.compression_index - ground.STATVAR.porosity .* ground.STATVAR.compression_index).^2);
            dPorosity_dWater = (ground.STATVAR.mineral + ground.STATVAR.organic) ./ ((ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic)).^2;
            
            ground.TEMP.dBearingCapacity_dWater = dBearingCapacity_dPorosity .* dPorosity_dWater;
               
        end
        
        %Richards equation excess ice, also handles fluxes of excess water
        %including formation of segregation ice - DISCONTINUED
        function ground = get_derivative_water_RichardsEq_Xice(ground) %adapts the fluxes automatically so that no checks are necessary when advancing the prognostic variable

            waterIce_saturation = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.XwaterIce);
            waterIce_saturation = max(0,min(1,waterIce_saturation));
            water_saturation = (ground.STATVAR.waterIce - ground.PARA.min_waterIce .* (ground.STATVAR.layerThick.*ground.STATVAR.area- ground.STATVAR.XwaterIce)) ./ ...
                (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.XwaterIce);
            water_saturation = max(0,min(1,water_saturation));
            water_saturation(water_saturation >= (1 - 1e-9)) = 1;
            
            guaranteed_flow = ground.TEMP.d_water_ET;  %add other external fluxes here
            guaranteed_flow_energy = ground.TEMP.d_water_ET_energy;
            
            k_eff = ground.STATVAR.hydraulicConductivity(1:end-1,1).*ground.STATVAR.hydraulicConductivity(2:end,1) .*...
                (ground.STATVAR.layerThick(2:end,1)  + ground.STATVAR.layerThick(1:end-1,1))./2 ./ ...
                (ground.STATVAR.hydraulicConductivity(1:end-1,1).*ground.STATVAR.layerThick(2:end,1)./2 + ...
                ground.STATVAR.hydraulicConductivity(2:end,1).*ground.STATVAR.layerThick(1:end-1,1)./2);
            %set conductivity to zero between each two cells that have Xice
            %k_eff = k_eff .* double((double(ground.STATVAR.Xice(1:end-1)>0) + double(ground.STATVAR.Xice(2:end)>0))<2);
            
            geopotential = -cumsum(ground.STATVAR.layerThick) + ground.STATVAR.layerThick./2; %midpoint of each cell - increases with with elevation
            
            %in m water equivalent - in principle has to take classses on
            %top into account as well, requires a class-specific
            %overburden_pressure function - decreases with elevation
            overburden_pressure = ((cumsum(ground.STATVAR.mineral) + ground.STATVAR.mineral./2) .* ground.CONST.rho_m ./ground.CONST.rho_w + ...
                (cumsum(ground.STATVAR.organic) + ground.STATVAR.organic./2) .* ground.CONST.rho_o ./ground.CONST.rho_w + ...
                (cumsum(ground.STATVAR.waterIce) + ground.STATVAR.waterIce./2 + cumsum(ground.STATVAR.XwaterIce) + ground.STATVAR.XwaterIce./2)) ./ ground.STATVAR.area; 
            
            total_potential = ground.STATVAR.waterPotential + geopotential + overburden_pressure .* double(ground.STATVAR.Xwater>0);
            fluxes = -k_eff .*(total_potential(1:end-1,1) - total_potential(2:end,1))./ ((ground.STATVAR.layerThick(2:end,1)  + ground.STATVAR.layerThick(1:end-1,1))./2) ...
                .* 0.5.* (ground.STATVAR.area(1:end-1,1) + ground.STATVAR.area(2:end,1));
            
            fluxes_uncorrected = fluxes;
            
            %account for guaranteed flow
            guaranteed_flow = ground.TEMP.d_water_ET;
            %1. reduce fluxes due to lack of water than can flow out
            guaranteed_inflow = guaranteed_flow.* double(guaranteed_flow > 0);
            flux_out_down = ground.STATVAR.hydraulicConductivity .* 0;
            flux_out_down(1:end-1,1) = -fluxes .* double(fluxes<0);
            flux_out_up = ground.STATVAR.hydraulicConductivity .* 0;
            flux_out_up(2:end,1) = fluxes .* double(fluxes>0);
            flux_out = flux_out_down + flux_out_up;
            
            flux_out = double(guaranteed_inflow >= flux_out) .* flux_out + double(guaranteed_inflow < flux_out) .* ...
                (guaranteed_inflow + (flux_out - guaranteed_inflow) .* reduction_factor_out(water_saturation, ground)); %this is positive when flowing out
            flux_out_down = flux_out ./ (flux_out_down + flux_out_up) .* flux_out_down;
            flux_out_down(isnan(flux_out_down)) = 0;
            flux_out_up = flux_out ./ (flux_out_down + flux_out_up) .* flux_out_up;
            flux_out_up(isnan(flux_out_up)) = 0;
            
            %2.inflow
            guaranteed_outflow = guaranteed_flow.* double(guaranteed_flow < 0);
            flux_in_from_above = ground.STATVAR.hydraulicConductivity .* 0;
            flux_in_from_above(2:end,1) = -fluxes .* double(fluxes<0);
            flux_in_from_below = ground.STATVAR.hydraulicConductivity .* 0;
            flux_in_from_below(1:end-1,1) = fluxes .* double(fluxes>0);
            flux_in = flux_in_from_above + flux_in_from_below;
            
            flux_in = double(-guaranteed_outflow >= flux_in) .* flux_in + double(-guaranteed_outflow < flux_in) .* ...
                (-guaranteed_outflow + (flux_in + guaranteed_outflow).* reduction_factor_in(waterIce_saturation, ground));
            flux_in_from_above = flux_in ./ (flux_in_from_above + flux_in_from_below) .* flux_in_from_above;
            flux_in_from_above(isnan(flux_in_from_above)) = 0;
            
            flux_in_from_below = flux_in ./ (flux_in_from_above + flux_in_from_below) .* flux_in_from_below;
            flux_in_from_below(isnan(flux_in_from_below)) = 0;
            
            fluxes = max(fluxes, -flux_out_down(1:end-1,1));
            fluxes = max(fluxes, -flux_in_from_above(2:end,1));
            fluxes = min(fluxes, flux_out_up(2:end,1));
            fluxes = min(fluxes, flux_in_from_below(1:end-1,1));

            %formation of excess ice
            %check for all the fluxes into saturated cells, if the pressure is high enough to overcome the overburden pressure
            %if yes, leave the flux as it is and mark these cells (needed for timestep)
            %if no, leave the flux as it is
            ground.TEMP.XwaterIce_formation = ground.STATVAR.waterIce .* 0;
            
            for i=1:size(fluxes,1)
                inflow_cell = i+double(sign(fluxes_uncorrected(i,1))==-1);
                if waterIce_saturation(inflow_cell,1)>=1-1e-9 && ground.STATVAR.Xwater(inflow_cell)<=0
                    total_potential(inflow_cell) = total_potential(inflow_cell) + overburden_pressure(inflow_cell);
                    flux_pot = -k_eff(i,1) .*(total_potential(i,1) - total_potential(i+1,1))./ ((ground.STATVAR.layerThick(i,1) + ground.STATVAR.layerThick(i+1,1))./2) ...
                        .* 0.5.* (ground.STATVAR.area(i,1) + ground.STATVAR.area(i+1,1));
                    if sign(fluxes_uncorrected(i,1))==sign(flux_pot)
                        fluxes(i,1)  = flux_pot;
                        ground.TEMP.XwaterIce_formation(inflow_cell,1) = 1;
                    end
                    total_potential(inflow_cell) = total_potential(inflow_cell) - overburden_pressure(inflow_cell);
                end
            end
  
            ground.TEMP.fluxes=fluxes;
            ground.TEMP.total_potential = total_potential;
            ground.TEMP.overburden_pressure = overburden_pressure;

            %same as for bucketW
            d_water_out = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_out(1:end-1,1) = -fluxes .* double(fluxes <0);
            d_water_out(2:end,1) = d_water_out(2:end,1)  + fluxes .* double(fluxes >0);
            
            d_water_in_from_above = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_in_from_above(2:end,1) = -fluxes .* double(fluxes<0);
            d_water_in_from_below = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_in_from_below(1:end-1,1) = fluxes .* double(fluxes>0);
            
            %energy advection
            d_water_out_energy = d_water_out .* (double(ground.STATVAR.T>=0) .* ground.CONST.c_w + double(ground.STATVAR.T<0) .* ground.CONST.c_i) .* ground.STATVAR.T;
            d_water_in_energy = d_water_out_energy.*0;
            d_water_in_energy(2:end,1) = d_water_in_energy(2:end,1) + d_water_in_from_above(2:end,1) .* ...
                (double(ground.STATVAR.T(1:end-1,1)>=0) .* ground.CONST.c_w + double(ground.STATVAR.T(1:end-1,1)<0) .* ground.CONST.c_i) .* ground.STATVAR.T(1:end-1,1);
            d_water_in_energy(1:end-1,1) = d_water_in_energy(1:end-1,1) + d_water_in_from_below(1:end-1,1) .* ...
                (double(ground.STATVAR.T(2:end,1)>=0) .* ground.CONST.c_w + double(ground.STATVAR.T(2:end,1)<0) .* ground.CONST.c_i) .* ground.STATVAR.T(2:end,1);
            
            
            %sum up               
            ground.TEMP.d_water = ground.TEMP.d_water + guaranteed_flow - d_water_out + d_water_in_from_above + d_water_in_from_below; 
            ground.TEMP.d_water_energy = ground.TEMP.d_water_energy + guaranteed_flow_energy - d_water_out_energy + d_water_in_energy;
            
            ground.TEMP.d_water_in_from_above = d_water_in_from_above; % at this stage nice-to-have variables, good for troubleshooting
            ground.TEMP.d_water_in_from_below = d_water_in_from_below;
            ground.TEMP.d_water_out = d_water_out;

        end
        
         %Richards equation excess ice
        function ground = get_derivative_water_RichardsEq_Xice2(ground) %adapts the fluxes automatically so that no checks are necessary when advancing the prognostic variable

            waterIce_saturation = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.XwaterIce);
            waterIce_saturation = max(0,min(1,waterIce_saturation));
            water_saturation = (ground.STATVAR.waterIce - ground.PARA.min_waterIce .* (ground.STATVAR.layerThick.*ground.STATVAR.area- ground.STATVAR.XwaterIce)) ./ ...
                (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.XwaterIce ...
                - ground.PARA.min_waterIce .* (ground.STATVAR.layerThick.*ground.STATVAR.area- ground.STATVAR.XwaterIce));
            water_saturation = max(0,min(1,water_saturation));
            water_saturation(water_saturation >= (1 - 1e-6 .* rand(size(water_saturation,1),1))) = 1;
            waterIce_saturation(waterIce_saturation >= (1 - 1e-6 .* rand(size(waterIce_saturation,1),1))) = 1;

            
%             guaranteed_flow = ground.TEMP.d_water_ET;  %add other external fluxes here
%             guaranteed_flow_energy = ground.TEMP.d_water_ET_energy;
           
            guaranteed_flow = ground.TEMP.d_water_ET + ground.TEMP.d_water;  %add fluxfrom first grid cell
            guaranteed_flow_energy = ground.TEMP.d_water_ET_energy + ground.TEMP.d_water_energy;
            
            
            k_eff = ground.STATVAR.hydraulicConductivity(1:end-1,1).*ground.STATVAR.hydraulicConductivity(2:end,1) ./ ...
                (ground.STATVAR.hydraulicConductivity(1:end-1,1).*ground.STATVAR.layerThick(2:end,1)./2 + ground.STATVAR.hydraulicConductivity(2:end,1).*ground.STATVAR.layerThick(1:end-1,1)./2);
            
            fluxes = -k_eff.* (ground.STATVAR.waterPotential(1:end-1,1) - ground.STATVAR.waterPotential(2:end,1) + (ground.STATVAR.layerThick(2:end,1)./2 + ground.STATVAR.layerThick(1:end-1,1)./2)) .* 0.5.* (ground.STATVAR.area(1:end-1,1) + ground.STATVAR.area(2:end,1));
            %minus means flux downwards
            
            %1. reduce fluxes due to lack of water than can flow out
            guaranteed_inflow = guaranteed_flow.* double(guaranteed_flow > 0); 
            flux_out_down = ground.STATVAR.hydraulicConductivity .* 0;
            flux_out_down(1:end-1,1) = -fluxes .* double(fluxes<0);
            flux_out_up = ground.STATVAR.hydraulicConductivity .* 0;
            flux_out_up(2:end,1) = fluxes .* double(fluxes>0);
            flux_out = flux_out_down + flux_out_up;
            
            flux_out = double(guaranteed_inflow >= flux_out) .* flux_out + double(guaranteed_inflow < flux_out) .* ...
                 (guaranteed_inflow + (flux_out - guaranteed_inflow) .* reduction_factor_out(water_saturation, ground)); %this is positive when flowing out
            flux_out_down = flux_out ./ (flux_out_down + flux_out_up) .* flux_out_down;
            flux_out_down(isnan(flux_out_down)) = 0;
            flux_out_up = flux_out ./ (flux_out_down + flux_out_up) .* flux_out_up;
            flux_out_up(isnan(flux_out_up)) = 0;
            
            %2.inflow             
            guaranteed_outflow = guaranteed_flow.* double(guaranteed_flow < 0); 
            flux_in_from_above = ground.STATVAR.hydraulicConductivity .* 0;
            flux_in_from_above(2:end,1) = -fluxes .* double(fluxes<0);
            flux_in_from_below = ground.STATVAR.hydraulicConductivity .* 0;
            flux_in_from_below(1:end-1,1) = fluxes .* double(fluxes>0);
            flux_in = flux_in_from_above + flux_in_from_below;
            
            flux_in = double(-guaranteed_outflow >= flux_in) .* flux_in + double(-guaranteed_outflow < flux_in) .* ...
                (-guaranteed_outflow + (flux_in + guaranteed_outflow).* reduction_factor_in(waterIce_saturation, ground));
            flux_in_from_above = flux_in ./ (flux_in_from_above + flux_in_from_below) .* flux_in_from_above;
            flux_in_from_above(isnan(flux_in_from_above)) = 0;
            
            flux_in_from_below = flux_in ./ (flux_in_from_above + flux_in_from_below) .* flux_in_from_below;
            flux_in_from_below(isnan(flux_in_from_below)) = 0;
            
            fluxes = max(fluxes, -flux_out_down(1:end-1,1));
            fluxes = max(fluxes, -flux_in_from_above(2:end,1));
            fluxes = min(fluxes, flux_out_up(2:end,1));
            fluxes = min(fluxes, flux_in_from_below(1:end-1,1));
            
            ground.TEMP.fluxes = fluxes;
            
            %same as for bucketW
            d_water_out = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_out(1:end-1,1) = -fluxes .* double(fluxes <0);
            d_water_out(2:end,1) = d_water_out(2:end,1)  + fluxes .* double(fluxes >0);
            
            d_water_in_from_above = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_in_from_above(2:end,1) = -fluxes .* double(fluxes<0);
            d_water_in_from_below = ground.STATVAR.hydraulicConductivity .* 0;
            d_water_in_from_below(1:end-1,1) = fluxes .* double(fluxes>0);
            
            %energy advection
            d_water_out_energy = d_water_out .* (double(ground.STATVAR.T>=0) .* ground.CONST.c_w + double(ground.STATVAR.T<0) .* ground.CONST.c_i) .* ground.STATVAR.T;
            d_water_in_energy = d_water_out_energy.*0;
            d_water_in_energy(2:end,1) = d_water_in_energy(2:end,1) + d_water_in_from_above(2:end,1) .* ...
                (double(ground.STATVAR.T(1:end-1,1)>=0) .* ground.CONST.c_w + double(ground.STATVAR.T(1:end-1,1)<0) .* ground.CONST.c_i) .* ground.STATVAR.T(1:end-1,1);
            d_water_in_energy(1:end-1,1) = d_water_in_energy(1:end-1,1) + d_water_in_from_below(1:end-1,1) .* ...
                (double(ground.STATVAR.T(2:end,1)>=0) .* ground.CONST.c_w + double(ground.STATVAR.T(2:end,1)<0) .* ground.CONST.c_i) .* ground.STATVAR.T(2:end,1);
            
            %sum up               
%             ground.TEMP.d_water = ground.TEMP.d_water + guaranteed_flow - d_water_out + d_water_in_from_above + d_water_in_from_below;
%             ground.TEMP.d_water_energy = ground.TEMP.d_water_energy + guaranteed_flow_energy - d_water_out_energy + d_water_in_energy;
            
            ground.TEMP.d_water = guaranteed_flow - d_water_out + d_water_in_from_above + d_water_in_from_below; %d_water part of guaranteed_flow
            ground.TEMP.d_water_energy = guaranteed_flow_energy - d_water_out_energy + d_water_in_energy;
            
            ground.TEMP.d_water_in_from_above = d_water_in_from_above; % at this stage nice-to-have variables, good for troubleshooting
            ground.TEMP.d_water_in_from_below = d_water_in_from_below;
            ground.TEMP.d_water_out = d_water_out;
      
        end
        
        function ground = gravitational_potential(ground) %calculates gravitational potential
            
            %Find all saturated cells
            saturated = ground.STATVAR.saturation > 1-1e-6;
            id_unsaturated = find(saturated == 0);
            
            %Calculate gravitational potential, assuming that all gridcells are unsaturated
            for i = size(saturated,1) : -1 : 1
                gravitationalPotential_unsaturated(i,1) = sum(ground.STATVAR.layerThick(i:size(saturated,1),1));
            end
            
            %Calculate gravitational potential
            gravitationalPotential = [];
            for i = 1 : size(saturated,1)
                if saturated(i,1) == 0
                    gravitationalPotential(i,1) = gravitationalPotential_unsaturated(i,1);
                    if i ~= size(saturated,1) && saturated(i+1,1) == 1 && ground.STATVAR.T(i+1,1) > 0 %&& abs(ground.STATVAR.waterPotential(i,1)-ground.STATVAR.waterPotential(i+1,1)) < gravitationalPotential_unsaturated(i,1)-gravitationalPotential_unsaturated(i+1,1)
                        %If gridcell below is saturated and unfrozen
                        %--> Water would be pressed into saturated zone --> should not be the case
                        %--> Add waterPotential to gravitationalPotential so that flux into saturated zone is prevented(d_head = 0)
                        gravitationalPotential(i+1,1) = gravitationalPotential_unsaturated(i,1) + ground.STATVAR.waterPotential(i,1);
                        %Should not be the case if the saturated cell below is just a waterfront after precipitation that moves down
                        %--> Water accumulates in saturated cell and porosity becomes unnaturally large
                        %--> Check if there is an unsaturated cell below
                        for j = 1 : size(id_unsaturated)
                            if id_unsaturated(j,1) > i
                                gravitationalPotential(i+1,1) = gravitationalPotential_unsaturated(i+1,1);
                                break
                            end
                        end    
                    end
                elseif saturated(i,1) == 1
                    if size(gravitationalPotential,1) == i
                        gravitationalPotential(i+1,1) = gravitationalPotential(i,1);
                    else
                        gravitationalPotential(i,1) = gravitationalPotential_unsaturated(i,1);
                        gravitationalPotential(i+1,1) = gravitationalPotential(i,1);
                    end
                end
            end
            if size(gravitationalPotential,1) > size(saturated,1)
                gravitationalPotential(end) = [];
            end
            
            ground.STATVAR.gravitationalPotential = gravitationalPotential;
  
        end
        
        
        function rf = reduction_factor_out(saturation, ground)  %part of get_derivative_water2(ground)
            %smoothness = 3e-2;
            smoothness = 3e-1;
            rf = (1-exp(-saturation./smoothness));
        end
        
        function rf = reduction_factor_in(saturation, ground)   %part of get_derivative_water2(ground)
            %smoothness = 3e-2;
            smoothness = 3e-1;
            rf = (1- exp((saturation-1)./smoothness));

%             rf = double(saturation <= 0.9) + double (saturation > 0.9) .* (1 - saturation) ./0.1;
        end
        
        
        function timestep = get_timestep_water(ground)
            %outflow + inflow
             %timestep = ( double(ground.TEMP.d_water <0 & ground.STATVAR.waterIce > ground.STATVAR.field_capacity .* ground.STATVAR.layerThick .* ground.STATVAR.area) .* (ground.STATVAR.waterIce - ground.STATVAR.field_capacity .* ground.STATVAR.layerThick .* ground.STATVAR.area) ./ -ground.TEMP.d_water + ...
             %    double(ground.TEMP.d_water > 0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce ) ./ ground.TEMP.d_water); %[m3 / (m3/sec) = sec]
             timestep =  double(ground.TEMP.d_water <0) .* ground.STATVAR.water ./ -ground.TEMP.d_water ./10 + ...
                 double(ground.TEMP.d_water > 0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce ) ./ ground.TEMP.d_water; %[m3 / (m3/sec) = sec]
           
rand_factor = 1e-6 .* (2.*rand(size(ground.STATVAR.waterIce,1),1) -1);
             timestep = double(ground.TEMP.d_water <0 & ground.STATVAR.waterIce + rand_factor > ground.STATVAR.field_capacity .* ...
                 ground.STATVAR.layerThick .* ground.STATVAR.area) .* (ground.STATVAR.waterIce - ground.STATVAR.field_capacity .* ground.STATVAR.layerThick .* ground.STATVAR.area) ./ -ground.TEMP.d_water + ...
                 double(ground.TEMP.d_water <0 & ground.STATVAR.waterIce +  rand_factor <= ground.STATVAR.field_capacity .* ground.STATVAR.layerThick .* ground.STATVAR.area) .* ground.STATVAR.water ./ -ground.TEMP.d_water./10 + ...
                 double(ground.TEMP.d_water > 0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce ) ./ ground.TEMP.d_water; %[m3 / (m3/sec) = sec]

             timestep(timestep<=0) = ground.PARA.dt_max;
             timestep=nanmin(timestep);
        end
        
        function timestep = get_timestep_water_Xice(ground)
            %outflow + inflow
%              timestep = ( double(ground.TEMP.d_water <0 & ground.STATVAR.waterIce > ground.STATVAR.field_capacity .* ...
%                  (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce)) .* ...
%                  (ground.STATVAR.waterIce - ground.STATVAR.field_capacity .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce)) ./ -ground.TEMP.d_water + ...
%                  double(ground.TEMP.d_water > 0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce - ground.STATVAR.XwaterIce) ...
%                  ./ ground.TEMP.d_water); %[m3 / (m3/sec) = sec]
             
             timestep = ( double(ground.TEMP.d_water <0 & ground.STATVAR.waterIce + 1e-9 .* (2.*rand(size(ground.STATVAR.waterIce,1),1) -1) > ground.STATVAR.field_capacity .* ...
                 (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce)) .* ...
                 (ground.STATVAR.waterIce - ground.STATVAR.field_capacity .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce)) ./ -ground.TEMP.d_water + ...
                 double(ground.TEMP.d_water > 0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce - ground.STATVAR.XwaterIce) ...
                 ./ ground.TEMP.d_water); %[m3 / (m3/sec) = sec]

             
             timestep(timestep<=0) = ground.PARA.dt_max;
             %[mini,pos] = nanmin(timestep)
             %[timestep,posi] = nanmin(timestep);
             timestep = nanmin(timestep);
           
             
             %if timestep < 1e-12
%                  volume_matrix = ground.STATVAR.layerThick(posi) .* ground.STATVAR.area(posi) - ground.STATVAR.XwaterIce(posi);
%                  saturation = (ground.STATVAR.waterIce(posi)  - ground.STATVAR.field_capacity(posi) .* volume_matrix)./...
%                      (volume_matrix - ground.STATVAR.mineral(posi) - ground.STATVAR.organic(posi) - ground.STATVAR.field_capacity(posi) .* volume_matrix);
%                 saturation(saturation >= (1 - 1e-9)) = 1;
%                 disp(timestep)
%                 disp(posi)
%                 disp(ground.TEMP.d_water(posi))
%                 disp(ground.STATVAR.layerThick(posi) .* ground.STATVAR.area(posi) - ground.STATVAR.mineral(posi) - ground.STATVAR.organic(posi) - ground.STATVAR.waterIce(posi) - ground.STATVAR.XwaterIce(posi))
%                 disp(1-saturation)
%                 disp(reduction_factor_in(saturation, ground))
%                 disp(ground.TEMP.F_ub_water)
%                 
%                 volume_matrix = ground.STATVAR.layerThick(1) .* ground.STATVAR.area(1) - ground.STATVAR.XwaterIce(1);
%                 saturation_first_cell = (ground.STATVAR.waterIce(1)  - ground.STATVAR.field_capacity(1) .* volume_matrix)./...
%                     (volume_matrix - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.field_capacity(1) .* volume_matrix);
%                 saturation_first_cell(saturation_first_cell >= (1 - 1e-9)) = 1;
%                 disp(1-saturation_first_cell)
%                 disp('next')
             %end
             
        end
        
        function timestep = get_timestep_water_RichardsEq(ground)
             %no negative values and no overtopping
             timestep = ( double(ground.TEMP.d_water <0)  .* ground.STATVAR.water./2 ./ -ground.TEMP.d_water + ...
                 double(ground.TEMP.d_water > 0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce ) ./ ground.TEMP.d_water); %[m3 / (m3/sec) = sec]
             timestep(timestep<=0) = ground.PARA.dt_max;
             timestep=nanmin(timestep);
             
             timestep2 = ground.PARA.dWater_max .* ground.STATVAR.layerThick .* ground.STATVAR.area ./ abs(ground.TEMP.d_water);
             timestep2 = nanmin(timestep2);
             timestep = min(timestep, timestep2);
        end
        
        function timestep = get_timestep_water_RichardsEq_pressure(ground)
            %no negative values
            range1 = find(ground.TEMP.d_water > 0 & ground.STATVAR.saturation <= 1-1e-6);
            range2 = find(ground.TEMP.d_water > 0 & ground.STATVAR.saturation > 1-1e-6);
            range3 = find(ground.TEMP.d_water < 0 & (ground.STATVAR.saturation <= 1-1e-6 | (ground.STATVAR.saturation > 1-1e-6 & ~ground.TEMP.no_air)));
            range4 = find(ground.TEMP.d_water < 0 & ground.STATVAR.saturation > 1-1e-6 & ground.TEMP.no_air); %no air available, cell contracts
            range6 = find(ground.STATVAR.saturation > 1-1e-6);
            
            timestep = zeros(6,1).*NaN;
            if ~isempty(range1)
                timestep(1) = min((ground.STATVAR.layerThick(range1) .* ground.STATVAR.area(range1) - ground.STATVAR.mineral(range1) - ground.STATVAR.organic(range1) - ground.STATVAR.waterIce(range1) ) ./ ground.TEMP.d_water(range1));
            end
            if ~isempty(range3)
                timestep(3) = min(ground.STATVAR.water(range3)./4 ./ -ground.TEMP.d_water(range3));
            end
            timestep(5) = nanmin(ground.PARA.dWater_max .* ground.STATVAR.layerThick .* ground.STATVAR.area ./ abs(ground.TEMP.d_water));

            timestep(timestep<=0) = ground.PARA.dt_max;
            timestep=nanmin(timestep);          
        end
        
        %Discontiued
        function timestep = get_timestep_water_RichardsEq_Xice(ground)
             %no negative values and no overtopping
             timestep =  double(ground.TEMP.d_water <0 & ground.STATVAR.Xwater./ground.STATVAR.area<=1e-6)  .* ground.STATVAR.water./2 ./ -ground.TEMP.d_water + ...
                 double(ground.TEMP.d_water <0 & ground.STATVAR.Xwater./ground.STATVAR.area > 1e-6)  .* ground.STATVAR.Xwater./ -ground.TEMP.d_water + ...
                 double(ground.TEMP.d_water > 0 & ground.TEMP.XwaterIce_formation==0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce ) ./ ground.TEMP.d_water; %[m3 / (m3/sec) = sec]
             timestep(timestep<=0) = ground.PARA.dt_max;
             

             timestep=nanmin(timestep);
             timestep2 = ground.PARA.dWater_max .* ground.STATVAR.layerThick .* ground.STATVAR.area ./ abs(ground.TEMP.d_water);
             timestep2 = nanmin(timestep2);
             timestep = min(timestep, timestep2);
           
        end
        
        function timestep = get_timestep_water_RichardsEq_Xice2(ground)
             %no negative values and no overtopping
             timestep = ( double(ground.TEMP.d_water <0)  .* ground.STATVAR.water./2 ./ -ground.TEMP.d_water + ...
                 double(ground.TEMP.d_water > 0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce - ground.STATVAR.XwaterIce) ./ ground.TEMP.d_water); %[m3 / (m3/sec) = sec]
             timestep(timestep<=0) = ground.PARA.dt_max;
             timestep=nanmin(timestep);
             
             timestep2 = ground.PARA.dWater_max .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.XwaterIce) ./ abs(ground.TEMP.d_water);
             timestep2 = nanmin(timestep2);
             timestep = min(timestep, timestep2);
           
        end
        
        %bucketW and Xice 
        function timestep = get_timestep_Xwater(ground) %not needed anymore, since Xwater fluxes are restricted by max timestep
            %only outflow
            timestep = double(ground.TEMP.d_water <0) .* ground.STATVAR.Xwater ./ -ground.TEMP.d_water;
             
             timestep(timestep<=0) = ground.PARA.dt_max;
             timestep=nanmin(timestep);

             %if negative, set to max_timestep
        end
        
        function timestep = get_timestep_water_SNOW(ground)
            %outflow + inflow
%             remaining_pore_space = ground.STATVAR.layerThick.* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.ice;
%             remaining_pore_space = max(remaining_pore_space,0);
            %timestep = ( double(ground.TEMP.d_water <0) .* (ground.STATVAR.waterIce - ground.PARA.field_capacity .* remaining_pore_space) ./ -ground.TEMP.d_water + ...
            %     double(ground.TEMP.d_water > 0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce ) ./ ground.TEMP.d_water); %[m3 / (m3/sec) = sec]
%             timestep = ( double(ground.TEMP.d_water <0 & ground.STATVAR.water > ground.PARA.field_capacity .* remaining_pore_space ) .* (ground.STATVAR.water - ground.PARA.field_capacity .* remaining_pore_space) ./ -ground.TEMP.d_water + ...
%                 double(ground.TEMP.d_water > 0) .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce ) ./ ground.TEMP.d_water) + ...
%                 double(ground.TEMP.d_water <0 & ground.STATVAR.water <= ground.PARA.field_capacity .* remaining_pore_space ) .*ground.STATVAR.water ./ -ground.TEMP.d_water ./10; %[m3 / (m3/sec) = sec]
%             
            timestep = double(ground.TEMP.d_water <0) .* ground.STATVAR.water ./ -ground.TEMP.d_water ./10 + ...
                double(ground.TEMP.d_water > 0) .* double(1-(ground.STATVAR.mineral + ground.STATVAR.organic + ground.STATVAR.waterIce)./(ground.STATVAR.layerThick .* ground.STATVAR.area)>1e-11) .* ...
                (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce ) ./ ground.TEMP.d_water;

            timestep(timestep<=0) = ground.PARA.dt_max;
            timestep=nanmin(timestep);
            
        end
        
        %Hydraulic conductivity
        
        %bucketW
        function ground = calculate_hydraulicConductivity(ground)  
            saturation = ground.STATVAR.water ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            saturation = max(0,min(1,saturation));
            ice_saturation = ground.STATVAR.ice ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            ice_saturation = max(0,min(1,ice_saturation));
            ground.STATVAR.hydraulicConductivity = ground.PARA.hydraulicConductivity .* saturation .* 10.^(-7.*ice_saturation); %final term from dall amico     
        end
        
        %bucketW
        function ground = calculate_hydraulicConductivity_Xice(ground) %hydraulic conductivity of the matrix part of the cell
            saturation = ground.STATVAR.water ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.XwaterIce);
            saturation = max(0,min(1,saturation));
            ice_saturation = ground.STATVAR.ice ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.XwaterIce);
            ice_saturation = max(0,min(1,ice_saturation)); %count both ice and excess ice
            n = ground.STATVAR.n;
            %ground.STATVAR.hydraulicConductivity = ground.PARA.hydraulicConductivity .* saturation .* 10.^(-7.*ice_saturation); %final term from dall amico   
            ground.STATVAR.hydraulicConductivity = ground.STATVAR.satHydraulicConductivity .* saturation .* 10.^(-7.*ice_saturation); %final term from dall amico
        end

        %Richards equation
        function ground = calculate_hydraulicConductivity_RichardsEq(ground) 
            saturation = ground.STATVAR.water ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            saturation = max(0,min(1,saturation));
            %ice_saturation = ground.STATVAR.ice ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            ice_saturation = ground.STATVAR.ice ./ ground.STATVAR.waterIce; %Changed Sebastian Hansen et al., 2004
            ice_saturation = max(0,min(1,ice_saturation));
            n = ground.STATVAR.n;
            %ground.STATVAR.hydraulicConductivity = ground.PARA.hydraulicConductivity .* saturation.^0.5 .* (1 - (1 - saturation.^(n./(n+1))).^(1-1./n)).^2 .* 10.^(-7.*ice_saturation); %dall amico 

            ground = calculate_viscosity_water(ground);
            hydr_cond = ground.STATVAR.permeability ./ ground.STATVAR.viscosity_water .* ground.CONST.rho_w .* ground.CONST.g; 
            ground.STATVAR.hydraulicConductivity = hydr_cond .* saturation.^0.5 .* (1 - (1 - saturation.^(n./(n+1))).^(1-1./n)).^2 .* 10.^(-7.*ice_saturation); %dall amico 

            ground.STATVAR.hydraulicConductivity = min(2e-6, ground.STATVAR.hydraulicConductivity);
      
            %SEBAS CHANGED
           % ground.STATVAR.hydraulicConductivity(ground.STATVAR.T<0) = 0;
        end
        
        function ground = calculate_hydraulicConductivity_RichardsEq_Xice(ground) %hydraulic conductivity of the matrix part of the cell
            saturation = ground.STATVAR.water ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.XwaterIce);
            saturation = max(0,min(1,saturation));
            ice_saturation = ground.STATVAR.ice ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.XwaterIce);
            ice_saturation = max(0,min(1,ice_saturation)); %count both ice and excess ice
            n = ground.STATVAR.n;
            %ground.STATVAR.hydraulicConductivity = ground.PARA.hydraulicConductivity .* saturation.^0.5 .* (1 - (1 - saturation.^(n./(n+1))).^(1-1./n)).^2 .* 10.^(-7.*ice_saturation); %dall amico            
            ground.STATVAR.hydraulicConductivity = ground.STATVAR.satHydraulicConductivity .* saturation.^0.5 .* (1 - (1 - saturation.^(n./(n+1))).^(1-1./n)).^2 .* 10.^(-7.*ice_saturation); %dall amico            
 
        end
        
        function ground = calculate_hydraulicConductivity_RichardsEq_Xice2(ground) %hydraulic conductivity of the matrix part of the cell
            saturation = (ground.STATVAR.water + ground.STATVAR.Xwater)./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            saturation = max(0,min(1,saturation));
            ice_saturation = (ground.STATVAR.ice + ground.STATVAR.Xice)./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            ice_saturation = max(0,min(1,ice_saturation)); %count both ice and excess ice
            n = ground.STATVAR.n;

            ground = calculate_viscosity_water(ground);
            hydr_cond = ground.STATVAR.permeability ./ ground.STATVAR.viscosity_water .* ground.CONST.rho_w .* ground.CONST.g; 
            ground.STATVAR.hydraulicConductivity = hydr_cond .* saturation.^0.5 .* (1 - (1 - saturation.^(n./(n+1))).^(1-1./n)).^2 .* 10.^(-7.*ice_saturation); %dall amico 

        end
        
        %SNOW classes
        function ground = calculate_hydraulicConductivity_SNOW(ground)
            ground.STATVAR.hydraulicConductivity = ground.PARA.hydraulicConductivity .* ground.STATVAR.water./ground.STATVAR.layerThick./ground.STATVAR.area;   
        end
        
        %viscosity
        function ground = calculate_viscosity_water(ground)
            
            T = max(0,ground.STATVAR.T)+273.15;
            
            %from Wikipedia, https://en.wikipedia.org/wiki/Temperature_dependence_of_viscosity#cite_note-Reid1987-13
            A = 1.856e-14;
            B = 4209;
            C = 0.04527;
            D = -3.376e-5;

            ground.STATVAR.viscosity_water = A.*exp(B./T + C.* T + D.* T.^2);
        end
        
        function ground = calculate_hydraulicConductivity_RichardsEq2(ground)
            saturation = ground.STATVAR.water ./ (ground.STATVAR.layerThick.*ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            saturation = max(0,min(1,saturation));
            ice_saturation = ground.STATVAR.ice ./ ground.STATVAR.waterIce; %Changed Sebastian Hansen et al., 2004
            ice_saturation = max(0,min(1,ice_saturation));
            n = ground.STATVAR.n;

            ground = calculate_viscosity_water(ground);
            
            ground.STATVAR.porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area; 
            
%             ground.STATVAR.diamater_pipe = 2./3 .*ground.STATVAR.porosity ./(1-ground.STATVAR.porosity) .* ground.STATVAR.grain_size;
%             ground.STATVAR.number_of_pipes = ground.STATVAR.porosity ./ (pi()./4 .* ground.STATVAR.diamater_pipe.^2 .* ground.PARA.tortuosity_water); 
%             ground.STATVAR.number_of_pipes(isnan(ground.STATVAR.number_of_pipes)) = 0; %if diameter is zero
            
          %  ground.STATVAR.satHydraulicConductivity = pi ./ 4.*ground.STATVAR.diamater_pipe.^4 ./ 8 ./ ground.STATVAR.viscosity_water .* ground.STATVAR.number_of_pipes .* ground.CONST.rho_w .* ground.CONST.g;
            
            permeability = ground.STATVAR.grain_size.^2 ./ 180 .* ground.STATVAR.porosity .^3 ./ (1 - ground.STATVAR.porosity).^2; %Carman Kozeny model
            ground.STATVAR.satHydraulicConductivity = permeability./ ground.STATVAR.viscosity_water .*ground.CONST.rho_w .* ground.CONST.g; %m
            ground.STATVAR.hydraulicConductivity = ground.STATVAR.satHydraulicConductivity .* saturation.^0.5 .* (1 - (1 - saturation.^(n./(n+1))).^(1-1./n)).^2 .* 10.^(-7.*ice_saturation); %dall amico
        end
        
        %--------------VEGETATION------------------
        function ground = get_boundary_condition_water_SNOW_canopy_m(ground, tile) % as function above, but for snow below canopy
            forcing = tile.FORCING;
            rainfall = ground.PARENT.PREVIOUS.TEMP.rain_thru .* ground.STATVAR.area(1);  
            
            %partition already here in infiltration and surface runoff,
            %considering ET losses and potentially external fluxes
            remaining_pore_space = ground.STATVAR.layerThick(1).* ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.ice(1);
            saturation_first_cell = (ground.STATVAR.waterIce(1) - ground.PARA.field_capacity .* remaining_pore_space) ./ ...
                (ground.STATVAR.layerThick(1).*ground.STATVAR.area(1) - remaining_pore_space); 
            saturation_first_cell = max(0,min(1,saturation_first_cell)); % 0 water at field capacity, 1: water at saturation
            %NEW SW
            saturation_first_cell(saturation_first_cell >= (1 - 1e-9)) = 1;
            
            ground.TEMP.F_ub_water = rainfall .* reduction_factor_in(saturation_first_cell, ground);
            ground.TEMP.surface_runoff = rainfall - ground.TEMP.F_ub_water;  %route this to surface pool
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
    end
end

