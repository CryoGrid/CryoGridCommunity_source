%========================================================================
% CryoGrid TIER1 INTERACTION (IA) class for functions related to water fluxes
% S. Westermann, October 2020
%========================================================================

classdef IA_WATER < IA_BASE
    
    methods
        %water fluxes between two normal GROUND classes with bucket water scheme
        function get_boundary_condition_BUCKET_m(ia_heat_water) %water fluxes between classes with bucket water scheme
            saturation_previous = (ia_heat_water.PREVIOUS.STATVAR.waterIce(end) - ia_heat_water.PREVIOUS.STATVAR.field_capacity(end) .* ia_heat_water.PREVIOUS.STATVAR.layerThick(end).*ia_heat_water.PREVIOUS.STATVAR.area(end))./ ...
                (ia_heat_water.PREVIOUS.STATVAR.layerThick(end).*ia_heat_water.PREVIOUS.STATVAR.area(end) - ia_heat_water.PREVIOUS.STATVAR.mineral(end) - ia_heat_water.PREVIOUS.STATVAR.organic(end) - ...
                ia_heat_water.PREVIOUS.STATVAR.field_capacity(end).*ia_heat_water.PREVIOUS.STATVAR.layerThick(end).*ia_heat_water.PREVIOUS.STATVAR.area(end));
            saturation_previous = max(0,min(1,saturation_previous)); % 0 water at field capacity, 1: water at saturation

            saturation_next = (ia_heat_water.NEXT.STATVAR.waterIce(1) - ia_heat_water.NEXT.STATVAR.field_capacity(1) .* ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1))./ ...
                (ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.mineral(1) - ia_heat_water.NEXT.STATVAR.organic(1) - ...
                ia_heat_water.NEXT.STATVAR.field_capacity(1).*ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1));
            saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            
            %outflow
            d_water_out = ia_heat_water.PREVIOUS.STATVAR.hydraulicConductivity(end) .* ia_heat_water.PREVIOUS.STATVAR.area(end); 
            d_water_out = d_water_out .* reduction_factor_out(saturation_previous, ia_heat_water); %this is positive when flowing out
            
            %inflow
            d_water_in = d_water_out .* reduction_factor_in(saturation_next, ia_heat_water);
               
            %readjust outflow
            d_water_out = d_water_in; %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* (double(ia_heat_water.PREVIOUS.STATVAR.T(end)>=0) .* ia_heat_water.PREVIOUS.CONST.c_w + ...
                double(ia_heat_water.PREVIOUS.STATVAR.T(end)<0) .* ia_heat_water.PREVIOUS.CONST.c_i).* ia_heat_water.PREVIOUS.STATVAR.T(end);
            d_water_in_energy = d_water_out_energy;            
            
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_wtaer_out is positive
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
            
            ia_heat_water.NEXT.TEMP.F_ub_water = d_water_in;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
            
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
            
            %no evaporation losses from the NEXT class
            %ia_heat_water.NEXT.TEMP.d_water_ET = 0;
            %ia_heat_water.NEXT.TEMP.d_water_ET_energy = 0;
        end
        
        %zero flux boundary condition
        function get_boundary_condition_ZEROFLUX_NEXT_m(ia_heat_water) %coupling between classes without (PREVIOUS) and with (NEXT) water balance
            ia_heat_water.NEXT.TEMP.F_ub_water = 0;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = 0;
            %ia_heat_water.NEXT.TEMP.d_water_ET = 0;
            %ia_heat_water.NEXT.TEMP.d_water_ET_energy = 0;
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
        end
        
        %zero flux boundary condition
        function get_boundary_condition_ZEROFLUX_PREVIOUS_m(ia_heat_water) %coupling between classes with (PREVIOUS) and without (NEXT) water balance
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = 0;
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = 0;
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
        end
        
        %water fluxes between LAKE and GROUND class with bucket water scheme
        function get_boundary_condition_BUCKET_LAKE_m(ia_heat_water)
            
            saturation_next = (ia_heat_water.NEXT.STATVAR.waterIce(1) - ia_heat_water.NEXT.STATVAR.field_capacity(1) .* ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1))./ ...
                (ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.mineral(1) - ia_heat_water.NEXT.STATVAR.organic(1) - ...
                ia_heat_water.NEXT.STATVAR.field_capacity(1).*ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1));
            saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            
            %outflow
            d_water_out = ia_heat_water.NEXT.STATVAR.hydraulicConductivity(end) .* ia_heat_water.PREVIOUS.STATVAR.area(end);
           
            %inflow
            d_water_in = d_water_out .* reduction_factor_in(saturation_next, ia_heat_water);
            
            %readjust outflow
            d_water_out = d_water_in; %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
            d_water_in_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_wtaer_out is positive
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
            
            ia_heat_water.NEXT.TEMP.F_ub_water = d_water_in;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
            
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
            
            %no evaporation losses from the NEXT class
            %ia_heat_water.NEXT.TEMP.d_water_ET = 0;
            %ia_heat_water.NEXT.TEMP.d_water_ET_energy = 0;
            
        end
        
        %water fluxes between SNOW and GROUND class with bucket water scheme
        function get_boundary_condition_BUCKET_SNOW_m(ia_heat_water) %water fluxes between classes with bucket water scheme
            remaining_pore_space = ia_heat_water.PREVIOUS.STATVAR.layerThick(end).* ia_heat_water.PREVIOUS.STATVAR.area(end) - ia_heat_water.PREVIOUS.STATVAR.mineral(end) - ia_heat_water.PREVIOUS.STATVAR.organic(end) - ia_heat_water.PREVIOUS.STATVAR.ice(end);
            saturation_previous = (ia_heat_water.PREVIOUS.STATVAR.water(end) - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space) ./ ...
                (remaining_pore_space - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space); 
            
            saturation_previous = max(0,min(1,saturation_previous)); % 0 water at field capacity, 1: water at saturation
            

            saturation_next = (ia_heat_water.NEXT.STATVAR.waterIce(1) - ia_heat_water.NEXT.STATVAR.field_capacity(1) .* ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1))./ ...
                (ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.mineral(1) - ia_heat_water.NEXT.STATVAR.organic(1) - ...
                ia_heat_water.NEXT.STATVAR.field_capacity(1).*ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1));
            saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            saturation_next(saturation_next >= (1 - 1e-9)) = 1;
            
            %outflow
            d_water_out = ia_heat_water.PREVIOUS.STATVAR.hydraulicConductivity(end) .* ia_heat_water.PREVIOUS.STATVAR.area(end); 
            d_water_out = d_water_out .* reduction_factor_out(saturation_previous, ia_heat_water); %this is positive when flowing out
            
            %inflow
            d_water_in = d_water_out .* reduction_factor_in(saturation_next, ia_heat_water);
               
            %readjust outflow
            d_water_out = d_water_in; %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
            d_water_in_energy = d_water_out_energy;            
            
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_wtaer_out is positive
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
            
            ia_heat_water.NEXT.TEMP.F_ub_water = d_water_in;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
            
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
            
            %no evaporation losses from the NEXT class
            %ia_heat_water.NEXT.TEMP.d_water_ET = 0;
            %ia_heat_water.NEXT.TEMP.d_water_ET_energy = 0;
        end
        
        %water fluxes between SNOW and excess ice GROUND class with bucket water scheme
        function get_boundary_condition_BUCKET_SNOW_XICE_m(ia_heat_water) %snow and classes with Xice scheme
            remaining_pore_space = ia_heat_water.PREVIOUS.STATVAR.layerThick(end).* ia_heat_water.PREVIOUS.STATVAR.area(end) - ia_heat_water.PREVIOUS.STATVAR.mineral(end) - ia_heat_water.PREVIOUS.STATVAR.organic(end) - ia_heat_water.PREVIOUS.STATVAR.ice(end);
            saturation_previous = (ia_heat_water.PREVIOUS.STATVAR.water(end) - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space) ./ ...
                (remaining_pore_space - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space); 
            saturation_previous = max(0,min(1,saturation_previous)); % 0 water at field capacity, 1: water at saturation
            

            volume_matrix = ia_heat_water.NEXT.STATVAR.layerThick(1) .* ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.XwaterIce(1);
            saturation_next = (ia_heat_water.NEXT.STATVAR.waterIce(1)  - ia_heat_water.NEXT.STATVAR.field_capacity(1) .* volume_matrix)./...
                (volume_matrix - ia_heat_water.NEXT.STATVAR.mineral(1) - ia_heat_water.NEXT.STATVAR.organic(1) - ia_heat_water.NEXT.STATVAR.field_capacity(1) .* volume_matrix);
            saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            saturation_next(saturation_next >= (1 - 1e-9)) = 1;
            
            %outflow
            d_water_out = ia_heat_water.PREVIOUS.STATVAR.hydraulicConductivity(end) .* ia_heat_water.PREVIOUS.STATVAR.area(end); 
            d_water_out = d_water_out .* reduction_factor_out(saturation_previous, ia_heat_water); %this is positive when flowing out
            
            %inflow
            d_water_in = d_water_out .* reduction_factor_in(saturation_next, ia_heat_water);
               
            %readjust outflow
            d_water_out = d_water_in; %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
            d_water_in_energy = d_water_out_energy;            
            
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_wtaer_out is positive
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
            
            ia_heat_water.NEXT.TEMP.F_ub_water = d_water_in;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
            
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
            
            %no evaporation losses from the NEXT class
            %ia_heat_water.NEXT.TEMP.d_water_ET = 0;
            %ia_heat_water.NEXT.TEMP.d_water_ET_energy = 0;
        end
        
        %water fluxes between LAKE and excess ice GROUND class with bucket water scheme
        function get_boundary_condition_BUCKET_LAKE_XICE_m(ia_heat_water) %LAKE to Xice classes, routes "normal" water down
            
            volume_matrix = ia_heat_water.NEXT.STATVAR.layerThick(1) .* ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.XwaterIce(1);
            saturation_next = (ia_heat_water.NEXT.STATVAR.waterIce(1)  - ia_heat_water.NEXT.STATVAR.field_capacity(1) .* volume_matrix)./...
                (volume_matrix - ia_heat_water.NEXT.STATVAR.mineral(1) - ia_heat_water.NEXT.STATVAR.organic(1) - ia_heat_water.NEXT.STATVAR.field_capacity(1) .* volume_matrix);
            saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            saturation_next(saturation_next >= (1 - 1e-9)) = 1;
            
            %outflow
            d_water_out = ia_heat_water.NEXT.STATVAR.hydraulicConductivity(1) .* ia_heat_water.PREVIOUS.STATVAR.area(end);
           
            %inflow
            d_water_in = d_water_out .* reduction_factor_in(saturation_next, ia_heat_water);
            
            %readjust outflow
            d_water_out = d_water_in; %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
            d_water_in_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_wtaer_out is positive
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
            
            ia_heat_water.NEXT.TEMP.F_ub_water = d_water_in;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
            
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
            
        end
        
        %upwards excess water fluxes between LAKE and excess ice GROUND class with bucket water scheme
        function get_boundary_condition_BUCKET_LAKE_XWATER_UP_m(ia_heat_water) %LAKE to Xice classes, routes Xwater flux up into lake

            d_Xwater_out = ia_heat_water.NEXT.STATVAR.hydraulicConductivity(1) .* ia_heat_water.NEXT.STATVAR.area(1);             
            d_Xwater_out = min(d_Xwater_out, ia_heat_water.NEXT.STATVAR.Xwater(1) ./ ia_heat_water.NEXT.PARA.dt_max); %makes explicit timestep check unnecessary
            d_Xwater_out_energy = d_Xwater_out .* ia_heat_water.NEXT.CONST.c_w .* ia_heat_water.NEXT.STATVAR.T(1);
            
            ia_heat_water.NEXT.TEMP.F_ub_Xwater = -d_Xwater_out;
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = d_Xwater_out; %this might overwrite the variable set in the "normal" water flux BC, but no problem 
            
            ia_heat_water.NEXT.TEMP.F_ub_Xwater_energy = -d_Xwater_out_energy;
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = d_Xwater_out_energy;
            
            ia_heat_water.NEXT.TEMP.d_Xwater(1) = ia_heat_water.NEXT.TEMP.d_Xwater(1) + ia_heat_water.NEXT.TEMP.F_ub_Xwater;
            ia_heat_water.NEXT.TEMP.d_Xwater_energy(1) = ia_heat_water.NEXT.TEMP.d_Xwater_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_Xwater_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;

        end
        
        %water fluxes between SNOW and LAKE class with bucket water scheme
        function get_boundary_condition_BUCKET_SNOW_LAKE_m(ia_heat_water) %water fluxes between classes with bucket water scheme
            remaining_pore_space = ia_heat_water.PREVIOUS.STATVAR.layerThick(end).* ia_heat_water.PREVIOUS.STATVAR.area(end) - ia_heat_water.PREVIOUS.STATVAR.mineral(end) - ia_heat_water.PREVIOUS.STATVAR.organic(end) - ia_heat_water.PREVIOUS.STATVAR.ice(end);
            saturation_previous = (ia_heat_water.PREVIOUS.STATVAR.water(end) - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space) ./ ...
                (remaining_pore_space - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space); 
            
            saturation_previous = max(0,min(1,saturation_previous)); % 0 water at field capacity, 1: water at saturation
            
            %check if water can penetrate through ice, and get first (almost) unfrozen cell
            i=1;
            first_water_cell_reached = 0;
            while i<= size(ia_heat_water.NEXT.STATVAR.layerThick,1) && ia_heat_water.NEXT.STATVAR.water(i,1) > 0 && ~first_water_cell_reached
                if ia_heat_water.NEXT.STATVAR.water(i,1) >= 0.9 .* ia_heat_water.NEXT.STATVAR.layerThick(i,1) %cell more or less melted
                    first_water_cell_reached =1;
                else
                    i=i+1;
                end
            end
            
            if first_water_cell_reached
                
                %outflow
                d_water_out = ia_heat_water.PREVIOUS.STATVAR.hydraulicConductivity(end) .* ia_heat_water.PREVIOUS.STATVAR.area(end);
                d_water_out = d_water_out .* reduction_factor_out(saturation_previous, ia_heat_water); %this is positive when flowing out
                
                
                %energy advection
                d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
                
                ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_water_out is positive
                ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
                
                ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
                ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
                
                % put water in cell i belwo the ice
                ia_heat_water.NEXT.TEMP.d_water(i) = ia_heat_water.NEXT.TEMP.d_water(i) + d_water_out;
                ia_heat_water.NEXT.TEMP.d_water_energy(i) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + d_water_out_energy;
                
            end
        end
           
        %water fluxes between LAKE and LAKE class with bucket water scheme
        function get_boundary_condition_BUCKET_LAKE_LAKE_m(ia_heat_water) %water fluxes between classes with bucket water scheme
            
            %check if water can penetrate through ice, and get first (almost) unfrozen cell
            i=1;
            first_water_cell_reached = 0;
            while i<= size(ia_heat_water.NEXT.STATVAR.layerThick,1) && ia_heat_water.NEXT.STATVAR.water(i,1) > 0 && ~first_water_cell_reached
                if ia_heat_water.NEXT.STATVAR.water(i,1) >= 0.9 .* ia_heat_water.NEXT.STATVAR.layerThick(i,1) %cell more or less melted
                    first_water_cell_reached =1;
                else
                    i=i+1;
                end
            end
            
            if first_water_cell_reached
                
                %outflow
                d_water_out = ia_heat_water.PREVIOUS.STATVAR.water(end)./ ia_heat_water.PREVIOUS.PARA.dt_max; % remove as much as possible
                
                %energy advection
                d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
                
                ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_water_out is positive
                ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
                
                ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
                ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
                
                % put water in cell i belwo the ice
                ia_heat_water.NEXT.TEMP.d_water(i) = ia_heat_water.NEXT.TEMP.d_water(i) + d_water_out;
                ia_heat_water.NEXT.TEMP.d_water_energy(i) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + d_water_out_energy;
                
            end
        end
        
        %----------Richards Equation----------------------
        
        %water fluxes between LAKE and GROUND class with Richards equation
        function get_boundary_condition_RichardsEq_LAKE_m(ia_heat_water) %should be done
            
%             saturation_next = (ia_heat_water.NEXT.STATVAR.waterIce(1) - ia_heat_water.NEXT.STATVAR.field_capacity(1) .* ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1))./ ...
%                 (ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.mineral(1) - ia_heat_water.NEXT.STATVAR.organic(1) - ...
%                 ia_heat_water.NEXT.STATVAR.field_capacity(1).*ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1));
%             saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            
            saturation_next = ia_heat_water.NEXT.STATVAR.waterIce(1,1) ./ (ia_heat_water.NEXT.STATVAR.layerThick(1,1).*ia_heat_water.NEXT.STATVAR.area(1,1) - ia_heat_water.NEXT.STATVAR.mineral(1,1) - ia_heat_water.NEXT.STATVAR.organic(1,1));
            saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            saturation_next(saturation_next >= (1 - 1e-9)) = 1;
            
            %outflow
            %d_water_out = ia_heat_water.NEXT.PARA.hydraulicConductivity .* ia_heat_water.PREVIOUS.STATVAR.water(end) ./ ia_heat_water.PREVIOUS.STATVAR.layerThick(end); % area cancels out; make this depended on both involved cells?
            
            d_water_out = ia_heat_water.NEXT.STATVAR.hydraulicConductivity(1,1) .* ground.STATVAR.area(1,1);
            %inflow
            d_water_in = d_water_out .* reduction_factor_in(saturation_next, ia_heat_water);
            
            %readjust outflow
            d_water_out = d_water_in; %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
            d_water_in_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_wtaer_out is positive
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
            
            ia_heat_water.NEXT.TEMP.F_ub_water = d_water_in;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
            
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
            
            %no evaporation losses from the NEXT class
            %ia_heat_water.NEXT.TEMP.d_water_ET = 0;
            %ia_heat_water.NEXT.TEMP.d_water_ET_energy = 0;
            
        end
        
        %water fluxes between SNOW class (bucket water) and GROUND class with Richards equation
        function get_boundary_condition_RichardsEq_SNOW_m(ia_heat_water) %water fluxes between classes with bucket water scheme
            remaining_pore_space = ia_heat_water.PREVIOUS.STATVAR.layerThick(end).* ia_heat_water.PREVIOUS.STATVAR.area(end) - ia_heat_water.PREVIOUS.STATVAR.mineral(end) - ia_heat_water.PREVIOUS.STATVAR.organic(end) - ia_heat_water.PREVIOUS.STATVAR.ice(end);
            saturation_previous = (ia_heat_water.PREVIOUS.STATVAR.water(end) - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space) ./ ...
                (remaining_pore_space - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space); 
            saturation_previous = max(0,min(1,saturation_previous)); % 0 water at field capacity, 1: water at saturation
            
            saturation_next = ia_heat_water.NEXT.STATVAR.waterIce(1,1) ./ (ia_heat_water.NEXT.STATVAR.layerThick(1,1).*ia_heat_water.NEXT.STATVAR.area(1,1) - ia_heat_water.NEXT.STATVAR.mineral(1,1) - ia_heat_water.NEXT.STATVAR.organic(1,1));
            saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            saturation_next(saturation_next >= (1 - 1e-9)) = 1;
            
            %outflow
            d_water_out = ia_heat_water.PREVIOUS.STATVAR.hydraulicConductivity(end) .* ia_heat_water.PREVIOUS.STATVAR.area(end); 
            d_water_out = d_water_out .* reduction_factor_out(saturation_previous, ia_heat_water); %this is positive when flowing out
            
            %inflow
            d_water_in = d_water_out .* reduction_factor_in(saturation_next, ia_heat_water);
               
            %readjust outflow
            d_water_out = d_water_in; %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
            d_water_in_energy = d_water_out_energy;            
            
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_wtaer_out is positive
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
            
            ia_heat_water.NEXT.TEMP.F_ub_water = d_water_in;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
            
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
            
        end
        
        %%water fluxes between LAKE and excess ice GROUND class with Richards equation
        function get_boundary_condition_RichardsEq_Xice_LAKE_m(ia_heat_water) 
            
%             saturation_next = (ia_heat_water.NEXT.STATVAR.waterIce(1) - ia_heat_water.NEXT.STATVAR.field_capacity(1) .* ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1))./ ...
%                 (ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.mineral(1) - ia_heat_water.NEXT.STATVAR.organic(1) - ...
%                 ia_heat_water.NEXT.STATVAR.field_capacity(1).*ia_heat_water.NEXT.STATVAR.layerThick(1).*ia_heat_water.NEXT.STATVAR.area(1));
%             saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            
            %downward flow
            saturation_next = ia_heat_water.NEXT.STATVAR.waterIce(1,1) ./ (ia_heat_water.NEXT.STATVAR.layerThick(1,1).*ia_heat_water.NEXT.STATVAR.area(1,1) ...
                - ia_heat_water.NEXT.STATVAR.mineral(1,1) - ia_heat_water.NEXT.STATVAR.organic(1,1) - ia_heat_water.NEXT.STATVAR.XwaterIce(1,1));
            saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            saturation_next(saturation_next >= (1 - 1e-9)) = 1;
            
            %outflow
            %d_water_out = ia_heat_water.NEXT.PARA.hydraulicConductivity .* ia_heat_water.PREVIOUS.STATVAR.water(end) ./ ia_heat_water.PREVIOUS.STATVAR.layerThick(end); % area cancels out; make this depended on both involved cells?
            d_water_out = ia_heat_water.NEXT.STATVAR.hydraulicConductivity(1,1) .* ia_heat_water.NEXT.STATVAR.area(1,1);
            %inflow
            d_water_in = d_water_out .* reduction_factor_in(saturation_next, ia_heat_water);
            
            %readjust outflow
            d_water_out = d_water_in; %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
            d_water_in_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_wtaer_out is positive
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
            
            ia_heat_water.NEXT.TEMP.F_ub_water = d_water_in;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
            
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
            
            if ia_heat_water.NEXT.STATVAR.Xwater(1) > 1e-6 .* ia_heat_water.NEXT.STATVAR.area(1) %Xwater moving up, 1e-6 corresponds to threshold in get_timestep.
                d_water_up = ia_heat_water.NEXT.STATVAR.hydraulicConductivity(1,1) .* ia_heat_water.NEXT.STATVAR.area(1,1); %does not include overburden pressure, but is right order of magnitude
                d_water_up_energy = d_water_up .* ia_heat_water.NEXT.CONST.c_w .* ia_heat_water.NEXT.STATVAR.T(1);
                
                ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + d_water_up;
                ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + d_water_up_energy;
                
                ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) - d_water_up;
                ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) - d_water_up_energy;
            end
            
        end
        
        %water fluxes between SNOW class (bucket water) and excess ice GROUND class with Richards equation
        function get_boundary_condition_RichardsEq_Xice_SNOW_m(ia_heat_water) %water fluxes between classes with bucket water scheme
            remaining_pore_space = ia_heat_water.PREVIOUS.STATVAR.layerThick(end).* ia_heat_water.PREVIOUS.STATVAR.area(end) - ia_heat_water.PREVIOUS.STATVAR.mineral(end) - ia_heat_water.PREVIOUS.STATVAR.organic(end) - ia_heat_water.PREVIOUS.STATVAR.ice(end);
            saturation_previous = (ia_heat_water.PREVIOUS.STATVAR.water(end) - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space) ./ ...
                (remaining_pore_space - ia_heat_water.PREVIOUS.PARA.field_capacity .* remaining_pore_space); 
            saturation_previous = max(0,min(1,saturation_previous)); % 0 water at field capacity, 1: water at saturation
            
            
            saturation_next = ia_heat_water.NEXT.STATVAR.waterIce(1,1) ./ (ia_heat_water.NEXT.STATVAR.layerThick(1,1).*ia_heat_water.NEXT.STATVAR.area(1,1) ...
                - ia_heat_water.NEXT.STATVAR.mineral(1,1) - ia_heat_water.NEXT.STATVAR.organic(1,1) - ia_heat_water.NEXT.STATVAR.XwaterIce(1,1));
            %saturation_next = ia_heat_water.NEXT.STATVAR.waterIce(1,1) ./ (ia_heat_water.NEXT.STATVAR.layerThick(1,1).*ia_heat_water.NEXT.STATVAR.area(1,1) - ia_heat_water.NEXT.STATVAR.mineral(1,1) - ia_heat_water.NEXT.STATVAR.organic(1,1));
            saturation_next = max(0,min(1,saturation_next)); % 0 water at field capacity, 1: water at saturation
            saturation_next(saturation_next >= (1 - 1e-9)) = 1;
            
            %outflow
            d_water_out = ia_heat_water.PREVIOUS.STATVAR.hydraulicConductivity(end) .* ia_heat_water.PREVIOUS.STATVAR.area(end);
            d_water_out = d_water_out .* reduction_factor_out(saturation_previous, ia_heat_water); %this is positive when flowing out
            
            %inflow
            d_water_in = d_water_out .* reduction_factor_in(saturation_next, ia_heat_water);
               
            %readjust outflow
            d_water_out = d_water_in; %reduce outflow if inflow is impossible
            
            %energy advection
            d_water_out_energy = d_water_out .* ia_heat_water.PREVIOUS.CONST.c_w .* ia_heat_water.PREVIOUS.STATVAR.T(end);
            d_water_in_energy = d_water_out_energy;            
            
            ia_heat_water.PREVIOUS.TEMP.F_lb_water = - d_water_out;  %negative as d_wtaer_out is positive
            ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy = - d_water_out_energy;
            
            ia_heat_water.NEXT.TEMP.F_ub_water = d_water_in;
            ia_heat_water.NEXT.TEMP.F_ub_water_energy = d_water_out_energy;
            
            ia_heat_water.PREVIOUS.TEMP.d_water(end) = ia_heat_water.PREVIOUS.TEMP.d_water(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water;
            ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) = ia_heat_water.PREVIOUS.TEMP.d_water_energy(end) + ia_heat_water.PREVIOUS.TEMP.F_lb_water_energy;
            
            ia_heat_water.NEXT.TEMP.d_water(1) = ia_heat_water.NEXT.TEMP.d_water(1) + ia_heat_water.NEXT.TEMP.F_ub_water;
            ia_heat_water.NEXT.TEMP.d_water_energy(1) = ia_heat_water.NEXT.TEMP.d_water_energy(1) + ia_heat_water.NEXT.TEMP.F_ub_water_energy;
            
        end
        
        %---------- Vegetation ----------------------
        
        function ia_seb_water = get_boundary_condition_RichardsEq_canopy_m(ia_seb_water, tile)
            % Equivalent to get_boundary_condition_u_RichardsEq(...) in WATER_FLUXES
            forcing = tile.FORCING;
            
            max_infiltration = max(0, ia_seb_water.NEXT.STATVAR.hydraulicConductivity(1,1).* ((0 - ia_seb_water.NEXT.STATVAR.waterPotential(1,1)) ./ (ia_seb_water.NEXT.STATVAR.layerThick(1,1) ./ 2) + 1) .* ia_seb_water.NEXT.STATVAR.area(1,1));
            
            rainfall = ia_seb_water.PREVIOUS.TEMP.rain_thru;  % throughfall from canopy

            %partition already here in infiltration and surface runoff,
            %considering ET losses and potentially external fluxes
            saturation_first_cell = ia_seb_water.NEXT.STATVAR.waterIce(1)./ (ia_seb_water.NEXT.STATVAR.layerThick(1).*ia_seb_water.NEXT.STATVAR.area(1) - ia_seb_water.NEXT.STATVAR.mineral(1) - ia_seb_water.NEXT.STATVAR.organic(1));
            saturation_first_cell = max(0,min(1,saturation_first_cell)); % 0 water at field capacity, 1: water at saturation
            
            evap = double(ia_seb_water.NEXT.TEMP.d_water_ET(1)<0).*ia_seb_water.NEXT.TEMP.d_water_ET(1);
            condensation = double(ia_seb_water.NEXT.TEMP.d_water_ET(1)>0).*ia_seb_water.NEXT.TEMP.d_water_ET(1);
            
            rainfall = rainfall + condensation; %add condensation to rainfall to avoid overflowing of grid cell
            excessRain = max(0, rainfall-max_infiltration);
            rainfall = min(rainfall, max_infiltration);
            
            ia_seb_water.NEXT.TEMP.d_water_ET(1) = evap; %evaporation (water loss) subrtacted in get_derivative
            
            ia_seb_water.NEXT.TEMP.F_ub_water = double(rainfall <= -evap) .* rainfall + ...
                double(rainfall > -evap) .* (-evap + (rainfall + evap) .* reduction_factor_in(saturation_first_cell, ia_seb_water));
            ia_seb_water.NEXT.TEMP.surface_runoff = rainfall - ia_seb_water.NEXT.TEMP.F_ub_water + excessRain;
            
            ia_seb_water.NEXT.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ia_seb_water.NEXT.TEMP.F_ub_water_energy = ia_seb_water.NEXT.TEMP.F_ub_water .* ia_seb_water.NEXT.CONST.c_w .* ia_seb_water.NEXT.TEMP.T_rainWater;
            
            ia_seb_water.NEXT.TEMP.d_water(1) = ia_seb_water.NEXT.TEMP.d_water(1) + ia_seb_water.NEXT.TEMP.F_ub_water;
            ia_seb_water.NEXT.TEMP.d_water_energy(1) = ia_seb_water.NEXT.TEMP.d_water_energy(1) + ia_seb_water.NEXT.TEMP.F_ub_water_energy;
            
        end
        
        function ia_seb_water = get_boundary_condition_water_canopy_SNOW_m(ia_seb_water, tile)
            % Equivalent to get_boundary_condition_u_water_SNOW(...) in WATER_FLUXES
            ground = ia_seb_water.NEXT;
            forcing = tile.FORCING;
            
            rainfall = ia_seb_water.PREVIOUS.TEMP.rain_thru;  % throughfall from canopy
            
            %partition already here in infiltration and surface runoff,
            %considering ET losses and potentially external fluxes
            remaining_pore_space = ground.STATVAR.layerThick(1).* ground.STATVAR.area(1) - ground.STATVAR.mineral(1) - ground.STATVAR.organic(1) - ground.STATVAR.ice(1);
            saturation_first_cell = (ground.STATVAR.waterIce(1) - ground.PARA.field_capacity .* remaining_pore_space) ./ ...
                (ground.STATVAR.layerThick(1).*ground.STATVAR.area(1) - remaining_pore_space); 
            saturation_first_cell = max(0,min(1,saturation_first_cell)); % 0 water at field capacity, 1: water at saturation
            saturation_first_cell(saturation_first_cell >= (1 - 1e-9)) = 1;
            
            ground.TEMP.F_ub_water = rainfall .* reduction_factor_in(saturation_first_cell, ground);
            ground.TEMP.surface_runoff = rainfall - ground.TEMP.F_ub_water;  %route this to surface pool
            
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
            ground.TEMP.F_ub_water_energy = ground.TEMP.F_ub_water .* ground.CONST.c_w .* ground.TEMP.T_rainWater;
            
            ground.TEMP.d_water(1) = ground.TEMP.d_water(1) + ground.TEMP.F_ub_water;
            ground.TEMP.d_water_energy(1) = ground.TEMP.d_water_energy(1) + ground.TEMP.F_ub_water_energy;
        end
        
        function ia_seb_water = get_water_transpiration(ia_seb_water)
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % soil
            transp_water = stratigraphy1.TEMP.transp;
            f_root = stratigraphy2.STATVAR.f_root;
            psi = stratigraphy2.STATVAR.waterPotential;
            psi_wilt = stratigraphy1.PARA.psi_wilt;
            
            water_out = transp_water .* f_root.*max(0,(psi_wilt-psi)./psi_wilt) ./sum(f_root.*max(0,(psi_wilt-psi)./psi_wilt));
            water_out(isnan(water_out)) = 0;
            water_out_energy = water_out .* stratigraphy2.STATVAR.T .* (double(stratigraphy2.STATVAR.T>=0).*stratigraphy2.CONST.c_w + double(stratigraphy2.STATVAR.T<0).*stratigraphy2.CONST.c_i);
            
            stratigraphy2.TEMP.d_water_ET = stratigraphy2.TEMP.d_water_ET - water_out;
            stratigraphy2.TEMP.d_water_ET_energy = stratigraphy2.TEMP.d_water_ET_energy - water_out_energy;
        end
        
        function ia_seb_water = canopy_drip(ia_seb_water, tile)
            stratigraphy1 = ia_seb_water.PREVIOUS; %canopy
            stratigraphy2 = ia_seb_water.NEXT; %ground
            
            water_capacity = stratigraphy1.PARA.Wmax*stratigraphy1.STATVAR.area*(stratigraphy1.STATVAR.LAI+stratigraphy1.STATVAR.SAI);
            if stratigraphy1.STATVAR.waterIce > water_capacity
                water_fraction = stratigraphy1.STATVAR.water./stratigraphy1.STATVAR.waterIce;
                ice_fraction = stratigraphy1.STATVAR.ice./stratigraphy1.STATVAR.waterIce;
                excess_waterIce = max(0,stratigraphy1.STATVAR.waterIce - water_capacity);
                excess_water = excess_waterIce.*water_fraction;
                excess_ice = excess_waterIce.*ice_fraction;
                excess_water_energy = excess_water.*stratigraphy1.CONST.c_w.*stratigraphy1.STATVAR.T(1);
                excess_ice_energy = excess_ice.*(stratigraphy1.CONST.c_i.*stratigraphy1.STATVAR.T(1)-stratigraphy2.CONST.L_f);
                
                stratigraphy1.STATVAR.waterIce = water_capacity;
                stratigraphy1.STATVAR.energy = stratigraphy1.STATVAR.energy - excess_water_energy - excess_ice_energy;
                
                available_pore_space = stratigraphy2.STATVAR.layerThick(1).*stratigraphy2.STATVAR.area(1) - stratigraphy2.STATVAR.waterIce(1) - stratigraphy2.STATVAR.mineral(1) - stratigraphy2.STATVAR.organic(1);
                if excess_water > 0 % only water can enter top cell for ground
                    stratigraphy2.STATVAR.waterIce(1) = stratigraphy2.STATVAR.waterIce(1) + min(excess_water, available_pore_space);
                    stratigraphy2.STATVAR.energy(1) = stratigraphy2.STATVAR.energy(1) + excess_water_energy.*min(excess_water, available_pore_space)./excess_water;
                end
                stratigraphy2.STATVAR.excessWater = stratigraphy2.STATVAR.excessWater + excess_ice + max(0,excess_water-available_pore_space); % excess ice is routed directly here

                stratigraphy1 = get_T_water_vegetation(stratigraphy1);
                stratigraphy2 = compute_diagnostic(stratigraphy2, tile);
            end
             
        end
        
        
        %---service functions-----------------
        %redce in and outflow close to field capacity and full saturation,
        %bucket water scheme
        %Note: these functions are duplicates of the respective functions in TIER1 WATER_FLUXES of the GROUND classes 
        function rf = reduction_factor_out(saturation, ia_heat_water)  %part of get_derivative_water2(ground)
            smoothness = 3e-2;
            rf = (1-exp(-saturation./smoothness));
        end
        
        function rf = reduction_factor_in(saturation, ia_heat_water)   %part of get_derivative_water2(ground)
            smoothness = 3e-2;
            rf = (1- exp((saturation-1)./smoothness));
        end
    end
end


           
        
      