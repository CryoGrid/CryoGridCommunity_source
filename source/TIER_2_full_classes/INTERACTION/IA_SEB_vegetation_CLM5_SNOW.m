%========================================================================
% CryoGrid INTERACTION (IA) class for surface energy balance of a SNOW
% class below a shading VEGETATION class
% R. B. Zwegiel, February 2022
%========================================================================

classdef IA_SEB_vegetation_CLM5_SNOW < IA_SEB & IA_WATER % Throughfall of water
    
    methods
        function  ia_seb_water = get_boundary_condition_m(ia_seb_water, tile) 
            % SEB of snow below canopy
            forcing = tile.FORCING;
            snow = ia_seb_water.NEXT;
            
            % 1. get snowfall
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600) .* snow.STATVAR.area(1,1); %snowfall is in mm/day -> [m3/sec]
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);  %[J/sec]
            % 2. Add rainfall
            ia_seb_water = get_boundary_condition_water_canopy_SNOW_m(ia_seb_water, tile);
            % 3. turbulent fluxes
            ia_seb_water = get_boundary_condition_Qh_CLM5_m(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_Qe_CLM5_SNOW_m(ia_seb_water, tile);
            % 4. add fluxes to uppermost cell (radiative fluxes are added by penetration)
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + (-snow.STATVAR.Qh - snow.STATVAR.Qe) .* snow.STATVAR.area(1);
            % 5. make new snow
            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing);
            % 6. store vars for use later
            snow.TEMP.wind = forcing.TEMP.wind;
            snow.TEMP.wind_surface = forcing.TEMP.wind;
        end
        
        function remove_excessWater_CHILD(ia_heat_water) %move excessWater from SNOW CHILD to water in canopy, from there it drips to the ground or snow below
            snow = ia_heat_water.PREVIOUS;
            canopy = ia_heat_water.NEXT;
            
            porosity = snow.STATVAR.area .* snow.STATVAR.layerThick - snow.STATVAR.ice;
            max_water = porosity .* snow.PARA.field_capacity;
            remove_water = max(0, snow.STATVAR.water - max_water);
            snow.STATVAR.waterIce = snow.STATVAR.waterIce - remove_water;
            snow.STATVAR.water = snow.STATVAR.water - remove_water;
            snow.STATVAR.excessWater = snow.STATVAR.excessWater + remove_water;
            canopy.STATVAR.waterIce = canopy.STATVAR.waterIce + snow.STATVAR.excessWater;
            snow.STATVAR.excessWater = 0;
            snow.STATVAR.excessWater_energy = 0;
            %no energy transfer, water in snow is at 0 degree
        end
        
        function ia_seb_water = canopy_snow_drip(ia_seb_water, tile) %move dripping snow to fully developed snow class below 
            vegetation = ia_seb_water.PREVIOUS;
            snow_below = ia_seb_water.NEXT; %snow below canopy
            snow = vegetation.CHILD;
            snow_fraction = snow.STATVAR.area ./ vegetation.STATVAR.area; %snow faction in canopy
            
            if snow_fraction > vegetation.PARA.max_snow_fraction %this could be much more complicated, make dependent on forcing,snow state, etc.
                remaining_area = vegetation.PARA.max_snow_fraction .* vegetation.STATVAR.area;
                removed_area = snow.STATVAR.area - remaining_area;
                remaining_fraction = remaining_area ./ snow.STATVAR.area;
                removed_fraction = 1 - remaining_fraction;
                %this is  coded for Crocus, make it a specific
                %split function later
                snow.STATVAR.volume = snow.STATVAR.layerThick .* snow.STATVAR.area;
                variable_list = {'energy'; 'waterIce'; 'water'; 'ice'; 'volume'; 'area'};
                for i=1:size(variable_list, 1)
                    drip_snow.STATVAR.(variable_list{i,1}) = removed_fraction .* snow.STATVAR.(variable_list{i,1});
                    snow.STATVAR.(variable_list{i,1}) = snow.STATVAR.(variable_list{i,1}) - drip_snow.STATVAR.(variable_list{i,1});
                end
                variable_list = {'d'; 's'; 'gs'; 'time_snowfall'; 'layerThick'; 'top_snow_date'; 'bottom_snow_date'; 'target_density'};
                for i=1:size(variable_list, 1)
                    drip_snow.STATVAR.(variable_list{i,1}) = snow.STATVAR.(variable_list{i,1});
                end
                drip_snow.STATVAR.layerThick = drip_snow.STATVAR.volume ./ snow.STATVAR.area(1);
                drip_snow.STATVAR.layerThick = snow.STATVAR.area(1);
                %merge dripping snow
                snow_below = merge_cells_intensive2(snow_below, 1, drip_snow, 1, {'d'; 's'; 'gs'; 'time_snowfall'}, 'ice');
                snow_below = merge_cells_extensive2(snow_below, 1, drip_snow, 1, {'waterIce'; 'energy'; 'layerThick'; 'ice'});
                snow_below = merge_cells_snowfall_times2(snow_below, 1, drip_snow, 1); %specific function merginging bottom and top snow dates
                
                snow_below = compute_diagnostic(snow_below, tile); %recomputes also for CHILD

            end
         end

        
        function q_g = get_humidity_surface(ia_seb_water, tile)
            q_g = get_humidity_surface_SNOW(ia_seb_water, tile);
        end
        
         function r_soil = ground_resistance_evap(ia_seb_water, tile)
             r_soil = 0; % No extra vapor resistance for snow surface
         end
         
         function ia_seb_water = canopy_drip(ia_seb_water, tile)
            ia_seb_water = canopy_drip@IA_WATER(ia_seb_water, tile);
        end
    end
    
end
