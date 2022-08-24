%========================================================================
% CryoGrid INTERACTION (IA) class for surface energy balance of a SNOW
% class below a shading VEGETATION class
% R. B. Zwegiel, February 2022
%========================================================================

classdef IA_SEB_vegetation_CLM5_SNOW < IA_SEB & IA_WATER % Throughfall of water
    
    methods
        function  ia_seb_water = get_boundary_condition_m(ia_seb_water, tile) 
            forcing = tile.FORCING;
            % SEB of snow below canopy
            snow = ia_seb_water.NEXT;
            % 1. get snowfall
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600) .* snow.STATVAR.area(1,1); %snowfall is in mm/day -> [m3/sec]
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);  %[J/sec]
            % 2. Add rainfall
            ia_seb_water = get_boundary_condition_water_canopy_SNOW_m(ia_seb_water, tile);
            % 3. turbulent fluxes
            ia_seb_water = get_boundary_condition_Qh_CLM5_m(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_Qe_CLM5_SNOW_m(ia_seb_water, tile);
            % 4. add fluxes to uppermost cell (ratiative fluxes are added by penetration)
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + (-snow.STATVAR.Qh - snow.STATVAR.Qe) .* snow.STATVAR.area(1);
            % 5. make new snow
            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing);
            % 6. store vars for use later
            snow.TEMP.wind = forcing.TEMP.wind;
            snow.TEMP.wind_surface = forcing.TEMP.wind;
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
