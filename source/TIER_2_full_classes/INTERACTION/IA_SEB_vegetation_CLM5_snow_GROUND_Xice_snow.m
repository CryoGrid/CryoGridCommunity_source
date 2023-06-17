%========================================================================
% CryoGrid INTERACTION (IA) class for surface energy balance of a GROUND
% class with snow CHILD below a shading VEGETATION class
% R. B. Zwegiel, February 2022
%========================================================================

classdef IA_SEB_vegetation_CLM5_snow_GROUND_Xice_snow < IA_SEB_vegetation_CLM5
   
    methods
        
        function ia_seb_water = get_boundary_condition_m(ia_seb_water, tile)
            % Roughly equivalent to
            % get_boundary_condition_u@GROUND_freezeC_RichardsEqW_seb_snow,
            % but with no radiation
            ground = ia_seb_water.NEXT; % for easy access
            vegetation = ground.PREVIOUS;
            
            if ground.CHILD == 0 % No CHILD exists
                ia_seb_water = get_boundary_condition_m_Xice(ia_seb_water, tile);
                
                if vegetation.TEMP.snow_thru > 0 % Create CHILD
                    %modify FORCING so that oroginal function can work
                    f_thru = vegetation.TEMP.snow_thru .* 1000 .* (24.*3600) ./tile.FORCING.TEMP.snowfall;
                    tile.FORCING.TEMP.snowfall = tile.FORCING.TEMP.snowfall .* f_thru;
                    
                    ground.CHILD = copy(tile.STORE.SNOW);
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground; 
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                    %would make sense to have a specifc function for making
                    %a CHILD which works both for CHILD in canopy and on
                    %ground below
                    ground.CHILD = get_boundary_condition_u_create_CHILD(ground.CHILD, tile);  %initialize with fresh snowfall
                    tile.FORCING.TEMP.snowfall = tile.FORCING.TEMP.snowfall ./ f_thru;
                end
            else % CHILD exists
                total.area = ground.STATVAR.area;
                total.waterIce  = ground.STATVAR.waterIce;
                total.water = ground.STATVAR.water;
                total.ice = ground.STATVAR.ice;
                total.mineral = ground.STATVAR.mineral;
                total.organic = ground.STATVAR.organic;
                total_XwaterIce = ground.STATVAR.XwaterIce;
                total_Xice = ground.STATVAR.Xice;
                total_Xwater = ground.STATVAR.Xwater;
                
                ground.STATVAR.area = ground.STATVAR.area - ground.CHILD.STATVAR.area(1); % use snow free area from now on
                reduction = ground.STATVAR.area(1) ./ total.area(1);
                ground.STATVAR.waterIce = ground.STATVAR.waterIce .* reduction;
                ground.STATVAR.water = ground.STATVAR.water .* reduction;
                ground.STATVAR.ice = ground.STATVAR.ice .* reduction;
                ground.STATVAR.mineral = ground.STATVAR.mineral .* reduction;
                ground.STATVAR.organic = ground.STATVAR.organic .* reduction;
                ground.STATVAR.XwaterIce = ground.STATVAR.XwaterIce .* reduction;
                ground.STATVAR.Xice = ground.STATVAR.Xice .* reduction;
                ground.STATVAR.Xwater = ground.STATVAR.Xwater .* reduction;
                
                %-------------
                
                ground.CHILD.STATVAR.Lstar = ground.STATVAR.Lstar;
                
                ia_seb_water = get_boundary_condition_m_vegetation_snow_CHILD_below(ia_seb_water, tile); %@IA_SEB; throughfall snow and rain on snow CHILD, Qh, Qe & d_energy
                
                %This is OK!
                ia_seb_water = get_boundary_condition_m_Xice(ia_seb_water, tile); % rain throughfall, Qh, Qe & d_energy
                
                get_IA_CHILD_boundary_condition_u(ground.IA_CHILD, tile); %call designated mandatory function for CHILD-PARENT interactions in the IA class governing IA between SNOW and GROUND
                
                % Mix the SEB fluxes from snow and ground
                ground.STATVAR.Qh = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qh + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qh) ./ total.area(1,1);
                ground.STATVAR.Qe = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qe + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qe) ./ total.area(1,1);
                
                %reassign the true totals of ground
                ground.STATVAR.area = total.area;
                ground.STATVAR.waterIce = total.waterIce;
                ground.STATVAR.water = total.water;
                ground.STATVAR.ice = total.ice;
                ground.STATVAR.mineral = total.mineral;
                ground.STATVAR.organic = total.organic;
                ground.STATVAR.XwaterIce = total_XwaterIce;
                ground.STATVAR.Xice = total_Xice;
                ground.STATVAR.Xwater = total_Xwater;
            end
        
        end
        
        function q_g = get_humidity_surface(ia_seb_water, tile)
            q_g = get_humidity_surface_GROUND(ia_seb_water, tile);
            if ia_seb_water.NEXT.CHILD ~= 0
                q_snow = get_humidity_surface_SNOW_CHILD(ia_seb_water, tile);
                snow_fraction = ia_seb_water.NEXT.CHILD.STATVAR.area./ia_seb_water.NEXT.STATVAR.area(1);
                q_g = q_g.*(1-snow_fraction) + q_snow.*snow_fraction;
            end
        end        
        
        %NEW SW, Dec 2022, must be separate function for Xice classes,
        %should also be changed for Snow_crocus2?
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
                
                stratigraphy2.STATVAR.XwaterIce(1,1) = stratigraphy2.STATVAR.XwaterIce(1,1) + excess_water + excess_ice;
                stratigraphy2.STATVAR.layerThick(1,1) = stratigraphy2.STATVAR.layerThick(1,1) + (excess_water + excess_ice)./stratigraphy2.STATVAR.area(1,1);
                stratigraphy2.STATVAR.energy(1,1) = stratigraphy2.STATVAR.energy(1,1) + excess_water_energy + excess_ice_energy;
                
                stratigraphy1 = get_T_water_vegetation(stratigraphy1);
                stratigraphy2 = compute_diagnostic(stratigraphy2, tile);
            end
        end
        
         function ia_seb_water = canopy_snow_drip(ia_seb_water, tile)
            vegetation = ia_seb_water.PREVIOUS;
            ground = ia_seb_water.NEXT; %ground with or without snow 
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
                %merge dripping snow
                if ground.CHILD ~= 0
                    ground.CHILD = merge_cells_intensive2(ground.CHILD, 1, drip_snow, 1, {'d'; 's'; 'gs'; 'time_snowfall'}, 'ice');
                    ground.CHILD = merge_cells_extensive2(ground.CHILD, 1, drip_snow, 1, {'waterIce'; 'energy'; 'volume'; 'ice'});
                    ground.CHILD = merge_cells_snowfall_times2(ground.CHILD, 1, drip_snow, 1); %specific function merginging bottom and top snow dates
                else %create new snow CHILD for ground, using the snow object from the vegetation
                    CHILD = copy(tile.STORE.SNOW);
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground; 
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                    
                    ground.CHILD = get_boundary_condition_u_create_CHILD(ground.CHILD, tile);  %fill with nonsense
                    variable_list = fieldnames(drip_snow.STATVAR);
                    for i=1:size(variable_list,1)
                        ground.CHILD.STATVAR.(variable_list{i,1}) = drip_snow.STATVAR.(variable_list{i,1});
                    end
                end
                ground = compute_diagnostic(ground, tile); %recomputes also for CHILD

            end
         end
        
    end 
    
end