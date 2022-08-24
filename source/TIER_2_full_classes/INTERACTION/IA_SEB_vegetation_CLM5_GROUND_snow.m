%========================================================================
% CryoGrid INTERACTION (IA) class for surface energy balance of a GROUND
% class with snow CHILD below a shading VEGETATION class
% R. B. Zwegiel, February 2022
%========================================================================

classdef IA_SEB_vegetation_CLM5_GROUND_snow < IA_SEB_vegetation_CLM5
   
    methods
        
        function ia_seb_water = get_boundary_condition_m(ia_seb_water, tile)
            % Roughly equivalent to
            % get_boundary_condition_u@GROUND_freezeC_RichardsEqW_seb_snow,
            % but with no radiation
            ground = ia_seb_water.NEXT; % for easy acces
            
            if ground.CHILD == 0 % No CHILD exists
                ia_seb_water = get_boundary_condition_m@IA_SEB_vegetation_CLM5(ia_seb_water, tile);
                
                if tile.FORCING.TEMP.snowfall > 0 % Create CHILD
                    ground.CHILD = copy(tile.STORE.SNOW);
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground; 
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                    
                    ground.CHILD = get_boundary_condition_u_create_CHILD(ground.CHILD, tile);  %initialize with fresh snowfall
                end
            else % CHILD exists
                total.area = ground.STATVAR.area;
                total.waterIce  = ground.STATVAR.waterIce;
                total.water = ground.STATVAR.water;
                total.ice = ground.STATVAR.ice;
                total.mineral = ground.STATVAR.mineral;
                total.organic = ground.STATVAR.organic;
                
                ground.STATVAR.area = ground.STATVAR.area - ground.CHILD.STATVAR.area(1); % use snow free area from now on
                reduction = ground.STATVAR.area(1) ./ total.area(1);
                ground.STATVAR.waterIce = ground.STATVAR.waterIce .* reduction;
                ground.STATVAR.water = ground.STATVAR.water .* reduction;
                ground.STATVAR.ice = ground.STATVAR.ice .* reduction;
                ground.STATVAR.mineral = ground.STATVAR.mineral .* reduction;
                ground.STATVAR.organic = ground.STATVAR.organic .* reduction;
                
                %-------------
                
                ground.CHILD.STATVAR.Lstar = ground.STATVAR.Lstar;
                
                ia_seb_water = get_boundary_condition_m_CHILD(ia_seb_water, tile); %@IA_SEB; snowfall, throughfall, Qh, Qe & d_energy
                
                ia_seb_water = get_boundary_condition_m@IA_SEB_vegetation_CLM5(ia_seb_water, tile); % throughfall, Qh, Qe & d_energy
                
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
            end
        
        end
        
        function q_g = get_humidity_surface(ia_seb_water, tile)
            q_g = get_humidity_surface@IA_SEB_vegetation_CLM5(ia_seb_water, tile);
            if ia_seb_water.NEXT.CHILD ~= 0
                q_snow = get_humidity_surface_SNOW_CHILD(ia_seb_water, tile);
                snow_fraction = ia_seb_water.NEXT.CHILD.STATVAR.area./ia_seb_water.NEXT.STATVAR.area(1);
                q_g = q_g*(1-snow_fraction) + q_snow*snow_fraction;
            end
        end        
    
    end 
    
end