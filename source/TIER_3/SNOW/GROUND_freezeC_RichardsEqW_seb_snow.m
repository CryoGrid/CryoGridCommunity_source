%========================================================================
% CryoGrid GROUND class GROUND_freezeC_RichardsEqW_seb_snow
% heat conduction, Richards equation water scheme, freeze curve based on
% freezing=drying assumption, surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freezeC_RichardsEqW_seb_snow < GROUND_freezeC_RichardsEqW_seb
    properties
        CHILD
        IA_CHILD
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

       function ground = provide_PARA(ground)  
            ground = provide_PARA@GROUND_freezeC_RichardsEqW_seb(ground);
       end
       
       function ground = provide_CONST(ground)  
           ground = provide_CONST@GROUND_freezeC_RichardsEqW_seb(ground);
       end
       
       function ground = provide_STATVAR(ground)  
           ground = provide_STATVAR@GROUND_freezeC_RichardsEqW_seb(ground);
       end
       
       function ground = finalize_init(ground, tile)
           ground = finalize_init@GROUND_freezeC_RichardsEqW_seb(ground, tile);
           ground.CHILD = 0; % no snow
           ground.IA_CHILD = 0;
       end

       
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            if ground.CHILD == 0  %CHILD does not exist
                ground = get_boundary_condition_u@GROUND_freezeC_RichardsEqW_seb(ground, tile); %call the native function for the ground class
                
                if forcing.TEMP.snowfall > 0  %create CHILD 
%                     CURRENT = ground.PREVIOUS;  %go to Top() and get the stored SNOW class
%                     while ~strcmp(class(CURRENT), 'Top')
%                         CURRENT = CURRENT.PREVIOUS;
%                     end
%                     ground.CHILD = copy(CURRENT.STORE.SNOW);
                    ground.CHILD = copy(tile.STORE.SNOW);
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground; 
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                    
                    ground.CHILD = get_boundary_condition_u_create_CHILD(ground.CHILD, tile);  %initialize with fresh snowfall
                end
            else %CHILD exists
                total_area = ground.STATVAR.area; %store the total area of the ground
                total_waterIce = ground.STATVAR.waterIce;
                total_water = ground.STATVAR.water;
                total_ice = ground.STATVAR.ice;
                total_mineral = ground.STATVAR.mineral;
                total_organic = ground.STATVAR.organic;
                
                %split up area in snow-covered (CHILD) and snow-free part (PARENT)
                ground.STATVAR.area = ground.STATVAR.area - ground.CHILD.STATVAR.area(1,1); %replace by snow-free area
                reduction = ground.STATVAR.area(1) ./ total_area(1);
                ground.STATVAR.waterIce = ground.STATVAR.waterIce .* reduction;
                ground.STATVAR.water = ground.STATVAR.water .* reduction;
                ground.STATVAR.ice = ground.STATVAR.ice .* reduction;
                ground.STATVAR.mineral = ground.STATVAR.mineral .* reduction;
                ground.STATVAR.organic = ground.STATVAR.organic .* reduction;
                
                %-------------
                
                ground.CHILD.STATVAR.Lstar = ground.STATVAR.Lstar;
                
                ground.CHILD = get_boundary_condition_u_CHILD(ground.CHILD, tile);
                ground = get_boundary_condition_u@GROUND_freezeC_RichardsEqW_seb(ground, tile);
                
                get_IA_CHILD_boundary_condition_u(ground.IA_CHILD, tile);
                %call designated mandatory function for CHILD-PARENT interactions in
                %the IA class governing IA between SNOW and GROUND
                
                ground.STATVAR.Lout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Lout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Lout) ./ total_area(1,1); %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Sout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Sout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Sout) ./ total_area(1,1);
                ground.STATVAR.Qh = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qh + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qh) ./ total_area(1,1);
                ground.STATVAR.Qe = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qe + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qe) ./ total_area(1,1);
                
                %reassign the true totals of ground
                ground.STATVAR.area = total_area;
                ground.STATVAR.waterIce = total_waterIce;
                ground.STATVAR.water = total_water;
                ground.STATVAR.ice = total_ice;
                ground.STATVAR.mineral = total_mineral;
                ground.STATVAR.organic = total_organic;
                
            end
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            if ground.CHILD ~= 0
                snow_fraction = ground.CHILD.STATVAR.area(1)./ground.STATVAR.area(1);
                ground_fraction = 1 - snow_fraction;
                [ground, S_up] = penetrate_SW@GROUND_freezeC_RichardsEqW_seb(ground, S_down.*ground_fraction);
                [ground.CHILD, S_up2] = penetrate_SW(ground.CHILD, S_down.*snow_fraction);
                S_up = S_up + sum(S_up2); % snow_crocus splits SW into spectral bands
            else
                [ground, S_up] = penetrate_SW@GROUND_freezeC_RichardsEqW_seb(ground, S_down);
            end
        end
        
        function [ground, S_up] = penetrate_SW_PARENT(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW_PARENT@GROUND_freezeC_RichardsEqW_seb(ground, S_down);
        end
        
        function [ground, L_up] = penetrate_LW(ground, L_down)  %mandatory function when used with class that features SW penetration
            % Lin is in W, not W/m2!
            
             if ground.CHILD ~= 0
                snow_fraction = ground.CHILD.STATVAR.area(1)./ground.STATVAR.area(1);
                ground_fraction = 1 - snow_fraction;
                L_up_ground = ground_fraction .*(1-ground.PARA.epsilon) .* L_down  + ground_fraction .* ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ ground.CONST.Tmfw).^4 .*ground.STATVAR.area(1); %area is full area here
                L_up_snow = snow_fraction .*(1-ground.CHILD.PARA.epsilon) .* L_down  + ground.CHILD.PARA.epsilon .* ground.CHILD.CONST.sigma .* (ground.CHILD.STATVAR.T+ ground.CHILD.CONST.Tmfw).^4 .*ground.CHILD.STATVAR.area(1);
                L_up = L_up_ground + L_up_snow;
                ground.STATVAR.Lout = L_up./ground.STATVAR.area(1);
                ground.STATVAR.Lin = L_down./ground.STATVAR.area(1);
                %ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) + L_down - L_up;
                ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) + L_down.* ground_fraction - L_up_ground;
                ground.CHILD.TEMP.d_energy(1,1) = ground.CHILD.TEMP.d_energy(1,1) + L_down.* snow_fraction - L_up_snow;
            else
                [ground, L_up] = penetrate_LW@GROUND_freezeC_RichardsEqW_seb(ground, L_down);
             end
        end
        
        function ground = get_boundary_condition_l(ground, tile)
              ground = get_boundary_condition_l@GROUND_freezeC_RichardsEqW_seb(ground, tile);
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            if ground.CHILD == 0  
                ground = get_derivatives_prognostic@GROUND_freezeC_RichardsEqW_seb(ground, tile); %call normal function
            else
                ground.CHILD = get_derivatives_prognostic_CHILD(ground.CHILD, tile);
                ground = get_derivatives_prognostic@GROUND_freezeC_RichardsEqW_seb(ground, tile); 
            end
        end
        
        function timestep = get_timestep(ground, tile) 
            if ground.CHILD == 0
                timestep = get_timestep@GROUND_freezeC_RichardsEqW_seb(ground, tile);
            else 
                timestep_snow = get_timestep_CHILD(ground.CHILD, tile);
                timestep_ground =  get_timestep@GROUND_freezeC_RichardsEqW_seb(ground, tile);
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
            end
        end
        
        function ground = advance_prognostic(ground, tile) 
            if ground.CHILD == 0
                ground =  advance_prognostic@GROUND_freezeC_RichardsEqW_seb(ground, tile);
            else                
                ground.CHILD = advance_prognostic_CHILD(ground.CHILD, tile);
                ground =  advance_prognostic@GROUND_freezeC_RichardsEqW_seb(ground, tile);
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, tile)
            if ground.CHILD == 0
                ground = compute_diagnostic@GROUND_freezeC_RichardsEqW_seb(ground, tile);
            else
                ground = compute_diagnostic@GROUND_freezeC_RichardsEqW_seb(ground, tile);
                ground.CHILD = compute_diagnostic_CHILD(ground.CHILD, tile);
                
            end
        end
        
        function ground = check_trigger(ground, tile)
            if ground.CHILD ~= 0
                %delete CHILD
                if ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) < 1e-6 %cutoff to get rid of remaining snow
                   ground.CHILD = 0;
                   ground.IA_CHILD = 0;
                %make SNOW CHILD full class   
                elseif ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) > 1 
                    %transforms dimensions and STAVAR
                    snow_volume = ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.layerThick;
                    ground.CHILD.STATVAR.area = ground.STATVAR.area(1,1);
                    ground.CHILD.STATVAR.layerThick = snow_volume ./ ground.CHILD.STATVAR.area;
                   
%                     %make snow a real class
%                     ground.CHILD.PARENT = 0;
%                     ground.CHILD.PREVIOUS = ground.PREVIOUS;
%                     ground.CHILD.NEXT = ground;
%                     ground.PREVIOUS.NEXT = ground.CHILD;
%                     ground.PREVIOUS = ground.CHILD;
%                     ground.CHILD = 0;
%                     ground.IA_PREVIOUS = ground.IA_CHILD; 
%                     ground.PREVIOUS.IA_NEXT = ground.IA_CHILD;
%                     ground.IA_CHILD = 0;
                    %make snow a real class
                    ground.CHILD.PARENT = 0;
                    ground.CHILD.PREVIOUS = ground.PREVIOUS;
                    ground.CHILD.NEXT = ground;
                    ground.PREVIOUS.NEXT = ground.CHILD;
                    ia_class = get_IA_class(class(ground.PREVIOUS), class(ground.CHILD));
                    ground.PREVIOUS.IA_NEXT = ia_class;
                    ground.CHILD.IA_PREVIOUS = ia_class;
                    ground.CHILD.IA_PREVIOUS.NEXT = ground.CHILD;
                    ground.CHILD.IA_PREVIOUS.PREVIOUS = ground.PREVIOUS;
                    finalize_init(ground.CHILD.IA_PREVIOUS, tile);
                    
                    ground.PREVIOUS = ground.CHILD;
                    ground.CHILD = 0;
                    ground.IA_PREVIOUS = ground.IA_CHILD; 
                    ground.PREVIOUS.IA_NEXT = ground.IA_CHILD;
                    ground.IA_CHILD = 0;
                    
                end
            end
        end
        
        
        function z0 = get_z0_surface(ground)
            if ground.CHILD ~= 0
                snow_fraction = ground.CHILD.STATVAR.area/ground.STATVAR.area(1);
                ground_fraction = 1 - snow_fraction;
                z0 = snow_fraction .* get_z0_surface(ground.CHILD) +  ground_fraction .* ground.PARA.z0;
            else
                z0 = ground.PARA.z0;
            end
        end
        
        function albedo = get_albedo(ground)
            if ground.CHILD ~= 0
                snow_fraction = ground.CHILD.STATVAR.area/ground.STATVAR.area(1);
                ground_fraction = 1 - snow_fraction;
                albedo = snow_fraction .* get_albedo(ground.CHILD) +  ground_fraction .* ground.PARA.albedo;
            else
                albedo = ground.PARA.albedo;
            end
        end
        
        function Ts = get_surface_T(ground, tile)
            if ground.CHILD ~= 0
                snow_fraction = ground.CHILD.STATVAR.area/ground.STATVAR.area(1);
                ground_fraction = 1 - snow_fraction;
                Ts = snow_fraction .* get_surface_T(ground.CHILD) +  ground_fraction .* ground.STATVAR.T(1);
            else
                Ts = ground.STATVAR.T(1);
            end
        end

                
        %----------
        %reset timestamp when changing TILES
        function ground = reset_timestamps(ground, tile)
            if ground.CHILD ~= 0
                ground.CHILD = reset_timestamps(ground.CHILD, tile);
            end
        end
        
    end
end