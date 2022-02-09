%========================================================================
% CryoGrid GROUND class GROUND_freezeC_bucketW_Xice_seb_snow
% heat conduction, bucket water scheme, freeze curve based on
% freezing=drying assumption, surface energy balance, excess ice
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freezeC_bucketW_Xice_seb_snow < GROUND_freezeC_bucketW_Xice_seb
    properties
        CHILD
        IA_CHILD
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = provide_PARA(ground)
            ground = provide_PARA@GROUND_freezeC_bucketW_Xice_seb(ground);
        end
        
        function ground = provide_CONST(ground)
            ground = provide_CONST@GROUND_freezeC_bucketW_Xice_seb(ground);
        end
        
        function ground = provide_STATVAR(ground)
            ground = provide_STATVAR@GROUND_freezeC_bucketW_Xice_seb(ground);
        end
        
        function ground = finalize_init(ground, tile)
            ground = finalize_init@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
            ground.CHILD = 0; % no snow
            ground.IA_CHILD = 0;
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            if ground.CHILD == 0  %CHILD does not exist
                ground = get_boundary_condition_u@GROUND_freezeC_bucketW_Xice_seb(ground, tile); %call the native function for the ground class
                
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
                total_area = ground.STATVAR.area; %store the totals of the GROUND class
                total_waterIce = ground.STATVAR.waterIce;
                total_water = ground.STATVAR.water;
                total_ice = ground.STATVAR.ice;
                total_mineral = ground.STATVAR.mineral;
                total_organic = ground.STATVAR.organic;
                total_XwaterIce = ground.STATVAR.XwaterIce;
                total_Xice = ground.STATVAR.Xice;
                total_Xwater = ground.STATVAR.Xwater;
                
                %split up area in snow-covered (CHILD) and snow-free part (PARENT)
                ground.STATVAR.area = ground.STATVAR.area - ground.CHILD.STATVAR.area(1,1); 
                reduction = ground.STATVAR.area(1) ./ total_area(1);
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
                
                ground.CHILD = get_boundary_condition_u_CHILD(ground.CHILD, tile);
                ground = get_boundary_condition_u@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
                
                get_IA_CHILD_boundary_condition_u(ground.IA_CHILD, tile);
                %call designated mandatory function for CHILD-PARENT interactions in
                %the IA class governing IA between SNOW and GROUND
                
                ground.STATVAR.Lout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Lout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Lout) ./ total_area(1,1); %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Sout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Sout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Sout) ./ total_area(1,1);
                ground.STATVAR.Qh = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qh + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qh) ./ total_area(1,1);
                ground.STATVAR.Qe = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qe + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qe) ./ total_area(1,1);
                
                %reassign the true area of ground
                ground.STATVAR.area = total_area; 
                ground.STATVAR.waterIce = total_waterIce;
                ground.STATVAR.water = total_water;
                ground.STATVAR.ice = total_ice;
                ground.STATVAR.mineral = total_mineral;
                ground.STATVAR.organic = total_organic;
                ground.STATVAR.XwaterIce = total_XwaterIce;
                ground.STATVAR.Xice = total_Xice;
                ground.STATVAR.Xwater = total_Xwater;
            end
            
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW@GROUND_freezeC_bucketW_Xice_seb(ground, S_down);
        end
        
        
        function ground = get_boundary_condition_l(ground, tile)
            ground = get_boundary_condition_l@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
        end
        
        
        function ground = get_derivatives_prognostic(ground, tile)
            if ground.CHILD == 0
                ground = get_derivatives_prognostic@GROUND_freezeC_bucketW_Xice_seb(ground, tile); %call normal function
            else
                ground.CHILD = get_derivatives_prognostic_CHILD(ground.CHILD, tile);
                ground = get_derivatives_prognostic@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
            end
        end
        
        function timestep = get_timestep(ground, tile)
            if ground.CHILD == 0
                timestep = get_timestep@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
            else
                timestep_snow = get_timestep_CHILD(ground.CHILD, tile);
                timestep_ground =  get_timestep@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
            end
        end
        
        function ground = advance_prognostic(ground, tile)
            if ground.CHILD == 0
                ground =  advance_prognostic@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
            else
                ground.CHILD = advance_prognostic_CHILD(ground.CHILD, tile);
                ground =  advance_prognostic@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, tile)
            if ground.CHILD == 0
                ground = compute_diagnostic@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
            else
                ground.CHILD = compute_diagnostic_CHILD(ground.CHILD, tile);
                ground = compute_diagnostic@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
            end
        end
        
        function ground = check_trigger(ground, tile)
            
            ground = check_trigger@GROUND_freezeC_bucketW_Xice_seb(ground, tile);
            
            if ground.CHILD ~= 0
                %delete CHILD
                if ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) < 1e-6 %cutoff to get rid of remaining snow
                    ground.CHILD = 0;
                    ground.IA_CHILD = 0;
                    %make SNOW CHILD full class
                elseif ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) > 1
                    %transform dimensions and STAVAR
                    snow_volume = ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.layerThick;
                    ground.CHILD.STATVAR.area = ground.STATVAR.area(1,1);
                    ground.CHILD.STATVAR.layerThick = snow_volume ./ ground.CHILD.STATVAR.area;
                    
                    %make snow a real class
                    ground.CHILD.PARENT = 0;
                    ground.CHILD.PREVIOUS = ground.PREVIOUS;
                    ground.CHILD.NEXT = ground;
                    ground.PREVIOUS.NEXT = ground.CHILD;
                    ground.PREVIOUS = ground.CHILD;
                    ground.CHILD = 0;
                    ground.IA_PREVIOUS = ground.IA_CHILD; 
                    ground.PREVIOUS.IA_NEXT = ground.IA_CHILD;
                    ground.IA_CHILD = 0;
                end
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