%========================================================================
% CryoGrid GROUND class GROUND_freezeC_ubT_snow forced by surface
% temperature
% heat conduction, constant water  +ice, freeze curve based on
% freezing=drying assumption, surface energy balance
% S. Westermann, October 2021
%========================================================================

classdef GROUND_freeW_ubT_snow < GROUND_freeW_ubT
    properties
        CHILD
        IA_CHILD
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
       function ground = provide_PARA(ground)  
            ground = provide_PARA@GROUND_freeW_ubT(ground);
       end
       
       function ground = provide_CONST(ground)  
           ground = provide_CONST@GROUND_freeW_ubT(ground);
       end
       
       function ground = provide_STATVAR(ground)  
           ground = provide_STATVAR@GROUND_freeW_ubT(ground);
       end


       function ground = finalize_init(ground, tile)
           ground = finalize_init@GROUND_freeW_ubT(ground, tile);
           ground.CHILD = 0; % no snow
           ground.IA_CHILD = 0;
       end

       
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            forcing = tile.FORCING;
            if ground.CHILD == 0  %CHILD does not exist
                ground = get_boundary_condition_u@GROUND_freeW_ubT(ground, tile); %call the native function for the ground class
                
                if forcing.TEMP.snow_depth > 0  %create CHILD 
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
                total_area = ground.STATVAR.area(1,1); %store the total area of the ground
                ground.STATVAR.area(1,1) = ground.STATVAR.area(1,1) - ground.CHILD.STATVAR.area(1,1); %replace by snow-free area
                
                ground.CHILD.STATVAR.Lstar = -100;
                
                ground.CHILD = get_boundary_condition_u_CHILD(ground.CHILD, tile);
                ground = get_boundary_condition_u@GROUND_freeW_ubT(ground, tile);
                
                get_IA_CHILD_boundary_condition_u(ground.IA_CHILD, tile);
                %call designated mandatory function for CHILD-PARENT interactions in
                %the IA class governing IA between SNOW and GROUND
                
                ground.STATVAR.area(1,1) = total_area; %reassign the true area of ground
            end
        end
        
%         function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
%             [ground, S_up] = penetrate_SW@GROUND_freeW_ubT(ground, S_down);
%         end
        
        function ground = get_boundary_condition_l(ground, tile)
              ground = get_boundary_condition_l@GROUND_freeW_ubT(ground, tile);

        end
        
        
        function ground = get_derivatives_prognostic(ground, tile)
            if ground.CHILD == 0  
                ground = get_derivatives_prognostic@GROUND_freeW_ubT(ground, tile); %call normal function
            else
                ground.CHILD = get_derivatives_prognostic_CHILD(ground.CHILD, tile);
                ground = get_derivatives_prognostic@GROUND_freeW_ubT(ground, tile); 
            end
        end
        
        function timestep = get_timestep(ground, tile) 
            if ground.CHILD == 0
                timestep = get_timestep@GROUND_freeW_ubT(ground, tile);
            else 
                timestep_snow = get_timestep_CHILD(ground.CHILD, tile);
                timestep_ground =  get_timestep@GROUND_freeW_ubT(ground, tile);
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
            end
        end
        
        function ground = advance_prognostic(ground, tile) 
            if ground.CHILD == 0
                ground =  advance_prognostic@GROUND_freeW_ubT(ground, tile);
            else                
                ground.CHILD = advance_prognostic_CHILD(ground.CHILD, tile);
                ground =  advance_prognostic@GROUND_freeW_ubT(ground, tile);
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)

        end
        
        function ground = compute_diagnostic(ground, tile)
            if ground.CHILD == 0
                ground = compute_diagnostic@GROUND_freeW_ubT(ground, tile);
            else
                ground = compute_diagnostic@GROUND_freeW_ubT(ground, tile);
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
                    %ground.CHILD = compute_diagnostic(ground.CHILD, forcing); %splits snow in 2 grid cells
                   
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
        
                
        %----------
        %reset timestamp when changing TILES
        function ground = reset_timestamps(ground, tile)
            if ground.CHILD ~= 0
                ground.CHILD = reset_timestamps(ground.CHILD, tile);
            end
        end
        
    end
end