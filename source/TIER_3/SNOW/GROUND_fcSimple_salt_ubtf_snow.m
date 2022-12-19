%========================================================================
% CryoGrid GROUND class GROUND_fcSimple_salt_UBTF_snow
% heat conduction, static soil water, salt diffusion, simple freeze curve including freezing point depression due to salt
% coupling to SNOW classes
% T. Ingeman-Nielsen, S. Westermann, J. Scheer, December 2021
%========================================================================

classdef GROUND_fcSimple_salt_ubtf_snow < GROUND_fcSimple_salt_ubtf
    properties
        CHILD
        IA_CHILD
    end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = finalize_init(ground, forcing)
            % call parent class finalize_init method
            ground = finalize_init@GROUND_fcSimple_salt_ubtf(ground, forcing);
            
            % then initialize child classes (for snow)
            ground.CHILD = 0; % no snow in the beginning
            ground.IA_CHILD = 0;
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            % Override parent class method in order to handle snow
            
            forcing = tile.FORCING;
            
            if ground.CHILD == 0  %SNOW CHILD does not exist
                % call the native method for the ground class to get upper
                % boundary condition
                ground = get_boundary_condition_u@GROUND_fcSimple_salt_ubtf(ground, tile); %call the native function for the ground class
                
                if forcing.TEMP.snowfall > 0  
                    % If snowfall in this timestep, create child class
                    
                    % go to Top() and get the stored SNOW class
                    CURRENT = ground.PREVIOUS;  
                    while ~strcmp(class(CURRENT), 'Top')
                        CURRENT = CURRENT.PREVIOUS;
                    end
                    
                    % set pointers between classes
                    ground.CHILD = copy(CURRENT.STORE.SNOW);
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground;
                    
                    % define interaction between snow class and ground class
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                    
                    % Now get upper boundary condition for child class
                    ground.CHILD = get_boundary_condition_u_create_CHILD(ground.CHILD, tile);  %initialize with fresh snowfall
                end
            else %SNOW CHILD exists
                % The child snow class already exist
                total_area = ground.STATVAR.area(1,1); %store the total area of the ground
                
                % temporarily replace area by snow-free area
                ground.STATVAR.area(1,1) = ground.STATVAR.area(1,1) - ground.CHILD.STATVAR.area(1,1); %replace by snow-free area
                
                % get upper boundary conditions for child and native ground class
                ground.CHILD = get_boundary_condition_u_CHILD(ground.CHILD, tile);
                ground = get_boundary_condition_u@GROUND_fcSimple_salt_ubtf(ground, tile);
                
                get_IA_CHILD_boundary_condition_u(ground.IA_CHILD, tile);
                %call designated mandatory function for CHILD-PARENT interactions in
                %the IA class governing IA between SNOW and GROUND
                
                %reassign the true area of ground
                ground.STATVAR.area(1,1) = total_area; %reassign the true area of ground
            end
        end
        
       
        function ground = get_derivatives_prognostic(ground, tile)
            if ground.CHILD == 0
                % when there is no child defined, call only method on
                % native ground class
                ground = get_derivatives_prognostic@GROUND_fcSimple_salt_seb(ground, tile); %call normal function
            else
                % when a child is defined, call method of child, then
                % method of native ground class                
                ground.CHILD = get_derivatives_prognostic_CHILD(ground.CHILD, tile);
                ground = get_derivatives_prognostic@GROUND_fcSimple_salt_seb(ground, tile);
            end
        end
        
        
        function timestep = get_timestep(ground, tile)
            if ground.CHILD == 0
                % when there is no child defined, get the timestep from the
                % native ground class
                timestep = get_timestep@GROUND_fcSimple_salt_seb(ground, tile);
            else
                % when a child is defined, get timestep of both child and
                % native ground class, and take the minimum
                % (except if timestep_snow is 0, it is neglected)
                timestep_snow = get_timestep_CHILD(ground.CHILD, tile);
                timestep_ground =  get_timestep@GROUND_fcSimple_salt_seb(ground, tile);
                
                % The following statement:
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
                % is equivalent to:
                % if timestep_snow > 0
                %     timestep = min(timestep_ground, timestep_snow);
                % else
                %     timestep = timestep_ground;
                % end
                
            end
        end
        
        
        function ground = advance_prognostic(ground, tile) 
            if ground.CHILD == 0
                % when there is no child defined, call only method on
                % native ground class
                ground =  advance_prognostic@GROUND_fcSimple_salt_seb(ground, tile);
            else
                % when a child is defined, call method of child, then
                % method of native ground class
                ground.CHILD = advance_prognostic_CHILD(ground.CHILD, tile);
                ground =  advance_prognostic@GROUND_fcSimple_salt_seb(ground, tile);
            end
        end
        
        
        function ground = compute_diagnostic(ground, tile)
            if ground.CHILD == 0
                % when there is no child defined, call only method on
                % native ground class
                ground = compute_diagnostic@GROUND_fcSimple_salt_seb(ground, tile);
            else
                % when a child is defined, call method of child, then
                % method of native ground class
                ground = compute_diagnostic@GROUND_fcSimple_salt_seb(ground, tile);
                ground.CHILD = compute_diagnostic_CHILD(ground.CHILD, tile);
            end
        end
        
        
        function ground = check_trigger(ground, tile)
            if ground.CHILD ~= 0 % CHILD class exists
                if ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) < 1e-6 % cutoff to get rid of remaining snow
                    % if the area of the snow class is less than threshold
                    % remove the snow class
                    ground.CHILD = 0;
                    ground.IA_CHILD = 0;
                
                elseif ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) > 1
                    % if area of SNOW CHILD exceeds the ground class area,
                    % convert snow to into a full class

                    %transforms dimensions and STAVAR
                    snow_volume = ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.layerThick;
                    ground.CHILD.STATVAR.area = ground.STATVAR.area(1,1);
                    ground.CHILD.STATVAR.layerThick = snow_volume ./ ground.CHILD.STATVAR.area;
                     
%                     % transfer the snow class from a CHILD of the ground 
%                     % to a part of the ordinary stratigraphy
%                     ground.CHILD.PARENT = 0;
%                     ground.CHILD.PREVIOUS = ground.PREVIOUS;
%                     ground.CHILD.NEXT = ground;
%                     ground.PREVIOUS.NEXT = ground.CHILD;
%                     ground.PREVIOUS = ground.CHILD;
%                     ground.CHILD = 0;
%                     
%                     % adjust interactions accordingly
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