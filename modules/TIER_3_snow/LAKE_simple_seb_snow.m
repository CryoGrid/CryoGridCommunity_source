

classdef LAKE_simple_seb_snow < LAKE_simple_seb
    properties
        CHILD
        IA_CHILD
    end
    
    
    methods
        
        %mandatory functions for each class
        
        function self = LAKE_simple_seb_snow(index, pprovider, cprovider, forcing)
            self@LAKE_simple_seb(index, pprovider, cprovider, forcing);
        end
        
%        function ground = initialize_from_LAKE_unfrozen(ground, LAKE_simple_unfrozen)
%             ground = initialize_from_LAKE_unfrozen@LAKE_simple_seb(ground, LAKE_simple_unfrozen);
%             ground.CHILD = 0; % no snow
%             ground.IA_CHILD = 0;
%        end
         
       function ground = initialize_from_LAKE_previous_season(ground, LAKE_simple_unfrozen)
            ground = initialize_from_LAKE_previous_season@LAKE_simple_seb(ground, LAKE_simple_unfrozen);
            ground.CHILD = 0; % no snow
            ground.IA_CHILD = 0;
         end

       function ground = provide_PARA(ground)  %initializes the subvariables as empty arrays
            ground = provide_PARA@LAKE_simple_seb(ground);
       end
       
       function ground = provide_CONST(ground)  %initializes the subvariables as empty arrays
           ground = provide_CONST@LAKE_simple_seb(ground);
       end
       
       function ground = provide_STATVAR(ground)  %initializes the subvariables as empty arrays
           ground = provide_STATVAR@LAKE_simple_seb(ground);
       end
       
       function ground = finalize_init(ground, forcing)
           ground = finalize_init@LAKE_simple_seb(ground, forcing);
           ground.CHILD = 0; % no snow
           ground.IA_CHILD = 0;
       end

       
       %-------------------
        
        function ground = get_boundary_condition_u(ground, forcing)
            
            if ground.CHILD == 0  %CHILD does not exist
                ground = get_boundary_condition_u@LAKE_simple_seb(ground, forcing); %call the native function for the ground class
                
                if forcing.TEMP.snowfall > 0  %create CHILD 
                    CURRENT = ground.PREVIOUS;  %go to Top() and get the stored SNOW class
                    while ~strcmp(class(CURRENT), 'Top')
                        CURRENT = CURRENT.PREVIOUS;
                    end
                    ground.CHILD = copy(CURRENT.STORE.SNOW);
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground; 
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                    
                    ground.CHILD = get_boundary_condition_u_create_CHILD(ground.CHILD, forcing);  %initialize with fresh snowfall
                end
            else %CHILD exists
                total_area = ground.STATVAR.area(1,1); %store the total area of the ground
                ground.STATVAR.area(1,1) = ground.STATVAR.area(1,1) - ground.CHILD.STATVAR.area(1,1); %replace by snow-free area
                
                ground.CHILD.STATVAR.Lstar = ground.STATVAR.Lstar;
                
                ground.CHILD = get_boundary_condition_u_CHILD(ground.CHILD, forcing);
                ground = get_boundary_condition_u@LAKE_simple_seb(ground, forcing);
                
                get_IA_CHILD_boundary_condition_u(ground.IA_CHILD);
                %call designated mandatory function for CHILD-PARENT interactions in
                %the IA class governing IA between SNOW and GROUND
                
                ground.STATVAR.Lout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Lout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Lout) ./ total_area; %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Sout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Sout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Sout) ./ total_area;
                ground.STATVAR.Qh = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qh + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qh) ./ total_area;
                ground.STATVAR.Qe = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qe + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qe) ./ total_area;
                
                ground.STATVAR.area(1,1) = total_area; %reassign the true area of ground
            end
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
              ground = get_boundary_condition_l@LAKE_simple_seb(ground, forcing);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW@LAKE_simple_seb(ground, S_down);
        end
        
        function ground = get_derivatives_prognostic(ground)
            if ground.CHILD == 0  
                ground = get_derivatives_prognostic@LAKE_simple_seb(ground); %call normal function
            else
                ground.CHILD = get_derivatives_prognostic_CHILD(ground.CHILD);
                ground = get_derivatives_prognostic@LAKE_simple_seb(ground); 
            end
        end
        
        function timestep = get_timestep(ground) 
            if ground.CHILD == 0
                timestep = get_timestep@LAKE_simple_seb(ground);
            else 
                timestep_snow = get_timestep_CHILD(ground.CHILD);
                timestep_ground =  get_timestep@LAKE_simple_seb(ground);
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
            end
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            if ground.CHILD == 0
                ground =  advance_prognostic@LAKE_simple_seb(ground, timestep);
            else                
                ground.CHILD = advance_prognostic_CHILD(ground.CHILD, timestep);
                ground =  advance_prognostic@LAKE_simple_seb(ground, timestep);
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing);
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, forcing)
            if ground.CHILD == 0
                ground = compute_diagnostic@LAKE_simple_seb(ground, forcing);
            else
                ground = compute_diagnostic@LAKE_simple_seb(ground, forcing);
                ground.CHILD = compute_diagnostic_CHILD(ground.CHILD, forcing);
                
            end
        end
        
        function ground = check_trigger(ground, forcing)


            %snow trigger
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
                   
                    %make snow a real class
                    ground.CHILD.PARENT = 0;
                    ground.CHILD.PREVIOUS = ground.PREVIOUS;
                    ground.CHILD.NEXT = ground;
                    ground.PREVIOUS.NEXT = ground.CHILD;
                    ground.PREVIOUS = ground.CHILD;
                    ground.CHILD = 0;
                    ground.IA_PREVIOUS = ground.IA_CHILD; %should already point right
                    ground.PREVIOUS.IA_NEXT = ground.IA_CHILD;
                    ground.IA_CHILD = 0;
                end
            end
           %lake trigger 
           dummy = check_trigger@LAKE_simple_seb(ground, forcing);
        end
        
    end
end