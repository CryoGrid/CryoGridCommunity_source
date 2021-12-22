%========================================================================
% CryoGrid LATERAL_IA class LAT3D_HEAT 
% simulates lateral conductive heat fluxes between different CryoGrid tiles. 
% S. Westermann, Oct 2020
%========================================================================

classdef LAT3D_HEAT < BASE_LATERAL

    
    methods
            
        %----mandatory functions---------------
        %----initialization--------------------

        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.ia_time_increment = []; %0.25; %must be a multiple of the time increment of the main lateral class
            %lateral.PARA.ia_time_next = [];
        end
        
        function lateral = provide_STATVAR(lateral)
            
        end

        function lateral = finalize_init(lateral, tile)
            %lateral.STATVAR.subsurface_run_off = 0;
        end
        
        %----time integration----------------- 
        
        function lateral = pull(lateral, tile)
            
            lateral.PARENT.STATVAR.depths_heat = [];
            lateral.PARENT.STATVAR.thermCond = [];
            lateral.PARENT.STATVAR.T_heat = [];

            CURRENT = lateral.PARENT.TOP.NEXT;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral3D_pull_heat(CURRENT, lateral);
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function lateral = get_derivatives(lateral, tile) %no need to loop through stratigraphy, all the information is in lateral.PARENT

            lateral.PARENT = get_overlap_cells(lateral.PARENT, 'depths_heat', 'overlap_heat');
            
            %calculate fluxes
            flux_heat = lateral.PARENT.STATVAR.thermCond .* 0;
 
            for j=1:size(lateral.PARENT.ENSEMBLE,1)
                if lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index)  %add if water is available at all
                    for i=1:size(lateral.PARENT.ENSEMBLE{j,1}.overlap_heat,1) %does not do anything when overlap is empty!
                        
                        tc1 = lateral.PARENT.STATVAR.thermCond(lateral.PARENT.ENSEMBLE{j,1}.overlap_heat(i,1),1);
                        tc2 = lateral.PARENT.ENSEMBLE{j,1}.thermCond(lateral.PARENT.ENSEMBLE{j,1}.overlap_heat(i,2),1);
                        distance = lateral.PARENT.PARA.distance(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                        thermCond = tc1 .* tc2 ./ (tc1 .* distance./2 + tc2 .* distance./2); %change to different distances laters
                        thermCond(isnan(thermCond)) = 0;
                        contact_length = lateral.PARENT.PARA.contact_length(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                        
                        flux_i = contact_length .* lateral.PARENT.ENSEMBLE{j,1}.overlap_heat(i,3) .* thermCond .* ...
                            -(lateral.PARENT.STATVAR.T_heat(lateral.PARENT.ENSEMBLE{j,1}.overlap_heat(i,1),1) - lateral.PARENT.ENSEMBLE{j,1}.T_heat(lateral.PARENT.ENSEMBLE{j,1}.overlap_heat(i,2),1));

                        flux_heat(lateral.PARENT.ENSEMBLE{j,1}.overlap_heat(i,1),1) = flux_heat(lateral.PARENT.ENSEMBLE{j,1}.overlap_heat(i,1),1) + flux_i;                        
                    end
                end
            end
            lateral.PARENT.STATVAR.heat_flux = flux_heat .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;            
        end

        
        function lateral = push(lateral, tile)
            
            CURRENT = lateral.PARENT.TOP.NEXT; %find correct stratigraphy class
            
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral3D_push_heat(CURRENT, lateral);
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 0;
            if t + lateral.PARENT.IA_TIME_INCREMENT >= lateral.PARA.ia_time_next -1e-7
                lateral.PARENT.ACTIVE(i,1) = 1;
                %lateral.PARA.ia_time_next = t + lateral.PARENT.IA_TIME_INCREMENT + lateral.PARA.ia_time_increment;
                lateral.PARA.ia_time_next = lateral.PARA.ia_time_next + lateral.PARA.ia_time_increment;
            end
        end
        
        function lateral = set_ia_time(lateral, t)
            lateral.PARA.ia_time_next = t;
        end
        
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of of LATERAL class timestep'};
        end
    end
    
end


