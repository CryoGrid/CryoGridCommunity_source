%========================================================================
% CryoGrid LATERAL_IA class LAT_HEAT 
% LAT_HEAT simulates lateral heat transfer with an external reservoir at
% constant temperature located a defined altitudinal interval

% S. Westermann, Aug 2021
%========================================================================


classdef LAT_HEAT < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        

        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = [];
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.reservoir_T = [];
            lateral.PARA.upperElevation = []; %Inf;
            lateral.PARA.lowerElevation = []; %20;
            lateral.PARA.distance_heatReservoir = []; %1; 
            lateral.PARA.heatReservoir_contact_length = []; %4;
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.heat_flux_accum = [];
        end

        function lateral = finalize_init(lateral, tile)

        end

        %----time integration------
        
        %only push function needed
        function lateral = push(lateral, tile)
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            lateral.TEMP.head = 0;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_heat(CURRENT, lateral);
                CURRENT = compute_diagnostic(CURRENT, tile);
                CURRENT = CURRENT.NEXT;
            end
            
        end
        
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 1;
        end
        
        function lateral = get_derivatives(lateral, tile)
            
        end
        
        function lateral = pull(lateral, tile)
            
        end
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];
            
            ground.PARA.default_value.reservoir_T = {''};
            ground.PARA.comment.reservoir_T = {'temperature of reservoir [degreeC]'};
            
            ground.PARA.default_value.upperElevation = {10000};
            ground.PARA.comment.upperElevation = {'upper elevation of reservoir [meter a.s.l.]'};
            
            ground.PARA.default_value.lowerElevation = {0};
            ground.PARA.comment.lowerElevation = {'lower elevation of reservoir [meter a.s.l.]'};
            
            ground.PARA.default_value.distance_heatReservoir = {1};
            ground.PARA.comment.distance_heatReservoir ={'distance to heat reservoir [m]'};
            
            ground.PARA.default_value.heatReservoir_contact_length = {1};
            ground.PARA.comment.heatReservoir_contact_length = {'lateral contact lengt to heat reservoir [m]'};
        end

        
    end
    
end


