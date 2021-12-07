%========================================================================
% CryoGrid LATERAL_IA class LAT_OVERLAND_FLOW 
% simulates overland flow accroding to Gauckler-Manning equation
% S. Westermann, June 2021
%========================================================================


classdef LAT_OVERLAND_FLOW < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        

        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = [];
            lateral.CONST.c_w = [];
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.gradient = [];
            lateral.PARA.GaMa_coefficient = []; %Gauckler-Manning coefficient, https://en.wikipedia.org/wiki/Manning_formula
            lateral.PARA.overland_flow_contact_length = []; %40;
            lateral.PARA.ia_time_increment = []; 
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.surface_flow = [];
        end

        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.surface_flow = 0;
        end

        %----time integration------
        
        %only push function needed
        function lateral = push(lateral, tile)
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            %lateral.STATVAR.water_depth = 0;
            %while ~(strcmp(class(CURRENT), 'Bottom'))
            CURRENT = lateral_push_remove_water_overland_flow(CURRENT, lateral);
            CURRENT = compute_diagnostic(CURRENT, tile);
            %CURRENT = CURRENT.NEXT;
            %end
            
        end
        
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 1;
        end
        
        function lateral = get_derivatives(lateral, tile)
            
        end
        
        function lateral = pull(lateral, tile)
            
        end
    end
    
end


