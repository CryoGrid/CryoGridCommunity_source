
classdef LAT_WATER_RESERVOIR < BASE_LATERAL

    
    methods
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = 24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.reservoir_elevation = 20;
            lateral.PARA.reservoir_temperature = 0.5; %only active for Xice classes - if empty, water added at the temperature of the respective grid cell
            lateral.PARA.hardBottom_cutoff = 0.03; %hard bottom if saturated and water content below
            lateral.PARA.distance_reservoir = 3; 
            lateral.PARA.reservoir_contact_length = 1;
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
        end
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function lateral = finalize_init(lateral)
            lateral.STATVAR.subsurface_run_off = 0;
        end
        
        function lateral = push(lateral, forcing)
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            lateral.TEMP.open_system = 1; %start with open system
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_water_reservoir(CURRENT, lateral);
                CURRENT = compute_diagnostic(CURRENT, forcing);
                CURRENT = CURRENT.NEXT;
            end
            
        end
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 1;
        end
        
        function lateral = get_derivatives(lateral)
            
        end
        
        function lateral = pull(lateral)
            
        end
    end
    
end


