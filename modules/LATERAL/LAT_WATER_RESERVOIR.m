
classdef LAT_WATER_RESERVOIR < BASE_LATERAL
    
    methods
        function self = LAT_WATER_RESERVOIR(index, pprovider, cprovider)  
            self@BASE_LATERAL(index, pprovider, cprovider);
        end
        
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = [];
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.reservoir_elevation = [];
            lateral.PARA.reservoir_temperature = []; %only active for Xice classes - if empty, water added at the temperature of the respective grid cell
            lateral.PARA.hardBottom_cutoff = []; %hard bottom if saturated and water content below
            lateral.PARA.distance_reservoir = []; 
            lateral.PARA.reservoir_contact_length = [];
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


