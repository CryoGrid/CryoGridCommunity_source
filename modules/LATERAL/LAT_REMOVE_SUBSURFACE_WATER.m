
classdef LAT_REMOVE_SUBSURFACE_WATER < BASE_LATERAL

    
    methods
        
        function lateral = provide_CONST(lateral)
            
        end
        
        function lateral = provide_PARA(lateral)
            
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
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_remove_subsurfaceWater(CURRENT, lateral);
                CURRENT = compute_diagnostic(CURRENT, forcing);
                CURRENT = CURRENT.NEXT;
            end
            
        end
    end
    
end


