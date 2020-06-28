classdef LAT_REMOVE_SURFACE_WATER < BASE_LATERAL
    properties
        
    end
    
    methods
        
        function lateral = provide_CONST(lateral)
            
        end
        
        function lateral = provide_PARA(lateral)
            
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.surface_run_off = [];
        end
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function lateral = finalize_init(lateral)
            lateral.STATVAR.surface_run_off = 0;
        end
        

        
        function lateral = push(lateral, forcing)
            %remove water from first class in stratigraphy only
            TOP.NEXT = lateral_push_remove_surfaceWater(lateral.PARENT.TOP.NEXT, lateral); 
            
        end

    end
end

