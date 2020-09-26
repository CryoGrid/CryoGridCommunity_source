classdef LAT_REMOVE_SURFACE_WATER < BASE_LATERAL
    
    methods
        function self = LAT_REMOVE_SURFACE_WATER(index, pprovider, cprovider)  
            self@BASE_LATERAL(index, pprovider, cprovider);
        end
        
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
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 1;
        end
        
        function lateral = get_derivatives(lateral)
            
        end
        
        function lateral = pull(lateral)
            
        end

    end
end

