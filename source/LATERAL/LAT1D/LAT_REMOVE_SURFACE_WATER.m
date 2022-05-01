%========================================================================
% CryoGrid LATERAL_IA class LAT_REMOVE_SURFACE_WATER
% LAT_REMOVE_SURFACE_WATER which removes all surface water overtopping 
% the first grid cell for GROUND classses
% This can for example prevent the triggering of a LAKE in some GROUND
% classes, or to analyze the water balance.
% S. Westermann, Oct 2020
%========================================================================

classdef LAT_REMOVE_SURFACE_WATER < BASE_LATERAL
    
    methods
                
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function lateral = provide_CONST(lateral)
            
        end
        
        function lateral = provide_PARA(lateral)
            
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.surface_run_off = [];
        end

        
        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.surface_run_off = 0;
        end
        
        %------time integration-------------
        
        %only push function needed
        function lateral = push(lateral, tile)
            %remove water from first class in stratigraphy only
            TOP.NEXT = lateral_push_remove_surfaceWater(lateral.PARENT.TOP.NEXT, lateral); 
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

