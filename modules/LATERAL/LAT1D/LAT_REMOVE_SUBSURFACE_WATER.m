%========================================================================
% CryoGrid LATERAL_IA class which removes all subsurface water for GROUND classses with water bucket scheme, 
% all water exceeding the field capacity is removed
% S. Westermann, Oct 2020
%========================================================================

classdef LAT_REMOVE_SUBSURFACE_WATER < BASE_LATERAL

    
    methods
        
                
        %----mandatory functions---------------
        %----initialization--------------------
        
        function lateral = LAT_REMOVE_SUBSURFACE_WATER(index, pprovider, cprovider)
            lateral@BASE_LATERAL(index, pprovider, cprovider);
        end
        
        function lateral = provide_CONST(lateral)
            
        end
        
        function lateral = provide_PARA(lateral)
            
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
        end

        function lateral = finalize_init(lateral)
            lateral.STATVAR.subsurface_run_off = 0;
        end
        
        %----time integration-----
        
        %only push function required
        function lateral = push(lateral, forcing)
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_remove_subsurfaceWater(CURRENT, lateral);
                CURRENT = compute_diagnostic(CURRENT, forcing);
                CURRENT = CURRENT.NEXT;
            end
            
        end
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 1;
        end
        
        function lateral = pull(lateral)
            
        end
    end
    
end


