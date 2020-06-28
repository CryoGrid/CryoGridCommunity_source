
classdef LAT_SEEPAGE_FACE_WATER < BASE_LATERAL

    
    methods
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = 24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.upperElevation = Inf;
            lateral.PARA.lowerElevation = 20;
            lateral.PARA.hardBottom_cutoff = 0.03; %hard bottom if saturated and water content below
            lateral.PARA.distance_seepageFace = 10; 
            lateral.PARA.seepage_contact_length = 1;
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
            lateral.TEMP.head = 0;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_remove_water_seepage(CURRENT, lateral);
                CURRENT = compute_diagnostic(CURRENT, forcing);
                CURRENT = CURRENT.NEXT;
            end
            
        end
    end
    
end


