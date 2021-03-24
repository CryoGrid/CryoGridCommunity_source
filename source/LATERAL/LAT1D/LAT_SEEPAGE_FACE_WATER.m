%========================================================================
% CryoGrid LATERAL_IA class LAT_SEEPAGE_FACE_WATER 
% simulates lateral water flow through a seepage face with defined upper
% and lower elevation (absolute elevation, not relative to the surface!), 
%  as well as contact lengths (i.e. width) and distance from the GROUND column.
% At the seepage face, air pressure is assumed.
% S. Westermann, Oct 2020
%========================================================================


classdef LAT_SEEPAGE_FACE_WATER < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
%         function lateral = LAT_SEEPAGE_FACE_WATER(index, pprovider, cprovider)
%             lateral@BASE_LATERAL(index, pprovider, cprovider);
%         end

        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = [];
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.upperElevation = []; %Inf;
            lateral.PARA.lowerElevation = []; %20;
            lateral.PARA.hardBottom_cutoff = []; %0.03; %hard bottom if saturated and water content below
            lateral.PARA.distance_seepageFace = []; %1; 
            lateral.PARA.seepage_contact_length = []; %4;
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
        end

        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.subsurface_run_off = 0;
        end

        %----time integration------
        
        %only push function needed
        function lateral = push(lateral, tile)
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            lateral.TEMP.head = 0;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_remove_water_seepage(CURRENT, lateral);
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
    end
    
end


