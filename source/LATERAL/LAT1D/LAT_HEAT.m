%========================================================================
% CryoGrid LATERAL_IA class LAT_HEAT 

% S. Westermann, Oct 2020
%========================================================================


classdef LAT_HEAT < BASE_LATERAL

    
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
            lateral.PARA.reservoir_T = [];
            lateral.PARA.upperElevation = []; %Inf;
            lateral.PARA.lowerElevation = []; %20;
            lateral.PARA.distance_heatReservoir = []; %1; 
            lateral.PARA.heatReservoir_contact_length = []; %4;
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.heat_flux_accum = [];
        end

        function lateral = finalize_init(lateral, tile)

        end

        %----time integration------
        
        %only push function needed
        function lateral = push(lateral, tile)
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            lateral.TEMP.head = 0;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_heat(CURRENT, lateral);
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


