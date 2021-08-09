%========================================================================
% CryoGrid LATERAL_IA class which removes all surface water overtopping the first grid cell for GROUND classses
% This can for example prevent the formation of a LAKE in some GROUND
% classes.
% S. Westermann, Oct 2020
%========================================================================

classdef LAT3D_REMOVE_SURFACE_WATER < BASE_LATERAL
    
    methods
                
        %----mandatory functions---------------
        %----initialization--------------------
        
%         function self = LAT_REMOVE_SURFACE_WATER(index, pprovider, cprovider)  
%             self@BASE_LATERAL(index, pprovider, cprovider);
%         end
        
        function lateral = provide_CONST(lateral)
            
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.ia_time_increment = [];
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
        
%         function lateral = set_ACTIVE(lateral, i, t)
%             lateral.PARENT.ACTIVE(i,1) = 1;
%         end
        
        function lateral = get_derivatives(lateral, tile)
            
        end
        
        function lateral = pull(lateral, tile)
            
        end
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 0;
            if t + lateral.PARENT.IA_TIME_INCREMENT >= lateral.PARA.ia_time_next -1e-9
                lateral.PARENT.ACTIVE(i,1) = 1;
                lateral.PARA.ia_time_next = t + lateral.PARENT.IA_TIME_INCREMENT + lateral.PARA.ia_time_increment;
            end
        end
        
        function lateral = set_ia_time(lateral, t)
            lateral.PARA.ia_time_next = t;
        end

    end
end

