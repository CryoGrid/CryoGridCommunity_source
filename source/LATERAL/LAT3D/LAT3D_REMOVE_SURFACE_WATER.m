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

        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of of LATERAL class timestep'};
        end
    end
end

