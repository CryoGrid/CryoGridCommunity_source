%defines a regular grid in geographical coordinates, with fixed resolution

classdef CLUSTER_SLOPE_ASPECT < matlab.mixin.Copyable

    properties
        SPATIAL
        PARA
        CONST
        STATVAR
    end
    
    methods
        function cluster = provide_PARA(cluster)

        end
        
        function cluster = provide_STATVAR(cluster)

        end
        
        function cluster = provide_CONST(cluster)
            
        end
        
        function cluster = finalize_init(cluster)
            
        end
        
        function data = provide_data(cluster)
            data = [-cosd(cluster.SPATIAL.STATVAR.aspect) sind(cluster.SPATIAL.STATVAR.aspect) sind(cluster.SPATIAL.STATVAR.slope_angle)];
        end
         
    end
end

