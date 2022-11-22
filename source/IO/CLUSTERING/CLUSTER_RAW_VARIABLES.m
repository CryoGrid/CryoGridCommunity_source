%defines a regular grid in geographical coordinates, with fixed resolution

classdef CLUSTER_RAW_VARIABLES < matlab.mixin.Copyable

    properties
        SPATIAL
        PARA
        CONST
        STATVAR
    end
    
    methods
        function cluster = provide_PARA(cluster)
            cluster.PARA.variables = [];
        end
        
        function cluster = provide_STATVAR(cluster)

        end
        
        function cluster = provide_CONST(cluster)
            
        end
        
        function cluster = finalize_init(cluster)
            
        end
        
        function data = provide_data(cluster)
            data = [];
            for i=1:size(cluster.PARA.variables,1)
                data = [data cluster.SPATIAL.STATVAR.(cluster.PARA.variables{i,1})];
            end
        end
         
    end
end

