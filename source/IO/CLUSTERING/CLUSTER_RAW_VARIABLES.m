%========================================================================
% service class providing data for CLUSTERING classes 
% CLUSTER_RAW_VARIABLES provides unchanged raw data 
%
% S. Westermann, Dec 2022
%========================================================================

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
         
        
        
        %-------------param file generation-----
        function data = param_file_info(data)
            data = provide_PARA(data);
            
            data.PARA.STATVAR = [];
            data.PARA.class_category = 'CLUSTERING DATA';
            data.PARA.default_value = [];
            
            data.PARA.comment.variables = {'list of variables to be used for clustering'};
            data.PARA.options.variables.name = 'H_LIST';
            data.PARA.options.variables.entries_x = {'latitude' 'longitude' 'altitude'};
                        
        end
    end
end

