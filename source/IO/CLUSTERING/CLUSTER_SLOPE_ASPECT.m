%========================================================================
% service class providing data for CLUSTERING classes 
% CLUSTER_SLOPE_ASPECT synthesizes slope and aspect data to that they can 
% be employed in the clustering 
%
% S. Westermann, Dec 2022
%========================================================================

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
         
        
        %-------------param file generation-----
        function data = param_file_info(data)
            data = provide_PARA(data);
            
            data.PARA.STATVAR = [];
            data.PARA.class_category = 'CLUSTERING DATA';
            data.PARA.default_value = [];
                        
        end
    end
end

