%========================================================================
% CryoGrid SPATIAL_REFERENCE class POINT_SIMPLE
% POINT class providing a minimum of information (geographic coordinates, 
% altitude, area) for a single target point. 
% Recommended for simple simulations. 
%
% S. Westermann, Dec 2022
%========================================================================

classdef POINT_SIMPLE < matlab.mixin.Copyable

    properties
        RUN_INFO
        PARA
        CONST
        STATVAR
        TEMP
        ACTION
    end
    
    methods
        function point = provide_PARA(point)
            point.PARA.latitude = [];
            point.PARA.longitude = [];
            point.PARA.altitude = [];
            point.PARA.area = [];
        end
        
        function point = provide_STATVAR(point)

        end
        
        function point = provide_CONST(point)
            
        end
        
        function point = finalize_init(point)
            point.STATVAR.latitude = point.PARA.latitude;
            point.STATVAR.longitude = point.PARA.longitude;
            point.STATVAR.altitude = point.PARA.altitude;
            point.STATVAR.area = point.PARA.area;
            
            point.STATVAR.slope_angle = 0;     %
            point.STATVAR.aspect = 0;     %
            point.STATVAR.skyview_factor = 0;     %
            point.STATVAR.horizon_angles = 0;   
            point.STATVAR.horizon_bins = 0;
        end
        

        %-------------param file generation-----
        function point = param_file_info(point)
            point = provide_PARA(point);
            
            point.PARA.STATVAR = [];
            point.PARA.class_category = 'SPATIAL_REFERENCE';
            point.PARA.options = [];
            
            point.PARA.comment.latitude = {'latitude in decimal degrees'};
            point.PARA.default_value.latitude = {78.9};
            
            point.PARA.comment.longitude = {'longitude in decimal degrees'};
            point.PARA.default_value.longitude = {11.1};
            
            point.PARA.comment.altitude = {'altitude in m a.s.l.'};
            point.PARA.default_value.altitude = {10};
            
            point.PARA.comment.area = {'area of target point in m2'};
            point.PARA.default_value.area = {1};
        end
        
    end
end

