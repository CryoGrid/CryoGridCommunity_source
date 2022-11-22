
classdef POINT_SLOPE < matlab.mixin.Copyable

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
            point.PARA.slope_angle = [];     %
            point.PARA.aspect = [];     %
            point.PARA.skyview_factor = [];     %
            point.PARA.horizon_angles = [];     
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
            point.STATVAR.slope_angle = point.PARA.slope_angle;     
            point.STATVAR.aspect = point.PARA.aspect;     
            
            point.PARA.horizon_angles = point.PARA.horizon_angles';
            
            if isempty(point.PARA.horizon_angles) || sum(isnan(point.PARA.horizon_angles(:)))>0
                point.STATVAR.skyview_factor = point.PARA.skyview_factor
                point.STATVAR.horizon_bins = 0;
                point.STATVAR.horizon_angles = point.PARA.skyview_factor .* 90;
            end
            if isempty(point.PARA.skyview_factor) || sum(isnan(point.PARA.skyview_factor))>0
                angle_slices = [point.PARA.horizon_angles(:,1); 360];
                angle_fraction = [angle_slices(2:end,1) - angle_slices(1:end-1,1)] ./ 360;
                point.STATVAR.skyview_factor = sum(angle_fraction .* point.PARA.horizon_angles(:,2)) ./90;
                point.STATVAR.horizon_bins = point.PARA.horizon_angles(:,1);
                point.STATVAR.horizon_angles = point.PARA.horizon_angles(:,2);
            end
            
            
        end
        
 
        
    end
end

