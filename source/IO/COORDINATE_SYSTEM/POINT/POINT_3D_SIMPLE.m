
classdef POINT_3D_SIMPLE < matlab.mixin.Copyable

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
            point.PARA.number_of_tiles = []; %3;
            point.PARA.latitude = [];
            point.PARA.longitude = [];
            point.PARA.altitude = []; 
            point.PARA.area = [];
            point.PARA.param_file_number = []; %[1;2;3];

            point.PARA.connected = [];
            point.PARA.contact_length = [];
            point.PARA.distance = [];
            
            %provide default values
            point.PARA.slope_angle = 0;     
            point.PARA.aspect = 0;     
            point.PARA.skyview_factor = 0;     
            point.PARA.horizon_bins = 0;
            point.PARA.horizon_angles = 0;
            
        end
        
        function point = provide_STATVAR(point)

        end
        
        function point = provide_CONST(point)
            
        end
        
        function point = finalize_init(point)
            if size(point.PARA.latitude,1) == 1
                point.STATVAR.latitude = repmat(point.PARA.latitude, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.latitude = point.PARA.latitude;
            end
            if size(point.PARA.longitude,1) == 1
                point.STATVAR.longitude = repmat(point.PARA.longitude, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.longitude = point.PARA.longitude;
            end
            if size(point.PARA.altitude,1) == 1
                point.STATVAR.altitude = repmat(point.PARA.altitude, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.altitude = point.PARA.altitude;
            end
            if size(point.PARA.area,1) == 1
                point.STATVAR.area = repmat(point.PARA.area, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.area = point.PARA.area;
            end
            if size(point.PARA.slope_angle,1) == 1
                point.STATVAR.slope_angle = repmat(point.PARA.slope_angle, point.PARA.number_of_tiles, 1);
            else
                point.PARA.slope_angle = point.PARA.slope_angle;
            end
            if size(point.PARA.slope_angle,1) == 1
                point.STATVAR.slope_angle = repmat(point.PARA.slope_angle, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.slope_angle = point.PARA.slope_angle;
            end
            if size(point.PARA.aspect,1) == 1
                point.STATVAR.aspect = repmat(point.PARA.aspect, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.aspect = point.PARA.aspect;
            end
            if size(point.PARA.skyview_factor,1) == 1
                point.STATVAR.skyview_factor = repmat(point.PARA.skyview_factor, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.skyview_factor = point.PARA.skyview_factor;
            end
            if size(point.PARA.horizon_angles,1) == 1
                point.STATVAR.horizon_angles = repmat(point.PARA.horizon_angles, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.horizon_angles = point.PARA.horizon_angles;
            end
            if size(point.PARA.horizon_bins,1) == 1
                point.STATVAR.horizon_bins = repmat(point.PARA.horizon_bins, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.horizon_bins = point.PARA.horizon_bins;
            end
        end
        
 
        
    end
end

