
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
            
%             point.PARA.slope_angle = [];     %
%             point.PARA.aspect = [];     %
%             point.PARA.skyview_factor = [];     %
%             point.PARA.horizon_angles = [];

            
        end
        
        function point = provide_STATVAR(point)

        end
        
        function point = provide_CONST(point)
            
        end
        
        function point = finalize_init(point)
            if size(point.STATVAR.latitude,1) == 1
                point.STATVAR.latitude = repmat(point.PARA.latitude, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.latitude = point.PARA.latitude;
            end
            if size(point.STATVAR.longitude,1) == 1
                point.STATVAR.longitude = repmat(point.PARA.longitude, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.longitude = point.PARA.longitude;
            end
            if size(point.STATVAR.altitude,1) == 1
                point.STATVAR.altitude = repmat(point.PARA.altitude, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.altitude = point.PARA.altitude;
            end
            if size(point.STATVAR.area,1) == 1
                point.STATVAR.area = repmat(point.PARA.area, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.area = point.PARA.area;
            end
%             point.STATVAR.slope_angle = point.PARA.slope_angle;     
%             point.STATVAR.aspect = point.PARA.aspect;                 
            

            
            
        end
        
 
        
    end
end

