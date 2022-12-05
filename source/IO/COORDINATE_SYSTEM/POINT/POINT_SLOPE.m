
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
            point.STATVAR.skyview_factor = point.PARA.skyview_factor;
            
            point.PARA.horizon_angles = point.PARA.horizon_angles';
            point.STATVAR.horizon_bins = point.PARA.horizon_angles(:,1);
            point.STATVAR.horizon_angles = point.PARA.horizon_angles(:,2);
            
            if isempty(point.PARA.horizon_angles) || sum(isnan(point.PARA.horizon_angles(:)))>0
                point.STATVAR.horizon_bins = 0;
                point.STATVAR.horizon_angles = 0;
            end
            if isempty(point.PARA.skyview_factor) || sum(isnan(point.PARA.skyview_factor))>0
                
                azmRadian = (pi/180).*point.STATVAR.horizon_bins;
                
                % convert output from horizon program to radians and translate to angle
                % from zenith
                H = (pi/180).*(90-point.STATVAR.horizon_angles(:));
                
                aspectRadian = (pi/180)*(point.STATVAR.aspect);
                % modify limits of integration for slopes facing away from horizons
                t = cosd(point.STATVAR.aspect-point.STATVAR.horizon_bins)<0;
                %Simplified trig, the original was H(t) = min(H(t),...
                %  acos(-cos(azmRadian(t)-aspectRadian)*sind(slopeDegrees)./...
                %  sqrt(cosd(slopeDegrees)^2+sind(slopeDegrees)^2*cos(azmRadian(t)-aspectRadian).^2)));
                % but same as
                H(t) = min(H(t), acos(sqrt(1-1./(1+tand(point.STATVAR.slope_angle)^2*cos(azmRadian(t)-aspectRadian).^2))));
                qIntegrand = (cosd(point.STATVAR.slope_angle)*sin(H).^2 + sind(point.STATVAR.slope_angle)*cos(aspectRadian-azmRadian).*(H-cos(H).*sin(H)))/2;
                
                % shouldn't be any negative, except perhaps rounding error, so just in case
                qIntegrand(qIntegrand<0) = 0;
                
                % integrate
                point.STATVAR.skyview_factor = trapz(azmRadian,qIntegrand)./pi;
                
            end
            
            
        end
        
 
        
    end
end

