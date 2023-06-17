%========================================================================
% CryoGrid SPATIAL_REFERENCE class POINT_SLOPE
% POINT class providing information for a single target point on a slope.
% Recommended for simple simulations in sloping terrain.
%
% S. Westermann, Dec 2022
%========================================================================

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
            %append 360 degree
            point.STATVAR.horizon_bins = [point.PARA.horizon_angles(:,1); 360];
            point.STATVAR.horizon_angles = [point.PARA.horizon_angles(:,2); point.PARA.horizon_angles(1,2)];
            
            if isempty(point.PARA.horizon_angles) || sum(isnan(point.PARA.horizon_angles(:)))>0
                point.STATVAR.horizon_bins = [0; 360];
                point.STATVAR.horizon_angles = [0; 0];
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
        
 
        %-------------param file generation-----
        function point = param_file_info(point)
            point = provide_PARA(point);
            
            point.PARA.STATVAR = [];
            point.PARA.class_category = 'SPATIAL_REFERENCE';
            
            point.PARA.comment.latitude = {'latitude in decimal degrees'};
            point.PARA.default_value.latitude = {78.9};
            
            point.PARA.comment.longitude = {'longitude in decimal degrees'};
            point.PARA.default_value.longitude = {11.1};
            
            point.PARA.comment.altitude = {'altitude in m a.s.l.'};
            point.PARA.default_value.altitude = {10};
            
            point.PARA.comment.area = {'area of target point in m2'};
            point.PARA.default_value.area = {1};
            
            point.PARA.comment.slope_angle = {'angle of the slope in degrees (0 degrees = flat)'};  
            point.PARA.default_value.slope_angle = {10};
            
            point.PARA.comment.aspect = {'aspect of the slope in degrees (0: S; 90: E; 180: N; 270: W)'};  ;  
            point.PARA.default_value.aspect = {180};
            
            point.PARA.comment.skyview_factor = {'skyview factor, if empty calculated from slope and horizon angles (recommeded in most situations)'};    
            
            point.PARA.comment.horizon_angles = {'horizon angles (0: flat, no terrain shading; 90: all direct raiation blocked) at selected angles, linear interolation between; the default corresponds to shading in S direction for sun angles < 30 degrees'};
            point.PARA.options.horizon_angles.name = 'MATRIX';
            point.PARA.options.horizon_angles.entries_matrix = {'0' '90' '180' '270'; '30' '0' '0' '0'}; 
            
        end
    end
end

