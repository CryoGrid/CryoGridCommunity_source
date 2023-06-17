%========================================================================
% CryoGrid SPATIAL_REFERENCE class POINT_DEM
% POINT class deriving information for a single target point from a digital 
% elevation model (DEM),including slope, aspect and terrain shading. The 
% class can ingest Copernicus 30m DEMs downloaded from 
% https://portal.opentopography.org/   
%
% S. Westermann, Dec 2022
%========================================================================


classdef POINT_DEM < DEM_BASE
    properties
        RUN_INFO
    end
    
    methods
        function point = provide_PARA(point)
            point.PARA.latitude = [];
            point.PARA.longitude = [];
            point.PARA.area = [];
            
            point.PARA.variables = []; % altitude or altitude, slope_angle, aspect or altitude, slope_angle, aspect, horizon_angles
            
            point.PARA.number_of_horizon_bins = []; %multiples of 4!
            point.PARA.DEM_folder = [];
            point.PARA.DEM_filename = [];
            point.PARA.reproject2utm = 1; %select 1 when using a geographic coordinate system and computing more than just altitude; select 0 to speed up altitde computation in big DEMs
        end
        
        function point = provide_STATVAR(point)

        end
        
        function point = provide_CONST(point)
            
        end
        
        function point = finalize_init(point)

            point.STATVAR.latitude = point.PARA.latitude;
            point.STATVAR.longitude = point.PARA.longitude;
            point.STATVAR.area = point.PARA.area;
            point.STATVAR.altitude = 0;
            point.STATVAR.slope_angle = 0;
            point.STATVAR.aspect = 0;
            point.STATVAR.skyview_factor = 0;
            point.STATVAR.horizon_bins = 0;
            point.STATVAR.horizon_angles = 0;
            
            point = read_DEM_raster(point);
            point = project_target_coordinates(point);
            point = compute_global_offset_from_north(point);
            
%             point.TEMP.Z = readgeoraster([point.PARA.DEM_folder point.PARA.DEM_filename]);
%             info = georasterinfo([point.PARA.DEM_folder point.PARA.DEM_filename]);
%             point.TEMP.info = info;
%             point.TEMP.Z(point.TEMP.Z==info.MissingDataIndicator) = NaN;
%             xlim = info.RasterReference.XWorldLimits; 
%             ylim = info.RasterReference.YWorldLimits;
%             dx = info.RasterReference.CellExtentInWorldX; 
%             dy = info.RasterReference.CellExtentInWorldY;
%             x=((xlim(1)+dx/2):dx:(xlim(end)-dx/2))';
%             y=flipud(((ylim(1)+dy/2):dy:(ylim(end)-dy/2))');
%             [point.TEMP.X,point.TEMP.Y] = meshgrid(x,y);  
%             [point.TEMP.X_target, point.TEMP.Y_target] = projfwd(info.CoordinateReferenceSystem, point.PARA.latitude, point.PARA.longitude);
%             %[LAT,LON]=projinv(info.CoordinateReferenceSystem,X,Y);
%             
%             
%             [X1, Y1] = projfwd(point.TEMP.info.CoordinateReferenceSystem, point.PARA.latitude-0.1, point.PARA.longitude);
%             [X2, Y2] = projfwd(point.TEMP.info.CoordinateReferenceSystem, point.PARA.latitude+0.1, point.PARA.longitude);
%             point.TEMP.offset_angle_trueNorth = atand((X1-X2)./(Y1-Y2)); %offset grid North and true North 
%             
            for i=1:size(point.PARA.variables,1)
                a = str2func(['get_' point.PARA.variables{i,1}]);
                point = a(point);
            end

        end
        
        function point = get_altitude(point)
            point = get_altitude_base(point);
            %point.STATVAR.altitude = interp2(point.TEMP.X, point.TEMP.Y, point.TEMP.Z, point.TEMP.X_target, point.TEMP.Y_target);
            
        end
       
        function point = get_slope_angle(point)            
            
            point = get_slope_angle_base(point);
            %Elevation derivatives. N.B. not computed for the boundary of the DEM.
%             dzdy=NaN.*point.TEMP.Z;  
%             dzdx=dzdy;
%             dzdy(2:end-1,2:end-1)=-1.*(point.TEMP.Z(3:end,2:end-1)-point.TEMP.Z(1:end-2,2:end-1))./(2.*point.TEMP.info.RasterReference.CellExtentInWorldY); % N.B. y decreases with i.
%             dzdx(2:end-1,2:end-1)=(point.TEMP.Z(2:end-1,3:end)-point.TEMP.Z(2:end-1,1:end-2))./(2.*point.TEMP.info.RasterReference.CellExtentInWorldX);
%             
%             slope_angle = atand(sqrt(dzdx.^2+dzdy.^2));
%             point.STATVAR.slope_angle = interp2(point.TEMP.X, point.TEMP.Y, slope_angle, point.TEMP.X_target, point.TEMP.Y_target);

        end
        
        function point = get_aspect(point)
            
            point = get_aspect_base(point);
%             %Elevation derivatives. N.B. not computed for the boundary of the DEM.
%             dzdy=NaN.*point.TEMP.Z;  
%             dzdx=dzdy;
%             dzdy(2:end-1,2:end-1)=-1.*(point.TEMP.Z(3:end,2:end-1)-point.TEMP.Z(1:end-2,2:end-1))./(2.*point.TEMP.info.RasterReference.CellExtentInWorldY); % N.B. y decreases with i.
%             dzdx(2:end-1,2:end-1)=(point.TEMP.Z(2:end-1,3:end)-point.TEMP.Z(2:end-1,1:end-2))./(2.*point.TEMP.info.RasterReference.CellExtentInWorldX);
%             
%             % Aspect definition: N,E,S,W = pi, pi/2, 0, -pi/2
%             aspect = atan2d(-1.*dzdx,dzdy);
%             
%             point.STATVAR.aspect = interp2(point.TEMP.X, point.TEMP.Y, aspect, point.TEMP.X_target, point.TEMP.Y_target) + point.TEMP.offset_angle_trueNorth; %check if + or -!!!
        end
        
        function point = get_horizon_angles(point)
            
            point = get_horizon_angles_single_point(point);
% %             [X1, Y1] = projfwd(point.TEMP.info.CoordinateReferenceSystem, point.PARA.latitude-0.1, point.PARA.longitude);
% %             [X2, Y2] = projfwd(point.TEMP.info.CoordinateReferenceSystem, point.PARA.latitude+0.1, point.PARA.longitude);
% %             rot_angle = atand((X1-X2)./(Y1-Y2)); %offset grid North and true North 
%             
% %             distance2target = (point.TEMP.X - point.TEMP.X_target).^2 + (point.TEMP.Y - point.TEMP.Y_target).^2;
% %             [ii,jj]=find(distance2target==min(distance2target(:)));
% %             
% %             grad_X = (point.TEMP.Z(ii,:) - point.TEMP.Z(ii,jj)) ./ abs(point.TEMP.X(ii,:) - point.TEMP.X(ii,jj));
% %             grad_Y = (point.TEMP.Z(:,jj) - point.TEMP.Z(ii,jj)) ./ abs(point.TEMP.Y(:,jj) - point.TEMP.Y(ii,jj));
% % 
% %             point.STATVAR.horizon_bins = [0; 90; 180; 270];
% %             point.STATVAR.horizon_angles = [atand(max(grad_Y(ii+1:end))); atand(max(grad_X(1:jj-1))); atand(max(grad_Y(1:ii-1))); atand(max(grad_X(jj+1:end)))];
% 
%             point.STATVAR.horizon_bins =[];
%             point.STATVAR.horizon_angles =[];
% 
%            % rotation_angle_increment = 90./2.^(log(point.PARA.number_of_horizon_bins)./log(2)-2);
%             rotation_angle_increment = 90./(point.PARA.number_of_horizon_bins./4);
%              
%             for rotation_angle = 0:rotation_angle_increment:90-1e-9
%                 X=imrotate(point.TEMP.X, rotation_angle+point.TEMP.offset_angle_trueNorth);
%                 Y=imrotate(point.TEMP.Y, rotation_angle+point.TEMP.offset_angle_trueNorth);
%                 Z=imrotate(point.TEMP.Z, rotation_angle+point.TEMP.offset_angle_trueNorth);
%                 
%                 distance2target = (X - point.TEMP.X_target).^2 + (Y - point.TEMP.Y_target).^2;
%                 [ii,jj]=find(distance2target==min(distance2target(:)));
%                 hor_angles=[];
%                 for i=1:size(ii,1)
%                     grad_X =  (Z(ii(i,1),:) - Z(ii(i,1),jj(i,1))) ./ sqrt((X(ii(i,1),:) - X(ii(i,1),jj(i,1))).^2 + (Y(ii(i,1),:) - Y(ii(i,1),jj(i,1))).^2);
%                     grad_Y = (Z(:,jj(i,1)) - Z(ii(i,1),jj(i,1))) ./ sqrt((Y(:,jj(i,1)) - Y(ii(i,1),jj(i,1))).^2 + (X(:,jj(i,1)) - X(ii(i,1),jj(i,1))).^2);
%                     hor_angles = [hor_angles [atand(max(grad_Y(ii+1:end))); atand(max(grad_X(1:jj-1))); atand(max(grad_Y(1:ii-1))); atand(max(grad_X(jj+1:end)))]];
% 
%                 end
%                 
%                 point.STATVAR.horizon_bins = [point.STATVAR.horizon_bins ; rotation_angle + [0; 90; 180; 270]];
%                 point.STATVAR.horizon_angles = [point.STATVAR.horizon_angles; mean(hor_angles,2)];
%                
%             end
% 
%             [point.STATVAR.horizon_bins, order] = sort(point.STATVAR.horizon_bins);
%             point.STATVAR.horizon_angles = point.STATVAR.horizon_angles(order,1);
%             point.STATVAR.horizon_angles(point.STATVAR.horizon_angles<0) = 0;
%         
%             %skyview factor
%             azmRadian = (pi/180).*point.STATVAR.horizon_bins;
%             
%             % convert output from horizon program to radians and translate to angle
%             % from zenith
%             H = (pi/180).*(90-point.STATVAR.horizon_angles(:));
% 
%             aspectRadian = (pi/180)*(point.STATVAR.aspect);
%             % modify limits of integration for slopes facing away from horizons
%             t = cosd(point.STATVAR.aspect-point.STATVAR.horizon_bins)<0;
%             %Simplified trig, the original was H(t) = min(H(t),...
%             %  acos(-cos(azmRadian(t)-aspectRadian)*sind(slopeDegrees)./...
%             %  sqrt(cosd(slopeDegrees)^2+sind(slopeDegrees)^2*cos(azmRadian(t)-aspectRadian).^2)));
%             % but same as
%             H(t) = min(H(t), acos(sqrt(1-1./(1+tand(point.STATVAR.slope_angle)^2*cos(azmRadian(t)-aspectRadian).^2))));
%             qIntegrand = (cosd(point.STATVAR.slope_angle)*sin(H).^2 + sind(point.STATVAR.slope_angle)*cos(aspectRadian-azmRadian).*(H-cos(H).*sin(H)))/2;
%             
%             % shouldn't be any negative, except perhaps rounding error, so just in case
%             qIntegrand(qIntegrand<0) = 0;
%             
%             % integrate
%             point.STATVAR.skyview_factor = trapz(azmRadian,qIntegrand)./pi;
%             
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
            
            point.PARA.comment.area = {'area of target point in m2'};
            point.PARA.default_value.area = {1};
            
            point.PARA.comment.variables = {'properties calculated from DEM: altitude OR altitude, slope_angle, aspect OR altitude, slope_angle, aspect, horizon_angles'};
            point.PARA.options.variables.name = 'H_LIST';
            point.PARA.options.variables.entries_x = {'altitude' 'slope_angle' 'aspect' 'horizon_angles'};
                        
            point.PARA.comment.number_of_horizon_bins = {'number of angular points for which horizon is calculated; must be multiple of 4'};  
            point.PARA.default_value.number_of_horizon_bins = {24};
            
            point.PARA.comment.DEM_folder = {'folder in which DEM file is located'}; 
            
            point.PARA.comment.DEM_filename = {'name of DEM file'}; 
            
            point.PARA.comment.reproject2utm = {'select 1 when using a DEM in geographic coordinates (or similar) and computing more than just altitude; select 0 to speed up altitde computation in big DEMs'};
            point.PARA.default_value.reproject2utm = {1};
        end
        
    end
end

