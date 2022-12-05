
classdef DEM_BASE < matlab.mixin.Copyable

    properties
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods
        function dem = provide_PARA(dem)

        end
        
        function dem = provide_STATVAR(dem)

        end
        
        function dem = provide_CONST(dem)
            
        end
        
        function dem = read_DEM_raster(dem)
            
            disp('read DEM')
            dem.TEMP.Z = readgeoraster([dem.PARA.DEM_folder dem.PARA.DEM_filename]);
            info = georasterinfo([dem.PARA.DEM_folder dem.PARA.DEM_filename]);
            dem.TEMP.info = info;
            dem.TEMP.Z(dem.TEMP.Z==info.MissingDataIndicator) = NaN;
            xlim = info.RasterReference.XWorldLimits; 
            ylim = info.RasterReference.YWorldLimits;
            dx = info.RasterReference.CellExtentInWorldX; 
            dy = info.RasterReference.CellExtentInWorldY;
            x=((xlim(1)+dx/2):dx:(xlim(end)-dx/2))';
            y=flipud(((ylim(1)+dy/2):dy:(ylim(end)-dy/2))');
            [dem.TEMP.X,dem.TEMP.Y] = meshgrid(x,y);  
            
        end
        
        function dem = project_target_coordinates(dem)
            
            [dem.TEMP.X_target, dem.TEMP.Y_target] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARA.latitude, dem.PARA.longitude);
            %[LAT,LON]=projinv(info.CoordinateReferenceSystem,X,Y);
        end
        
        function dem = compute_global_offset_from_north(dem)
            
            [X1, Y1] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARA.latitude-0.1, dem.PARA.longitude);
            [X2, Y2] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARA.latitude+0.1, dem.PARA.longitude);
            dem.TEMP.offset_angle_trueNorth = atand((X1-X2)./(Y1-Y2)); %offset grid North and true North 

        end
        
        function dem = get_altitude_base(dem)
            disp('get altitude')
            dem.STATVAR.altitude = double(interp2(dem.TEMP.X, dem.TEMP.Y, dem.TEMP.Z, dem.TEMP.X_target, dem.TEMP.Y_target));
            
        end
       
        function dem = get_slope_angle_base(dem)
            disp('get slope')
            %Elevation derivatives. N.B. not computed for the boundary of the DEM.
            dzdy=NaN.*dem.TEMP.Z;  
            dzdx=dzdy;
            dzdy(2:end-1,2:end-1)=-1.*(dem.TEMP.Z(3:end,2:end-1)-dem.TEMP.Z(1:end-2,2:end-1))./(2.*dem.TEMP.info.RasterReference.CellExtentInWorldY); % N.B. y decreases with i.
            dzdx(2:end-1,2:end-1)=(dem.TEMP.Z(2:end-1,3:end)-dem.TEMP.Z(2:end-1,1:end-2))./(2.*dem.TEMP.info.RasterReference.CellExtentInWorldX);
            
            slope_angle = atand(sqrt(dzdx.^2+dzdy.^2));
            dem.STATVAR.slope_angle = double(interp2(dem.TEMP.X, dem.TEMP.Y, slope_angle, dem.TEMP.X_target, dem.TEMP.Y_target));

        end
        
        function dem = get_aspect_base(dem)
            
            disp('get aspect')
            %Elevation derivatives. N.B. not computed for the boundary of the DEM.
            dzdy=NaN.*dem.TEMP.Z;  
            dzdx=dzdy;
            dzdy(2:end-1,2:end-1)=-1.*(dem.TEMP.Z(3:end,2:end-1)-dem.TEMP.Z(1:end-2,2:end-1))./(2.*dem.TEMP.info.RasterReference.CellExtentInWorldY); % N.B. y decreases with i.
            dzdx(2:end-1,2:end-1)=(dem.TEMP.Z(2:end-1,3:end)-dem.TEMP.Z(2:end-1,1:end-2))./(2.*dem.TEMP.info.RasterReference.CellExtentInWorldX);
            
            % Aspect definition: N,E,S,W = pi, pi/2, 0, -pi/2
            % Aspect definition: N,E,S,W = 180, 90, 0, 270
%             aspect = atan2d(-1.*dzdx,dzdy);
%             aspect(aspect<0)=aspect(aspect<0)+360;
%             
% %             imagesc(dem.TEMP.X(1,:),dem.TEMP.Y(:,1), aspect)
% %             hold on
% %             plot(dem.TEMP.X_target, dem.TEMP.Y_target, '+')
%             
%             dem.STATVAR.aspect = interp2(dem.TEMP.X, dem.TEMP.Y, aspect, dem.TEMP.X_target, dem.TEMP.Y_target) + dem.TEMP.offset_angle_trueNorth; %check if + or -!!!
% %             dem.STATVAR.aspect(dem.STATVAR.aspect<0) = dem.STATVAR.aspect(dem.STATVAR.aspect<0) + 360;

            %interpolate first ebfore taking atangens, otherwise problems
            %interpolating between 0 and 360 can create artificial slopes
            dzdy2 = interp2(dem.TEMP.X, dem.TEMP.Y, dzdy, dem.TEMP.X_target, dem.TEMP.Y_target); %check if + or -!!!
            dzdx2 = interp2(dem.TEMP.X, dem.TEMP.Y, dzdx, dem.TEMP.X_target, dem.TEMP.Y_target); %check if + or -!!!
            aspect = atan2d(-1.*dzdx2,dzdy2);
            aspect(aspect<0)=aspect(aspect<0)+360;
            dem.STATVAR.aspect = double(aspect +  dem.TEMP.offset_angle_trueNorth);
        end
        
        function dem = get_horizon_angles_single_point(dem)
            disp('get horizon angles')
            dem.STATVAR.horizon_bins =[];
            dem.STATVAR.horizon_angles =[];

           % rotation_angle_increment = 90./2.^(log(dem.PARA.number_of_horizon_bins)./log(2)-2);
            rotation_angle_increment = 90./(dem.PARA.number_of_horizon_bins./4);
             
            for rotation_angle = 0:rotation_angle_increment:90-1e-9 %rotation_angle = 0:rotation_angle_increment:90-1e-9
                X=imrotate(dem.TEMP.X, -rotation_angle+dem.TEMP.offset_angle_trueNorth);
                Y=imrotate(dem.TEMP.Y, -rotation_angle+dem.TEMP.offset_angle_trueNorth);
                Z=imrotate(dem.TEMP.Z, -rotation_angle+dem.TEMP.offset_angle_trueNorth);

                
                distance2target = (X - dem.TEMP.X_target).^2 + (Y - dem.TEMP.Y_target).^2;
                [ii,jj]=find(distance2target==min(distance2target(:)));
                hor_angles=[];
                for i=1:size(ii,1)
                    grad_X =  (Z(ii(i,1),:) - Z(ii(i,1),jj(i,1))) ./ sqrt((X(ii(i,1),:) - X(ii(i,1),jj(i,1))).^2 + (Y(ii(i,1),:) - Y(ii(i,1),jj(i,1))).^2);
                    grad_Y = (Z(:,jj(i,1)) - Z(ii(i,1),jj(i,1))) ./ sqrt((Y(:,jj(i,1)) - Y(ii(i,1),jj(i,1))).^2 + (X(:,jj(i,1)) - X(ii(i,1),jj(i,1))).^2);
                    hor_angles = [hor_angles [atand(max(grad_Y(ii+1:end))); atand(max(grad_X(1:jj-1))); atand(max(grad_Y(1:ii-1))); atand(max(grad_X(jj+1:end)))]];
                    
                end
                
                dem.STATVAR.horizon_bins = [dem.STATVAR.horizon_bins ; rotation_angle + [0; 270; 180; 90]]; 
%                 dem.STATVAR.horizon_bins = [dem.STATVAR.horizon_bins ; rotation_angle + [0; 90; 180; 270]]; 
                dem.STATVAR.horizon_angles = double([dem.STATVAR.horizon_angles; mean(hor_angles,2)]);
               
            end

            [dem.STATVAR.horizon_bins, order] = sort(dem.STATVAR.horizon_bins);
            dem.STATVAR.horizon_angles = dem.STATVAR.horizon_angles(order,1);
            dem.STATVAR.horizon_angles(dem.STATVAR.horizon_angles<0) = 0;
        
            %skyview factor
            azmRadian = (pi/180).*dem.STATVAR.horizon_bins;
            
            % convert output from horizon program to radians and translate to angle
            % from zenith
            H = (pi/180).*(90-dem.STATVAR.horizon_angles(:));

            aspectRadian = (pi/180)*(dem.STATVAR.aspect);
            % modify limits of integration for slopes facing away from horizons
            t = cosd(dem.STATVAR.aspect-dem.STATVAR.horizon_bins)<0;
            %Simplified trig, the original was H(t) = min(H(t),...
            %  acos(-cos(azmRadian(t)-aspectRadian)*sind(slopeDegrees)./...
            %  sqrt(cosd(slopeDegrees)^2+sind(slopeDegrees)^2*cos(azmRadian(t)-aspectRadian).^2)));
            % but same as
            H(t) = min(H(t), acos(sqrt(1-1./(1+tand(dem.STATVAR.slope_angle)^2*cos(azmRadian(t)-aspectRadian).^2))));
            qIntegrand = (cosd(dem.STATVAR.slope_angle)*sin(H).^2 + sind(dem.STATVAR.slope_angle)*cos(aspectRadian-azmRadian).*(H-cos(H).*sin(H)))/2;
            
            % shouldn't be any negative, except perhaps rounding error, so just in case
            qIntegrand(qIntegrand<0) = 0;
            
            % integrate
            dem.STATVAR.skyview_factor = trapz(azmRadian,qIntegrand)./pi;
            
            %add 360 degree angle
            dem.STATVAR.horizon_bins = [dem.STATVAR.horizon_bins; 360];
            dem.STATVAR.horizon_angles = [dem.STATVAR.horizon_angles; dem.STATVAR.horizon_angles(1,:)];
            
        end
        
        function dem = get_horizon_angles_multiple_point(dem)
            %to be done, as in Kris original script
        end
        
        function dem = get_potential_solar_radiation(dem)
           %to be done, use slope, aspect, etc.
        end
        
    end
end

