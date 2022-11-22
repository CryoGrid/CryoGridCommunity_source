
classdef POINT_DEM < matlab.mixin.Copyable

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
            point.PARA.area = [];
            
            point.PARA.variables = []; % altitude or altitude, slope_angle, aspect or altitude, slope_angle, aspect, horizon_angles
                        
            point.PARA.DEM_folder = [];
            point.PARA.DEM_filename = [];
            
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
            
            point.TEMP.Z = readgeoraster([point.PARA.DEM_folder point.PARA.DEM_filename]);
            info = georasterinfo([point.PARA.DEM_folder point.PARA.DEM_filename]);
            point.TEMP.info = info;
            point.TEMP.Z(point.TEMP.Z==info.MissingDataIndicator) = NaN;
            xlim = info.RasterReference.XWorldLimits; 
            ylim = info.RasterReference.YWorldLimits;
            dx = info.RasterReference.CellExtentInWorldX; 
            dy = info.RasterReference.CellExtentInWorldY;
            x=((xlim(1)+dx/2):dx:(xlim(end)-dx/2))';
            y=flipud(((ylim(1)+dy/2):dy:(ylim(end)-dy/2))');
            [point.TEMP.X,point.TEMP.Y] = meshgrid(x,y);  
            [point.TEMP.X_target, point.TEMP.Y_target] = projfwd(info.CoordinateReferenceSystem, point.PARA.latitude, point.PARA.longitude);
            %[LAT,LON]=projinv(info.CoordinateReferenceSystem,X,Y);
            
            for i=1:size(point.PARA.variables,1)
                a = str2func(['get_' point.PARA.variables{i,1}]);
                point = a(point);
            end
            
            
        end
        
        function point = get_altitude(point)
            point.STATVAR.altitude = interp2(point.TEMP.X, point.TEMP.Y, point.TEMP.Z, point.TEMP.X_target, point.TEMP.Y_target);
            
        end
       
        function point = get_slope_angle(point)            
            %Elevation derivatives. N.B. not computed for the boundary of the DEM.
            dzdy=NaN.*point.TEMP.Z;  
            dzdx=dzdy;
            dzdy(2:end-1,2:end-1)=-1.*(point.TEMP.Z(3:end,2:end-1)-point.TEMP.Z(1:end-2,2:end-1))./(2.*point.TEMP.info.RasterReference.CellExtentInWorldY); % N.B. y decreases with i.
            dzdx(2:end-1,2:end-1)=(point.TEMP.Z(2:end-1,3:end)-point.TEMP.Z(2:end-1,1:end-2))./(2.*point.TEMP.info.RasterReference.CellExtentInWorldX);
            
            slope_angle = atand(sqrt(dzdx.^2+dzdy.^2));
            point.STATVAR.slope_angle = interp2(point.TEMP.X, point.TEMP.Y, slope_angle, point.TEMP.X_target, point.TEMP.Y_target);
            
        
            
        end
        
        function point = get_aspect(point)
            %Elevation derivatives. N.B. not computed for the boundary of the DEM.
            dzdy=NaN.*point.TEMP.Z;  
            dzdx=dzdy;
            dzdy(2:end-1,2:end-1)=-1.*(point.TEMP.Z(3:end,2:end-1)-point.TEMP.Z(1:end-2,2:end-1))./(2.*point.TEMP.info.RasterReference.CellExtentInWorldY); % N.B. y decreases with i.
            dzdx(2:end-1,2:end-1)=(point.TEMP.Z(2:end-1,3:end)-point.TEMP.Z(2:end-1,1:end-2))./(2.*point.TEMP.info.RasterReference.CellExtentInWorldX);
            
            % Aspect definition: N,E,S,W = pi, pi/2, 0, -pi/2
            aspect = atan2d(-1.*dzdx,dzdy);
            point.STATVAR.aspect = interp2(point.TEMP.X, point.TEMP.Y, aspect, point.TEMP.X_target, point.TEMP.Y_target);
        end
        
        function point = get_horizon_angles(point)
            %talk to Kris
        end
        
        
    end
end

