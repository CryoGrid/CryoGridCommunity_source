
classdef DEM < matlab.mixin.Copyable

    properties
        PARENT
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods
        function dem = provide_PARA(dem)            
            dem.PARA.variables = []; % altitude or altitude, slope_angle, aspect or altitude, slope_angle, aspect, horizon_angles
                        
            dem.PARA.DEM_folder = [];
            dem.PARA.DEM_filename = [];
            
        end
        
        function dem = provide_STATVAR(dem)

        end
        
        function dem = provide_CONST(dem)
            
        end
        
        function dem = finalize_init(dem)
            
            dem.STATVAR.Z = readgeoraster([dem.PARA.DEM_folder dem.PARA.DEM_filename]);
            info = georasterinfo([dem.PARA.DEM_folder dem.PARA.DEM_filename]);
            dem.STATVAR.info = info;
            dem.STATVAR.Z(dem.STATVAR.Z==single(info.MissingDataIndicator)) = NaN;
            xlim = info.RasterReference.XWorldLimits; 
            ylim = info.RasterReference.YWorldLimits;
            dx = info.RasterReference.CellExtentInWorldX; 
            dy = info.RasterReference.CellExtentInWorldY;
            x=((xlim(1)+dx/2):dx:(xlim(end)-dx/2))';
            y=flipud(((ylim(1)+dy/2):dy:(ylim(end)-dy/2))');
            [dem.STATVAR.X,dem.STATVAR.Y] = meshgrid(x,y);  
            
        end
        
        function dem = load_data(dem)
            [dem.STATVAR.X_target, dem.STATVAR.Y_target] = projfwd(dem.STATVAR.info.CoordinateReferenceSystem, dem.PARENT.STATVAR.latitude, dem.PARENT.STATVAR.longitude);
            for i=1:size(dem.PARA.variables,1)
                a = str2func(['get_' dem.PARA.variables{i,1}]);
                dem = a(dem);
            end
        end
        
        function dem = get_altitude(dem)
            dem.PARENT.STATVAR.altitude = interp2(dem.STATVAR.X, dem.STATVAR.Y, dem.STATVAR.Z, dem.STATVAR.X_target, dem.STATVAR.Y_target);
            
        end
       
        function dem = get_slope_angle(dem)            
            %Elevation derivatives. N.B. not computed for the boundary of the DEM.
            dzdy=NaN.*dem.STATVAR.Z;  
            dzdx=dzdy;
            dzdy(2:end-1,2:end-1)=-1.*(dem.STATVAR.Z(3:end,2:end-1)-dem.STATVAR.Z(1:end-2,2:end-1))./(2.*dem.STATVAR.info.RasterReference.CellExtentInWorldY); % N.B. y decreases with i.
            dzdx(2:end-1,2:end-1)=(dem.STATVAR.Z(2:end-1,3:end)-dem.STATVAR.Z(2:end-1,1:end-2))./(2.*dem.STATVAR.info.RasterReference.CellExtentInWorldX);
            
            slope_angle = atand(sqrt(dzdx.^2+dzdy.^2));
            dem.PARENT.STATVAR.slope_angle = interp2(dem.STATVAR.X, dem.STATVAR.Y, slope_angle, dem.STATVAR.X_target, dem.STATVAR.Y_target);

        end
        
        function dem = get_aspect(dem)
            %Elevation derivatives. N.B. not computed for the boundary of the DEM.
            dzdy=NaN.*dem.STATVAR.Z;  
            dzdx=dzdy;
            dzdy(2:end-1,2:end-1)=-1.*(dem.STATVAR.Z(3:end,2:end-1)-dem.STATVAR.Z(1:end-2,2:end-1))./(2.*dem.STATVAR.info.RasterReference.CellExtentInWorldY); % N.B. y decreases with i.
            dzdx(2:end-1,2:end-1)=(dem.STATVAR.Z(2:end-1,3:end)-dem.STATVAR.Z(2:end-1,1:end-2))./(2.*dem.STATVAR.info.RasterReference.CellExtentInWorldX);
            
            % Aspect definition: N,E,S,W = pi, pi/2, 0, -pi/2
            aspect = atan2d(-1.*dzdx,dzdy);
            dem.PARENT.STATVAR.aspect = interp2(dem.STATVAR.X, dem.STATVAR.Y, aspect, dem.STATVAR.X_target, dem.STATVAR.Y_target);
        end
        
        function dem = get_horizon_angles(dem)
            %talk to Kris
        end
        
        
    end
end

