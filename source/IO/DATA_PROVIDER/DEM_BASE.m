
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
            if ~isempty(info.MissingDataIndicator)
                dem.TEMP.Z(dem.TEMP.Z==info.MissingDataIndicator) = NaN;
            end
            if strcmp(dem.TEMP.info.RasterReference.CoordinateSystemType, 'planar')
                
                xlim = info.RasterReference.XWorldLimits;
                ylim = info.RasterReference.YWorldLimits;
                dx = info.RasterReference.CellExtentInWorldX;
                dy = info.RasterReference.CellExtentInWorldY;
                x=((xlim(1)+dx/2):dx:(xlim(end)-dx/2))';
                y=flipud(((ylim(1)+dy/2):dy:(ylim(end)-dy/2))');
                [dem.TEMP.X,dem.TEMP.Y] = meshgrid(x,y);
                dem.TEMP.cell_size_X = dem.TEMP.info.RasterReference.CellExtentInWorldX;
                dem.TEMP.cell_size_Y = dem.TEMP.info.RasterReference.CellExtentInWorldY;
            elseif strcmp(dem.TEMP.info.RasterReference.CoordinateSystemType, 'geographic')
                
                R=dem.TEMP.info.RasterReference;
                
                lat=R.LatitudeLimits(2):-R.CellExtentInLatitude:R.LatitudeLimits(1);
                lon=R.LongitudeLimits(1):R.CellExtentInLongitude:R.LongitudeLimits(2);
                
                lat2=(lat(1:end-1)+lat(2:end))/2;
                lon2=(lon(1:end-1)+lon(2:end))/2;
                
                [lon2, lat2] = meshgrid(lon2, lat2');
                if dem.PARA.reproject2utm
                    
                    [X,Y, utm_zone] = ll2utm(dem, lat2, lon2, []);
                    disp(utm_zone)
                    dem.TEMP.utm_zone = utm_zone;
                    cell_size = 30;
                    dem.TEMP.cell_size_X = cell_size;
                    dem.TEMP.cell_size_Y = cell_size;
                    dem.TEMP.X =[round(min(X(:))):cell_size:round(max(X(:)))]';
                    dem.TEMP.Y =[round(max(Y(:))):-cell_size:round(min(Y(:)))]';
                    dem.TEMP.Z = griddata(X,Y,double(dem.TEMP.Z), dem.TEMP.X', dem.TEMP.Y);
                    [dem.TEMP.X,dem.TEMP.Y] = meshgrid(dem.TEMP.X,dem.TEMP.Y);
                else
                    dem.TEMP.X  = lon2;
                    dem.TEMP.Y = lat2;
                    
                end
            end
            
        end
        
        function dem = project_target_coordinates(dem)
            
            if ~strcmp(dem.TEMP.info.RasterReference.CoordinateSystemType, 'geographic') 
                [dem.TEMP.X_target, dem.TEMP.Y_target] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARA.latitude, dem.PARA.longitude);
            else
                if dem.PARA.reproject2utm
                    [dem.TEMP.X_target,dem.TEMP.Y_target] = ll2utm(dem, dem.PARA.latitude, dem.PARA.longitude, dem.TEMP.utm_zone);
                else
                    dem.TEMP.X_target = dem.PARA.longitude ;
                    dem.TEMP.Y_target = dem.PARA.latitude;
                    %[LAT,LON]=projinv(info.CoordinateReferenceSystem,X,Y);
                end
            end
        end
        
        function dem = compute_global_offset_from_north(dem)
            if ~strcmp(dem.TEMP.info.RasterReference.CoordinateSystemType, 'geographic') 
                [X1, Y1] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARA.latitude-0.1, dem.PARA.longitude);
                [X2, Y2] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARA.latitude+0.1, dem.PARA.longitude);
                dem.TEMP.offset_angle_trueNorth = atand((X1-X2)./(Y1-Y2)); %offset grid North and true North
            else
                if dem.PARA.reproject2utm
                    [X1,Y1] = ll2utm(dem, dem.PARA.latitude-0.1, dem.PARA.longitude, dem.TEMP.utm_zone);
                    [X2,Y2] = ll2utm(dem, dem.PARA.latitude+0.1, dem.PARA.longitude, dem.TEMP.utm_zone);
                    dem.TEMP.offset_angle_trueNorth = atand((X1-X2)./(Y1-Y2)); %offset grid North and true North
                else
                    dem.TEMP.offset_angle_trueNorth = 0;
                end
            end
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
            dzdy(2:end-1,2:end-1)=-1.*(dem.TEMP.Z(3:end,2:end-1)-dem.TEMP.Z(1:end-2,2:end-1))./(2.*dem.TEMP.cell_size_Y); % N.B. y decreases with i.
            dzdx(2:end-1,2:end-1)=(dem.TEMP.Z(2:end-1,3:end)-dem.TEMP.Z(2:end-1,1:end-2))./(2.*dem.TEMP.cell_size_X);
            
            slope_angle = atand(sqrt(dzdx.^2+dzdy.^2));
            dem.STATVAR.slope_angle = double(interp2(dem.TEMP.X, dem.TEMP.Y, slope_angle, dem.TEMP.X_target, dem.TEMP.Y_target));
            
        end
        
        function dem = get_aspect_base(dem)
            
            disp('get aspect')
            %Elevation derivatives. N.B. not computed for the boundary of the DEM.
            dzdy=NaN.*dem.TEMP.Z;
            dzdx=dzdy;
            dzdy(2:end-1,2:end-1)=-1.*(dem.TEMP.Z(3:end,2:end-1)-dem.TEMP.Z(1:end-2,2:end-1))./(2.*dem.TEMP.cell_size_Y); % N.B. y decreases with i.
            dzdx(2:end-1,2:end-1)=(dem.TEMP.Z(2:end-1,3:end)-dem.TEMP.Z(2:end-1,1:end-2))./(2.*dem.TEMP.cell_size_X);
            
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
        
        function [x,y,f]=ll2utm(dem, lat, lon, zone)
            
            %LL2UTM Lat/Lon to UTM coordinates precise conversion.
            %	[X,Y]=LL2UTM2(LAT,LON) or LL2UTM([LAT,LON]) converts coordinates
            %	LAT,LON (in degrees) to UTM X and Y (in meters). Default datum is WGS84.
            %
            %	LAT and LON can be scalars, vectors or matrix. Outputs X and Y will
            %	have the same size as inputs.
            %
            %	LL2UTM(...,DATUM) uses specific DATUM for conversion. DATUM can be one
            %	of the following char strings:
            %		'wgs84': World Geodetic System 1984 (default)
            %		'nad27': North American Datum 1927
            %		'clk66': Clarke 1866
            %		'nad83': North American Datum 1983
            %		'grs80': Geodetic Reference System 1980
            %		'int24': International 1924 / Hayford 1909
            %	or DATUM can be a 2-element vector [A,F] where A is semimajor axis (in
            %	meters)	and F is flattening of the user-defined ellipsoid.
            %
            %	LL2UTM(...,ZONE) forces the UTM ZONE (scalar integer or same size as
            %   LAT and LON) instead of automatic set.
            %
            %	[X,Y,ZONE]=LL2UTM(...) returns also the computed UTM ZONE (negative
            %	value for southern hemisphere points).
            %
            %
            %	XY=LL2UTM(...) or without any output argument returns a 2-column
            %	matrix [X,Y].
            %
            %	Note:
            %		- LL2UTM does not perform cross-datum conversion.
            %		- precision is near a millimeter.
            %
            %
            %	Reference:
            %		I.G.N., Projection cartographique Mercator Transverse: Algorithmes,
            %		   Notes Techniques NT/G 76, janvier 1995.
            %
            %	Acknowledgments: Mathieu, Frederic Christen.
            %
            %
            %	Author: Francois Beauducel, <beauducel@ipgp.fr>
            %	Created: 2003-12-02
            %	Updated: 2019-05-29
            %	Copyright (c) 2001-2019, François Beauducel, covered by BSD License.
            %	All rights reserved.
            %
            %	Redistribution and use in source and binary forms, with or without
            %	modification, are permitted provided that the following conditions are
            %	met:
            %
            %	   * Redistributions of source code must retain the above copyright
            %	     notice, this list of conditions and the following disclaimer.
            %	   * Redistributions in binary form must reproduce the above copyright
            %	     notice, this list of conditions and the following disclaimer in
            %	     the documentation and/or other materials provided with the distribution
            %
            %	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
            %	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
            %	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
            %	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
            %	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
            %	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
            %	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
            %	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
            %	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
            %	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
            %	POSSIBILITY OF SUCH DAMAGE.
            % Available datums
            datums = [ ...
                { 'wgs84', 6378137.0, 298.257223563 };
                { 'nad83', 6378137.0, 298.257222101 };
                { 'grs80', 6378137.0, 298.257222101 };
                { 'nad27', 6378206.4, 294.978698214 };
                { 'int24', 6378388.0, 297.000000000 };
                { 'clk66', 6378206.4, 294.978698214 };
                ];
            % constants
            D0 = 180/pi;	% conversion rad to deg
            K0 = 0.9996;	% UTM scale factor
            X0 = 500000;	% UTM false East (m)
            % defaults
            datum = 'wgs84'; 
            if ischar(datum)
                % LL2UTM(...,DATUM) with DATUM as char
                if ~any(strcmpi(datum,datums(:,1)))
                    error('Unkown DATUM name "%s"',datum);
                end
                k = find(strcmpi(datum,datums(:,1)));
                A1 = datums{k,2};
                F1 = datums{k,3};
            else
                % LL2UTM(...,DATUM) with DATUM as [A,F] user-defined
                A1 = datum(1);
                F1 = datum(2);
            end
            p1 = lat/D0;			% Phi = Latitude (rad)
            l1 = lon/D0;			% Lambda = Longitude (rad)
            % UTM zone automatic setting
            if isempty(zone)
                F0 = round((l1*D0 + 183)/6);
            else
                F0 = abs(zone);
            end
            B1 = A1*(1 - 1/F1);
            E1 = sqrt((A1*A1 - B1*B1)/(A1*A1));
            P0 = 0/D0;
            L0 = (6*F0 - 183)/D0;	% UTM origin longitude (rad)
            Y0 = 1e7*(p1 < 0);		% UTM false northern (m)
            N = K0*A1;
            C = coef(E1,0);
            B = C(1)*P0 + C(2)*sin(2*P0) + C(3)*sin(4*P0) + C(4)*sin(6*P0) + C(5)*sin(8*P0);
            YS = Y0 - N*B;
            C = coef(E1,2);
            L = log(tan(pi/4 + p1/2).*(((1 - E1*sin(p1))./(1 + E1*sin(p1))).^(E1/2)));
            z = complex(atan(sinh(L)./cos(l1 - L0)),log(tan(pi/4 + asin(sin(l1 - L0)./cosh(L))/2)));
            Z = N.*C(1).*z + N.*(C(2)*sin(2*z) + C(3)*sin(4*z) + C(4)*sin(6*z) + C(5)*sin(8*z));
            xs = imag(Z) + X0;
            ys = real(Z) + YS;
            % outputs zone if needed: scalar value if unique, or vector/matrix of the
            % same size as x/y in case of crossed zones
            if nargout > 2
                f = F0.*sign(lat);
                fu = unique(f);
                if isscalar(fu)
                    f = fu;
                end
            end
            if nargout < 2
                x = [xs(:),ys(:)];
            else
                x = xs;
                y = ys;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function c = coef(e,m)
                %COEF Projection coefficients
                %	COEF(E,M) returns a vector of 5 coefficients from:
                %		E = first ellipsoid excentricity
                %		M = 0 for transverse mercator
                %		M = 1 for transverse mercator reverse coefficients
                %		M = 2 for merdian arc
                if nargin < 2
                    m = 0;
                end
                switch m
                    case 0
                        c0 = [-175/16384, 0,   -5/256, 0,  -3/64, 0, -1/4, 0, 1;
                            -105/4096, 0, -45/1024, 0,  -3/32, 0, -3/8, 0, 0;
                            525/16384, 0,  45/1024, 0, 15/256, 0,    0, 0, 0;
                            -175/12288, 0, -35/3072, 0,      0, 0,    0, 0, 0;
                            315/131072, 0,        0, 0,      0, 0,    0, 0, 0];
                        
                    case 1
                        c0 = [-175/16384, 0,   -5/256, 0,  -3/64, 0, -1/4, 0, 1;
                            1/61440, 0,   7/2048, 0,   1/48, 0,  1/8, 0, 0;
                            559/368640, 0,   3/1280, 0,  1/768, 0,    0, 0, 0;
                            283/430080, 0, 17/30720, 0,      0, 0,    0, 0, 0;
                            4397/41287680, 0,        0, 0,      0, 0,    0, 0, 0];
                    case 2
                        c0 = [-175/16384, 0,   -5/256, 0,  -3/64, 0, -1/4, 0, 1;
                            -901/184320, 0,  -9/1024, 0,  -1/96, 0,  1/8, 0, 0;
                            -311/737280, 0,  17/5120, 0, 13/768, 0,    0, 0, 0;
                            899/430080, 0, 61/15360, 0,      0, 0,    0, 0, 0;
                            49561/41287680, 0,        0, 0,      0, 0,    0, 0, 0];
                        
                end
                c = zeros(size(c0,1),1);
                for i = 1:size(c0,1)
                    c(i) = polyval(c0(i,:),e);
                end
            end
        end
        
    end
end

