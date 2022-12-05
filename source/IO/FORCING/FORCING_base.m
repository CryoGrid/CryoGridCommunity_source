%========================================================================
% Base class for FORCING, containing functions that can be called from all
% inheriting, specific FORCING classes. Built with functions from classes
% available in previous FORCING class structure.
% R. B. Zweigel, September 2022
%========================================================================

classdef FORCING_base < matlab.mixin.Copyable
    
    properties
        forcing_index
        STATUS
        PARA
        CONST
        TEMP
        STATVAR
        DATA            % forcing data time series
    end
    
    
    methods
        
        function forcing = provide_PARA(forcing)
            % empty; PARAs are provided in the relevant forcing classes
        end
        
        function forcing = provide_CONST(forcing)
            % empty; CONSTs are provided in the relevant forcing classes
        end
        
        function forcing = provide_STATVAR(forcing)
            % empty; STATVARs are provided in the relevant forcing classes
        end
        
        %%%%%% ----- finalize init functions ----- %%%%%%
        
        
        function forcing = initialize_TEMP(forcing)
            % The TEMP variables required for all classes
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Lin=0;
            forcing.TEMP.Sin=0;
            forcing.TEMP.Tair=0;
            forcing.TEMP.wind=0;
            forcing.TEMP.q=0;
            forcing.TEMP.p=0;
        end
        
        function forcing = initialize_TEMP_slope(forcing)
            % additional TEMPs required for runs on slopes
            forcing.TEMP.S_TOA=0;
            %             forcing.TEMP.albedo_foot=0;
            forcing.TEMP.Sin_dif = 0;
            forcing.TEMP.Sin_dir = 0;
            forcing.TEMP.sunElevation = 0;
            forcing.TEMP.azimuth = 0;
        end
        
%         function forcing = initialize_terrain(forcing)
%             % load parameters required for terrain shielding
%             % All angles in degrees and clockwise from N (maybe change)
% 
%             % Added by THIN, 2022-10-26
%             error('Error. Terrain parameters should be requested from TERRAIN class!')
%             % This funciton should not be called.
%             % Terrain parameters should instead be obtained from the 
%             % TERRAIN class which can be accessed at tile.TERRAIN.
% 
%             i = forcing.PARA.tp_number;
%             path = forcing.PARA.forcing_path;
%             tp_file = forcing.PARA.terrain_file;
%             load([path tp_file])
%             
%             forcing.TEMP.azimuth = 0;
%             forcing.PARA.h = rad2deg(tp.h(i,:));
%             forcing.PARA.hbins = 180 - 180/pi*tp.hbins;
%             forcing.PARA.slope_angle = rad2deg(tp.slp(i));
%             forcing.PARA.aspect = 180 - 180/pi*tp.asp(i);
%             forcing.PARA.sky_view_factor = tp.svf(i);
%             
%         end
        
        function forcing = check_and_correct(forcing) % corrects typical forcing data errors
            
            % Check for consistent timesteps
            if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1))~=0
                disp('timestamp of forcing data is not in regular intervals -> check, fix and restart')
                forcing.STATUS=0;
                return
            else
                forcing.STATUS=1;
            end
            
            % Correct known isues
            if isfield(forcing.DATA, 'wind')
                forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence
            end

            if isfield(forcing.DATA, 'Lin') && isfield(forcing.DATA, 'Tair')
                forcing.DATA.Lin(find(forcing.DATA.Lin==0)) = 5.67e-8 .* (forcing.DATA.Tair(find(forcing.DATA.Lin==0))+273.15).^4;
            end
        end
        
        function forcing = set_start_and_end_time(forcing)
            % Assign forcing start and end times for current run
            % -> if not specified in params file, use forcing data length
            
            if isempty(forcing.PARA.start_time) || isnan(forcing.PARA.start_time(1,1))
                forcing.PARA.start_time = forcing.DATA.timeForcing(1,1);
            else
                forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            end
            
            if isempty(forcing.PARA.end_time) || isnan(forcing.PARA.end_time(1,1))
                forcing.PARA.end_time = floor(forcing.DATA.timeForcing(end,1));
            else
                forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));
            end
            
        end
        
        function forcing = split_precip_Tair(forcing)
            % Split total precip into snow and rain depending on air temp.
            Tair = forcing.DATA.Tair;
            T_all_rain = forcing.PARA.all_rain_T;
            T_all_snow = forcing.PARA.all_snow_T;
            
            forcing.DATA.snowfall = forcing.DATA.precip .* (double(Tair <= T_all_snow)  + ...
                double(Tair > T_all_snow & Tair < T_all_rain) .* (1- (Tair - T_all_snow) ./ max(1e-12, (T_all_rain - T_all_snow))));
            forcing.DATA.rainfall = forcing.DATA.precip .* (double(Tair >= T_all_rain)  + ...
                double(Tair > T_all_snow & Tair < T_all_rain) .* (Tair - T_all_snow) ./ max(1e-12, (T_all_rain - T_all_snow)));
            
        end
        
        function forcing = reduce_precip_slope(forcing, tile)
            slope = tile.PARA.slope_angle;
            % scale precip due to local slope and user definer factor
            forcing.DATA.rainfall = forcing.DATA.rainfall .*cosd(slope);
            forcing.DATA.snowfall = forcing.DATA.snowfall .*cosd(slope);
        end
        
        function forcing = split_Sin(forcing) % Consider making different functions for different parameterizations
            % Split Sin into direct and diffuse parts, see Fiddes(2014)
            Sin = forcing.DATA.Sin;
            Sin_TOA = forcing.DATA.S_TOA;
            
            kt = Sin ./ Sin_TOA; %clearness index
            kt(isnan(kt))=0;
            kd = 0.952 - 1.041.*exp(-exp(2.300 - 4.702 .* kt)); %diffuse fraction
            kd(isnan(kd))=0;
            kd = max(0,kd);
            
            forcing.DATA.Sin_dif = kd .* forcing.DATA.Sin; % full hemisphere diffuse Sin
            forcing.DATA.Sin_dir = (1-kd) .* forcing.DATA.Sin; % direct Sin on horizontal surface
            
        end
        
        function forcing = terrain_corr_Sin_dif(forcing, tile)
            % include effect of terrain on diffuse by removin the fraction
            % of the hemisphere covered by the horizon, and adding
            % reflected Sin from surrounding terrain
            
            Sin = forcing.DATA.Sin; % Total Sin (horizontal)
            Sin_dif = forcing.DATA.Sin_dif;% Diffuse Sin (horizontal)
            alpha = forcing.PARA.albedo_surrounding_terrain; %Albedo at the foot of the slope
            svf = tile.PARA.skyview_factor; % hemispheric fraction of sky not occluded by terrain
            
            forcing.DATA.Sin_dif = Sin_dif.*svf + Sin.*alpha.*(1-svf);
        end
        
        function forcing = reproject_Sin_dir(forcing, tile)
            Sin_dir = forcing.DATA.Sin_dir;
            aspect = tile.PARA.aspect;
            slope = tile.PARA.slope_angle;
            
            surf_norm_vec = repmat([0.0,0.0,1.0], size(forcing.DATA.Sin,1), 1); %Unit vector on the horizontal
            
            alpha = aspect.*pi./180; %Degree to radians of exposition of the slope
            beta = (90-slope).*pi./180; %Degree to radians of inclination of the slope
            face_vec = repmat([sin(alpha).*cos(beta) cos(alpha).*cos(beta) sin(beta)], size(forcing.DATA.Sin,1), 1); %Unit vector of the slope
            
            % Calculation the solar azimuth and elevation angle relative to
            % the coordinates of the site (revised after Darin C. Koblick)
            % ->CHANGED to Kris script
            forcing = SolarAzEl(forcing, tile);
            
            alpha = forcing.DATA.azimuth.*pi/180; %Degree to radians of the azimuth
            beta = forcing.DATA.sunElevation.*pi/180; %Degree to radians of elevation
            sun_vec = [sin(alpha).*cos(beta) cos(alpha).*cos(beta) sin(beta)]; %Unit vector of the radiation
            
            delta_angle_surf_norm = acos(dot(surf_norm_vec' ,sun_vec')').*180./pi; %angle between the radiation and the normal on the horizontal in degrees
            delta_angle_face=acos(dot(face_vec', sun_vec')').*180./pi; %angle between the radiation and the normal on the slope in degrees
            
            Sin_sun_direction = Sin_dir./cos(delta_angle_surf_norm.*pi./180); %direct Sin in direction of the sun
            Sin_face_direction = double(delta_angle_surf_norm < 85 & delta_angle_face < 85) .* Sin_sun_direction.*cos(delta_angle_face.*pi/180); % Direct Sin in direction of surface face
            
            forcing.DATA.Sin_dir = Sin_face_direction;
        end
        
        function forcing = terrain_corr_Lin(forcing, tile)
            % Reduce Lin from atmosphere to sky view fraction, add Lin from
            % surroundings assuming Tair everywhere
            sigma = forcing.CONST.sigma; %Stefan-Boltzman constant
            Lin = forcing.DATA.Lin;
            Tair = forcing.DATA.Tair;
            svf = tile.PARA.skyview_factor;
            Tmfw = forcing.CONST.Tmfw;
            
            forcing.DATA.Lin = svf.*Lin + (1-svf).*sigma.*(Tair+Tmfw).^4;
        end
        
        function forcing = interpolate_forcing(forcing, tile)
            t = tile.t;
            
            posit=floor((t-forcing.DATA.timeForcing(1,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)))+1;
            
            variables = fieldnames(forcing.TEMP);
            t_weight = (t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)); % current distance from last timestep (0,1)
            
            for i = 1:length(variables)
                if ~strcmp(variables{i},'t') && isfield(forcing.DATA, variables{i})
                    forcing.TEMP.(variables{i}) = forcing.DATA.(variables{i})(posit,1)+(forcing.DATA.(variables{i})(posit+1,1)-forcing.DATA.(variables{i})(posit,1)).*t_weight;
                end
            end
            forcing.TEMP.t = t;
        end

        
        function forcing = terrain_shade(forcing, tile)
            az = forcing.DATA.azimuth;
            el = forcing.DATA.sunElevation;
            hbins = tile.PARA.horizon_bins;
            h = tile.PARA.horizon_angles;
            
            I = knnsearch(hbins,az); % hbin containing current solar azimuth -> INTERPOLATE INSTEAD??? in the DEM analysis, this is points along lines, not bins!!
            forcing.DATA.Sin_dir(h(I)>el) = 0; % remove direct Sin if sun is below horizon
            %forcing.DATA.Sin(h(I)>el) = forcing.DATA.Sin_dif; %needs to be
            %done later!!
        end
        
        %%%%% ----- Support functions ----- %%%%%
        
        function p = satPresIce(forcing, T)
            Tmfw = forcing.CONST.Tmfw;
            p = 6.112.* 100.* exp(22.46.*(T-Tmfw)./(272.61+T-Tmfw));
        end
        
        function p = satPresWater(forcing, T)
            Tmfw = forcing.CONST.Tmfw;
            p = 6.112 .* 100 .* exp(17.62.*(T-Tmfw)./(243.12+T-Tmfw));
        end
        
        %Juditah's script, not used anymore - this is just a random choice, not
        %sure whichone is better
        
%         function forcing = SolarAzEl(forcing, tile)
%             % SolarAzEl will ingest a Universal Time, and specific site location on earth
%             % it will then output the solar Azimuth and Elevation angles relative to
%             % that site.
%             %
%             % Programed by Darin C. Koblick 2/17/2009
%             %              Darin C. Koblick 4/16/2013 Vectorized for Speed
%             %                             Allow for MATLAB Datevec input in
%             %                             addition to a UTC string.
%             %                             Cleaned up comments and code to
%             %                             avoid warnings in MATLAB editor.
%             %              Solar Position obtained from: http://stjarnhimlen.se/comp/tutorial.html#5
%             Lat = tile.PARA.latitude;
%             Lon = tile.PARA.longitude;
%             Alt = tile.PARA.altitude ./ 1000;
%             t_span = forcing.DATA.timeForcing;
%             UTC = datestr(t_span);
%             %i_end = size(UTC);
%             
%             %Loop through all timesteps
%             %for i = 1 : i_end(1,1)
%             
%             %compute JD from UTC
%             [year, month, day, hour, min, sec] = datevec(UTC); %datevec(UTC(i,:));
%             idx = (month <= 2);
%             year(idx) = year(idx)-1;
%             month(idx) = month(idx)+12;
%             jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
%                 floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + ...
%                 (hour + min/60 + sec/3600)/24;
%             d = jd-2451543.5;
%             
%             % Keplerian Elements for the Sun (geocentric)
%             w = 282.9404+4.70935e-5*d;     %(longitude of perihelion degrees)
%             e = 0.016709-1.151e-9.*d;       %(eccentricity)
%             M = mod(356.0470+0.9856002585.*d,360);   %(mean anomaly degrees)
%             L = w + M;                     %(Sun's mean longitude degrees)
%             oblecl = 23.4393-3.563e-7.*d;  %(Sun's obliquity of the ecliptic)
%             
%             %auxiliary angle
%             E = M+(180/pi).*e.*sin(M.*(pi/180)).*(1+e.*cos(M.*(pi/180)));
%             
%             %rectangular coordinates in the plane of the ecliptic (x axis toward perhilion)
%             x = cos(E.*(pi/180))-e;
%             y = sin(E.*(pi/180)).*sqrt(1-e.^2);
%             
%             %find the distance and true anomaly
%             r = sqrt(x.^2 + y.^2);
%             v = atan2(y,x).*(180/pi);
%             
%             %find the longitude of the sun
%             lon = v + w;
%             
%             %compute the ecliptic rectangular coordinates
%             xeclip = r.*cos(lon.*(pi/180));
%             yeclip = r.*sin(lon.*(pi/180));
%             zeclip = 0.0;
%             
%             %rotate these coordinates to equitorial rectangular coordinates
%             xequat = xeclip;
%             yequat = yeclip.*cos(oblecl.*(pi/180))+zeclip.*sin(oblecl.*(pi/180));
%             zequat = yeclip.*sin(23.4406.*(pi/180))+zeclip.*cos(oblecl.*(pi/180));
%             
%             %convert equatorial rectangular coordinates to RA and Decl:
%             r = sqrt(xequat.^2 + yequat.^2 + zequat.^2)-(Alt./149598000); %roll up the altitude correction
%             RA = atan2(yequat,xequat).*(180/pi);
%             delta = asin(zequat./r).*(180/pi);
%             
%             %Following the RA DEC to Az Alt conversion sequence explained here: http://www.stargazing.net/kepler/altaz.html
%             
%             %Find the J2000 value
%             J2000 = jd - 2451545.0;
%             %                 hourvec = datevec(UTC);  %(i,:)
%             %                 UTH = hourvec(:,4) + hourvec(:,5)/60 + hourvec(:,6)/3600;
%             UTH = hour + min./60 + sec./3600;
%             
%             %Calculate local siderial time
%             GMST0=mod(L+180,360)./15;
%             SIDTIME = GMST0 + UTH + Lon./15;
%             
%             %Replace RA with hour angle HA
%             HA = (SIDTIME.*15 - RA);
%             
%             %convert to rectangular coordinate system
%             x = cos(HA.*(pi/180)).*cos(delta.*(pi/180));
%             y = sin(HA.*(pi/180)).*cos(delta.*(pi/180));
%             z = sin(delta.*(pi/180));
%             
%             %rotate this along an axis going east-west.
%             xhor = x.*cos((90-Lat).*(pi/180))-z.*sin((90-Lat).*(pi/180));
%             yhor = y;
%             zhor = x.*sin((90-Lat).*(pi/180))+z.*cos((90-Lat).*(pi/180));
%             
%             %Find the h and AZ
%             forcing.DATA.azimuth = atan2(yhor,xhor).*(180./pi) + 180;
%             forcing.DATA.sunElevation = asin(zhor).*(180./pi); %in degrees
%             
%         end
        
        %wrapper for Kris function
        function forcing = SolarAzEl(forcing, tile)
            [forcing.DATA.azimuth,forcing.DATA.sunElevation] = solargeom(forcing, forcing.DATA.timeForcing ,tile.PARA.latitude,tile.PARA.longitude);
            forcing.DATA.azimuth = rad2deg(forcing.DATA.azimuth);
            forcing.DATA.azimuth(forcing.DATA.azimuth<0) = forcing.DATA.azimuth(forcing.DATA.azimuth<0) + 360;
            forcing.DATA.sunElevation = 90-rad2deg(forcing.DATA.sunElevation);
        end
        
        %Kris script, used now to be consistent with TopoScale
        function [solar_azimuth,solar_zenith]=solargeom(forcing, Time, Latitude,Longitude)
            %% [saz,szen]=solargeom(time,latitude,longitude)
            % Adopted from the Sandia National Labs PVL Toolbox ephemeris routine.
            % Inputs:
            %   time = Time stamp vector (matlab datenum format) assumed to be in UTC
            %   latitude = Latitude
            %   longitude = Longitude
            % Outputs:
            %   saz = Solar azimuth angle [radians, anticlockwise from south]
            %   szen = Solar zentih angle [radians].
            % Link to the original toolbox:
            % https://pvpmc.sandia.gov/applications/pv_lib-toolbox/
            % References:
            % Stein et al. (2012), doi:10.1109/PVSC.2012.6318225 [MATLAB version]
            % Holmgren et al. (2018), doi:10.21105/joss.00884 [Python version]
            
            
            TZone=0;
            Longitude=-Longitude;
            % tv=datevec(Time);
            Year=year(Time);
            % v0=zeros(size(Year)); v1=ones(size(Year));
            DayOfYear=floor(Time-datenum(Year,1, 1))+1;
            DecHours=(Time - floor(Time)) .* 24;
            RadtoDeg=180/pi;
            DegtoRad=pi/180;
            Abber = 20/3600;
            LatR = Latitude * DegtoRad;
            UnivDate = DayOfYear + floor((DecHours + TZone)/24);
            UnivHr = mod((DecHours + TZone), 24);
            Yr = Year-1900;
            YrBegin = 365 * Yr + floor((Yr-1)/4)-0.5;
            Ezero = YrBegin + UnivDate;
            T = Ezero / 36525;
            GMST0 = 6/24 +38/1440 + (45.836 + 8640184.542 * T + 0.0929 * T.^2)/86400;
            GMST0 = 360 * (GMST0 - floor(GMST0));
            GMSTi = mod(GMST0 + 360*(1.0027379093 * UnivHr / 24),360);
            LocAST = mod((360 + GMSTi - Longitude), 360);
            EpochDate = Ezero + UnivHr / 24;
            T1 = EpochDate / 36525;
            ObliquityR = DegtoRad * (23.452294 - 0.0130125 * T1 - 0.00000164 * T1.^2 ...
                + 0.000000503 * T1.^3);
            MlPerigee = 281.22083 + 0.0000470684 * EpochDate + 0.000453 * T1 .^ 2 + ...
                0.000003 * T1 .^ 3;
            MeanAnom = mod((358.47583 + 0.985600267 * EpochDate - 0.00015 * T1 .^ 2 - ...
                0.000003 * T1 .^ 3), 360);
            Eccen = 0.01675104 - 0.0000418 * T1 - 0.000000126 * T1 .^ 2;
            EccenAnom = MeanAnom;
            E=0;
            while max(abs(EccenAnom - E)) > 0.0001
                E = EccenAnom;
                EccenAnom = MeanAnom + RadtoDeg .* Eccen .* sin(DegtoRad .* E);
            end
            TrueAnom = 2 * mod(RadtoDeg * atan2(((1 + Eccen) ./ (1 - Eccen)).^ 0.5 .* tan(DegtoRad * EccenAnom / 2), 1), 360) ;
            EcLon = mod(MlPerigee + TrueAnom, 360) - Abber ;
            EcLonR = DegtoRad * EcLon;
            DecR = asin(sin(ObliquityR) .* sin(EcLonR));
            %Dec = RadtoDeg * DecR;
            RtAscen = RadtoDeg * atan2(cos(ObliquityR).*(sin(EcLonR)),cos(EcLonR));
            HrAngle = LocAST - RtAscen ;
            HrAngleR = DegtoRad .* HrAngle ;
            %HrAngle = HrAngle - (360 .* sign(HrAngle) .* (abs(HrAngle) > 180));
            SunAz = RadtoDeg .* atan2(-1 * sin(HrAngleR), cos(LatR) .* tan(DecR) - sin(LatR) .* cos(HrAngleR));
            SunAz = SunAz + (SunAz < 0) * 360; %shift from range of [-180,180] to [0,360]
            SunEl = asind(cos(LatR) .* cos(DecR) .* cos(HrAngleR) + sin(LatR) .* sin(DecR));
            
            % Convert solar azimuth angle from [N,E,S,W]=[0,90,180,270] to [180, 90, 0
            % -90], i.e. the same as the aspect and horizon angle system.
            solar_azimuth=deg2rad(SunAz);
            solar_azimuth=(5*pi/2)-solar_azimuth;
            solar_azimuth=solar_azimuth-2.*pi.*(solar_azimuth>2.*pi);
            solar_azimuth=solar_azimuth+pi/2;
            solar_azimuth=solar_azimuth-2.*pi.*(solar_azimuth>pi);
            
            % Calculate solar zenith angle from solar elevation angle
            SunEl=deg2rad(SunEl);
            solar_zenith=(pi/2)-SunEl;
        end
        
    end

end