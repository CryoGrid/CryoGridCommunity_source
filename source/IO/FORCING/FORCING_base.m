%========================================================================
% Base class for FORCING, containing functions that can be called from all
% inheriting, specific FORCING classes. Built with functions from classes
% available in previous FORCING class structure.
% R. B. Zweigel, September 2022
%========================================================================

classdef FORCING_base < matlab.mixin.Copyable
    
    properties
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
        
        function forcing = readNC(forcing, tile) % Consider making seperate functions dependant on how data is stored
            % Read forcing data (ERA5 format) from NetCDF files, and
            % populate forcing.DATA
            
            variables = {'t2m'; 'd2m'; 'u10'; 'v10'; 'ssrd'; 'strd'; 'tp'; 'sp'; 'tisr'};
            for i=1:size(variables,1)
                temp.(variables{i,1}) = squeeze(ncread([forcing.PARA.forcing_path variables{i,1} '.nc'], variables{i,1}));
            end
            temp.time = ncread([forcing.PARA.forcing_path 't2m.nc'], 'time');
            
            temp.info = ncinfo([forcing.PARA.forcing_path 't2m.nc']);
            Index1 = find(contains({temp.info.Variables.Name},'time'));
            Index2 = find(contains({temp.info.Variables(Index1).Attributes.Name},'units'));
            reference = split(temp.info.Variables(Index1).Attributes(Index2).Value);
            reference = reference(3);
            reference_date = split(reference,"-");
            reference_date = str2double(reference_date);
            forcing.DATA.timeForcing = datenum(reference_date(1),reference_date(2),reference_date(3)) + double(temp.time)./24;
            
            forcing.DATA.Tair = temp.t2m - forcing.CONST.Tmfw;
            forcing.DATA.wind = sqrt(temp.u10.^2 + temp.u10.^2);
            forcing.DATA.Sin = temp.ssrd./3600;
            forcing.DATA.Sin = [0;0;forcing.DATA.Sin];
            forcing.DATA.Lin = temp.strd./3600;
            forcing.DATA.Lin = [0;0;forcing.DATA.Lin];
            forcing.DATA.p = temp.sp; %in [Pa]!900.*100 + forcing.DATA.Tair.*0;
            forcing.DATA.q = (double(forcing.DATA.Tair<0).*satPresIce(forcing, temp.d2m) + double(forcing.DATA.Tair>=0).*satPresWater(forcing, temp.d2m))./ forcing.DATA.p;
            forcing.DATA.precip = [0;0; temp.tp]; %append missing first two timesteps
            
            % RBZ Sep22: Might remove Sin TOA and calculate this using ERA5
            % parameterizations -> consistent solar geometry as in source
            forcing.DATA.S_TOA = temp.tisr ./ 3600;
            forcing.DATA.S_TOA = [0;0;forcing.DATA.S_TOA];
            
        end
        
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
        end
        
        function forcing = initialize_terrain(forcing)
            % load parameters required for terrain shielding
            % All angles in degrees and clockwise from N (maybe change)
            i = forcing.PARA.tp_number;
            path = forcing.PARA.forcing_path;
            tp_file = forcing.PARA.terrain_file;
            load([path tp_file])
            
            forcing.TEMP.azimuth = 0;
            forcing.PARA.h = rad2deg(tp.h(i,:));
            forcing.PARA.hbins = 180 - 180/pi*tp.hbins;
            forcing.PARA.slope_angle = rad2deg(tp.slp(i));
            forcing.PARA.aspect = 180 - 180/pi*tp.asp(i);
            forcing.PARA.sky_view_factor = tp.svf(i);
            
        end
        
        function forcing = checkAndCorrect(forcing) % corrects typical forcing data errors
            
            % Check for consistent timesteps
            if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1))~=0
                disp('timestamp of forcing data is not in regular intervals -> check, fix and restart')
                return
            end
            
            % Correct known isues
            forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence
            forcing.DATA.Lin(find(forcing.DATA.Lin==0)) = 5.67e-8 .* (forcing.DATA.Tair(find(forcing.DATA.Lin==0))+273.15).^4;
            
        end
        
        function forcing = get_time_forcing(forcing)
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
            
            forcing.DATA.snowfall = forcing.DATA.precip .*24.*1000 .* (double(Tair <= T_all_snow)  + ...
                double(Tair > T_all_snow & Tair <= T_all_rain) .* (Tair - T_all_snow) ./ max(1e-12, (T_all_rain - T_all_snow)));
            forcing.DATA.rainfall = forcing.DATA.precip .*24.*1000 .* (double(Tair >= T_all_rain)  + ...
                double(Tair > T_all_snow & Tair < T_all_rain) .* (1 - (Tair - T_all_snow) ./ max(1e-12, (T_all_rain - T_all_snow))));
            
        end
        
        function forcing = reduce_precip_slope(forcing)
            % scale precip due to local slope and user definer factor
            forcing.DATA.rainfall = forcing.DATA.rainfall .* forcing.PARA.rain_fraction.*cosd(forcing.PARA.slope_angle);
            forcing.DATA.snowfall = forcing.DATA.snowfall .* forcing.PARA.snow_fraction.*cosd(forcing.PARA.slope_angle);
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
        
        function forcing = terrain_corr_Sin_dif(forcing)
            % include effect of terrain on diffuse by removin the fraction
            % of the hemisphere covered by the horizon, and adding
            % reflected Sin from surrounding terrain
            
            Sin = forcing.DATA.Sin; % Total Sin (horizontal)
            Sin_dif = forcing.DATA.Sin_dif;% Diffuse Sin (horizontal)
            alpha = forcing.PARA.albedo_surrounding_terrain; %Albedo at the foot of the slope
            svf = forcing.PARA.sky_view_factor; % hemispheric fraction of sky not occluded by terrain
            
            forcing.DATA.Sin_dif = Sin_dif.*svf + Sin.*alpha.*(1-svf);
        end
        
        function forcing = reproject_Sin_dir(forcing, tile)
            Sin_dir = forcing.DATA.Sin_dir;
            aspect = forcing.PARA.aspect;
            slope = forcing.PARA.slope_angle;
            
            surf_norm_vec = repmat([0.0,0.0,1.0], size(forcing.DATA.Sin,1), 1); %Unit vector on the horizontal
            
            alpha = aspect.*pi./180; %Degree to radians of exposition of the slope
            beta = (90-slope).*pi./180; %Degree to radians of inclination of the slope
            face_vec = repmat([sin(alpha).*cos(beta) cos(alpha).*cos(beta) sin(beta)], size(forcing.DATA.Sin,1), 1); %Unit vector of the slope
            
            % Calculation the solar azimuth and elevation angle relative to the coordinates of the site (revised after Darin C. Koblick)
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
        
        function forcing = terrain_corr_Lin(forcing)
            % Reduce Lin from atmosphere to sky view fraction, add Lin from
            % surroundings assuming Tair everywhere
            sigma = forcing.CONST.sigma; %Stefan-Boltzman constant
            Lin = forcing.DATA.Lin;
            Tair = forcing.DATA.Tair;
            svf = forcing.PARA.sky_view_factor;
            Tmfw = forcing.CONST.Tmfw;
            
            forcing.DATA.Lin = svf.*Lin + (1-svf).*sigma.*(Tair+Tmfw).^4;
        end
        
        function forcing = interpolate_forcing(forcing, tile)
            t = tile.t;
            
            posit=floor((t-forcing.DATA.timeForcing(1,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)))+1;
            
            variables = fieldnames(forcing.TEMP);
            t_weight = (t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)); % current distance from last timestep (0,1)
            
            for i = 1:length(variables)
                if ~strcmp(variables{i},'t')
                    forcing.TEMP.(variables{i}) = forcing.DATA.(variables{i})(posit,1)+(forcing.DATA.(variables{i})(posit+1,1)-forcing.DATA.(variables{i})(posit,1)).*t_weight;
                end
            end
            forcing.TEMP.t = t;
        end
        
        function forcing = terrain_shade(forcing)
            az = forcing.TEMP.azimuth;
            el = forcing.TEMP.sunElevation;
            hbins = forcing.PARA.hbins;
            h = forcing.PARA.h;
            
            I = knnsearch(hbins,az); % hbin containing current solar azimuth
            forcing.TEMP.Sin_dir(h(I)>el) = 0; % remove direct Sin if sun is below horizon
            forcing.TEMP.Sin(h(I)>el) = forcing.TEMP.Sin_dif;
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
        
        function forcing = SolarAzEl(forcing, tile)
            % SolarAzEl will ingest a Universal Time, and specific site location on earth
            % it will then output the solar Azimuth and Elevation angles relative to
            % that site.
            %
            % Programed by Darin C. Koblick 2/17/2009
            %              Darin C. Koblick 4/16/2013 Vectorized for Speed
            %                             Allow for MATLAB Datevec input in
            %                             addition to a UTC string.
            %                             Cleaned up comments and code to
            %                             avoid warnings in MATLAB editor.
            %              Solar Position obtained from: http://stjarnhimlen.se/comp/tutorial.html#5
            Lat = tile.PARA.latitude;
            Lon = tile.PARA.longitude;
            Alt = tile.PARA.altitude ./ 1000;
            t_span = forcing.DATA.timeForcing;
            UTC = datestr(t_span);
            %i_end = size(UTC);
            
            %Loop through all timesteps
            %for i = 1 : i_end(1,1)
            
            %compute JD from UTC
            [year month day hour min sec] = datevec(UTC); %datevec(UTC(i,:));
            idx = (month <= 2);
            year(idx) = year(idx)-1;
            month(idx) = month(idx)+12;
            jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
                floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + ...
                (hour + min/60 + sec/3600)/24;
            d = jd-2451543.5;
            
            % Keplerian Elements for the Sun (geocentric)
            w = 282.9404+4.70935e-5*d;     %(longitude of perihelion degrees)
            e = 0.016709-1.151e-9.*d;       %(eccentricity)
            M = mod(356.0470+0.9856002585.*d,360);   %(mean anomaly degrees)
            L = w + M;                     %(Sun's mean longitude degrees)
            oblecl = 23.4393-3.563e-7.*d;  %(Sun's obliquity of the ecliptic)
            
            %auxiliary angle
            E = M+(180/pi).*e.*sin(M.*(pi/180)).*(1+e.*cos(M.*(pi/180)));
            
            %rectangular coordinates in the plane of the ecliptic (x axis toward perhilion)
            x = cos(E.*(pi/180))-e;
            y = sin(E.*(pi/180)).*sqrt(1-e.^2);
            
            %find the distance and true anomaly
            r = sqrt(x.^2 + y.^2);
            v = atan2(y,x).*(180/pi);
            
            %find the longitude of the sun
            lon = v + w;
            
            %compute the ecliptic rectangular coordinates
            xeclip = r.*cos(lon.*(pi/180));
            yeclip = r.*sin(lon.*(pi/180));
            zeclip = 0.0;
            
            %rotate these coordinates to equitorial rectangular coordinates
            xequat = xeclip;
            yequat = yeclip.*cos(oblecl.*(pi/180))+zeclip.*sin(oblecl.*(pi/180));
            zequat = yeclip.*sin(23.4406.*(pi/180))+zeclip.*cos(oblecl.*(pi/180));
            
            %convert equatorial rectangular coordinates to RA and Decl:
            r = sqrt(xequat.^2 + yequat.^2 + zequat.^2)-(Alt./149598000); %roll up the altitude correction
            RA = atan2(yequat,xequat).*(180/pi);
            delta = asin(zequat./r).*(180/pi);
            
            %Following the RA DEC to Az Alt conversion sequence explained here: http://www.stargazing.net/kepler/altaz.html
            
            %Find the J2000 value
            J2000 = jd - 2451545.0;
            %                 hourvec = datevec(UTC);  %(i,:)
            %                 UTH = hourvec(:,4) + hourvec(:,5)/60 + hourvec(:,6)/3600;
            UTH = hour + min./60 + sec./3600;
            
            %Calculate local siderial time
            GMST0=mod(L+180,360)./15;
            SIDTIME = GMST0 + UTH + Lon./15;
            
            %Replace RA with hour angle HA
            HA = (SIDTIME.*15 - RA);
            
            %convert to rectangular coordinate system
            x = cos(HA.*(pi/180)).*cos(delta.*(pi/180));
            y = sin(HA.*(pi/180)).*cos(delta.*(pi/180));
            z = sin(delta.*(pi/180));
            
            %rotate this along an axis going east-west.
            xhor = x.*cos((90-Lat).*(pi/180))-z.*sin((90-Lat).*(pi/180));
            yhor = y;
            zhor = x.*sin((90-Lat).*(pi/180))+z.*cos((90-Lat).*(pi/180));
            
            %Find the h and AZ
            forcing.DATA.azimuth = atan2(yhor,xhor).*(180./pi) + 180;
            forcing.DATA.sunElevation = asin(zhor).*(180./pi); %in degrees
            
        end
        
    end
end