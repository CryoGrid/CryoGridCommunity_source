%========================================================================
% Same as FORCING_slope_seb_readNC, with small adaptations to fit with 
% Vegetation classes.
% R. B. Zweigel, August 2021
%========================================================================

classdef FORCING_slope_seb_readNc2 < SEB %matlab.mixin.Copyable
    
    properties
        DATA            % forcing data time series
        STATUS        
    end
    
    
    methods
        
        %mandatory functions
        
        function forcing = provide_PARA(forcing)         
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  

            forcing.PARA.forcing_path = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)        
            forcing.PARA.all_rain_T = [];
            forcing.PARA.all_snow_T = [];
            forcing.PARA.slope_angle = []; %slope angle in degrees
            forcing.PARA.aspect = []; %aspect of the slope in degrees
            forcing.PARA.albedo_surrounding_terrain = [];
            forcing.PARA.sky_view_factor = []; %sky view factor (0.5 for vertical rock walls)
            %forcing.PARA.albedo = []; %Albedo of the slope
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
        end
        
        function forcing = provide_CONST(forcing)
            forcing.CONST.Tmfw = [];
            forcing.CONST.sigma = [];
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end

        
        function forcing = finalize_init(forcing, tile)
            
          
            %variables = {'t2m'; 'd2m'; 'u10'; 'v10'; 'ssrd'; 'strd'; 'tp'};
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
            
            %forcing.DATA.timeForcing = datenum(1900,1,1) + double(temp.time)./24;
            forcing.DATA.Tair = temp.t2m - forcing.CONST.Tmfw;
            forcing.DATA.wind = sqrt(temp.u10.^2 + temp.u10.^2);
            forcing.DATA.Sin = temp.ssrd./3600;
            forcing.DATA.Sin = [0;0;forcing.DATA.Sin];
            forcing.DATA.Lin = temp.strd./3600;
            forcing.DATA.Lin = [0;0;forcing.DATA.Lin];
            forcing.DATA.p = temp.sp; %in [Pa]!900.*100 + forcing.DATA.Tair.*0;
            forcing.DATA.q = (double(forcing.DATA.Tair<0).*satPresIce(forcing, temp.d2m) + double(forcing.DATA.Tair>=0).*satPresWater(forcing, temp.d2m))./ forcing.DATA.p;
            
            forcing.DATA.S_TOA = temp.tisr ./ 3600;
            forcing.DATA.S_TOA = [0;0;forcing.DATA.S_TOA];
            
%             size(temp.tp)
%             size(forcing.DATA.Tair)
            
            temp.tp = [0;0; temp.tp]; %append missing first two timesteps
            forcing.DATA.snowfall = temp.tp .*24.*1000 .* (double(forcing.DATA.Tair <= forcing.PARA.all_snow_T)  + ...
                double(forcing.DATA.Tair > forcing.PARA.all_snow_T & forcing.DATA.Tair <= forcing.PARA.all_rain_T) .* ...
                (forcing.DATA.Tair - forcing.PARA.all_snow_T) ./ max(1e-12, (forcing.PARA.all_rain_T - forcing.PARA.all_snow_T)));
            forcing.DATA.rainfall = temp.tp .*24.*1000 .* (double(forcing.DATA.Tair >= forcing.PARA.all_rain_T)  + ...
                double(forcing.DATA.Tair > forcing.PARA.all_snow_T & forcing.DATA.Tair < forcing.PARA.all_rain_T) .* ...
                (1 - (forcing.DATA.Tair - forcing.PARA.all_snow_T) ./ max(1e-12, (forcing.PARA.all_rain_T - forcing.PARA.all_snow_T))));
            forcing.DATA.rainfall = forcing.DATA.rainfall .* forcing.PARA.rain_fraction.*cosd(forcing.PARA.slope_angle);
            forcing.DATA.snowfall = forcing.DATA.snowfall .* forcing.PARA.snow_fraction.*cosd(forcing.PARA.slope_angle);
            
            forcing.DATA.albedo_foot = forcing.PARA.albedo_surrounding_terrain; %Albedo at the foot of the slope

            
%             % Additional forcing data for slopes:
%             forcing.DATA.S_TOA = temp.FORCING.data.S_TOA; %short-wave radiation at the top of the atmosphere

%             
%             % Non-mandatory forcing data for slopes:
%             if isfield(temp.FORCING.data,'seaT') == 1
%                 forcing.DATA.seaT = temp.FORCING.data.seaT; %seawater temperature
%             end
%             if isfield(temp.FORCING.data,'seaIce') == 1
%                 forcing.DATA.seaIce = temp.FORCING.data.seaIce; %time steps with (1) or without (0) sea ice
%             end
            
            if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1))~=0
                disp('timestamp of forcing data is not in regular intervals -> check, fix and restart')
                forcing.STATUS=0;
                return
            else
                forcing.STATUS=1;
            end
            
            %here, consistency checks, RH->q calculation, set threhsolds for wind, etc. could be placed
            
            forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence
            forcing.DATA.Lin(find(forcing.DATA.Lin==0)) = 5.67e-8 .* (forcing.DATA.Tair(find(forcing.DATA.Lin==0))+273.15).^4;
            
            
            if isempty(forcing.PARA.start_time) || isnan(forcing.PARA.start_time(1,1)) % || ~ischar(forcing.PARA.start_time)
                forcing.PARA.start_time = forcing.DATA.timeForcing(1,1);
            else
                %forcing.PARA.start_time = datenum(forcing.PARA.start_time, 'dd.mm.yyyy');
                forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            end
            if isempty(forcing.PARA.end_time) || isnan(forcing.PARA.end_time(1,1)) %|| ~ischar(forcing.PARA.end_time)
                forcing.PARA.end_time = floor(forcing.DATA.timeForcing(end,1));
            else
                %forcing.PARA.end_time = datenum(forcing.PARA.end_time, 'dd.mm.yyyy');
                forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));
            end
            
            %initialize TEMP
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Lin=0;
            forcing.TEMP.Sin=0;
            forcing.TEMP.Tair=0;
            forcing.TEMP.wind=0;
            forcing.TEMP.q=0;
            forcing.TEMP.p=0;
            
            % Additional forcing data for slopes:
            forcing.TEMP.S_TOA=0;
            forcing.TEMP.albedo_foot=0;
            forcing.TEMP.Sin_dif = 0;
            forcing.TEMP.Sin_dir = 0;
            forcing.TEMP.sunElevation = 0;
            
            % Non-mandatory forcing data for slopes:
%             if isfield(temp.FORCING.data,'seaT') == 1
%                 forcing.TEMP.seaT=0;
%             end
%             if isfield(temp.FORCING.data,'seaIce') == 1
%                 forcing.TEMP.seaIce=0;
%             end
            
             forcing = scale_radiation(tile.PARA, forcing);
        end
            

        function forcing = interpolate_forcing(forcing, tile)
            t = tile.t;
            
            posit=floor((t-forcing.DATA.timeForcing(1,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)))+1;
            
            forcing.TEMP.snowfall=forcing.DATA.snowfall(posit,1)+(forcing.DATA.snowfall(posit+1,1)-forcing.DATA.snowfall(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.rainfall=forcing.DATA.rainfall(posit,1)+(forcing.DATA.rainfall(posit+1,1)-forcing.DATA.rainfall(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.Lin=forcing.DATA.Lin(posit,1)+(forcing.DATA.Lin(posit+1,1)-forcing.DATA.Lin(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.Sin=forcing.DATA.Sin(posit,1)+(forcing.DATA.Sin(posit+1,1)-forcing.DATA.Sin(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.Tair=forcing.DATA.Tair(posit,1)+(forcing.DATA.Tair(posit+1,1)-forcing.DATA.Tair(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.wind=forcing.DATA.wind(posit,1)+(forcing.DATA.wind(posit+1,1)-forcing.DATA.wind(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.q=forcing.DATA.q(posit,1)+(forcing.DATA.q(posit+1,1)-forcing.DATA.q(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.p=forcing.DATA.p(posit,1)+(forcing.DATA.p(posit+1,1)-forcing.DATA.p(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
            forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
            
            forcing.TEMP.Sin_dif=forcing.DATA.Sin_dif(posit,1)+(forcing.DATA.Sin_dif(posit+1,1)-forcing.DATA.Sin_dif(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            forcing.TEMP.Sin_dir=forcing.DATA.Sin_dir(posit,1)+(forcing.DATA.Sin_dir(posit+1,1)-forcing.DATA.Sin_dir(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            forcing.TEMP.sunElevation=forcing.DATA.sunElevation(posit,1)+(forcing.DATA.sunElevation(posit+1,1)-forcing.DATA.sunElevation(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));

            forcing.TEMP.t = t;
            
        end
        

        
        % non-mandatory functions
        
        function forcing = scale_radiation(PARA, forcing)
            
            % Projection of the forcing data to the slope
            % This function includes:
            % 1. Calculation of direct Sin and projection of direct Sin from a
            % from a horizontal surface to an inclined slope
            % 2. Calculation of diffuse Sin and reduction by a sky view factor
            % 3. Calculation of the reflected Sin on the foot of the slope
            % 4. Reduction of Lin by a sky view factor
            % 5. Calculation of long-wave emission of the close environment
        
            lat = PARA.latitude;
            lon = PARA.longitude;
            alt = PARA.altitude / 1000;
            aspect = forcing.PARA.aspect;
            slope_angle = forcing.PARA.slope_angle;
            %albedo = forcing.PARA.albedo;
            sky_view_factor = forcing.PARA.sky_view_factor;
            
            t_span = forcing.DATA.timeForcing;
            %UTC = datestr(t_span);
            %i_end = size(UTC); %LOESCHEN??
            
            %Diffuse and direct Sin
            
            kt = forcing.DATA.Sin ./ forcing.DATA.S_TOA; %clearness index, see Fiddes(2014)
            kt(isnan(kt))=0;
            kd = 0.952 - 1.041.*exp(-exp(2.300 - 4.702 .* kt)); %diffuse fraction, see Fiddes(2014)
            kd(isnan(kd))=0; 
            kd = max(0,kd);
            
            SW_diff_total = kd .* forcing.DATA.Sin; %Amount of total diffuse SW
            SW_dir_total = (1 - kd) .* forcing.DATA.Sin; %Amount of total direct SW
            
            %Calculation of diffuse SW

            SW_diff = SW_diff_total .* sky_view_factor; %Reduction of diffuse SW by sky view factor

            % Calculation of reflected SW

            SW_refl = (forcing.DATA.Sin .* forcing.DATA.albedo_foot) .* (1 - sky_view_factor); % corrected RBZ Aug-21
            
            %Calculation of reprojected direct SW

            surf_norm_vec = repmat([0.0,0.0,1.0], size(forcing.DATA.Sin,1), 1); %Unit vector on the horizontal
            sigma = forcing.CONST.sigma; %5.670373e-08; %Stefan-Boltzman constant

            alpha = aspect.*pi./180; %Degree to radians of exposition of the slope
            beta = (90-slope_angle).*pi./180; %Degree to radians of inclination of the slope
            face_vec = repmat([sin(alpha).*cos(beta) cos(alpha).*cos(beta) sin(beta)], size(forcing.DATA.Sin,1), 1); %Unit vector of the slope

            % Calculation the solar azimuth and elevation angle relative to
            % the coordinates of the site (revised after Darin C. Koblick)
            [Az,El] = SolarAzEl(PARA,forcing); %Calculation of the azimuth and elevation of the sun
%             i_end = size(Az);
% 
%             for i = 1 : i_end(1,1)
                
                alpha = Az.*pi/180; %Degree to radians of the azimuth
                beta = El.*pi/180; %Degree to radians of elevation
                sun_vec = [sin(alpha).*cos(beta) cos(alpha).*cos(beta) sin(beta)]; %Unit vector of the radiation
                
%                 size(surf_norm_vec)
%                 size(sun_vec)
                
                delta_angle_surf_norm = acos(dot(surf_norm_vec' ,sun_vec')').*180./pi; %angle between the radiation and the normal on the horizontal in degrees
                delta_angle_face=acos(dot(face_vec', sun_vec')').*180./pi; %angle between the radiation and the normal on the slope in degrees
                
                Sin_sun_direction = SW_dir_total./cos(delta_angle_surf_norm.*pi./180); %Sin in direction of the sun
                
%                 if delta_angle_surf_norm(i,1)<85 && delta_angle_face(i,1)<85
%                     %Sun has to be at least 5 degrees above the horizon and has to shine at
%                     %least with 85 degrees on the wall (otherwise values get too high)
%                     Sin_face_direction(i,1)=Sin_sun_direction(i,1)*cos(delta_angle_face(i,1)*pi/180); %Sin normal to the slope
%                 else %Sun is under the horizon
%                     Sin_face_direction(i,1)=0.0; %No short-wave radiation if sun is under the horizon or behind the wall
%                 end
                Sin_face_direction = double(delta_angle_surf_norm < 85 & delta_angle_face < 85) .* Sin_sun_direction.*cos(delta_angle_face.*pi/180); 
            
            % Calculation of Lin

%             if isfield(forcing.DATA,'seaT') == 1 %water at the foot of the slope
%                 if forcing.DATA.seaIce(i,1) == 0 %no sea ice --> take sea temperature
%                 Lin(i,1) = sky_view_factor * forcing.DATA.Lin(i,1) + (1 - sky_view_factor) * sigma * (forcing.DATA.seaT(i,1) + 273.15)^4; %T of ocean for calculation
%                 elseif forcing.DATA.seaIce(i,1) == 1 %sea ice --> take air temperature
%                 Lin(i,1) = sky_view_factor * forcing.DATA.Lin(i,1) + (1 - sky_view_factor) * sigma * (forcing.DATA.Tair(i,1) + 273.15)^4; %Tair for calculation
%                 end
%             else %no water at the foot of the slope
            Lin = sky_view_factor .* forcing.DATA.Lin + (1 - sky_view_factor) .* sigma .* (forcing.DATA.Tair + forcing.CONST.Tmfw).^4; %Tair for calculation
%             end
         
            %end
            
            forcing.DATA.Sin = Sin_face_direction + SW_diff + SW_refl;
            forcing.DATA.Lin = Lin;
            
            forcing.DATA.Sin_dif = SW_diff + SW_refl;
            forcing.DATA.Sin_dir = Sin_face_direction;
            
            
            forcing.DATA.azimuth = Az;
            forcing.DATA.sunElevation = El;
            
        end
            

        function [Az,El] = SolarAzEl(PARA,forcing)
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
            Lat = PARA.latitude;
            Lon = PARA.longitude;
            Alt = PARA.altitude ./ 1000;
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
            Az = atan2(yhor,xhor).*(180./pi) + 180; %Az(i,1)
            El = asin(zhor).*(180./pi); %El(i,1)  %in degrees

                %clearvars -except Az El Lat Lon Alt UTC path_out name forcing
           % end
        end

                
        end
end