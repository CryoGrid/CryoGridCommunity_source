classdef SINUSOIDAL_MODIS_multi_tile < matlab.mixin.Copyable

    properties
        RUN_INFO
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods
        function proj = provide_PARA(proj)
            proj.PARA.horizontal = [];
            proj.PARA.vertical = [];
            proj.PARA.number_of_pixels = [];
            proj.PARA.MODIS_path = [];
            
            proj.PARA.mask_class = [];
            proj.PARA.mask_class_index = [];
        end
        
        function proj = provide_STATVAR(proj)
            proj.STATVAR.time = [];
            proj.STATVAR.LST = [];
        end
        
        function proj = provide_CONST(proj)
            
        end
        
        function proj = finalize_init(proj)
            %store temprorary variables so that exact structure of
            %SINUSOIDAL_MODIS_single_tile can be used
            h_list = proj.PARA.horizontal;
            v_list = proj.PARA.vertical;
            proj.PARA.hv_string = [];
            
            
            lat_final = [];
            lon_final = [];
            key_final = [];
            list_of_MODIS_tiles = [];
            start_index = 1;
            
            for index=1:size(h_list,1)
                proj.PARA.horizontal = h_list(index,1);
                proj.PARA.vertical = v_list(index,1);
            
                proj = sin2llMOD(proj);
                
                %apply masks
                proj.STATVAR.mask = logical(zeros(proj.PARA.number_of_pixels, proj.PARA.number_of_pixels));
                for i=1:size(proj.PARA.mask_class_index,1)
                    mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.mask_class{i,1}){proj.PARA.mask_class_index(i,1),1});
                    mask_class.PARENT = proj;
                    mask_class = finalize_init(mask_class);
                    mask_class = apply_mask(mask_class); %can be additive or subtractive
                end
                
                
                proj.STATVAR.latitude = proj.STATVAR.latitude(:);
                proj.STATVAR.longitude = proj.STATVAR.longitude(:);
                proj.STATVAR.key = [1:size(proj.STATVAR.latitude,1)]';
                
                proj.STATVAR.latitude(~proj.STATVAR.mask)=[];
                proj.STATVAR.longitude(~proj.STATVAR.mask)=[];
                proj.STATVAR.key(~proj.STATVAR.mask)=[];
                if~isempty(proj.STATVAR.key)
                    
                    lat_final = [lat_final; proj.STATVAR.latitude];
                    lon_final = [lon_final; proj.STATVAR.longitude];
                    key_final = [key_final; proj.STATVAR.key];

                    list_of_MODIS_tiles = [list_of_MODIS_tiles; [proj.PARA.horizontal proj.PARA.vertical start_index start_index - 1 + size(proj.STATVAR.key,1)]]; %use to loop over several tiles, make class MODIS_multi_tile
                    
                    start_index = start_index + size(proj.STATVAR.key,1);
                    
                    
                    h = h_list(index,1);
                    h=h+100;
                    h=num2str(h);
                    h=h(2:3);
                    v = v_list(index,1);
                    v=v+100;
                    v=num2str(v);
                    v=v(2:3);

                    proj.PARA.hv_string = [proj.PARA.hv_string; ['h' h 'v' v]];
                end
            end
            
            proj.STATVAR.latitude = lat_final;
            proj.STATVAR.longitude = lon_final;
            proj.STATVAR.key = key_final;
            proj.STATVAR.list_of_MODIS_tiles = list_of_MODIS_tiles;
            
            proj.PARA.horizontal = h_list;
            proj.PARA.vertical = v_list;
            
        end
        
%         function proj = get_data(proj, start_time, number_of_days)
%             
%             key = proj.RUN_INFO.STATVAR.key;
%             
%             proj.STATVAR.timestamp = [];
%             proj.STATVAR.LST = [];
%             for index= 1:size(proj.RUN_INFO.STATVAR.list_of_MODIS_tiles,1)
%                 
%                 timestamp =[];
%                 LST=[];
%                 
%                 for time = start_time:start_time + number_of_days-1
%                     
%                     %computes the day f year
%                     f=zeros(size(time));
%                     [y,m,d,~,~,~]=datevec(time);
%                     doy = datenum(y,m,d,f,f,f)-datenum(y,f,f,f,f,f);
%                     doy_str=doy+1000;
%                     doy_str=num2str(doy_str);
%                     doy_str=doy_str(2:4);
% 
%                     range = [proj.RUN_INFO.STATVAR.list_of_MODIS_tiles(index,3):proj.RUN_INFO.STATVAR.list_of_MODIS_tiles(index,4)]';
%                     key = proj.RUN_INFO.STATVAR.key(range,1);
%                     
%                     ID = dir([proj.PARA.MODIS_path proj.PARA.hv_string(index,:) '/M*11A1.A'  datestr(time,'yyyy') doy_str '.' proj.PARA.hv_string(index,:) '*.hdf']);
%                     ID = struct2cell(ID);
%                     ID = ID';
%                     ID = ID(:,1);
%                     ID = cell2mat(ID);
%                     %ID = find(proj.TEMP.year_list(:,1) == str2num(datestr(time,'yyyy')) & proj.TEMP.doy_list(:,1) == doy );
%                     if ~isempty(ID)
%                         for j=1:size(ID,1)  %loops over Terra and Aqua
%                             %try
%                             Day_time = hdfread([proj.PARA.MODIS_path  proj.PARA.hv_string(index,:) '/' ID(j,:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'Day_view_time');
%                             Day_time = double(Day_time(key));
%                             
%                             Day_time(Day_time==255)=NaN;
%                             Day_time = time + Day_time.*0.1./24; %MODIS time in local solar noon
%                             Day_time = SolarTime2LocalTime(proj, Day_time, doy, proj.RUN_INFO.STATVAR.longitude(range), 0); %convert to UTC using longitudes; IMPORTANT: use last argument 0 to convert to UTC
%                             timestamp = [timestamp Day_time(:)];
%                             
%                             LST_Day = hdfread([proj.PARA.MODIS_path proj.PARA.hv_string(index,:) '/' ID(j,:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'LST_Day_1km');
%                             LST_Day = LST_Day(key);
%                             LST=[LST double(LST_Day(:)).*0.02];
%                             
%                             
%                             Night_time = hdfread([proj.PARA.MODIS_path  proj.PARA.hv_string(index,:) '/' ID(j,:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'Night_view_time');
%                             Night_time = double(Night_time(key));
%                             Night_time(Night_time==255)=NaN;
%                             Night_time = time + Night_time.*0.1./24; %MODIS time in local solar noon
%                             Night_time = SolarTime2LocalTime(proj, Night_time, doy, proj.RUN_INFO.STATVAR.longitude(range), 0); %convert to UTC using longitudes; IMPORTANT: use last argument 0 to convert to UTC
%                             timestamp=[timestamp Night_time(:)];
%                             %timestamp_LST=[timestamp_LST i];
%                             
%                             LST_Night = hdfread([proj.PARA.MODIS_path  proj.PARA.hv_string(index,:) '/' ID(j,:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'LST_Night_1km');
%                             LST_Night = LST_Night(key);
%                             LST=[LST double(LST_Night(:)).*0.02];
%                             %end
%                         end
%                     end
%                 end
%                 proj.STATVAR.timestamp = [proj.STATVAR.timestamp; timestamp];
%                 proj.STATVAR.LST = [proj.STATVAR.LST; LST];
%             end
%             proj.STATVAR.LST = single(proj.STATVAR.LST);
%             
%             if size(proj.STATVAR.timestamp,1) == 0
%                 proj.STATVAR.timestamp = key.*NaN;
%             end
%             
%             if size(proj.STATVAR.LST,1) == 0
%                 proj.STATVAR.LST = key.*NaN;
%             end
%             
%         end
%         
        function proj = get_data(proj, start_time, number_of_days)
            
            key = proj.RUN_INFO.STATVAR.key;
            
            proj.STATVAR.timestamp = zeros(size(key,1), number_of_days.*4) .* NaN ;
            proj.STATVAR.LST = zeros(size(key,1), number_of_days.*4) ;

            for index= 1:size(proj.RUN_INFO.STATVAR.list_of_MODIS_tiles,1)

                i= 1;
                for time = start_time:start_time + number_of_days-1

                    
                    %computes the day f year
                    f=zeros(size(time));
                    [y,m,d,~,~,~]=datevec(time);
                    doy = datenum(y,m,d,f,f,f)-datenum(y,f,f,f,f,f);
                    doy_str=doy+1000;
                    doy_str=num2str(doy_str);
                    doy_str=doy_str(2:4);

                    range = [proj.RUN_INFO.STATVAR.list_of_MODIS_tiles(index,3):proj.RUN_INFO.STATVAR.list_of_MODIS_tiles(index,4)]';
                    key = proj.RUN_INFO.STATVAR.key(range,1);
                    
                    ID = dir([proj.PARA.MODIS_path proj.PARA.hv_string(index,:) '/M*11A1.A'  datestr(time,'yyyy') doy_str '.' proj.PARA.hv_string(index,:) '*.hdf']);
                    ID = struct2cell(ID);
                    ID = ID';
                    ID = ID(:,1);
                    ID = cell2mat(ID);
                    %ID = find(proj.TEMP.year_list(:,1) == str2num(datestr(time,'yyyy')) & proj.TEMP.doy_list(:,1) == doy );
                    if ~isempty(ID)
                        for j=1:size(ID,1)  %loops over Terra and Aqua
                            %try
                            Day_time = hdfread([proj.PARA.MODIS_path  proj.PARA.hv_string(index,:) '/' ID(j,:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'Day_view_time');
                            Day_time = double(Day_time(key));
                            
                            Day_time(Day_time==255)=NaN;
                            Day_time = time + Day_time.*0.1./24; %MODIS time in local solar noon
                            Day_time = SolarTime2LocalTime(proj, Day_time, doy, proj.RUN_INFO.STATVAR.longitude(range), 0); %convert to UTC using longitudes; IMPORTANT: use last argument 0 to convert to UTC
                            %timestamp = [timestamp Day_time(:)];
                            proj.STATVAR.timestamp(range,i) = Day_time(:);
                            
                            LST_Day = hdfread([proj.PARA.MODIS_path proj.PARA.hv_string(index,:) '/' ID(j,:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'LST_Day_1km');
                            LST_Day = LST_Day(key);
                            %LST=[LST double(LST_Day(:)).*0.02];
                            proj.STATVAR.LST(range,i) = single(LST_Day(:)).*0.02;
                            
                            Night_time = hdfread([proj.PARA.MODIS_path  proj.PARA.hv_string(index,:) '/' ID(j,:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'Night_view_time');
                            Night_time = double(Night_time(key));
                            Night_time(Night_time==255)=NaN;
                            Night_time = time + Night_time.*0.1./24; %MODIS time in local solar noon
                            Night_time = SolarTime2LocalTime(proj, Night_time, doy, proj.RUN_INFO.STATVAR.longitude(range), 0); %convert to UTC using longitudes; IMPORTANT: use last argument 0 to convert to UTC
                            %timestamp=[timestamp Night_time(:)];
                            proj.STATVAR.timestamp(range,i+1) = Night_time(:);
                            
                            LST_Night = hdfread([proj.PARA.MODIS_path  proj.PARA.hv_string(index,:) '/' ID(j,:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'LST_Night_1km');
                            LST_Night = LST_Night(key);
                            %LST=[LST double(LST_Night(:)).*0.02];
                            proj.STATVAR.LST(range,i+1) = single(LST_Night(:)).*0.02;
                            
                            i=i+2;
                            
                        end
                    end
                end
%                 proj.STATVAR.timestamp = [proj.STATVAR.timestamp; timestamp];
%                 proj.STATVAR.LST = [proj.STATVAR.LST; LST];
            end
            
            proj.STATVAR.LST = single(proj.STATVAR.LST);
            
            if size(proj.STATVAR.timestamp,1) == 0
                proj.STATVAR.timestamp = key.*NaN;
            end
            
            if size(proj.STATVAR.LST,1) == 0
                proj.STATVAR.LST = key.*NaN;
            end
            
        end
        
    
        
        function plot_MODIS_multi_tile(proj, variable)
            worldmap([min(proj.RUN_INFO.STATVAR.latitude) max(proj.RUN_INFO.STATVAR.latitude)], [min(proj.RUN_INFO.STATVAR.longitude) max(proj.RUN_INFO.STATVAR.longitude)])
            
            PARA2 = proj.PARA;
            
            for index = 1:size(proj.RUN_INFO.STATVAR.list_of_MODIS_tiles,1)
                mask = zeros(proj.PARA.number_of_pixels, proj.PARA.number_of_pixels) .* NaN;
                range = [proj.RUN_INFO.STATVAR.list_of_MODIS_tiles(index,3):proj.RUN_INFO.STATVAR.list_of_MODIS_tiles(index,4)]';
%                 size(range)
%                 proj.RUN_INFO.STATVAR.key(range)
                mask(proj.RUN_INFO.STATVAR.key(range)) = variable(range);
                
                proj.PARA.horizontal = PARA2.horizontal(index,1);
                proj.PARA.vertical = PARA2.vertical(index,1);
                
                proj = sin2llMOD(proj);
                
                surfm(proj.STATVAR.latitude, proj.STATVAR.longitude, mask)
            end
            proj.PARA = PARA2;
            
        end
        
        
        
        function proj = sin2llMOD(proj)
            % number_of_pixels: number of pixels in east direction
            % horizontal: horz-tile
            % vertical: vert-tile
            % convert sinusoidal projection like MODIS products to long-lat.
            %
            % from:
            % MODIS Collection 5 Burned Area Product - MCD45
            % User?s Guide Version 2.0, November 2009
            % Luigi Boschetti, University of Maryland
            % David Roy, South Dakota State University
            % Anja A. Hoffmann, LM University of Munich
            %
            % PAGE 30-31
            % The MODIS data are re-projected using an equiareal sinusoidal projection, defined on a
            % sphere of radius
            %rho=6371007.181 m,
            %and with the Greenwich meridian as the central meridian of the projection.
            %Defining (x, y) as the East and North coordinate in meters of a point in the map space, and
            %its latitude and longitude in degrees, the direct formulas is:
            %
            % Torbjørn, nov 2011
         

            rho     = 6371007.181;  % sphere radius [m]
            t       = 1111950;      % m
            ny      = proj.PARA.number_of_pixels;           % number of pixels in north
            sz_pix  = t./proj.PARA.number_of_pixels;         % pixel size
            
            % array
            i=0:proj.PARA.number_of_pixels-1;
            j=0:proj.PARA.number_of_pixels-1;
            
            % zero in london
            
            h = proj.PARA.horizontal-18;
            
            % mesh, since x,y interact in 2longlat
            [i,j]= meshgrid(i,j);
            
            % projection in map space
            x = (i +.5).*sz_pix + h.*t;
            y = (9-proj.PARA.vertical).*t - (j+.5).*sz_pix;
            
            
            % to long lat
            proj.STATVAR.latitude     = y .* (180./pi)./rho;
            proj.STATVAR.longitude    = x./(rho.*cosd((y.*180)./(rho.*pi))) .* 180./pi;
            
            %conversion between geographic and sinusodial coordinates using minvtran
            %inialisation of coordinate system struct
            % sinmstruct = defaultm('sinusoid');
            % sinmstruct.origin = [0 0 0];
            % % %using wgs84 elipsoid
            % sinmstruct.geoid = almanac('earth','wgs84','meters');
            % %using spheroid
            % %sinmstruct.geoid = [6371007.181 0];
            % sinmstruct = defaultm(sinmstruct);
        end
        
        function LT = SolarTime2LocalTime(proj, LST_time, day_of_year, longitude, GMTdiff)
            % LT = SolarTime2LocalTime(LST,longitude,GMTdiff)
            %
            % LSTM: Local Standard Time Meridian (LSTM)
            % LST : local solar time
            % LT  : local time
            % EoT : equation of time
            % GMT diff: local time zone diff from GMT
            %
            % AWS1 @ ASF
            % convert MODIS local solar time to AWS UTC+2 time
            % longitude = 22.4;
            % GMTdiff   = 2; (summer time)
            
            
            % LSTM
            LSTM = 15 * GMTdiff; %time diff from GMT (zone)
            
            % EoT: equation of time (minutes)
            % The equation of time (EoT) (in minutes) is an empirical equation that corrects for the eccentricity of the Earth's orbit and the Earth's axial tilt.
            B   = (360/365) * (day_of_year-81);
            EoT = 9.87*sind(2*B) - 7.53*cosd(B) + 1.5*sind(B);
            
            % Time correction factor (minutes)
            %The net Time Correction Factor (in minutes) accounts for the variation of the Local Solar Time (LST) within a given time zone due to the longitude variations within the time zone and also incorporates the EoT above.
            
            TC = 4*(LSTM-longitude) + EoT;
            
            
            % LT (daynums) Reversed sign since we want local time
            LT = LST_time + TC/(60*24);
        end
        
    end
end

