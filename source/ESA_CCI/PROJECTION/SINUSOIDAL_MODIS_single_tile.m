classdef SINUSOIDAL_MODIS_single_tile < matlab.mixin.Copyable

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
            proj.STATVAR.list_of_MODIS_tiles = [proj.PARA.horizontal proj.PARA.vertical 1 size(proj.STATVAR.key,1)]; %use to loop over several tiles, make class MODIS_multi_tile

            file_info = dir([proj.PARA.MODIS_path '*.hdf']);
            file_info = struct2cell(file_info);
            file_info = file_info';
            file_names = file_info(:,1);
            proj.TEMP.file_names = cell2mat(file_names);
            proj.TEMP.year_list = str2num(proj.TEMP.file_names(:,10:13));
            proj.TEMP.doy_list = str2num(proj.TEMP.file_names(:,14:16));
        end
        
        function proj = get_data(proj, time, key)
            
            proj.STATVAR.time = [];
            proj.STATVAR.LST = [];
            
            %computes the day f year
            f=zeros(size(time));
            [y,m,d,~,~,~]=datevec(time);
            doy = datenum(y,m,d,f,f,f)-datenum(y,f,f,f,f,f);
            
            
            ID = find(proj.TEMP.year_list(:,1) == str2num(datestr(time,'yyyy')) & proj.TEMP.doy_list(:,1) == doy );
            if ~isempty(ID)
                
                for j=1:size(ID,1)  %loops over Terra and Aqua
                    %try
                        Day_time = hdfread([proj.PARA.MODIS_path  proj.TEMP.file_names(ID(j),:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'Day_view_time');
                        Day_time = double(Day_time(key));
                        
                        Day_time(Day_time==255)=NaN;
                        Day_time = time + Day_time.*0.1./24; %MODIS time in local solar noon
                        Day_time = SolarTime2LocalTime(proj, Day_time, doy, proj.STATVAR.longitude(key), 0); %convert to UTC using longitudes; IMPORTANT: use last argument 0 to convert to UTC
                        proj.STATVAR.time = [proj.STATVAR.time Day_time(:)];
                        
                        LST_Day = hdfread([proj.PARA.MODIS_path proj.TEMP.file_names(ID(j),:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'LST_Day_1km');
                        LST_Day = LST_Day(key);
                        proj.STATVAR.LST=[proj.STATVAR.LST double(LST_Day(:)).*0.02];
                        
                        
                        Night_time = hdfread([proj.PARA.MODIS_path  proj.TEMP.file_names(ID(j),:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'Night_view_time');
                        Night_time = double(Night_time(key));
                        Night_time(Night_time==255)=NaN;
                        Night_time = time + Night_time.*0.1./24; %MODIS time in local solar noon
                        Night_time = SolarTime2LocalTime(proj, Night_time, doy, proj.STATVAR.longitude(key), 0); %convert to UTC using longitudes; IMPORTANT: use last argument 0 to convert to UTC
                        proj.STATVAR.time=[proj.STATVAR.time Night_time(:)];
                        %timestamp_LST=[timestamp_LST i];
                        
                        LST_Night = hdfread([proj.PARA.MODIS_path  proj.TEMP.file_names(ID(j),:)], 'MODIS_Grid_Daily_1km_LST', 'Fields', 'LST_Night_1km');
                        LST_Night = LST_Night(key);
                        proj.STATVAR.LST=[proj.STATVAR.LST double(LST_Night(:)).*0.02];
                    %end
                end
            end
                       
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

