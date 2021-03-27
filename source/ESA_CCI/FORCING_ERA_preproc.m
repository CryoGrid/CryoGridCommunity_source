
%========================================================================

classdef FORCING_ERA_preproc < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        TEMP
    end
    
    
    methods
        
        
        function forcing = provide_PARA(forcing)
            forcing.PARA.ERA_path = [];
            forcing.PARA.filename_SURF_airT = [];
            forcing.PARA.filename_SURF_geopotential = [];
            forcing.PARA.filename_SURF_precip = [];
            forcing.PARA.filename_PLEV_airT = [];
            forcing.PARA.filename_PLEV_geopotential = [];
            forcing.PARA.number_of_values_per_day = [];
            forcing.PARA.start_time = [];
            forcing.PARA.end_time = [];
            
        end
        
        
        
        function forcing = provide_CONST(forcing)
           forcing.CONST.g =[]; 
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        function forcing = finalize_init(forcing, tile)
            
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1), forcing.PARA.end_time(3,1));
            
            test_year = 2019; %choose a year that exists to read the time-constant data sets
            filename = forcing.PARA.filename_SURF_airT;
            filename(21:24) = num2str(test_year);

            ERA_lat = ncread([forcing.PARA.ERA_path filename], 'latitude');
            ERA_lon = ncread([forcing.PARA.ERA_path filename], 'longitude');
            ERA_lon =[ERA_lon; single(360)]; %add the 0-degree-longitude to the end as 360-degree-longitude
            
            max_target_lat = min(round(max(tile.RUN_INFO.STATVAR.latitude) *4) /4 + 0.5, 90) ;  %ERA5 in 0.25 degree resolution
            min_target_lat= max(round(min(tile.RUN_INFO.STATVAR.latitude)*4) /4 - 0.5, -90);
            
            %not needed?
%             target_lat = INFO.lat;
%             target_alt = INFO.altitude;
            
            target_lon = tile.RUN_INFO.STATVAR.longitude;
            target_lon(target_lon < 0) = target_lon(target_lon < 0) + 360; %change from -180 to +180 to 0 to 360
            max_target_lon = min(round(max(target_lon)*4) /4 + 0.5, 360);
            min_target_lon = max(round(min(target_lon)*4) / 4 - 0.5, 0);
            
            forcing.TEMP.lat_index_start = find(ERA_lat(:,1) == max_target_lat);  
            forcing.TEMP.lat_index_end = find(ERA_lat(:,1) == min_target_lat);
            forcing.TEMP.lon_index_start = find(ERA_lon(:,1) == min_target_lon);
            forcing.TEMP.lon_index_end = find(ERA_lon(:,1) == max_target_lon);
            
            forcing.TEMP.max_target_lon = max_target_lon;
            forcing.TEMP.target_lon = target_lon;
            
            filename = forcing.PARA.filename_SURF_geopotential;
            filename(19:22) = num2str(test_year);
            
            forcing.TEMP.z_surface = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename], 'z', 1, 1) ./ forcing.CONST.g;
            
            [forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh] = meshgrid(ERA_lon(forcing.TEMP.lon_index_start:forcing.TEMP.lon_index_end), ERA_lat(forcing.TEMP.lat_index_start:forcing.TEMP.lat_index_end));
            
            forcing.TEMP.interpolated_ERA_orography = single(interp2(forcing.TEMP.ERA_lon_mesh , forcing.TEMP.ERA_lat_mesh, forcing.TEMP.z_surface', target_lon, tile.RUN_INFO.STATVAR.latitude, 'cubic'));

        end
        
        function forcing = interpolate_forcing(forcing, tile)
            %disp('busy writing ERA5 files')

            forcing.TEMP.ERA_time=[];
            forcing.TEMP.ERA_T_downscaled=[];
            forcing.TEMP.ERA_precip_downcaled=[];
            
            for time = tile.t:tile.t + tile.timestep %9 days
                year = str2num(datestr(tile.t,'yyyy'));
                doy = tile.t - datenum(year,1,1) + 1;
                time_index_start = (doy-1) .* 4 + 1;
                
                filename_SURF_precip = forcing.PARA.filename_SURF_precip;
                filename_SURF_precip(26:29) = num2str(year);
                
                filename_SURF_airT = forcing.PARA.filename_SURF_airT;
                filename_SURF_airT(21:24) = num2str(year);
                
                filename_PLEV_airT = forcing.PARA.filename_PLEV_airT;
                filename_PLEV_airT(18:21) = num2str(year);
                
                filename_PLEV_geopotential = forcing.PARA.filename_PLEV_geopotential;
                filename_PLEV_geopotential(19:22) = num2str(year);
                

                if exist([forcing.PARA.ERA_path filename_SURF_airT])==2
                    T2m = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename_SURF_airT], 't2m', time_index_start, 4);
                    
                    tot_prec = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename_SURF_precip], 'tp', time_index_start, 4); %hourly resolution!!
                    
                    T_700 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_airT], 't', time_index_start, 4, 2);
                    
                    z_700 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_geopotential], 'z', time_index_start, 4,2) ./ forcing.CONST.g;
                    
                    T_500 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_airT], 't', time_index_start, 4, 3);
                    
                    z_500 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_geopotential], 'z', time_index_start, 4, 3) ./ forcing.CONST.g;
                    
                    T_200 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_airT], 't', time_index_start, 4, 4);
                    
                    z_200 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_geopotential], 'z', time_index_start, 4, 4) ./ forcing.CONST.g;
                    
                    diff700 = z_700 - forcing.TEMP.z_surface;
                    diff500 = z_500 - forcing.TEMP.z_surface;
                    diff200 = z_200 - forcing.TEMP.z_surface;
                    
                    factor700 = double(diff700>500) .* normpdf(diff700, 3000, 800); %this gives a half-way smooth transition, the 200mbar level is not optimal, more than 10km, should be a bit lower, like 300 or 350mbar
                    factor500 = double(diff500>500) .* normpdf(diff500, 3000, 800);
                    factor200 = double(diff200>500) .* normpdf(diff200, 5000, 1600);
                    
                    sum_factors = factor700+factor500+factor200; %normalize
                    factor700 = factor700./sum_factors;
                    factor500 = factor500./sum_factors;
                    factor200 = factor200./sum_factors;
                    
                    lapse_rate = factor700 .* (T_700-T2m) ./ (z_700 - forcing.TEMP.z_surface) + factor500 .* (T_500-T2m) ./ (z_500 - forcing.TEMP.z_surface) + factor200 .* (T_200-T2m) ./ (z_200 - forcing.TEMP.z_surface);

                    for i=1:4
                        T2m_slice=squeeze(T2m(:,:,i));
                        interpolated_T2m = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, T2m_slice', forcing.TEMP.target_lon, tile.RUN_INFO.STATVAR.latitude, 'cubic');
                        lr_slice=squeeze(lapse_rate(:,:,i)); 
                        interpolated_lapse_rate = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, lr_slice', forcing.TEMP.target_lon, tile.RUN_INFO.STATVAR.latitude, 'cubic');
                        
                        T_downscaled =  interpolated_T2m + interpolated_lapse_rate .* (tile.RUN_INFO.STATVAR.altitude - forcing.TEMP.interpolated_ERA_orography);
                        
                        %append values to file
                        
                        forcing.TEMP.ERA_time=[forcing.TEMP.ERA_time tile.t + (i-1).*0.25];
                        forcing.TEMP.ERA_T_downscaled=[forcing.TEMP.ERA_T_downscaled T_downscaled];
                        
                    end
                    
                    for jjj=1:4  %jjj=1:6:19
                        
                        tot_prec_slice = tot_prec(:,:,jjj) .* 24 ./ forcing.PARA.number_of_values_per_day; % 6;  %accumulate over 6h slice, assuming precip is constant
                        %   tot_prec_slice = sum(tot_prec(:,:,jjj:jjj+5),3);  %accumulate over 6h slice
                        interpolated_tot_prec = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, tot_prec_slice', forcing.TEMP.target_lon, tile.RUN_INFO.STATVAR.latitude, 'cubic');
                        
                        tot_prec_downscaled = interpolated_tot_prec .* 1.04.^((tile.RUN_INFO.STATVAR.altitude - forcing.TEMP.interpolated_ERA_orography)./100);  %5 percent increase per 100m elevation increase -> change to Jaros formula later
                        tot_prec_downscaled(tot_prec_downscaled<0) = 0;
                        %append values to file
                        forcing.TEMP.ERA_precip_downcaled =[forcing.TEMP.ERA_precip_downcaled tot_prec_downscaled.*1000];
                    end
                else
                    forcing.TEMP.ERA_time=[forcing.TEMP.ERA_time tile.t + [0 0.25 0.5 0.75]];
                    forcing.TEMP.ERA_T_downscaled=[forcing.TEMP.ERA_T_downscaled repmat(273.15-10,size(tile.RUN_INFO.STATVAR.key,1), 4)];
                    forcing.TEMP.ERA_precip_downcaled =[forcing.TEMP.ERA_precip_downcaled repmat(0,size(tile.RUN_INFO.STATVAR.key,1), 4)];
                end
                
                
            end
            
            
        end
        
        
        function output = ncread_with_360_degree(forcing, fname, variable, time_index_start, time_number_of_values)
            
            if forcing.TEMP.max_target_lon < 360 %"normal" case, no need to change anything
                output = ncread(fname, variable, [ forcing.TEMP.lon_index_start forcing.TEMP.lat_index_start time_index_start], [forcing.TEMP.lon_index_end - forcing.TEMP.lon_index_start+1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 time_number_of_values], [1 1 1]);
                
            else  %this part is not yet tested, it only becomes effective at the 0 degree longitude -> the only affected PF seems to be in the Pyrenees, which is fairly irrelevant
                forcing.TEMP.lon_index_end=forcing.TEMP.lon_index_end-1; %read to lon = 359.75 degree
                output = ncread(fname, variable, [ forcing.TEMP.lon_index_start forcing.TEMP.lat_index_start time_index_start], [forcing.TEMP.lon_index_end - forcing.TEMP.lon_index_start+1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 time_number_of_values], [1 1 1]);
                output0degree = ncread(fname, variable, [ 1 forcing.TEMP.lat_index_start time_index_start], [1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 time_number_of_values], [1 1 1]); %read first "line" of longitude = 0 degree
                output = cat(1, output, output0degree); %concat in dimension1 (longitude)
            end
            
        end
        
        function output = ncread_with_360_degree_pl(forcing, fname, variable, time_index_start, time_number_of_values, pl)
            
            %pl = 1 1000; pl = 2 700; pl = 3 500; pl = 4 200
            
            if forcing.TEMP.max_target_lon < 360 %"normal" case, no need to change anything
                output = ncread(fname, variable, [ forcing.TEMP.lon_index_start forcing.TEMP.lat_index_start pl time_index_start], [forcing.TEMP.lon_index_end - forcing.TEMP.lon_index_start+1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 1 time_number_of_values], [1 1 1 1]);
            else  %this part is not yet tested, it only becomes effective at the 0 degree longitude -> the only affected PF seems to be in the Pyrenees, which is fairly irrelevant
                forcing.TEMP.lon_index_end=forcing.TEMP.lon_index_end-1; %read to lon = 359.75 degree
                output = ncread(fname, variable, [ forcing.TEMP.lon_index_start forcing.TEMP.lat_index_start pl time_index_start], [forcing.TEMP.lon_index_end - forcing.TEMP.lon_index_start+1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 1 time_number_of_values], [1 1 1 1]);
                output0degree = ncread(fname, variable, [ 1 forcing.TEMP.lat_index_start pl time_index_start], [1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 1 time_number_of_values], [1 1 1 1]); %read first "line" of longitude
                output = cat(1, output, output0degree); %concat in dimension1 (longitude)
            end
        end
        
    end
end