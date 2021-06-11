
%========================================================================

classdef FORCING_ERA_preproc_apply_anomaly_SEB < FORCING_ERA_preproc_SEB
    

    
    
    methods
        
        
        function forcing = provide_PARA(forcing)
            
            forcing.PARA.ERA_path = [];
            forcing.PARA.filename_SURF_airT = [];
            forcing.PARA.filename_SURF_geopotential = [];
            forcing.PARA.filename_SURF_dewpointT = [];
            forcing.PARA.filename_SURF_solar_radiation = [];
            forcing.PARA.filename_SURF_thermal_radiation = [];
            forcing.PARA.filename_SURF_TOA_radiation = [];
            forcing.PARA.filename_SURF_precip = [];
            forcing.PARA.filename_PLEV_airT = [];
            forcing.PARA.filename_PLEV_geopotential = [];
            forcing.PARA.filename_PLEV_RH = [];
            
            forcing.PARA.anomaly_file_path =[];
            forcing.PARA.anomaly_file = [];
            
            forcing.PARA.number_of_values_per_day = [];
            forcing.PARA.start_time = [];
            forcing.PARA.end_time = [];
            
        end
        

        
        function forcing = interpolate_forcing(forcing, tile)

                forcing = interpolate_forcing@FORCING_ERA_preproc_SEB(forcing, tile);
                forcing = apply_anomaly(forcing, tile);
                
        end
            
           
        function forcing = apply_anomaly(forcing, tile)
            
            %                 [anomaly_temperature, anomaly_precip] = get_anomaly(target_lat, target_lon, jj)
            
            target_lat = tile.RUN_INFO.STATVAR.latitude;
            target_lon = tile.RUN_INFO.STATVAR.longitude;
            
            jj = (tile.t - datenum(str2num(datestr(tile.t, 'yyyy')), 1,1)) ./ tile.timestep + 1;
            
            load([forcing.PARA.anomaly_file_path forcing.PARA.anomaly_file], 'anomaly_airT_interp')
            anomaly_airT_interp=squeeze(double(anomaly_airT_interp(:,:,jj)));
            load([forcing.PARA.anomaly_file_path forcing.PARA.anomaly_file], 'anomaly_precip_interp')
            anomaly_precip_interp=squeeze(double(anomaly_precip_interp(:,:,jj)));
            load([forcing.PARA.anomaly_file_path forcing.PARA.anomaly_file], 'anomaly_LW_interp')
            anomaly_LW_interp=squeeze(double(anomaly_LW_interp(:,:,jj)));
            load([forcing.PARA.anomaly_file_path forcing.PARA.anomaly_file], 'anomaly_SW_interp')
            anomaly_SW_interp=squeeze(double(anomaly_SW_interp(:,:,jj)));
            load([forcing.PARA.anomaly_file_path forcing.PARA.anomaly_file], 'ERA_lat')
            load([forcing.PARA.anomaly_file_path forcing.PARA.anomaly_file], 'ERA_lon')
            
            max_target_lat = min(round(max(target_lat) /2) *2 + 2, 90) ;  %ERA5 in 0.25 degree resolution
            min_target_lat= max(round(min(target_lat)/2) *2 - 2, -90);
            
            
            target_lon(target_lon < 0) = target_lon(target_lon < 0) + 360; %change from -180 to +180 to 0 to 360
            max_target_lon = min(round(max(target_lon)/2) *2 + 2, 360);
            min_target_lon = max(round(min(target_lon)/2) *2 - 2, 0);
            
            lat_index_start = find(ERA_lat(:,1) == max_target_lat);
            lat_index_end = find(ERA_lat(:,1) == min_target_lat);
            lon_index_start = find(ERA_lon(:,1) == min_target_lon);
            lon_index_end = find(ERA_lon(:,1) == max_target_lon);
            
            
            [ERA_lon_mesh, ERA_lat_mesh] = meshgrid(ERA_lon(lon_index_start:lon_index_end), ERA_lat(lat_index_start:lat_index_end));
            
            
            anomaly_temperature = interp2(ERA_lon_mesh ,ERA_lat_mesh, anomaly_airT_interp(lon_index_start:lon_index_end, lat_index_start:lat_index_end)', target_lon, target_lat, 'cubic');
            anomaly_precip = interp2(ERA_lon_mesh ,ERA_lat_mesh, anomaly_precip_interp(lon_index_start:lon_index_end, lat_index_start:lat_index_end)', target_lon, target_lat, 'cubic');
            anomaly_LW = interp2(ERA_lon_mesh ,ERA_lat_mesh, anomaly_LW_interp(lon_index_start:lon_index_end, lat_index_start:lat_index_end)', target_lon, target_lat, 'cubic');
            anomaly_SW = interp2(ERA_lon_mesh ,ERA_lat_mesh, anomaly_SW_interp(lon_index_start:lon_index_end, lat_index_start:lat_index_end)', target_lon, target_lat, 'cubic');
            

            forcing.TEMP.ERA_T_downscaled = forcing.TEMP.ERA_T_downscaled + anomaly_temperature;
            forcing.TEMP.ERA_precip_downcaled = forcing.TEMP.ERA_precip_downcaled .* anomaly_precip;
            
            forcing.TEMP.ERA_Lin_downscaled = forcing.TEMP.ERA_Lin_downscaled + anomaly_LW;
            forcing.TEMP.ERA_Sin_downscaled = forcing.TEMP.ERA_Sin_downscaled .* anomaly_SW;
        end
        
        
        
        
        
    end
end