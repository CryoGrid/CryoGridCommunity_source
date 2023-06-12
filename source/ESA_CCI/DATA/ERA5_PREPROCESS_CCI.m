%========================================================================
% CryoGrid DATA_PROVIDER class DEM
% DEM class deriving information from a digital 
% elevation model (DEM),including slope, aspect (terrain shading still 
% lacking). The class can ingest Copernicus 30m DEMs downloaded from 
% https://portal.opentopography.org/   
%
% S. Westermann, Dec 2022
%========================================================================

classdef ERA5_PREPROCESS_CCI < FORCING_seb_CCI_full

    properties
        PARENT
    end
    
    methods
        function era = provide_PARA(era)    
            era.PARA.ERA5_path = [];
            era.PARA.target_path = [];
            era.PARA.add_run_name_folder = [];
            
            era.PARA.filename_SURF_airT = [];
            era.PARA.filename_SURF_geopotential = [];
            era.PARA.filename_SURF_dewpointT = [];
            era.PARA.filename_SURF_solar_radiation = [];
            era.PARA.filename_SURF_thermal_radiation = [];
            era.PARA.filename_SURF_TOA_radiation = [];
            era.PARA.filename_SURF_precip = [];
            era.PARA.filename_PLEV_airT = [];
            era.PARA.filename_PLEV_geopotential = [];
            era.PARA.filename_PLEV_RH = [];
            
            era.PARA.start_year = [];
            era.PARA.end_year = [];
        end
        
        function era = provide_STATVAR(era)

        end
        
        function era = provide_CONST(era)
            
        end
        
        function era = finalize_init(era)
            %make a list of all jobs to use for parallelization
            %9 variables x (era.PARA.end_year - era.PARA.start_year + 1)
            variable_list = [1:9]';
            year_list = [era.PARA.start_year:era.PARA.end_year]';
            [variable_list, year_list] = meshgrid(variable_list, year_list);
            variable_list = variable_list(:);
            year_list = year_list(:);
            
            
            if era.PARA.add_run_name_folder
                if ~(exist([era.PARA.target_path era.PARENT.RUN_INFO.PARA.run_name])==7)
                    mkdir([era.PARA.target_path era.PARENT.RUN_INFO.PARA.run_name])
                end
                era.PARA.target_path = [era.PARA.target_path era.PARENT.RUN_INFO.PARA.run_name '/'];
            end
            test_year = era.PARA.start_year;
            filename = era.PARA.filename_SURF_airT;
            filename(21:24) = num2str(test_year);
            
            ERA_lat = ncread([era.PARA.ERA5_path filename], 'latitude');
            ERA_lon = ncread([era.PARA.ERA5_path filename], 'longitude');
            ERA_lon =[ERA_lon; single(360)]; %add the 0-degree-longitude to the end as 360-degree-longitude
            
            max_target_lat = min(round(max(era.PARENT.STATVAR.latitude) *4) /4 + 0.25, 90) ;  %ERA5 in 0.25 degree resolution
            min_target_lat= max(round(min(era.PARENT.STATVAR.latitude)*4) /4 - 0.25, -90);
            
            target_lon = era.PARENT.STATVAR.longitude;
            target_lon(target_lon < 0) = target_lon(target_lon < 0) + 360; %change from -180 to +180 to 0 to 360
            max_target_lon = min(round(max(target_lon)*4) /4 + 0.25, 360);
            min_target_lon = max(round(min(target_lon)*4) / 4 - 0.25, 0);
            
            lat_index_start = find(ERA_lat(:,1) == max_target_lat);
            lat_index_end = find(ERA_lat(:,1) == min_target_lat);
            lon_index_start = find(ERA_lon(:,1) == min_target_lon);
            lon_index_end = find(ERA_lon(:,1) == max_target_lon);
            
            era.TEMP.lat_index_start = find(ERA_lat(:,1) == max_target_lat);
            era.TEMP.lat_index_end = find(ERA_lat(:,1) == min_target_lat);
            era.TEMP.lon_index_start = find(ERA_lon(:,1) == min_target_lon);
            era.TEMP.lon_index_end = find(ERA_lon(:,1) == max_target_lon);
            
            era.TEMP.max_target_lon = max_target_lon;
            era.TEMP.target_lon = target_lon;
            
            if era.PARENT.RUN_INFO.PARA.worker_number == 1
                %geopotential
                filename = era.PARA.filename_SURF_geopotential;
                filename(19:22) = num2str(test_year);
                
                z_surface = ncread_with_360_degree(era, [era.PARA.ERA5_path filename], 'z', 1, 1);
                
                nccreate([era.PARA.target_path filename], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                ncwrite([era.PARA.target_path filename], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                nccreate([era.PARA.target_path filename], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                ncwrite([era.PARA.target_path filename], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                nccreate([era.PARA.target_path filename], 'z', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(z_surface,1),'Y',size(z_surface,2)})
                ncwrite([era.PARA.target_path filename], 'z', single(z_surface))
            end
            
            for i = era.PARENT.RUN_INFO.PARA.worker_number:era.PARENT.RUN_INFO.PARA.number_of_cores:size(year_list,1)
                current_year = year_list(i,1);
                
                filename_SURF_precip = era.PARA.filename_SURF_precip;
                filename_SURF_precip(26:29) = num2str(current_year);
                
                filename_SURF_airT = era.PARA.filename_SURF_airT;
                filename_SURF_airT(21:24) = num2str(current_year);
                
                filename_SURF_dewpointT = era.PARA.filename_SURF_dewpointT;
                filename_SURF_dewpointT(30:33) = num2str(current_year);
                
                filename_SURF_solar_radiation = era.PARA.filename_SURF_solar_radiation;
                filename_SURF_solar_radiation(40:43) = num2str(current_year);
                
                filename_SURF_thermal_radiation = era.PARA.filename_SURF_thermal_radiation;
                filename_SURF_thermal_radiation(42:45) = num2str(current_year);
                
                filename_SURF_TOA_radiation = era.PARA.filename_SURF_TOA_radiation;
                filename_SURF_TOA_radiation(35:38) = num2str(current_year);
                
                filename_PLEV_airT = era.PARA.filename_PLEV_airT;
                filename_PLEV_airT(18:21) = num2str(current_year);
                
                filename_PLEV_geopotential = era.PARA.filename_PLEV_geopotential;
                filename_PLEV_geopotential(19:22) = num2str(current_year);
                
                filename_PLEV_RH = era.PARA.filename_PLEV_RH;
                filename_PLEV_RH(24:27) = num2str(current_year);
                
                if variable_list(i,1) == 1
                    t2m = ncread_with_360_degree(era, [era.PARA.ERA5_path filename_SURF_airT], 't2m', 1, Inf);
                    nccreate([era.PARA.target_path filename_SURF_airT], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_airT], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_airT], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_airT], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_airT], 't2m', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(t2m,1),'Y',size(t2m,2), 'time', size(t2m,3)})
                    ncwrite([era.PARA.target_path filename_SURF_airT], 't2m', single(t2m))
                    t2m = [];
                end
                if variable_list(i,1) == 2
                    tp = ncread_with_360_degree(era, [era.PARA.ERA5_path filename_SURF_precip], 'tp', 1, Inf);
                    nccreate([era.PARA.target_path filename_SURF_precip], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_precip], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_precip], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_precip], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_precip], 'tp', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(tp,1),'Y',size(tp,2), 'time', size(tp,3)})
                    ncwrite([era.PARA.target_path filename_SURF_precip], 'tp', single(tp))
                    tp=[];
                end
                if variable_list(i,1) == 3
                    t = ncread_with_360_degree_pl(era, [era.PARA.ERA5_path filename_PLEV_airT], 't', 1, Inf);
                    nccreate([era.PARA.target_path filename_PLEV_airT], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                    ncwrite([era.PARA.target_path filename_PLEV_airT], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                    nccreate([era.PARA.target_path filename_PLEV_airT], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                    ncwrite([era.PARA.target_path filename_PLEV_airT], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                    nccreate([era.PARA.target_path filename_PLEV_airT], 't', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(t,1),'Y',size(t,2), 'pl', size(t,3), 'time', size(t,4)})
                    ncwrite([era.PARA.target_path filename_PLEV_airT], 't', single(t))
                    t = [];
                end
                if variable_list(i,1) == 4
                    z = ncread_with_360_degree_pl(era, [era.PARA.ERA5_path filename_PLEV_geopotential], 'z', 1, Inf);
                    nccreate([era.PARA.target_path filename_PLEV_geopotential], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                    ncwrite([era.PARA.target_path filename_PLEV_geopotential], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                    nccreate([era.PARA.target_path filename_PLEV_geopotential], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                    ncwrite([era.PARA.target_path filename_PLEV_geopotential], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                    nccreate([era.PARA.target_path filename_PLEV_geopotential], 'z', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(z,1),'Y',size(z,2), 'pl', size(z,3), 'time', size(z,4)})
                    ncwrite([era.PARA.target_path filename_PLEV_geopotential], 'z', single(z))
                    z=[];
                end
                if variable_list(i,1) == 5
                    d2m = ncread_with_360_degree(era, [era.PARA.ERA5_path filename_SURF_dewpointT], 'd2m', 1, Inf);
                    nccreate([era.PARA.target_path filename_SURF_dewpointT], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_dewpointT], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_dewpointT], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_dewpointT], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_dewpointT], 'd2m', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(d2m,1),'Y',size(d2m,2), 'time', size(d2m,3)})
                    ncwrite([era.PARA.target_path filename_SURF_dewpointT], 'd2m', single(d2m))
                    d2m = [];
                end
                if variable_list(i,1) == 6
                    ssrd = ncread_with_360_degree(era, [era.PARA.ERA5_path filename_SURF_solar_radiation], 'ssrd', 1, Inf);
                    nccreate([era.PARA.target_path filename_SURF_solar_radiation], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_solar_radiation], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_solar_radiation], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_solar_radiation], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_solar_radiation], 'ssrd', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(ssrd,1),'Y',size(ssrd,2), 'time', size(ssrd,3)})
                    ncwrite([era.PARA.target_path filename_SURF_solar_radiation], 'ssrd', single(ssrd))
                    ssrd = [];
                end
                if variable_list(i,1) == 7
                    strd = ncread_with_360_degree(era, [era.PARA.ERA5_path filename_SURF_thermal_radiation], 'strd', 1, Inf);
                    nccreate([era.PARA.target_path filename_SURF_thermal_radiation], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_thermal_radiation], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_thermal_radiation], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_thermal_radiation], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_thermal_radiation], 'strd', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(strd,1),'Y',size(strd,2), 'time', size(strd,3)})
                    ncwrite([era.PARA.target_path filename_SURF_thermal_radiation], 'strd', single(strd))
                    strd=[];
                end
                if variable_list(i,1) == 8
                    tisr = ncread_with_360_degree(era, [era.PARA.ERA5_path filename_SURF_TOA_radiation], 'tisr', 1, Inf);
                    nccreate([era.PARA.target_path filename_SURF_TOA_radiation], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_TOA_radiation], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_TOA_radiation], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                    ncwrite([era.PARA.target_path filename_SURF_TOA_radiation], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                    nccreate([era.PARA.target_path filename_SURF_TOA_radiation], 'tisr', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(tisr,1),'Y',size(tisr,2), 'time', size(tisr,3)})
                    ncwrite([era.PARA.target_path filename_SURF_TOA_radiation], 'tisr', single(tisr))
                    tisr = [];
                end
                if variable_list(i,1) == 9
                    r = ncread_with_360_degree_pl(era, [era.PARA.ERA5_path filename_PLEV_RH], 'r', 1, Inf);
                    nccreate([era.PARA.target_path filename_PLEV_RH], 'latitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'Y',lat_index_end-lat_index_start+1})
                    ncwrite([era.PARA.target_path filename_PLEV_RH], 'latitude', ERA_lat(lat_index_start:lat_index_end,1))
                    nccreate([era.PARA.target_path filename_PLEV_RH], 'longitude', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',lon_index_end-lon_index_start+1})
                    ncwrite([era.PARA.target_path filename_PLEV_RH], 'longitude', ERA_lon(lon_index_start:lon_index_end,1))
                    nccreate([era.PARA.target_path filename_PLEV_RH], 'r', 'Format', 'netcdf4', 'Datatype', 'single', 'Dimensions',{'X',size(r,1),'Y',size(r,2), 'pl', size(r,3), 'time', size(r,4)})
                    ncwrite([era.PARA.target_path filename_PLEV_RH], 'r', single(r))
                    r=[];
                end
            end
            if era.PARENT.RUN_INFO.PARA.parallelized
                NMPI_Barrier();
            end
        end

        
        
        function era = load_data(era)

            
        end
        

       
        
  
                
        %-----service functions-----------------
        
        
        function output = ncread_with_360_degree(era, fname, variable, time_index_start, time_number_of_values)
            
            if era.TEMP.max_target_lon < 360 %"normal" case, no need to change anything
                output = ncread(fname, variable, [ era.TEMP.lon_index_start era.TEMP.lat_index_start time_index_start], [era.TEMP.lon_index_end - era.TEMP.lon_index_start+1 era.TEMP.lat_index_end - era.TEMP.lat_index_start+1 time_number_of_values], [1 1 1]);
                
            else  %this part is not yet tested, it only becomes effective at the 0 degree longitude -> the only affected PF seems to be in the Pyrenees, which is fairly irrelevant
                era.TEMP.lon_index_end=era.TEMP.lon_index_end-1; %read to lon = 359.75 degree
                output = ncread(fname, variable, [ era.TEMP.lon_index_start era.TEMP.lat_index_start time_index_start], [era.TEMP.lon_index_end - era.TEMP.lon_index_start+1 era.TEMP.lat_index_end - era.TEMP.lat_index_start+1 time_number_of_values], [1 1 1]);
                output0degree = ncread(fname, variable, [ 1 era.TEMP.lat_index_start time_index_start], [1 era.TEMP.lat_index_end - era.TEMP.lat_index_start+1 time_number_of_values], [1 1 1]); %read first "line" of longitude = 0 degree
                output = cat(1, output, output0degree); %concat in dimension1 (longitude)
		era.TEMP.lon_index_end=era.TEMP.lon_index_end+1;
            end
            
        end
        
        function output = ncread_with_360_degree_pl(era, fname, variable, time_index_start, time_number_of_values)
            
            %pl = 1 1000; pl = 2 700; pl = 3 500; pl = 4 200
            
            if era.TEMP.max_target_lon < 360 %"normal" case, no need to change anything
                output = ncread(fname, variable, [ era.TEMP.lon_index_start era.TEMP.lat_index_start 1 time_index_start], [era.TEMP.lon_index_end - era.TEMP.lon_index_start+1 era.TEMP.lat_index_end - era.TEMP.lat_index_start+1 Inf time_number_of_values], [1 1 1 1]);
            else  %this part is not yet tested, it only becomes effective at the 0 degree longitude -> the only affected PF seems to be in the Pyrenees, which is fairly irrelevant
                era.TEMP.lon_index_end=era.TEMP.lon_index_end-1; %read to lon = 359.75 degree
                output = ncread(fname, variable, [ era.TEMP.lon_index_start era.TEMP.lat_index_start 1 time_index_start], [era.TEMP.lon_index_end - era.TEMP.lon_index_start+1 era.TEMP.lat_index_end - era.TEMP.lat_index_start+1 Inf time_number_of_values], [1 1 1 1]);
                output0degree = ncread(fname, variable, [ 1 era.TEMP.lat_index_start 1 time_index_start], [1 era.TEMP.lat_index_end - era.TEMP.lat_index_start+1 Inf time_number_of_values], [1 1 1 1]); %read first "line" of longitude
                output = cat(1, output, output0degree); %concat in dimension1 (longitude)
		era.TEMP.lon_index_end=era.TEMP.lon_index_end+1;
            end
        end
        
    end
end

