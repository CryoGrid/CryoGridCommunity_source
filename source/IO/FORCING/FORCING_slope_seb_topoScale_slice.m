%========================================================================
% CryoGrid FORCING class FORCING_slope_seb_topoScale_slice
%
% model forcing for GROUND classes computing the surface energy 
% balance (keyword "seb"). 
% 
% The data is obtained by interpolating surface level and pressure level 
% data (generally from ERA reanalysis) to the coordinate of the target 
% location and applying the TopoScale downscaling routine (Fiddes et al., 
% 2014). The input data are nc-files containing all variables with 
% either monthy, quarterly or yearly data per file. The files can be 
% downloaded by the the python scripts "request_sl.py" and "request_gp.py".
% 
% The generated forcing variables are:
%
% Tair:      Air temperature (in degree Celsius)
% Lin:       incoming long-wave radiation (in W/m2)
% Sin:       incoming short-wave radiation (in W/m2)
% Sin_dir:   direct incoming short-wave radiation (in W/m2)
% Sin_dif:   diffuse incoming short-wave radiation (in W/m2)
% rainfall:  Rainfall (in mm/day)
% snowfall:  Snowfall (in mm/day)
% q:         absolute humidity (in kg water vapor / kg air)
% p:         air pressure (OPTIONAL, in hPa)
% wind:       wind speed (in m/sec)
% 
% The forcing data are read sequentially, thus keeping only a small chunk 
% of forcing data in memory. This class is recommended for large input data 
% sets, in particular if hourly timestamps are used. 
%
%
% Authors:
% K. Aalstad, S. Westermann, December 2022
%
%========================================================================

classdef FORCING_slope_seb_topoScale_slice < FORCING_base & READ_FORCING_NC & PROCESS_FORCING_topoScale
    
    properties
        
    end
    
    methods
        function forcing = provide_PARA(forcing)
            
            forcing.PARA.nc_folder = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_path = [];
            forcing.PARA.time_resolution_input = [];
            forcing.PARA.top_pressure_level = [];
            forcing.PARA.bottom_pressure_level = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.all_rain_T = [];     % Temperature above which all precipitation is considered as rain
            forcing.PARA.all_snow_T = [];     % Temperature below which all precipitation is considered as snow
            forcing.PARA.albedo_surrounding_terrain = [];
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
            
            forcing.PARA.post_proc_class = [];  %optional post-processing classes
            forcing.PARA.post_proc_class_index = [];
        end
        
        
        function forcing = provide_CONST(forcing)
            forcing.CONST.Tmfw = [];
            forcing.CONST.sigma = [];
        end
        
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        
        function forcing = finalize_init(forcing, tile)
            
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));
            forcing.TEMP.current_month = month(forcing.PARA.start_time);
            forcing.TEMP.current_quarter = floor((forcing.TEMP.current_month-1)./3)+1;
            forcing.TEMP.current_year = year(forcing.PARA.start_time);
            
            forcing = load_ERA_slice(forcing);

            forcing = process_topoScale(forcing, tile);
            forcing.TEMP.era = []; 
            
            forcing = split_precip_Tair(forcing); % distinguish snow-/rainfall

            forcing.DATA.rainfall =  forcing.DATA.rainfall.*forcing.PARA.rain_fraction;
            forcing.DATA.snowfall =  forcing.DATA.snowfall.*forcing.PARA.snow_fraction;
            
            forcing = check_and_correct(forcing); % Remove known errors
%             forcing = set_start_and_end_time(forcing); % assign start/end time
            forcing = initialize_TEMP(forcing);
            forcing = initialize_TEMP_slope(forcing);

            forcing = reduce_precip_slope(forcing, tile);
            
            forcing = convert_accumulated2instantaneous_Sin(forcing, tile);

            forcing = SolarAzEl(forcing, tile);
            forcing = terrain_corr_Sin_dif(forcing, tile);
            forcing = reproject_Sin_dir(forcing, tile);
            forcing = terrain_shade(forcing, tile);
            forcing = terrain_corr_Lin(forcing, tile);
            forcing.DATA.Sin = forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
            
            %optional post-processing with dedicated classes
            if ~isempty(forcing.PARA.post_proc_class) && sum(isnan(forcing.PARA.post_proc_class_index)) == 0
                for i=1:size(forcing.PARA.post_proc_class,1)
                    post_proc_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.post_proc_class{i,1}){forcing.PARA.post_proc_class_index(i,1),1});
                    post_proc_class = finalize_init(post_proc_class, tile);
                    forcing = post_process(post_proc_class, forcing, tile);
                end
            end
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing = interpolate_forcing@FORCING_base(forcing, tile);
            
            forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
            forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
            
            %load new chunk
            if tile.t > forcing.DATA.timeForcing(end,1)-2 && forcing.DATA.timeForcing(end,1) < forcing.PARA.end_time
                
                old_slice = forcing.DATA;
                
                if strcmp(forcing.PARA.time_resolution_input, 'quarter')
                    forcing.TEMP.current_quarter = forcing.TEMP.current_quarter + 1;
                    if forcing.TEMP.current_quarter == 5
                        forcing.TEMP.current_quarter =1;
                        forcing.TEMP.current_year = forcing.TEMP.current_year +1;
                    end
                elseif strcmp(forcing.PARA.time_resolution_input, 'month')
                    forcing.TEMP.current_month = forcing.TEMP.current_month + 1;
                    if forcing.TEMP.current_month == 13
                        forcing.TEMP.current_month =1;
                        forcing.TEMP.current_year = forcing.TEMP.current_year +1;
                    end
                elseif strcmp(forcing.PARA.time_resolution_input, 'year')
                        forcing.TEMP.current_year = forcing.TEMP.current_year +1;
                end
                
                forcing = load_ERA_slice(forcing); %overwrites forcing.DATA
                forcing = process_topoScale(forcing, tile);
                forcing.TEMP.era = [];
                
                forcing = split_precip_Tair(forcing); % distinguish snow-/rainfall
                
                forcing.DATA.rainfall =  forcing.DATA.rainfall.*forcing.PARA.rain_fraction;
                forcing.DATA.snowfall =  forcing.DATA.snowfall.*forcing.PARA.snow_fraction;
                
                forcing = check_and_correct(forcing); % Remove known errors
                forcing = reduce_precip_slope(forcing, tile);
                
                forcing = convert_accumulated2instantaneous_Sin(forcing, tile);

                forcing = SolarAzEl(forcing, tile);
                forcing = terrain_corr_Sin_dif(forcing, tile);
                forcing = reproject_Sin_dir(forcing, tile);
                forcing = terrain_shade(forcing, tile);
                forcing = terrain_corr_Lin(forcing, tile);
                forcing.DATA.Sin = forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
                
                %optional post-processing with dedicated classes
                if ~isempty(forcing.PARA.post_proc_class) && sum(isnan(forcing.PARA.post_proc_class_index)) == 0
                    for i=1:size(forcing.PARA.post_proc_class,1)
                        post_proc_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.post_proc_class{i,1}){forcing.PARA.post_proc_class_index(i,1),1});
                        post_proc_class = finalize_init(post_proc_class, tile);
                        forcing = post_process(post_proc_class, forcing, tile);
                    end
                end
                %overwrite with new data completed
                
                %cut away the old chunk except for the last 4 days (2 days
                %buffer)
                variables = fieldnames(old_slice);
                offset_index = min(size(old_slice.timeForcing,1)+1, floor(4./(old_slice.timeForcing(2)-old_slice.timeForcing(1))));
                for i=1:size(variables,1)
                    data_chunk = old_slice.(variables{i,1});
                    forcing.DATA.(variables{i,1}) = [data_chunk(end-offset_index:end,1); forcing.DATA.(variables{i,1})];
                end
            end
            
        end
        
        %-------------param file generation-----
        function forcing = param_file_info(forcing)
            forcing = provide_PARA(forcing);

            forcing.PARA.STATVAR = [];
            forcing.PARA.class_category = 'FORCING';
            
            forcing.PARA.comment.nc_folder = {'folder containing nc-files with forcing data'};
            
            forcing.PARA.default_value.forcing_path = {'../CryoGridCommunity_forcing/'};
            forcing.PARA.comment.forcing_path = {'path where forcing data file is located'};
            
            forcing.PARA.default_value.time_resolution_input = {'month'};
            forcing.PARA.comment.time_resolution_input = {'time resolution of nc-files; three options: "month", "quarter" or "year"'};
            
            forcing.PARA.default_value.top_pressure_level = {800};
            forcing.PARA.comment.top_pressure_level = {'highest pressure level that is read (in hPa) - should be above the highest target location'};
            
            forcing.PARA.default_value.bottom_pressure_level = {1000};
            forcing.PARA.comment.bottom_pressure_level = {'lowest pressure level that is read (in hPa), surface level data are assigned below'};
            
            forcing.PARA.comment.start_time = {'start time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.start_time.name =  'H_LIST';
            forcing.PARA.options.start_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.comment.end_time = {'end_time time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.end_time.name =  'H_LIST'; % 
            forcing.PARA.options.end_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.default_value.rain_fraction = {1};  
            forcing.PARA.comment.rain_fraction = {'rainfall fraction assumed in simulations (rainfall from the forcing data file is multiplied by this parameter)'};
            
            forcing.PARA.default_value.snow_fraction = {1};  
            forcing.PARA.comment.snow_fraction = {'snowfall fraction assumed in simulations (rainfall from the forcing data file is multiplied by this parameter)'};

            forcing.PARA.default_value.all_rain_T = {0.5};    
            forcing.PARA.comment.all_rain_T = {'Temperature above which all precipitation is considered as rain'};

            forcing.PARA.default_value.all_snow_T =  {-0.5};    
            forcing.PARA.comment.all_snow_T = {'Temperature below which all precipitation is considered as snow'};
             
            forcing.PARA.default_value.albedo_surrounding_terrain = {0.2};
            forcing.PARA.comment.albedo_surrounding_terrain = {'albedo of terrain in the field of view of the target location, reflecting short-wave radiation; considered static throughout the year'}; 
            
            forcing.PARA.default_value.heatFlux_lb = {0.05};
            forcing.PARA.comment.heatFlux_lb = {'heat flux at the lower boundary [W/m2] - positive values correspond to energy gain'};
            
            forcing.PARA.default_value.airT_height = {2};  
            forcing.PARA.comment.airT_height = {'height above ground surface where air temperature from forcing data is applied'};
            
            forcing.PARA.comment.post_proc_class = {'list of postprocessing classes to modify forcing data in user-defined ways; no post-processing applied when empty'};
            forcing.PARA.options.post_proc_class.name = 'H_LIST';
            
            forcing.PARA.comment.post_proc_class_index = {''};
            forcing.PARA.options.post_proc_class_index.name = 'H_LIST';
        end
        
    end
    
end