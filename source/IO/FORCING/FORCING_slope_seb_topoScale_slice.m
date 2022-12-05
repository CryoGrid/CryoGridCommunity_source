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
            
            forcing = reduce_precip_slope(forcing, tile);
            
            forcing = terrain_corr_Sin_dif(forcing, tile);
            forcing = reproject_Sin_dir(forcing, tile);
            forcing = terrain_shade(forcing, tile);
            forcing = terrain_corr_Lin(forcing, tile);
            forcing.DATA.Sin = forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
            
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
                
                forcing = terrain_corr_Sin_dif(forcing, tile);
                forcing = reproject_Sin_dir(forcing, tile);
                forcing = terrain_shade(forcing, tile);
                forcing = terrain_corr_Lin(forcing, tile);
                forcing.DATA.Sin = forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
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
        
    end
    
end