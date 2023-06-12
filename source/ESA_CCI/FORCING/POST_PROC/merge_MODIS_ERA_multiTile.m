%========================================================================
% CryoGrid FORCING post-processing 
%
%
% Authors:
% S. Westermann, January 2023
%
%========================================================================

classdef merge_MODIS_ERA_multiTile < merge_MODIS_ERA 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.MODIS_path = [];
            post_proc.PARA.deg_tile_list_file = [];
            post_proc.PARA.deg_tile_list_folder = [];
            
            post_proc.PARA.threshold_diff_MODIS_ERA = [];
            
            post_proc.PARA.averaging_period = [];  
            post_proc.PARA.year_range = [];
            
            post_proc.PARA.annual = [];

        end
        
        
        function post_proc = provide_CONST(post_proc)

        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)
            load([post_proc.PARA.deg_tile_list_folder post_proc.PARA.deg_tile_list_file], 'MODIS_deg_list');
            post_proc.TEMP.MODIS_deg_list = MODIS_deg_list;
            post_proc.TEMP.offset_years = 0;
        end
        
        
        %         function forcing = post_process(post_proc, forcing, tile)
        %
        %             if forcing.TEMP.current_year >= post_proc.PARA.year_range(1) && forcing.TEMP.current_year <= post_proc.PARA.year_range(end)
        %                 pos_year = forcing.TEMP.current_year - (forcing.PARA.ERA_data_years(1)- forcing.PARA.number_of_spin_up_years) + 1 - post_proc.TEMP.offset_years;
        %                 post_proc = load_MODIS(post_proc, forcing, tile);
        %                 forcing = merge_MODIS_ERA_and_average(post_proc, forcing, tile);
        %             end
        %         end
        
        function post_proc = load_MODIS(post_proc, forcing, tile)
            
            post_proc.TEMP.MODIS_LST = zeros(size(tile.PARA.range,1), 4.*(datenum(forcing.TEMP.current_year+1,1,1)-datenum(forcing.TEMP.current_year,1,1))) .*NaN; % [MOD_LST_Day_1km MOD_LST_Night_1km MYD_LST_Day_1km MYD_LST_Night_1km];
            post_proc.TEMP.MODIS_LST_time = zeros(size(tile.PARA.range,1), 4*(datenum(forcing.TEMP.current_year+1,1,1)-datenum(forcing.TEMP.current_year,1,1))) .* NaN; %[MOD_Day_view_time MOD_Night_view_time MYD_Day_view_time MYD_Night_view_time];
            
            for i=1:size(post_proc.TEMP.MODIS_deg_list,1)
                
                range_in_tile = find(tile.PARA.latitude(:,1) > post_proc.TEMP.MODIS_deg_list(i,1) & ...
                    tile.PARA.latitude(:,1) < post_proc.TEMP.MODIS_deg_list(i,2) & tile.PARA.longitude(:,1) > post_proc.TEMP.MODIS_deg_list(i,3) & ...
                    tile.PARA.longitude(:,1) < post_proc.TEMP.MODIS_deg_list(i,4));
                if ~isempty(range_in_tile)
                    MODIS_file = ['MODIS_' num2str(post_proc.TEMP.MODIS_deg_list(i,1)) '_' num2str(post_proc.TEMP.MODIS_deg_list(i,2)) '_' ...
                        num2str(post_proc.TEMP.MODIS_deg_list(i,3)) '_' num2str(post_proc.TEMP.MODIS_deg_list(i,4))];
                    fname = [post_proc.PARA.MODIS_path MODIS_file '_' num2str(forcing.TEMP.current_year) '.nc'];
                    MODIS_lat = ncread(fname, 'latitude', 1, Inf);
                    MODIS_lon = ncread(fname, 'longitude', 1, Inf);
                    MODIS_delta_lat = ncread(fname, 'delta_lat', 1, Inf);
                    MODIS_delta_lon = ncread(fname, 'delta_lon', 1, Inf);
                    for j=1:size(range_in_tile,1)
                        
                        %find closest pixel in MODIS file
                        [dist, closest_pixel_id] = min((tile.PARA.longitude(range_in_tile(j),1) - MODIS_lon).^2 + (tile.PARA.latitude(range_in_tile(j),1) - MODIS_lat).^2);
                        if 1 %dist < (MODIS_delta_lat).^2 + (MODIS_delta_lon).^2 %pixel actually available in MODIS file
                            
                            MOD_LST_Day_1km = single(ncread(fname, 'MOD_LST_Day_1km', [closest_pixel_id 1], [1 Inf])).*0.02;
                            MOD_LST_Night_1km = single(ncread(fname, 'MOD_LST_Night_1km', [closest_pixel_id 1], [1 Inf])).*0.02;
                            MYD_LST_Day_1km = single(ncread(fname, 'MYD_LST_Day_1km', [closest_pixel_id 1], [1 Inf])).*0.02;
                            MYD_LST_Night_1km = single(ncread(fname, 'MYD_LST_Night_1km', [closest_pixel_id(1) 1], [1 Inf])).*0.02;
                            
                            MOD_Day_view_time = double(ncread(fname, 'MOD_Day_view_time', [closest_pixel_id 1], [1 Inf]));
                            MOD_Day_view_time(MOD_Day_view_time==255) = NaN;
                            MOD_Night_view_time = double(ncread(fname, 'MOD_Night_view_time', [closest_pixel_id 1], [1 Inf]));
                            MOD_Night_view_time(MOD_Night_view_time == 255) = NaN;
                            MYD_Day_view_time = double(ncread(fname, 'MYD_Day_view_time', [closest_pixel_id 1], [1 Inf]));
                            MYD_Day_view_time(MYD_Day_view_time==255) = NaN;
                            MYD_Night_view_time = double(ncread(fname, 'MYD_Night_view_time', [closest_pixel_id 1], [1 Inf]));
                            MYD_Night_view_time(MYD_Night_view_time==255) = NaN;
                            
                            doy = floor([1:size(MOD_LST_Day_1km,2)]); %1:365/366
                            daily_timestamp = datenum(forcing.TEMP.current_year,1,1) + floor([1:size(MOD_LST_Day_1km,2)])-1;
                            
                            MOD_Day_view_time = MOD_Day_view_time.*0.1./24 + daily_timestamp;
                            MOD_Night_view_time = MOD_Night_view_time.*0.1./24 + daily_timestamp;
                            MYD_Day_view_time = MYD_Day_view_time.*0.1./24 + daily_timestamp;
                            MYD_Night_view_time = MYD_Night_view_time.*0.1./24 + daily_timestamp;
                            
                            MOD_Day_view_time = SolarTime2LocalTime(post_proc, MOD_Day_view_time, doy, tile.PARA.longitude(range_in_tile(j),:), 0);
                            MOD_Night_view_time = SolarTime2LocalTime(post_proc, MOD_Night_view_time, doy, tile.PARA.longitude(range_in_tile(j),:), 0);
                            MYD_Day_view_time = SolarTime2LocalTime(post_proc, MYD_Day_view_time, doy, tile.PARA.longitude(range_in_tile(j),:), 0);
                            MYD_Night_view_time = SolarTime2LocalTime(post_proc, MYD_Night_view_time, doy, tile.PARA.longitude(range_in_tile(j),:), 0);
                            
                            post_proc.TEMP.MODIS_LST(range_in_tile(j),:) = [MOD_LST_Day_1km MOD_LST_Night_1km MYD_LST_Day_1km MYD_LST_Night_1km];
                            post_proc.TEMP.MODIS_LST_time(range_in_tile(j),:) = [MOD_Day_view_time MOD_Night_view_time MYD_Day_view_time MYD_Night_view_time];
                        end
                    end
                end
            end
        end
        
%         function forcing = merge_MODIS_ERA_and_average(post_proc, forcing, tile)  %merge MODIS and ERA
%             
%             %             compute_slice = 10000;
%             %             preproc.STATVAR.final_av_T = preproc.STATVAR.ERA_T_downscaled.*0;
%             %             preproc.STATVAR.final_MODIS_weight = preproc.STATVAR.ERA_T_downscaled.*0;
%             %
%             %             for count = 1:compute_slice:size(post_proc.MODIS_CLASS.STATVAR.timestamp,1)
%             %                 range = [count:min(count+compute_slice-1, size(post_proc.MODIS_CLASS.STATVAR.timestamp,1))];
%             
%             %1st dim is the different locations , second dimension is
%             %time
%             MODIS_LST_time = post_proc.TEMP.MODIS_LST_time;
%             MODIS_LST = post_proc.TEMP.MODIS_LST;
%             
%             %append one value at both sides (needed so that LST
%             %interpolation works at the very start and end of time series works)
%             ERA_time = [forcing.TEMP.ERA_time(1,1)-0.25 forcing.TEMP.ERA_time forcing.TEMP.ERA_time(1,end)+0.25];
%             ERA_T2m_downscaled = [forcing.TEMP.ERA_T_downscaled(:,4) forcing.TEMP.ERA_T_downscaled forcing.TEMP.ERA_T_downscaled(:,end-3)];
%             
%             %make the 4 dependent on
%             %tile.FORCING.PARA.number_of_values_per_day
%             
%             %upper lower is the index of the ERA value before and after
%             %the MODIS value
%             lower = (floor(MODIS_LST_time.*4)./4 - ERA_time(1,1)).*4 + 1;
%             upper= (floor(MODIS_LST_time.*4)./4 + 0.25-ERA_time(1,1)).*4 + 1;
%             
%             %weights are the weights for replacing the ERA value to the
%             %left and right
%             weight_left = 1 - (MODIS_LST_time - floor(MODIS_LST_time.*4)./4 ) ./0.25;
%             weight_right= 1 - weight_left;
%             
%             lower(isnan(lower)) = 1;
%             upper(isnan(upper)) = 2;
%             
%             %new
%             %                 weight_left(upper<=1) = 0;
%             %                 weight_right(upper<=1) = 1;
%             %                 lower(upper<=1) = 1;
%             %                 upper(upper<=1) = 1;
%             %                 weight_left(lower >= size(ERA_time,2)) = 1;
%             %                 weight_right(lower >= size(ERA_time,2)) = 0;
%             %                 lower(lower >= size(ERA_time,2)) = size(ERA_time,2);
%             %                 upper(lower >= size(ERA_time,2)) = size(ERA_time,2);
%             %end new
%             
%             weight_left(lower<1 | upper>size(ERA_time,2)) = NaN;
%             weight_right(lower<1 | upper>size(ERA_time,2)) = NaN;
%             
%             %MODIS_LST(lower<1)=0;
%             upper(lower<1)=2;
%             lower(lower<1)=1;
%             %MODIS_LST(upper>size(ERA_time,2))=0;
%             lower(upper>size(ERA_time,2)) = size(ERA_time,2) - 1;
%             upper(upper>size(ERA_time,2)) = size(ERA_time,2);
%             
%             
%             
%             MODIS_LST_interp=MODIS_LST.*0;  %interpolate ERA to MODIS LST timestamp
%             
%             matrix = repmat([1:size(ERA_T2m_downscaled,2)], size(ERA_T2m_downscaled,1),1);
%             for i=1:size(MODIS_LST,2)
%                 %MODIS_LST_interp(:,i) = weight_left(:,i).*diag(ERA_T2m_downscaled(:,lower(:,i))) + weight_right(:,i) .* diag(ERA_T2m_downscaled(:,upper(:,i)));
%                 %faster version
%                 upper_matrix = double(matrix == repmat(upper(:,i), 1, size(ERA_T2m_downscaled,2)));
%                 lower_matrix = double(matrix == repmat(lower(:,i), 1, size(ERA_T2m_downscaled,2)));
%                 
%                 MODIS_LST_interp(:,i) = weight_left(:,i) .* sum(ERA_T2m_downscaled .* lower_matrix,2) + weight_right(:,i) .* sum(ERA_T2m_downscaled .* upper_matrix,2);
%                 
%             end
%             MODIS_LST_interp(isnan(MODIS_LST_interp))=0;
%             
%             
%             %post_proc.PARA.threshold_diff = 10; %MODIS LST deleted if more than 10K colder
%             %weight_left(MODIS_LST<MODIS_LST_interp-10) = NaN;
%             %weight_right(MODIS_LST<MODIS_LST_interp-10) = NaN;
%             %MODIS_LST(MODIS_LST<MODIS_LST_interp-threshold_diff) = 0;
%             %             post_proc.PARA.threshold_diff = 10;
%             
%             %MODIS LST deleted if more than 10K colder
%             weight_left(MODIS_LST<MODIS_LST_interp - post_proc.PARA.threshold_diff_MODIS_ERA) = NaN;
%             weight_right(MODIS_LST<MODIS_LST_interp - post_proc.PARA.threshold_diff_MODIS_ERA) = NaN;
%             MODIS_LST(MODIS_LST<MODIS_LST_interp - post_proc.PARA.threshold_diff_MODIS_ERA) = 0;
%             
%             %MODIS_LST_mean_reduced = sum(MODIS_LST,2) ./sum(double(MODIS_LST~=0),2)-273.15;
%             
%             weight_left(isnan(weight_left))=0;
%             weight_right(isnan(weight_right))=0;
%             
%             %get the total weight(=sum of all weights from surrounding
%             %MODIS LST values) for each ERA value
%             total_weights=ERA_T2m_downscaled.*0;
%             for j=1:size(MODIS_LST,2)
%                 for i=1:size(MODIS_LST,1)
%                     total_weights(i, lower(i,j)) = total_weights(i, lower(i,j)) + weight_left(i,j);
%                     total_weights(i, upper(i,j)) = total_weights(i, upper(i,j)) + weight_right(i,j);
%                 end
%             end
%             
%             %if weight is > 1, the ERA weight is zero anyway
%             total_weights_reduced = total_weights ./ max(1,total_weights);
%             ERA_weights = 1-total_weights_reduced;
%             %can be taken away?
%             ERA_weights(:,end)=0;
%             ERA_weights(:,1) = 0;
%             
%             reduction_of_MODIS_weights = total_weights_reduced ./ max(1e-20, total_weights);
%             
%             MODIS_LST_weight = MODIS_LST.*0;
%             
%             matrix = repmat([1:size(reduction_of_MODIS_weights,2)], size(reduction_of_MODIS_weights,1),1);
%             for i=1:size(MODIS_LST,2)
%                 %MODIS_LST_weight(:,i) = weight_left(:,i) .* diag(reduction_of_MODIS_weights(:,lower(:,i))) + weight_right(:,i) .* diag(reduction_of_MODIS_weights(:,upper(:,i)));
%                 %faster version
%                 upper_matrix = double(matrix == repmat(upper(:,i), 1, size(reduction_of_MODIS_weights,2)));
%                 lower_matrix = double(matrix == repmat(lower(:,i), 1, size(reduction_of_MODIS_weights,2)));
%                 MODIS_LST_weight(:,i) = weight_left(:,i) .* sum(reduction_of_MODIS_weights .* lower_matrix,2) + weight_right(:,i) .* sum(reduction_of_MODIS_weights.* upper_matrix,2);
%             end
%             
%             ERA_weights = ERA_weights(:,2:end-1);
%             ERA_time = ERA_time(:,2:end-1);
%             ERA_T2m_downscaled = ERA_T2m_downscaled(:,2:end-1);
%             
%             for i=1:floor(365./post_proc.PARA.averaging_period)+1
%                 pos_year = forcing.TEMP.current_year - (forcing.PARA.ERA_data_years(1) - forcing.PARA.number_of_spin_up_years) + 1 - post_proc.TEMP.offset_years;
%                 
%                 start_index = (i-1).*post_proc.PARA.averaging_period.*4+1;
%                 end_index = min(i*post_proc.PARA.averaging_period*4, size(forcing.TEMP.ERA_T_downscaled,2));
%                 range_MODIS = double(repmat(ERA_time(1, start_index), size(MODIS_LST,1), 1)<= MODIS_LST_time & repmat(ERA_time(1, end_index)+0.25, size(MODIS_LST,1), 1) > MODIS_LST_time);
%                 
%                 forcing.DATA.final_av_T(:,i,pos_year) = (sum((MODIS_LST_weight .* range_MODIS .* MODIS_LST),2) + sum(ERA_weights(:, start_index:end_index) .* ERA_T2m_downscaled(:, start_index:end_index) ,2)) ./(sum(MODIS_LST_weight .* range_MODIS,2) + sum(ERA_weights(:,start_index:end_index),2))-273.15;
%                 
%                 forcing.DATA.final_MODIS_weight(:,i,pos_year) = sum(MODIS_LST_weight .* range_MODIS,2) ./ (sum(MODIS_LST_weight.* range_MODIS,2) + sum(ERA_weights(:,start_index:end_index),2));
%             end
%         end
%         
%         
%         
%         function LT = SolarTime2LocalTime(post_proc, LST_time, day_of_year, longitude, GMTdiff)
%             % LT = SolarTime2LocalTime(LST,longitude,GMTdiff)
%             %
%             % LSTM: Local Standard Time Meridian (LSTM)
%             % LST : local solar time
%             % LT  : local time
%             % EoT : equation of time
%             % GMT diff: local time zone diff from GMT
%             %
%             % AWS1 @ ASF
%             % convert MODIS local solar time to AWS UTC+2 time
%             % longitude = 22.4;
%             % GMTdiff   = 2; (summer time)
%             
%             
%             % LSTM
%             LSTM = 15 * GMTdiff; %time diff from GMT (zone)
%             
%             % EoT: equation of time (minutes)
%             % The equation of time (EoT) (in minutes) is an empirical equation that corrects for the eccentricity of the Earth's orbit and the Earth's axial tilt.
%             B   = (360/365) * (day_of_year-81);
%             EoT = 9.87*sind(2*B) - 7.53*cosd(B) + 1.5*sind(B);
%             
%             % Time correction factor (minutes)
%             %The net Time Correction Factor (in minutes) accounts for the variation of the Local Solar Time (LST) within a given time zone due to the longitude variations within the time zone and also incorporates the EoT above.
%             
%             TC = 4*(LSTM-longitude) + EoT;
%             
%             
%             % LT (daynums) Reversed sign since we want local time
%             LT = LST_time + TC/(60*24);
%         end
        
%                 %-------------param file generation-----
%         function post_proc = param_file_info(post_proc)
%             post_proc = provide_PARA(post_proc);
% 
%             post_proc.PARA.STATVAR = [];
%             post_proc.PARA.class_category = 'FORCING POST_PROCESSING';
%             post_proc.PARA.options = [];
%             
%             post_proc.PARA.eliminate_fraction = [];
%             post_proc.PARA.survive_fraction = [];
%                         
%             post_proc.PARA.default_value.window_size = {7};
%             post_proc.PARA.comment.window_size = {'window size in days within which precipitation is reallocated'};
%             
%             post_proc.PARA.default_value.eliminate_fraction = {0.5};
%             post_proc.PARA.comment.eliminate_fraction = {'fraction of smallest precipitation events (= timestamps with precipitation) that is reallocated to larger events'};
%             
%             post_proc.PARA.default_value.survive_fraction = {0.5};  
%             post_proc.PARA.comment.survive_fraction = {'fraction of largest precipitation events (= timestamps with precipitation) that the small events are reallocated to'};
%             
%         end
        
    end
    
end