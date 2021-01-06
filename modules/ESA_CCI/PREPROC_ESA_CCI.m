
classdef PREPROC_ESA_CCI < BASE
    
    properties
        MODIS_CLASS
    end
    
    
    methods
        
        %-----initialize-----------------
        
        function preproc = provide_PARA(preproc)
            preproc.PARA.timestep = []; %[days]
            preproc.PARA.threshold_diff_MODIS_ERA = []; %[degree C]
            preproc.PARA.threshold_T_snowfall = []; %[degree C]
            preproc.PARA.threshold_T_snowmelt = [];

        end
        
        
        function preproc = provide_CONST(preproc)
            preproc.CONST.L_f = []; %3.34e8;
 
        end
       
        
        function preproc = provide_STATVAR(preproc)

     
        end
        
        function preproc = finalize_init(preproc, tile)
            preproc.MODIS_CLASS = tile.RUN_INFO.STATVAR.spatial;

        end
        
        %-----mandatory functions------------------------
  
        function preproc = get_boundary_condition_u(preproc, tile) %get_MODIS
            
            preproc.MODIS_CLASS = get_data(preproc.MODIS_CLASS, tile.t, tile.timestep);
        
        end
      
        function preproc = get_derivatives_prognostic(preproc, tile)  %calculate snowfall and melt
                        
            preproc.STATVAR.ERA_snowfall_downscaled = tile.FORCING.TEMP.ERA_precip_downcaled(:,1:end-4) .*double(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4)-273.15 <= preproc.PARA.threshold_T_snowfall) .* double(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4)~=0);
            preproc.STATVAR.ERA_snowfall_downscaled = sum(preproc.STATVAR.ERA_snowfall_downscaled,2);
            
            %preproc.STATVAR.ERA_snowfall_downscaled = preproc.STATVAR.ERA_snowfall_downscaled(:,1:end-4); %remove ninth day
            
            
            %degree_day_factor =  DDFsFromDayLenAndSunAngle(datenum(year, 1, 1) + (jj-1).*8 + 3.5, tile.RUN_INFO.STATVAR.latitude); %mm/day/degreeC -> using Jaros's function taking latitudes into account
            
            preproc = DDFsFromDayLenAndSunAngle(preproc, tile);
            %average T for the 8d period
            preproc.STATVAR.ERA_melt_degree_days = sum(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4) - 273.15,2) ./ sum(double(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4) > 0),2); %average T removing all NaN values
            preproc.STATVAR.ERA_melt_degree_days = 8 .* ( preproc.STATVAR.ERA_melt_degree_days - preproc.PARA.threshold_T_snowmelt) .* double(preproc.STATVAR.ERA_melt_degree_days >= preproc.PARA.threshold_T_snowmelt);
            preproc.STATVAR.ERA_melt_T_index = preproc.STATVAR.ERA_melt_degree_days .*preproc.STATVAR.degreeDayFactor; %in mm
            
            preproc.STATVAR.ERA_T_downscaled = nanmean(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4) - 273.15, 2);
            
        end
        
        function preproc = compute_diagnostic(preproc, tile)  %merge MODIS and ERA 
            %[final_av_T, final_MODIS_weight] = merge_MODIS_ERA(MODIS_LST, MODIS_LST_time, ERA_T2m_downscaled, ERA_time)
            
            MODIS_LST_time = preproc.MODIS_CLASS.STATVAR.timestamp;
            MODIS_LST = preproc.MODIS_CLASS.STATVAR.LST;
            
            ERA_time = tile.FORCING.TEMP.ERA_time;
            ERA_T2m_downscaled = tile.FORCING.TEMP.ERA_T_downscaled;
            
            %make the 4 dependent on
            %tile.FORCING.PARA.number_of_values_per_day
            
            lower = (floor(MODIS_LST_time.*4)./4 - ERA_time(1,1)).*4 + 1;
            upper= (floor(MODIS_LST_time.*4)./4 + 0.25-ERA_time(1,1)).*4 + 1;
            
            
            
            weight_left = 1 - (MODIS_LST_time - floor(MODIS_LST_time.*4)./4 ) ./0.25;
            weight_right= 1 - weight_left;
            
            lower(isnan(lower)) = 1;
            upper(isnan(upper)) = 2;
            
            weight_left(lower<1 | upper>size(ERA_time,2)) = NaN;
            weight_right(lower<1 | upper>size(ERA_time,1)) = NaN;
            
            %MODIS_LST(lower<1)=0;
            upper(lower<1)=2;
            lower(lower<1)=1;
            %MODIS_LST(upper>size(ERA_time,2))=0;
            lower(upper>size(ERA_time,2)) = size(ERA_time,2) - 1;
            upper(upper>size(ERA_time,2)) = size(ERA_time,2);
            
            
            
            MODIS_LST_interp=MODIS_LST.*0;  %interpolate ERA to MODIS LST timestamp
            
            matrix = repmat([1:size(ERA_T2m_downscaled,2)], size(ERA_T2m_downscaled,1),1);            
            for i=1:size(MODIS_LST,2)
                %MODIS_LST_interp(:,i) = weight_left(:,i).*diag(ERA_T2m_downscaled(:,lower(:,i))) + weight_right(:,i) .* diag(ERA_T2m_downscaled(:,upper(:,i)));
                %faster version
                upper_matrix = double(matrix == repmat(upper(:,i), 1, size(ERA_T2m_downscaled,2)));
                lower_matrix = double(matrix == repmat(lower(:,i), 1, size(ERA_T2m_downscaled,2)));
                
                MODIS_LST_interp(:,i) = weight_left(:,i) .* sum(ERA_T2m_downscaled .* lower_matrix,2) + weight_right(:,i) .* sum(ERA_T2m_downscaled .* upper_matrix,2);
                
            end
            MODIS_LST_interp(isnan(MODIS_LST_interp))=0;
            
            
            %preproc.PARA.threshold_diff = 10; %MODIS LST deleted if more than 10K colder
            %weight_left(MODIS_LST<MODIS_LST_interp-10) = NaN;
            %weight_right(MODIS_LST<MODIS_LST_interp-10) = NaN;
            %MODIS_LST(MODIS_LST<MODIS_LST_interp-threshold_diff) = 0;
            
           %             preproc.PARA.threshold_diff = 10; %MODIS LST deleted if more than 10K colder
            weight_left(MODIS_LST<MODIS_LST_interp-preproc.PARA.threshold_diff_MODIS_ERA) = NaN;
            weight_right(MODIS_LST<MODIS_LST_interp-preproc.PARA.threshold_diff_MODIS_ERA) = NaN;
            MODIS_LST(MODIS_LST<MODIS_LST_interp-preproc.PARA.threshold_diff_MODIS_ERA) = 0;
            
            %MODIS_LST_mean_reduced = sum(MODIS_LST,2) ./sum(double(MODIS_LST~=0),2)-273.15;
            
            weight_left(isnan(weight_left))=0;
            weight_right(isnan(weight_right))=0;
            
            total_weights=ERA_T2m_downscaled.*0;
            for j=1:size(MODIS_LST,2)
                for i=1:size(MODIS_LST,1)
                    total_weights(i, lower(i,j)) = total_weights(i, lower(i,j)) + weight_left(i,j);
                    total_weights(i, upper(i,j)) = total_weights(i, upper(i,j)) + weight_right(i,j);
                end
            end
            
            total_weights_reduced = total_weights ./ max(1,total_weights);
            ERA_weights = 1-total_weights_reduced;
            ERA_weights(:,end)=0;
            
            reduction_of_MODIS_weights = total_weights_reduced ./ max(1e-20, total_weights);
            
            MODIS_LST_weight = MODIS_LST.*0;

            matrix = repmat([1:size(reduction_of_MODIS_weights,2)], size(reduction_of_MODIS_weights,1),1);     
            for i=1:size(MODIS_LST,2)
                %MODIS_LST_weight(:,i) = weight_left(:,i) .* diag(reduction_of_MODIS_weights(:,lower(:,i))) + weight_right(:,i) .* diag(reduction_of_MODIS_weights(:,upper(:,i)));
                %faster version
                upper_matrix = double(matrix == repmat(upper(:,i), 1, size(reduction_of_MODIS_weights,2)));
                lower_matrix = double(matrix == repmat(lower(:,i), 1, size(reduction_of_MODIS_weights,2)));
                MODIS_LST_weight(:,i) = weight_left(:,i) .* sum(reduction_of_MODIS_weights .* lower_matrix,2) + weight_right(:,i) .* sum(reduction_of_MODIS_weights.* upper_matrix,2);
                
            end
            
            % size(MODIS_LST_weight)
            % size(MODIS_LST)
            % size(ERA_weights)
            % size(ERA_T2m_downscaled)
            
            preproc.STATVAR.final_av_T = (sum((MODIS_LST_weight .* MODIS_LST),2) + sum(ERA_weights(:,1:end-1) .* ERA_T2m_downscaled(:,1:end-1) ,2)) ./(sum(MODIS_LST_weight,2) + sum(ERA_weights(:,1:end-1),2))-273.15;
            
            preproc.STATVAR.final_MODIS_weight = sum(MODIS_LST_weight,2) ./ (sum(MODIS_LST_weight,2) + sum(ERA_weights(:,1:end-1),2));
            
        end
        
        
        %non-madatory functions
        
        function preproc = DDFsFromDayLenAndSunAngle(preproc, tile)   %mm/day/degreeC -> using Jaros's function taking latitudes into account
            %timestamp,latitude)
            
            %Function calculates average degree day factors (DDFs) for a specific day
            %based on sun culmination angle and day length. For this, it uses "day_length" and
            %"solarCulminationAngle" functions. The DDFs are linearly scaled between a
            %product of culmination angle and day length.
            
            %maximum and minimum DDFs to which DDFs are scaled
            ddfmin = 2;
            ddfmax = 12;
            
            %maximum product between culmination angle and day length/2 in degrees.
            maxProduct = 10169.237;
            
            %---------------------------------
            %solar culmination angle
            %curcul = solarCulminationAngle(timestamp,latitude);
            
            %day of the year
            doy = tile.t - datenum(str2num(datestr(tile.t, 'yyyy')),1,1) + 1 + 3.5;  %timestamp in the middle of 8d interval

            %declination angle formula
            declination_angle = 23.45*sind(360/365*(284+doy));
            culmination_angle = 90 - tile.RUN_INFO.STATVAR.latitude + declination_angle;

            %only positive
            culmination_angle(culmination_angle < 0) = 0;
            
            %day of the year according to winter solstice
            dlDoy = doy + 11;
            
            %day lengths in hours
            %replaces
            %curlen = day_length(dlDoy,latitude);
            
            Axis=23.439*pi./180;

            j=pi./182.625;%constant (radians)
            m=1 - tan(tile.RUN_INFO.STATVAR.latitude.*pi/180).*tan(Axis.*cos(j.*dlDoy));
            
            m(m>2)=2;%saturate value for artic
            m(m<0)=0;
            
            b=acos(1-m)/pi;%fraction of the day the sun is up
            curlen=b*24;%hours of sunlight
            curlen(curlen < 0) = 0;
            %scaling hours to degrees
            curlen_deg = (curlen/24)*360;
            
            %multiplying culmination angle with half of the daylength degrees. The
            %product equals approx. area of sun path above horizont.
            daylculm = culmination_angle.*(curlen_deg/2);
            %converting to fraction based on maximum šroduct
            daylculm_frac = daylculm/maxProduct;
            
            %scaling DDFs to maximum and minimum DDF
            preproc.STATVAR.degreeDayFactor= daylculm_frac*(ddfmax-ddfmin)+ddfmin;
            
        end
        
    end
    
end

