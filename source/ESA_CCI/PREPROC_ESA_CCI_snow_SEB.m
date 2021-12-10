
classdef PREPROC_ESA_CCI_snow_SEB < BASE
    
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
            
            preproc.PARA.Tr=276.15; % Threshold temperature for rainfall.
            preproc.PARA.Ts=272.15; % Threshold tempearture for snowfall.

            preproc.PARA.emissivity_snow = 0.99; % Snow emissivity (assumed known).
            
            preproc.PARA.taus=0.01; % Threshold snowfall for resetting to maximum [m w.e.].
            preproc.PARA.taua=0.008; % Time constant for snow albedo change in non-melting conditions [/day].
            preproc.PARA.tauf=0.24; % Time constant for snow albedo change in melting conditions [/day].
            
            preproc.PARA.albsmax_bare=0.85; % Maximum snow albedo.
            preproc.PARA.albsmin_bare=0.5; % Minimum snow albedo.
            
            preproc.PARA.albsmax_forest=0.6; % Maximum snow albedo.
            preproc.PARA.albsmin_forest=0.3; % Minimum snow albedo.
            
            preproc.PARA.ar=0.2/1e2; % Restricted degree day factor (m*degC/day , value from Burbaker et al. 1996)
            
            preproc.PARA.canopy_transmissivity = 0.3;
            preproc.PARA.emissivity_canopy = 0.96;
            
        end
        
        
        function preproc = provide_CONST(preproc)
            preproc.CONST.L_f = []; %3.34e8;
            preproc.CONST.sigma = [];
            preproc.CONST.day_sec = [];
            preproc.CONST.Tmfw = [];
         end
       
        
        function preproc = provide_STATVAR(preproc)

     
        end
        
        function preproc = finalize_init(preproc, tile)
            preproc.MODIS_CLASS = tile.RUN_INFO.STATVAR.spatial;

            preproc.PARA.Lupwelling = preproc.PARA.emissivity_snow.*preproc.CONST.sigma.*preproc.CONST.Tmfw.^4; % upwelling longwave radiation for T=273.15K
            
            
            preproc.STATVAR.albedo_bare = repmat((preproc.PARA.albsmax_bare + preproc.PARA.albsmin_bare)./2, size(tile.RUN_INFO.STATVAR.key,1), 1);
            preproc.STATVAR.albedo_forest = repmat((preproc.PARA.albsmax_forest + preproc.PARA.albsmin_forest)./2, size(tile.RUN_INFO.STATVAR.key,1), 1);
        end
        
        %-----mandatory functions------------------------
  
        function preproc = get_boundary_condition_u(preproc, tile) %get_MODIS
            
            preproc.MODIS_CLASS = get_data(preproc.MODIS_CLASS, tile.t, tile.timestep);
        
        end
      
        function preproc = get_derivatives_prognostic(preproc, tile)  %calculate snowfall and melt
%             tile.FORCING.TEMP.ERA_Sin_downscaled(1,:)
%             tile.FORCING.TEMP.ERA_Lin_downscaled(1,:)
%             tile.FORCING.TEMP.ERA_time

            preproc = get_snowfall(preproc, tile);
            preproc = get_melt_SEB(preproc, tile);
            
            preproc.STATVAR.ERA_snowfall_downscaled = mean(preproc.STATVAR.ERA_snowfall_downscaled, 2); %in mm/day
            
            preproc.STATVAR.ERA_T_downscaled = nanmean(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4) - 273.15, 2);

            
%             preproc.STATVAR.ERA_snowfall_downscaled = tile.FORCING.TEMP.ERA_precip_downcaled(:,1:end-4) .*double(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4)-273.15 <= preproc.PARA.threshold_T_snowfall) .* double(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4)~=0);
%             preproc.STATVAR.ERA_snowfall_downscaled = sum(preproc.STATVAR.ERA_snowfall_downscaled,2) ./ 8; %in mm/day
%             
%             %preproc.STATVAR.ERA_snowfall_downscaled = preproc.STATVAR.ERA_snowfall_downscaled(:,1:end-4); %remove ninth day
%             
%             
%             %degree_day_factor =  DDFsFromDayLenAndSunAngle(datenum(year, 1, 1) + (jj-1).*8 + 3.5, tile.RUN_INFO.STATVAR.latitude); %mm/day/degreeC -> using Jaros's function taking latitudes into account
%             
%             preproc = DDFsFromDayLenAndSunAngle(preproc, tile);
%             
%             %average T for the 8d period
%             preproc.STATVAR.ERA_melt_degree_days = sum(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4) - 273.15,2) ./ sum(double(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4) > 0),2); %average T removing all NaN values
%             preproc.STATVAR.ERA_melt_degree_days = 8 .* ( preproc.STATVAR.ERA_melt_degree_days - preproc.PARA.threshold_T_snowmelt) .* double(preproc.STATVAR.ERA_melt_degree_days >= preproc.PARA.threshold_T_snowmelt);
%             preproc.STATVAR.ERA_melt_T_index = preproc.STATVAR.ERA_melt_degree_days .*preproc.STATVAR.degreeDayFactor ./ 8; %in mm/day
%             
%             preproc.STATVAR.ERA_T_downscaled = nanmean(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4) - 273.15, 2);
            
        end
        
        function preproc = compute_diagnostic(preproc, tile)  %merge MODIS and ERA 
            %[final_av_T, final_MODIS_weight] = merge_MODIS_ERA(MODIS_LST, MODIS_LST_time, ERA_T2m_downscaled, ERA_time)
            
            compute_slice = 10000;
            preproc.STATVAR.final_av_T = preproc.STATVAR.ERA_T_downscaled.*0;
            preproc.STATVAR.final_MODIS_weight = preproc.STATVAR.ERA_T_downscaled.*0;
            
            for count = 1:compute_slice:size(preproc.MODIS_CLASS.STATVAR.timestamp,1)
                range = [count:min(count+compute_slice-1, size(preproc.MODIS_CLASS.STATVAR.timestamp,1))];
                                
                MODIS_LST_time = preproc.MODIS_CLASS.STATVAR.timestamp(range,:);
                MODIS_LST = preproc.MODIS_CLASS.STATVAR.LST(range,:);
                
                ERA_time = tile.FORCING.TEMP.ERA_time;
                ERA_T2m_downscaled = tile.FORCING.TEMP.ERA_T_downscaled(range,:);
                
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
                
                
                preproc.STATVAR.final_av_T(range,1) = (sum((MODIS_LST_weight .* MODIS_LST),2) + sum(ERA_weights(:,1:end-1) .* ERA_T2m_downscaled(:,1:end-1) ,2)) ./(sum(MODIS_LST_weight,2) + sum(ERA_weights(:,1:end-1),2))-273.15;
                
                preproc.STATVAR.final_MODIS_weight(range,1) = sum(MODIS_LST_weight,2) ./ (sum(MODIS_LST_weight,2) + sum(ERA_weights(:,1:end-1),2));
                
            end
%             preproc.STATVAR.final_av_T = preproc.STATVAR.ERA_T_downscaled;
%            
%             preproc.STATVAR.final_MODIS_weight = preproc.STATVAR.ERA_T_downscaled;

        end
        
        
        %non-madatory functions
        
        
        function preproc = get_snowfall(preproc, tile)

            frain = (tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4)-preproc.PARA.Ts)./(preproc.PARA.Tr-preproc.PARA.Ts);  % Rain fraction.
            frain(frain>1)=1; 
            frain(frain<0)=0; % Must be between 0 and 1.
            preproc.STATVAR.ERA_snowfall_downscaled = (1-frain) .* tile.FORCING.TEMP.ERA_precip_downcaled(:,1:end-4) ; % mm/day
            
            preproc.STATVAR.ERA_snowfall_downscaled = preproc.STATVAR.ERA_snowfall_downscaled .* double(tile.FORCING.TEMP.ERA_T_downscaled(:,1:end-4)~=0);

        end
        
         function preproc = get_melt_SEB(preproc, tile)
            
            % Time loop.
            melt_depth_bare = 0;
            melt_depth_forest = 0;
            for i= 1:8
                daily_melt_depth_bare =0;
                daily_melt_depth_forest =0;
                for n=1:4
                    
                    % Ablation term
                    Ldnet_bare = preproc.PARA.emissivity_snow.*tile.FORCING.TEMP.ERA_Lin_downscaled(:, (i-1).*4+n) - preproc.PARA.Lupwelling; % Net downwelling longwave
                    Ldnet_forest = preproc.PARA.emissivity_snow.*tile.FORCING.TEMP.ERA_Lin_downscaled(:, (i-1).*4+n) .* preproc.PARA.canopy_transmissivity + ...
                        (1-preproc.PARA.canopy_transmissivity) .* preproc.PARA.emissivity_canopy .* preproc.CONST.sigma .* (tile.FORCING.TEMP.ERA_T_downscaled(:, (i-1).*4+n)) .^4 ...
                        - preproc.PARA.Lupwelling;
                    
                    Sdnet_bare = (1-preproc.STATVAR.albedo_bare).*tile.FORCING.TEMP.ERA_Sin_downscaled(:, (i-1).*4+n); % Net downwelling shortwave
                    Sdnet_forest = (1-preproc.STATVAR.albedo_forest).*tile.FORCING.TEMP.ERA_Sin_downscaled(:, (i-1).*4+n) .* preproc.PARA.canopy_transmissivity; % Net downwelling shortwave

                    Tnet=preproc.PARA.ar.*(tile.FORCING.TEMP.ERA_T_downscaled(:, (i-1).*4+n) - preproc.CONST.Tmfw); % Warming through turbulent heat fluxes, parametrized using a restricted degree day approach.
                    
                    daily_melt_depth_bare = daily_melt_depth_bare + Sdnet_bare + Ldnet_bare + Tnet; % Melt depth over the time step.
                    daily_melt_depth_forest = daily_melt_depth_forest + Sdnet_forest + Ldnet_forest + Tnet; % Melt depth over the time step.
                    
                    % The above can be negative if -Rad>Tnet, it's important to count this
                    % in diurnal accumulated melt depths, since this implicitly accounts for cold
                    % content on the subdaily time-scale.

                end
                daily_melt_depth_bare = daily_melt_depth_bare ./ 4 .* preproc.CONST.day_sec ./ preproc.CONST.L_f .* 1000;  %in mm/day
                daily_melt_depth_forest = daily_melt_depth_forest ./ 4 .* preproc.CONST.day_sec ./ preproc.CONST.L_f .* 1000;  %in mm/day
                
                melt_depth_bare = melt_depth_bare + daily_melt_depth_bare;
                melt_depth_forest = melt_depth_forest + daily_melt_depth_forest;
                
                % Update snow albedo for next step.
                % Latest ECMWF "continuous reset" snow albedo scheme (Dutra et al. 2010)
                new_snow = mean(preproc.STATVAR.ERA_snowfall_downscaled(:, (i-1)*4+1:i*4), 2); %in mm/day
                
                %bare
                net_acc = new_snow - daily_melt_depth_bare; % Net accumulation for one day time-step.
                accumulation = net_acc>0;
                preproc.STATVAR.albedo_bare(accumulation,1) = preproc.STATVAR.albedo_bare(accumulation,1) + min(1,net_acc(accumulation,1)./preproc.PARA.taus) .* ...
                    (preproc.PARA.albsmax_bare - preproc.STATVAR.albedo_bare(accumulation,1));

                no_melting = (net_acc==0); % "Steady" case (linear decay)
                preproc.STATVAR.albedo_bare(no_melting,1) = preproc.STATVAR.albedo_bare(no_melting,1) - preproc.PARA.taua;
                
                melting = net_acc<0; % Ablating case (exponential decay)
                preproc.STATVAR.albedo_bare(melting,1) = (preproc.STATVAR.albedo_bare(melting,1) - preproc.PARA.albsmin_bare) .* exp(-preproc.PARA.tauf) + preproc.PARA.albsmin_bare;
                preproc.STATVAR.albedo_bare(preproc.STATVAR.albedo_bare < preproc.PARA.albsmin_bare) = preproc.PARA.albsmin_bare;
                
                %forest
                net_acc = new_snow - daily_melt_depth_forest; % Net accumulation for one day time-step.
                accumulation = net_acc>0;
                preproc.STATVAR.albedo_forest(accumulation,1) = preproc.STATVAR.albedo_forest(accumulation,1) + min(1,net_acc(accumulation,1)./preproc.PARA.taus) .* ...
                    (preproc.PARA.albsmax_forest - preproc.STATVAR.albedo_forest(accumulation,1));
                
                no_melting = (net_acc==0); % "Steady" case (linear decay)
                preproc.STATVAR.albedo_forest(no_melting,1) = preproc.STATVAR.albedo_forest(no_melting,1) - preproc.PARA.taua;
                
                melting = net_acc<0; % Ablating case (exponential decay)
                preproc.STATVAR.albedo_forest(melting,1) = (preproc.STATVAR.albedo_forest(melting,1) - preproc.PARA.albsmin_forest) .* exp(-preproc.PARA.tauf) + preproc.PARA.albsmin_forest;
                preproc.STATVAR.albedo_forest(preproc.STATVAR.albedo_forest < preproc.PARA.albsmin_forest) = preproc.PARA.albsmin_forest;
                
            end
            preproc.STATVAR.ERA_melt_bare = single(max(0, melt_depth_bare ./ 8)); %in mm/day
            preproc.STATVAR.ERA_melt_forest = single(max(0, melt_depth_forest ./ 8)); %in mm/day

        end
                
%         function preproc = get_melt_SEB(preproc, tile)
%             
%             % Time loop.
%             melt_depth_bare = 0;
%             melt_depth_forest = 0;
%             for i= 1:8
%                 daily_melt_depth_bare =0;
%                 daily_melt_depth_forest =0;
%                 for n=1:4
%                     
%                     % Ablation term
%                     Ldnet = preproc.PARA.emissivity_snow.*tile.FORCING.TEMP.ERA_Lin_downscaled(:, (i-1).*4+n) - preproc.PARA.Lupwelling; % Net downwelling longwave
%                     
%                     Sdnet_bare = (1-preproc.STATVAR.albedo_bare).*tile.FORCING.TEMP.ERA_Sin_downscaled(:, (i-1).*4+n); % Net downwelling shortwave
%                     Sdnet_forest = (1-preproc.STATVAR.albedo_forest).*tile.FORCING.TEMP.ERA_Sin_downscaled(:, (i-1).*4+n); % Net downwelling shortwave
% 
%                     Tnet=preproc.PARA.ar.*(tile.FORCING.TEMP.ERA_T_downscaled(:, (i-1).*4+n) - preproc.CONST.Tmfw); % Warming through turbulent heat fluxes, parametrized using a restricted degree day approach.
%                     
%                     daily_melt_depth_bare = daily_melt_depth_bare + Sdnet_bare + Ldnet + Tnet; % Melt depth over the time step.
%                     daily_melt_depth_forest = daily_melt_depth_forest + Sdnet_forest + Ldnet + Tnet; % Melt depth over the time step.
%                     
%                     % The above can be negative if -Rad>Tnet, it's important to count this
%                     % in diurnal accumulated melt depths, since this implicitly accounts for cold
%                     % content on the subdaily time-scale.
% 
%                 end
%                 daily_melt_depth_bare = daily_melt_depth_bare ./ 4 .* preproc.CONST.day_sec ./ preproc.CONST.L_f .* 1000;  %in mm/day
%                 daily_melt_depth_forest = daily_melt_depth_forest ./ 4 .* preproc.CONST.day_sec ./ preproc.CONST.L_f .* 1000;  %in mm/day
%                 
%                 melt_depth_bare = melt_depth_bare + daily_melt_depth_bare;
%                 melt_depth_forest = melt_depth_forest + daily_melt_depth_forest;
%                 
%                 % Update snow albedo for next step.
%                 % Latest ECMWF "continuous reset" snow albedo scheme (Dutra et al. 2010)
%                 new_snow = mean(preproc.STATVAR.ERA_snowfall_downscaled(:, (i-1)*4+1:i*4), 2); %in mm/day
%                 
%                 %bare
%                 net_acc = new_snow - daily_melt_depth_bare; % Net accumulation for one day time-step.
%                 accumulation = net_acc>0;
%                 preproc.STATVAR.albedo_bare(accumulation,1) = preproc.STATVAR.albedo_bare(accumulation,1) + min(1,net_acc(accumulation,1)./preproc.PARA.taus) .* ...
%                     (preproc.PARA.albsmax_bare - preproc.STATVAR.albedo_bare(accumulation,1));
% 
%                 no_melting = (net_acc==0); % "Steady" case (linear decay)
%                 preproc.STATVAR.albedo_bare(no_melting,1) = preproc.STATVAR.albedo_bare(no_melting,1) - preproc.PARA.taua;
%                 
%                 melting = net_acc<0; % Ablating case (exponential decay)
%                 preproc.STATVAR.albedo_bare(melting,1) = (preproc.STATVAR.albedo_bare(melting,1) - preproc.PARA.albsmin_bare) .* exp(-preproc.PARA.tauf) + preproc.PARA.albsmin_bare;
%                 preproc.STATVAR.albedo_bare(preproc.STATVAR.albedo_bare < preproc.PARA.albsmin_bare) = preproc.PARA.albsmin_bare;
%                 
%                 %forest
%                 net_acc = new_snow - daily_melt_depth_forest; % Net accumulation for one day time-step.
%                 accumulation = net_acc>0;
%                 preproc.STATVAR.albedo_forest(accumulation,1) = preproc.STATVAR.albedo_forest(accumulation,1) + min(1,net_acc(accumulation,1)./preproc.PARA.taus) .* ...
%                     (preproc.PARA.albsmax_forest - preproc.STATVAR.albedo_forest(accumulation,1));
%                 
%                 no_melting = (net_acc==0); % "Steady" case (linear decay)
%                 preproc.STATVAR.albedo_forest(no_melting,1) = preproc.STATVAR.albedo_forest(no_melting,1) - preproc.PARA.taua;
%                 
%                 melting = net_acc<0; % Ablating case (exponential decay)
%                 preproc.STATVAR.albedo_forest(melting,1) = (preproc.STATVAR.albedo_forest(melting,1) - preproc.PARA.albsmin_forest) .* exp(-preproc.PARA.tauf) + preproc.PARA.albsmin_forest;
%                 preproc.STATVAR.albedo_forest(preproc.STATVAR.albedo_forest < preproc.PARA.albsmin_forest) = preproc.PARA.albsmin_forest;
%                 
%             end
%             preproc.STATVAR.ERA_melt_bare = single(max(0, melt_depth_bare ./ 8)); %in mm/day
%             preproc.STATVAR.ERA_melt_forest = single(max(0, melt_depth_forest ./ 8)); %in mm/day
% 
%         end
        
        
        
        
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

