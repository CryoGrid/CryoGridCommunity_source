%========================================================================
% CryoGrid FORCING post-processing 
%
%
% Authors:
% S. Westermann, January 2023
%
%========================================================================

classdef get_snowfall_melt < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.averaging_period = [];  
            post_proc.PARA.year_range = [];
            
            post_proc.PARA.annual = [];
            
%            post_proc.PARA.threshold_T_snowmelt = [];
            
            post_proc.PARA.Tr=274.15; % Threshold temperature for rainfall.
            post_proc.PARA.Ts=272.15; % Threshold tempearture for snowfall.

            post_proc.PARA.emissivity_snow = 0.99; % Snow emissivity (assumed known).
            
            post_proc.PARA.taus=0.0025; % Threshold snowfall for resetting to maximum [m w.e.].
            post_proc.PARA.taua=0.008; % Time constant for snow albedo change in non-melting conditions [/day].
            post_proc.PARA.tauf=0.24; % Time constant for snow albedo change in melting conditions [/day].
            
            post_proc.PARA.albsmax_bare=0.85; % Maximum snow albedo.
            post_proc.PARA.albsmin_bare=0.5; % Minimum snow albedo.
            
            post_proc.PARA.albsmax_forest=0.6; % Maximum snow albedo.
            post_proc.PARA.albsmin_forest=0.3; % Minimum snow albedo.
            
            post_proc.PARA.ar=0.2/1e2; % Restricted degree day factor (m*degC/day , value from Burbaker et al. 1996)
            
            post_proc.PARA.canopy_transmissivity = 0.3;
            post_proc.PARA.emissivity_canopy = 0.96;
            
        end
        
        
        function post_proc = provide_CONST(post_proc)
            post_proc.CONST.L_f = []; 
            post_proc.CONST.sigma = [];
            post_proc.CONST.day_sec = [];
            post_proc.CONST.Tmfw = [];
        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)

            post_proc.PARA.Lupwelling = post_proc.PARA.emissivity_snow.*post_proc.CONST.sigma.*post_proc.CONST.Tmfw.^4; % upwelling longwave radiation for melting snow, T=273.15K
            post_proc.STATVAR.albedo_bare = repmat(post_proc.PARA.albsmax_bare, size(tile.PARA.latitude,1), 1); %initialize 1st albedo values 
            post_proc.STATVAR.albedo_forest = repmat(post_proc.PARA.albsmax_forest, size(tile.PARA.latitude,1), 1);
            
            post_proc.TEMP.offset_years = 0;
        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
            if forcing.TEMP.current_year >= post_proc.PARA.year_range(1) && forcing.TEMP.current_year <= post_proc.PARA.year_range(end)
                pos_year = forcing.TEMP.current_year - (forcing.PARA.ERA_data_years(1)- forcing.PARA.number_of_spin_up_years) + 1 - post_proc.TEMP.offset_years;
                snowfall = get_snowfall(post_proc, tile);
                for i=1:floor(365./post_proc.PARA.averaging_period)+1
                    forcing.DATA.ERA_snowfall_downscaled(:,i,pos_year) = nanmean(snowfall(:, (i-1).*post_proc.PARA.averaging_period.*4+1:min(i*post_proc.PARA.averaging_period*4, size(snowfall,2))),2);
                end
                post_proc.STATVAR.ERA_snowfall_downscaled = snowfall; %needed to get albedo for melt
                [melt_depth_bare, melt_depth_forest] = get_melt_SEB(post_proc, tile);
                for i=1:floor(365./post_proc.PARA.averaging_period)+1
                    forcing.DATA.ERA_melt_bare(:,i,pos_year) = nanmean(melt_depth_bare(:, (i-1).*post_proc.PARA.averaging_period+1:min(i*post_proc.PARA.averaging_period, size(melt_depth_bare,2))),2);
                    forcing.DATA.ERA_melt_forest(:,i,pos_year) = nanmean(melt_depth_forest(:, (i-1).*post_proc.PARA.averaging_period+1:min(i*post_proc.PARA.averaging_period, size(melt_depth_forest,2))),2);
                end
                forcing.DATA.ERA_melt_bare(forcing.DATA.ERA_melt_bare<0) = 0;
                forcing.DATA.ERA_melt_forest(forcing.DATA.ERA_melt_forest<0) = 0;
                snowfall = [];
                melt_depth_bare = [];
                melt_depth_forest = [];
            end
        end
        
        
        
        %service functions
        function snowfall = get_snowfall(post_proc, tile)

            frain = (tile.FORCING.TEMP.ERA_T_downscaled-post_proc.PARA.Ts)./(post_proc.PARA.Tr-post_proc.PARA.Ts);  % Rain fraction.
            frain(frain>1)=1; 
            frain(frain<0)=0; % Must be between 0 and 1.
            snowfall = (1-frain) .* tile.FORCING.TEMP.ERA_precip_downcaled; % mm/day
            snowfall = snowfall .* double(tile.FORCING.TEMP.ERA_T_downscaled~=0);

        end
        
        
        function [melt_depth_bare, melt_depth_forest] = get_melt_SEB(post_proc, tile)
            melt_depth_bare = tile.FORCING.TEMP.ERA_Lin_downscaled(:, 1:size(tile.FORCING.TEMP.ERA_Lin_downscaled,2)/4) .* 0;
            melt_depth_forest = melt_depth_bare;
            
%             melt_depth_bare = 0; 
%             melt_depth_forest = 0;
            for i = 1:size(tile.FORCING.TEMP.ERA_Lin_downscaled, 2)/4
                daily_melt_depth_bare = 0;
                daily_melt_depth_forest = 0;
                for n=1:4
                    % Ablation term
                    Ldnet_bare = post_proc.PARA.emissivity_snow.*tile.FORCING.TEMP.ERA_Lin_downscaled(:, (i-1).*4+n) - post_proc.PARA.Lupwelling; % Net downwelling longwave
                    Ldnet_forest = post_proc.PARA.emissivity_snow.*tile.FORCING.TEMP.ERA_Lin_downscaled(:, (i-1).*4+n) .* post_proc.PARA.canopy_transmissivity + ...
                        (1-post_proc.PARA.canopy_transmissivity) .* post_proc.PARA.emissivity_canopy .* post_proc.CONST.sigma .* (tile.FORCING.TEMP.ERA_T_downscaled(:, (i-1).*4+n)) .^4 ...
                        - post_proc.PARA.Lupwelling;
                    
                    Sdnet_bare = (1-post_proc.STATVAR.albedo_bare).*tile.FORCING.TEMP.ERA_Sin_downscaled(:, (i-1).*4+n); % Net downwelling shortwave
                    Sdnet_forest = (1-post_proc.STATVAR.albedo_forest).*tile.FORCING.TEMP.ERA_Sin_downscaled(:, (i-1).*4+n) .* post_proc.PARA.canopy_transmissivity; % Net downwelling shortwave

                    Tnet=post_proc.PARA.ar.*(tile.FORCING.TEMP.ERA_T_downscaled(:, (i-1).*4+n) - post_proc.CONST.Tmfw); % Warming through turbulent heat fluxes, parametrized using a restricted degree day approach.
                    
                    daily_melt_depth_bare = daily_melt_depth_bare + Sdnet_bare + Ldnet_bare + Tnet; % Melt depth over the time step.
                    daily_melt_depth_forest = daily_melt_depth_forest + Sdnet_forest + Ldnet_forest + Tnet; % Melt depth over the time step.
                    
                    % The above can be negative if -Rad>Tnet, it's important to count this
                    % in diurnal accumulated melt depths, since this implicitly accounts for cold
                    % content on the subdaily time-scale.

                end
                daily_melt_depth_bare = daily_melt_depth_bare ./ 4 .* post_proc.CONST.day_sec ./ post_proc.CONST.L_f .* 1000;  %in mm/day
                daily_melt_depth_forest = daily_melt_depth_forest ./ 4 .* post_proc.CONST.day_sec ./ post_proc.CONST.L_f .* 1000;  %in mm/day
                
                melt_depth_bare(:,i) = daily_melt_depth_bare;
                melt_depth_forest(:,i) = daily_melt_depth_forest;
%                 melt_depth_bare = melt_depth_bare + daily_melt_depth_bare;
%                 melt_depth_forest = melt_depth_forest + daily_melt_depth_forest;
                
                % Update snow albedo for next step.
                % Latest ECMWF "continuous reset" snow albedo scheme (Dutra et al. 2010)
                new_snow = mean(post_proc.STATVAR.ERA_snowfall_downscaled(:, (i-1)*4+1:i*4), 2); %in mm/day
                
                %bare
                net_acc = new_snow - max(0,daily_melt_depth_bare); % Net accumulation for one day time-step.
                accumulation = net_acc>0;
                post_proc.STATVAR.albedo_bare(accumulation,1) = post_proc.STATVAR.albedo_bare(accumulation,1) + min(1,net_acc(accumulation,1)./(post_proc.PARA.taus .* 1000)) .* ...
                    (post_proc.PARA.albsmax_bare - post_proc.STATVAR.albedo_bare(accumulation,1));

                no_melting = (net_acc==0); % "Steady" case (linear decay)
                post_proc.STATVAR.albedo_bare(no_melting,1) = post_proc.STATVAR.albedo_bare(no_melting,1) - post_proc.PARA.taua;
                
                melting = net_acc<0; % Ablating case (exponential decay)
                post_proc.STATVAR.albedo_bare(melting,1) = (post_proc.STATVAR.albedo_bare(melting,1) - post_proc.PARA.albsmin_bare) .* exp(-post_proc.PARA.tauf) + post_proc.PARA.albsmin_bare;
                post_proc.STATVAR.albedo_bare(post_proc.STATVAR.albedo_bare < post_proc.PARA.albsmin_bare) = post_proc.PARA.albsmin_bare;
                
                %forest
                net_acc = new_snow - max(0,daily_melt_depth_forest); % Net accumulation for one day time-step.
                accumulation = net_acc>0;
                post_proc.STATVAR.albedo_forest(accumulation,1) = post_proc.STATVAR.albedo_forest(accumulation,1) + min(1,net_acc(accumulation,1)./(post_proc.PARA.taus.*1000)) .* ...
                    (post_proc.PARA.albsmax_forest - post_proc.STATVAR.albedo_forest(accumulation,1));
                
                no_melting = (net_acc==0); % "Steady" case (linear decay)
                post_proc.STATVAR.albedo_forest(no_melting,1) = post_proc.STATVAR.albedo_forest(no_melting,1) - post_proc.PARA.taua;
                
                melting = net_acc<0; % Ablating case (exponential decay)
                post_proc.STATVAR.albedo_forest(melting,1) = (post_proc.STATVAR.albedo_forest(melting,1) - post_proc.PARA.albsmin_forest) .* exp(-post_proc.PARA.tauf) + post_proc.PARA.albsmin_forest;
                post_proc.STATVAR.albedo_forest(post_proc.STATVAR.albedo_forest < post_proc.PARA.albsmin_forest) = post_proc.PARA.albsmin_forest;
                
            end
%             post_proc.STATVAR.ERA_melt_bare = single(max(0, melt_depth_bare ./ 8)); %in mm/day
%             post_proc.STATVAR.ERA_melt_forest = single(max(0, melt_depth_forest ./ 8)); %in mm/day
        end
        
        
        
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