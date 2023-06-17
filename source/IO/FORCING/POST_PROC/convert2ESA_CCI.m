%========================================================================
% CryoGrid FORCING post-processing class condense_precip
%
% The class changes the time distribution of precipitation (both rain- and
% snowfall) by moving the precipitation from small events to large events.
% 
% It is recommended to compare the resulting precipitation statistics to
% measurements of other data sources.
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef convert2ESA_CCI < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.averaging_period = [];  %in days

            
            post_proc.PARA.emissivity_snow = 0.99; % Snow emissivity (assumed known).
            
            post_proc.PARA.taus=0.0025; % Threshold snowfall for resetting to maximum [m w.e.].
            post_proc.PARA.taua=0.008; % Time constant for snow albedo change in non-melting conditions [/day].
            post_proc.PARA.tauf=0.24; % Time constant for snow albedo change in melting conditions [/day].
            
            post_proc.PARA.albsmax=0.85; % Maximum snow albedo.
            post_proc.PARA.albsmin=0.5; % Minimum snow albedo.
            
            post_proc.PARA.ar=0.2/1e2; % Restricted degree day factor (m*degC/day , value from Burbaker et al. 1996)
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
            post_proc.STATVAR.Lupwelling = post_proc.PARA.emissivity_snow.*post_proc.CONST.sigma.*post_proc.CONST.Tmfw.^4; % upwelling longwave radiation for melting snow, T=273.15K
            post_proc.STATVAR.albedo = post_proc.PARA.albsmax; %initialize 1st albedo values 
        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
            data_full = forcing.DATA;
            forcing.DATA = [];
            forcing.DATA.snowfall = [];
            forcing.DATA.melt = [];
            forcing.DATA.surfT = [];
            forcing.DATA.timeForcing = [];
            
            for i = data_full.timeForcing(1,1):post_proc.PARA.averaging_period:data_full.timeForcing(end,1)-post_proc.PARA.averaging_period
                range = data_full.timeForcing>=i & data_full.timeForcing < min(data_full.timeForcing(end,1), i+ post_proc.PARA.averaging_period);
                forcing.DATA.timeForcing = [forcing.DATA.timeForcing; mean(data_full.timeForcing(range,1))];
                forcing.DATA.surfT = [forcing.DATA.surfT; mean(data_full.Tair(range,1))];
                forcing.DATA.snowfall = [forcing.DATA.snowfall; mean(data_full.snowfall(range,1))];
                
                melt_depth = 0;
                for j = 0:post_proc.PARA.averaging_period-1 %loop over individual days
                    range = data_full.timeForcing>=i+j & data_full.timeForcing < min(data_full.timeForcing(end,1), i+j+1);
                    
                    % Ablation term
                    LW_net = mean(post_proc.PARA.emissivity_snow.*data_full.Lin(range,1) - post_proc.STATVAR.Lupwelling); % Net  longwave
                    SW_net = mean((1-post_proc.STATVAR.albedo).*data_full.Sin(range,1)); % Net shortwave
                    SH_net = mean(post_proc.PARA.ar .* data_full.Tair(range,1)); % Warming through turbulent heat fluxes, parametrized using a restricted degree day approach.
                    
                    daily_melt_depth = (LW_net + SW_net + SH_net)  .* post_proc.CONST.day_sec ./ post_proc.CONST.L_f .* 1000;
                    melt_depth = melt_depth + daily_melt_depth; % Melt depth over the time step.

                    
                    % Update snow albedo for next step.
                    % Latest ECMWF "continuous reset" snow albedo scheme (Dutra et al. 2010)
                    new_snow = mean(data_full.snowfall(range,1)); %in mm/day

                    net_acc = new_snow - max(0,daily_melt_depth); % Net accumulation for one day time-step.
                    if net_acc>0
                        post_proc.STATVAR.albedo = post_proc.STATVAR.albedo + min(1,net_acc./(post_proc.PARA.taus .* 1000)) .* (post_proc.PARA.albsmax - post_proc.STATVAR.albedo);
                    elseif net_acc==0 %"Steady" case (linear decay)
                        post_proc.STATVAR.albedo = post_proc.STATVAR.albedo - post_proc.PARA.taua;
                    else
                        post_proc.STATVAR.albedo = (post_proc.STATVAR.albedo - post_proc.PARA.albsmin) .* exp(-post_proc.PARA.tauf) + post_proc.PARA.albsmin;
                        post_proc.STATVAR.albedo(post_proc.STATVAR.albedo < post_proc.PARA.albsmin) = post_proc.PARA.albsmin;
                    end

                end
                melt_depth(melt_depth <0) = 0;
                forcing.DATA.melt = [forcing.DATA.melt; melt_depth ./ post_proc.PARA.averaging_period];  %in mm/day
                
            end
            
            tile.PARA.geothermal = forcing.PARA.heatFlux_lb;
            
            %overwrite target variables in TEMP in FORCING
            forcing.TEMP = [];
            forcing.TEMP.snowfall=0;
            forcing.TEMP.melt = 0;
            forcing.TEMP.surfT = 0;
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