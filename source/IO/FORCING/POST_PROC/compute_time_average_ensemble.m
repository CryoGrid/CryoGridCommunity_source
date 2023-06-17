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

classdef compute_time_average_ensemble < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.averaging_period = [];  %in days
            post_proc.PARA.all_snow_T = [];
            post_proc.PARA.all_rain_T = [];
            
            %these can be written by ensemble class
            post_proc.PARA.ensemble_size = 1;
            post_proc.PARA.absolute_change_Tair = 0;
            post_proc.PARA.snow_fraction = 1;
            post_proc.PARA.rain_fraction = 1;
            post_proc.PARA.relative_change_Sin = 1;                 
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
            
            if size(post_proc.PARA.absolute_change_Tair,2)==1 
                post_proc.PARA.absolute_change_Tair = repmat(post_proc.PARA.absolute_change_Tair,1, post_proc.PARA.ensemble_size);
            end   
            if size(post_proc.PARA.snow_fraction,2)==1 
                post_proc.PARA.snow_fraction = repmat(post_proc.PARA.snow_fraction,1, post_proc.PARA.ensemble_size);
            end
            if size(post_proc.PARA.rain_fraction,2)==1
                post_proc.PARA.rain_fraction = repmat(post_proc.PARA.rain_fraction,1, post_proc.PARA.ensemble_size);
            end
            if size(post_proc.PARA.relative_change_Sin,2)==1
                post_proc.PARA.relative_change_Sin = repmat(post_proc.PARA.relative_change_Sin,1, post_proc.PARA.ensemble_size);
            end
        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
            data_full = forcing.DATA;
            forcing.DATA = [];
            forcing.DATA.snowfall = [];
            forcing.DATA.rainfall = [];
            forcing.DATA.Tair = [];
            forcing.DATA.wind = [];
            forcing.DATA.q = [];
            forcing.DATA.p = [];
            forcing.DATA.Lin = [];
            forcing.DATA.Sin = [];
            forcing.DATA.timeForcing = [];
            
            for i = data_full.timeForcing(1,1):post_proc.PARA.averaging_period:data_full.timeForcing(end,1)-post_proc.PARA.averaging_period
                range = find(data_full.timeForcing>=i & data_full.timeForcing < min(data_full.timeForcing(end,1), i + post_proc.PARA.averaging_period));
                forcing.DATA.timeForcing = [forcing.DATA.timeForcing; mean(data_full.timeForcing(range,1))];
                forcing.DATA.Tair = [forcing.DATA.Tair; mean(data_full.Tair(range,1)) + post_proc.PARA.absolute_change_Tair];
                forcing.DATA.wind = [forcing.DATA.wind; mean(data_full.wind(range,1))];
                forcing.DATA.q = [forcing.DATA.q; mean(data_full.q(range,1))];
                forcing.DATA.p = [forcing.DATA.p; mean(data_full.p(range,1))];

                sky_emissivity = data_full.Lin(range,1) ./ (data_full.Tair(range,1)+273.15).^4 ./ post_proc.CONST.sigma;
                forcing.DATA.Lin = [forcing.DATA.Lin; mean(sky_emissivity .* post_proc.CONST.sigma .* (data_full.Tair(range,1) + 273.15 + post_proc.PARA.absolute_change_Tair).^4)];
                forcing.DATA.Sin = [forcing.DATA.Sin; mean(data_full.Sin(range,1) .*  (1+post_proc.PARA.relative_change_Sin))];
                
                sf = 0;
                rf = 0;
                for j=1:size(range,1)
                    precip = data_full.snowfall(range(j),1) + data_full.rainfall(range(j),1);
                    factor = max(0, min(1, (data_full.Tair(range(j),1) + post_proc.PARA.absolute_change_Tair - post_proc.PARA.all_snow_T) ./ max(1e-12, (post_proc.PARA.all_rain_T - post_proc.PARA.all_snow_T))));
                    sf = sf + precip.*(1 - factor);
                    rf = rf + precip.*factor;
                end
                forcing.DATA.snowfall = [forcing.DATA.snowfall; sf./size(range,1) .* post_proc.PARA.snow_fraction];
                forcing.DATA.rainfall = [forcing.DATA.rainfall; rf./size(range,1) .* post_proc.PARA.rain_fraction];
                                
            end
            
            if size(forcing.PARA.heatFlux_lb,2) == 1
                tile.PARA.geothermal = repmat(forcing.PARA.heatFlux_lb, 1, post_proc.PARA.ensemble_size);
            end
            
            %overwrite target variables in TEMP in FORCING
            forcing.TEMP = [];
            forcing.TEMP.Lin = 0;
            forcing.TEMP.Sin = 0;
            forcing.TEMP.q = 0;
            forcing.TEMP.p = 0;
            forcing.TEMP.Tair = 0;
            forcing.TEMP.wind = 0;
            forcing.TEMP.rainfall = 0;
            forcing.TEMP.snowfall = 0;
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