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

classdef remove_melt_event < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.start_year = [];  
            post_proc.PARA.end_year = [];
            post_proc.PARA.start_month = [];
            post_proc.PARA.end_month = [];
            post_proc.PARA.max_airT = [];
            post_proc.PARA.all_rain_T = [];
            post_proc.PARA.all_snow_T = [];

        end
        
        
        function post_proc = provide_CONST(post_proc)

        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)

        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
            range = find(forcing.DATA.timeForcing >= datenum(post_proc.PARA.start_year, post_proc.PARA.start_month, 1) & ...
                forcing.DATA.timeForcing < datenum(post_proc.PARA.end_year, post_proc.PARA.end_month+1, 1) & ...
                forcing.DATA.Tair > post_proc.PARA.max_airT);
            forcing.DATA.Tair(range,1) = post_proc.PARA.max_airT;
            precip = forcing.DATA.snowfall(range,1) + forcing.DATA.rainfall(range,1);
            forcing.DATA.snowfall(range,1) = precip .* (double(forcing.DATA.Tair(range,1) <= post_proc.PARA.all_snow_T)  + ...
                double(forcing.DATA.Tair(range,1) > post_proc.PARA.all_snow_T & forcing.DATA.Tair(range,1) < post_proc.PARA.all_rain_T) .* ...
                (1- (forcing.DATA.Tair(range,1) - post_proc.PARA.all_snow_T) ./ max(1e-12, (post_proc.PARA.all_rain_T - post_proc.PARA.all_snow_T))));
            forcing.DATA.rainfall(range,1) = precip .* (double(forcing.DATA.Tair(range,1) >= post_proc.PARA.all_rain_T)  + ...
                double(forcing.DATA.Tair(range,1) > post_proc.PARA.all_snow_T & forcing.DATA.Tair(range,1) < post_proc.PARA.all_rain_T) .* ...
                (forcing.DATA.Tair(range,1) - post_proc.PARA.all_snow_T) ./ max(1e-12, (post_proc.PARA.all_rain_T - post_proc.PARA.all_snow_T)));
        
        end
        
        
                %-------------param file generation-----
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