%========================================================================
% CryoGrid FORCING post-processing 
%
%
% Authors:
% S. Westermann, January 2023
%
%========================================================================

classdef time_average_temperature < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.averaging_period = [];  
            post_proc.PARA.year_range = [];
            
            post_proc.PARA.annual = [];
            
        end
        
        
        function post_proc = provide_CONST(post_proc)

        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)
            post_proc.TEMP.offset_years = 0;
        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
            if forcing.TEMP.current_year >= post_proc.PARA.year_range(1) && forcing.TEMP.current_year <= post_proc.PARA.year_range(end)
                pos_year = forcing.TEMP.current_year - (forcing.PARA.ERA_data_years(1) - forcing.PARA.number_of_spin_up_years) + 1 - post_proc.TEMP.offset_years;
                for i=1:floor(365./post_proc.PARA.averaging_period)+1
                    forcing.DATA.ERA_T_downscaled(:,i,pos_year) = nanmean(forcing.TEMP.ERA_T_downscaled(:, (i-1).*post_proc.PARA.averaging_period.*4+1:min(i*post_proc.PARA.averaging_period*4, size(forcing.TEMP.ERA_T_downscaled,2))),2)-273.15;
                end
            end
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