%========================================================================
% CryoGrid FORCING post-processing 
%
%
% Authors:
% S. Westermann, January 2023
%
%========================================================================

classdef statistical_downscaling_MODIS < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.overlap_years = [];  
            post_proc.PARA.downscaling_years = [];
            
            post_proc.PARA.annual = [];
            
        end
        
        
        function post_proc = provide_CONST(post_proc)

        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)

        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
            %statistical downscaling
            
            sum_xy = zeros(size(forcing.DATA.final_av_T,1), size(forcing.DATA.final_av_T,2));
            sum_xx = zeros(size(forcing.DATA.final_av_T,1), size(forcing.DATA.final_av_T,2));
            sum_x = zeros(size(forcing.DATA.final_av_T,1), size(forcing.DATA.final_av_T,2));
            sum_y = zeros(size(forcing.DATA.final_av_T,1), size(forcing.DATA.final_av_T,2));
            n = zeros(size(forcing.DATA.final_av_T,1), size(forcing.DATA.final_av_T,2));
            fit_window = 5; % 5*8 = 40 days
            
            
            for ii=post_proc.PARA.overlap_years(1)-forcing.DATA.years(1)+1:post_proc.PARA.overlap_years(2)-forcing.DATA.years(1)+1 %loop over the years
                for i=1:46
                    for jj=1:fit_window  % fit routine ERA vs T_final here
                        j = mod(i+jj - 4,46)+1;
                        sum_xy(:,j) = sum_xy(:,j) + forcing.DATA.final_av_T(:, i, ii) .* forcing.DATA.ERA_T_downscaled(:, i, ii);
                        sum_xx(:,j) = sum_xx(:,j) + forcing.DATA.ERA_T_downscaled(:, i, ii).^2;
                        sum_x(:,j) = sum_x(:,j) + forcing.DATA.ERA_T_downscaled(:, i, ii);
                        sum_y(:,j) = sum_y(:,j) + forcing.DATA.final_av_T(:,i, ii);
                        n(:,j) =  n(:,j) + 1;
                    end
                end
            end
            
            slope = (sum_xy - sum_x .* sum_y ./ n) ./ ( sum_xx - sum_x.^2 ./ n);
            intercept = sum_y ./ n - slope .* sum_x ./ n;
            
            forcing.DATA.slope = slope;
            forcing.DATA.intercept = intercept;
            
            for ii=post_proc.PARA.downscaling_years(1)-forcing.DATA.years(1)+1:post_proc.PARA.downscaling_years(2)-forcing.DATA.years(1)+1 %loop over the years
                for i=1:46
                    forcing.DATA.final_av_T(:,i,ii) = intercept(:,i) + slope(:,i) .* forcing.DATA.ERA_T_downscaled(:,i,ii);
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