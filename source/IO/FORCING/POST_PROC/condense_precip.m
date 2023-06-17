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

classdef condense_precip < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.window_size = [];  
            post_proc.PARA.eliminate_fraction = [];
            post_proc.PARA.survive_fraction = [];
        end
        
        
        function post_proc = provide_CONST(post_proc)

        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)
            post_proc.PARA.survive_fraction = min(post_proc.PARA.survive_fraction, 1-post_proc.PARA.eliminate_fraction);
        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
            window_size =  floor(post_proc.PARA.window_size ./ (forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1))); 
            vars={'rainfall'; 'snowfall'};
            for k=1:size(vars,1)
                precip = forcing.DATA.(vars{k,1});
                precip2=precip.*0;
                for i=1:window_size:size(precip,1)
                    [a,b]= sort(precip(i:min(length(precip),i+window_size-1)));
                    a(a==0)=[];
                    for j=1:floor(post_proc.PARA.eliminate_fraction .* size(a,1))
                        ind = floor(rand().*size(a,1).* post_proc.PARA.survive_fraction);
                        a(end-ind)=a(end-ind)+ a(j);
                        a(j) = 0;
                    end
                    a=[zeros(size(b,1)-size(a,1),1); a];
                    
                    c=a*0;
                    c(b) = a;
                    precip2(i:min(length(precip),i+window_size-1)) = c;
                end
                forcing.DATA.(vars{k,1}) = precip2;
            end
        end
        
        
                %-------------param file generation-----
        function post_proc = param_file_info(post_proc)
            post_proc = provide_PARA(post_proc);

            post_proc.PARA.STATVAR = [];
            post_proc.PARA.class_category = 'FORCING POST_PROCESSING';
            post_proc.PARA.options = [];
            
            post_proc.PARA.eliminate_fraction = [];
            post_proc.PARA.survive_fraction = [];
                        
            post_proc.PARA.default_value.window_size = {7};
            post_proc.PARA.comment.window_size = {'window size in days within which precipitation is reallocated'};
            
            post_proc.PARA.default_value.eliminate_fraction = {0.5};
            post_proc.PARA.comment.eliminate_fraction = {'fraction of smallest precipitation events (= timestamps with precipitation) that is reallocated to larger events'};
            
            post_proc.PARA.default_value.survive_fraction = {0.5};  
            post_proc.PARA.comment.survive_fraction = {'fraction of largest precipitation events (= timestamps with precipitation) that the small events are reallocated to'};
            
        end
        
    end
    
end