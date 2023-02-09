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

classdef bias_correct_BASE < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
        end
        
        
        function post_proc = provide_CONST(post_proc)

        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)

        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
        end
        
        function overlap_pairs = find_overlap_pairs(downscale, forcing_time, forcing, reference_time, reference)
            overlap_pairs=[];
            start_search = find(forcing_time(:,1) - reference_time(1,1) < 0, 1,'last'); %last entry of forcing for which first entry of ref time is larger
            if isempty(start_search)
                start_search = 1;
            end
            end_search = find(forcing_time(:,1) - reference_time(end,1) > 0, 1); %first entry of foricng for which last entry of reference time is smaller
            if isempty(end_search)
                end_search = size(forcing_time,1);
            end
            for i = start_search:end_search %376937:403233
                index=find(abs(reference_time(:,1)-forcing_time(i,1)) < downscale.PARA.tolerance);
                if ~isempty(index)
                    doy = floor(forcing_time(i,1) - datenum(year(forcing_time(i,1)), 1, 1))+1;
                    overlap_pairs=[overlap_pairs; [doy reference(index(1),1) forcing(i,1)]];
                end
            end
            overlap_pairs(isnan(overlap_pairs(:,2)) | isnan(overlap_pairs(:,3)) | abs(overlap_pairs(:,2)-overlap_pairs(:,3)) > downscale.PARA.invalid_threshold, :) = [];
        end
        

        function forcing_corrected = linear_regression_doy(downscale, overlap_pairs, forcing_time, forcing)
            
            forcing_corrected = forcing;
            slope = [];
            intercept = [];
            
            for i=1:365
                doy_corrected = overlap_pairs(:,1);
                doy_corrected(doy_corrected > i+downscale.PARA.window+1) = doy_corrected(doy_corrected > i+downscale.PARA.window+1) - 365;
                doy_corrected(doy_corrected < i-downscale.PARA.window-1) = doy_corrected(doy_corrected < i-downscale.PARA.window-1) + 365;
                range = find(doy_corrected >= i-downscale.PARA.window & doy_corrected <= i+downscale.PARA.window);
                
                
                P = polyfit(overlap_pairs(range,3), overlap_pairs(range,2), 1);
                
%                 P(1) = min(max(P(1), 0.5), 2);
                
%                 disp([i  P(1) P(2)])
                
                range = find(floor(forcing_time - datenum(year(forcing_time), 1, 1))+1 == i);
                forcing_corrected(range,1) = P(2) + P(1) .* forcing_corrected(range,1);
            end
            i = 366;
            range = find(floor(forcing_time - datenum(year(forcing_time), 1, 1))+1 == i); %same as i=365
            forcing_corrected(range,1) = P(2) + P(1) .* forcing_corrected(range,1);
        end
        
        
        function forcing_corrected = linear_regression_zero_intercept_doy(downscale, overlap_pairs, forcing_time, forcing)
            
            warning off all
            
            overlap_pairs(overlap_pairs(:,2)<200 | overlap_pairs(:,3)<200, :) = [];
            
            forcing_corrected = forcing;
            slope = [];
            intercept = [];
            
            for i=1:365
                doy_corrected = overlap_pairs(:,1);
                doy_corrected(doy_corrected > i+downscale.PARA.window+1) = doy_corrected(doy_corrected > i+downscale.PARA.window+1) - 365;
                doy_corrected(doy_corrected < i-downscale.PARA.window-1) = doy_corrected(doy_corrected < i-downscale.PARA.window-1) + 365;
                range = find(doy_corrected >= i-downscale.PARA.window & doy_corrected <= i+downscale.PARA.window);

                if length(range)>20
                    
                    P = nlinfit(overlap_pairs(range,3), overlap_pairs(range,2), @(a,x) a.*x, 1);
                    P = min(max(P, 0.5), 2);
                    
                else
                    P=1;
                end
                
                range = find(floor(forcing_time - datenum(year(forcing_time), 1, 1))+1 == i & forcing>200);
%                 forcing_corrected(range,1) = P(2) + P(1) .* forcing_corrected(range,1);
                forcing_corrected(range,1) = P .* forcing_corrected(range,1);
            end
            i = 366;
            range = find(floor(forcing_time - datenum(year(forcing_time), 1, 1))+1 == i & forcing>200); %same as i=365
%             forcing_corrected(range,1) = P(2) + P(1) .* forcing_corrected(range,1);
            forcing_corrected(range,1) =  P .* forcing_corrected(range,1);
        end
        
        function q_new = correct_q_constant_RH(downscale, q_old, forcing_T_old, forcing_T_new)
            q_new = q_old;
            forcing_T_old = forcing_T_old + 273.15;
            forcing_T_new = forcing_T_new + 273.15;
            
            range = find(forcing_T_old>=273.15);
            q_new(range) = q_new(range) .* exp(17.62.*(forcing_T_new(range)-273.15)./(243.12-273.15+forcing_T_new(range))) ./ exp(17.62.*(forcing_T_old(range)-273.15)./(243.12-273.15+forcing_T_old(range)));
            range = find(forcing_T_old<273.15);
            q_new(range) = q_new(range) .*  exp(22.46.*(forcing_T_new(range)-273.15)./(272.61-273.15+forcing_T_new(range))) ./ exp(22.46.*(forcing_T_old(range)-273.15)./(272.61-273.15+forcing_T_old(range)));

        end
        
        function [q_new, RH] = correct_q_constant_RH_dewpointT(downscale, forcing_Td_old, forcing_T_old, forcing_T_new, forcing_p)
            q_new = forcing_Td_old.*0;
            RH = forcing_Td_old.*0;
            forcing_T_old = forcing_T_old + 273.15;
            forcing_T_new = forcing_T_new + 273.15;
            forcing_Td_old = forcing_Td_old + 273.15;
            
            range = find(forcing_Td_old>=273.15);
            RH(range) = exp(17.62.*(forcing_Td_old(range)-273.15)./(243.12-273.15+forcing_Td_old(range))) ./ exp(17.62.*(forcing_T_old(range)-273.15)./(243.12-273.15+forcing_T_old(range)));
            q_new(range) = RH(range) .*0.622.* 6.112 .* 100 .* exp(17.62.*(forcing_T_new(range)-273.15)./(243.12-273.15+forcing_T_new(range))) ./ forcing_p(range);
            range = find(forcing_Td_old<273.5);
            RH(range) = exp(22.46.*(forcing_Td_old(range)-273.15)./(272.61-273.15+forcing_Td_old(range))) ./ exp(22.46.*(forcing_T_old(range)-273.15)./(272.61-273.15+forcing_T_old(range)));
            q_new(range) = RH(range) .* 0.622.* 6.112 .* 100 .* exp(22.46.*(forcing_T_new(range)-273.15)./(272.61-273.15+forcing_T_new(range)))./ forcing_p(range);

        end
        
        
    end
    
end