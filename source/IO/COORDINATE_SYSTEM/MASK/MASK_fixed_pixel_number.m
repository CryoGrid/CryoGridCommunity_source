%========================================================================
% CryoGrid MASK class MASK_fraction
% splits the target points in a number of sections. This provides an 
% additional level of parallelization for very large data sets, which can
% be useful on clusters with defined walltime.
%
% S. Westermann, Dec 2022
%========================================================================

classdef MASK_fixed_pixel_number < matlab.mixin.Copyable

    
    properties
        PARENT
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function mask = provide_PARA(mask)
            mask.PARA.number_of_pixels = [];
        end

        function mask = provide_STATVAR(mask)

        end
        
        function mask = provide_CONST(mask)
            
        end
        
        function mask = finalize_init(mask)
            
        end
        

        function mask = apply_mask(mask)
                
                index = find(mask.PARENT.STATVAR.mask==1);
                
                number_of_existing = size(index,1);
                yes = randperm(number_of_existing, min(mask.PARA.number_of_pixels, number_of_existing));
                index = index(yes,1);
       
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask .*0;
                mask.PARENT.STATVAR.mask(index) = 1;

        end
        
        
        
        %-------------param file generation-----
        function mask = param_file_info(mask)
            mask = provide_PARA(mask);
            
            mask.PARA.STATVAR = [];
            mask.PARA.class_category = 'MASK';
            mask.PARA.default_value = [];
            mask.PARA.options = [];
            
            mask.PARA.comment.section = {'number of currently modeled section [1 ... number_of_sections]'};
            
            mask.PARA.comment.number_of_sections = {'number of sections in which the data set is split'};
        end
 
    end
end

