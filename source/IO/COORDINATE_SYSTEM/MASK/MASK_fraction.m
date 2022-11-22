classdef MASK_fraction < matlab.mixin.Copyable

    
    properties
        PARENT
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function mask = provide_PARA(mask)
            mask.PARA.section = [];
            mask.PARA.number_of_sections = [];
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
                range = round((mask.PARA.section -1) .* number_of_existing./mask.PARA.number_of_sections)+1:round(mask.PARA.section .* number_of_existing./mask.PARA.number_of_sections);
                yes = logical(zeros(number_of_existing,1));
                yes(range) = 1;
                index(~yes) = [];
       
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask .*0;
                mask.PARENT.STATVAR.mask(index) = 1;

        end
        
 
    end
end

