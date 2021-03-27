classdef MASK_GlobPermafrost < matlab.mixin.Copyable

    
    properties
        PARENT
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function mask = provide_PARA(mask)
            mask.PARA.GlobPermafrost_path = [];
            mask.PARA.filename = [];
            mask.PARA.threshold_T = [];
            mask.PARA.fraction_modeled = [];
            mask.PARA.additive = [];
        end

        function mask = provide_STATVAR(mask)

        end
        
        function mask = provide_CONST(mask)
            
        end
        
        function mask = finalize_init(mask)
            
        end
        

        function mask = apply_mask(mask)
            h = mask.PARENT.PARA.horizontal;
            h=h+100;
            h=num2str(h);
            h=h(2:3);
            v = mask.PARENT.PARA.vertical;
            v=v+100;
            v=num2str(v);
            v=v(2:3);
            mask.PARA.filename(11:12) = h;
            mask.PARA.filename(14:15) = v;
            TTOP = dlmread([mask.PARA.GlobPermafrost_path mask.PARA.filename]);
            
            mask_temp = mask.PARENT.STATVAR.mask .* 0;
            for i=size(mask.PARA.fraction_modeled,1):-1:1
                a = min(1,round(rand(size(mask.PARENT.STATVAR.mask,1),size(mask.PARENT.STATVAR.mask,2)) .* 0.5 ./(1-mask.PARA.fraction_modeled(i,1))));
                b = TTOP<=mask.PARA.threshold_T(i,1);
                mask_temp = mask_temp | (a & b);
            end
            
            if mask.PARA.additive
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask | mask_temp;
            else
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask & mask_temp;
            end
            
        end
            
            
    end
end

