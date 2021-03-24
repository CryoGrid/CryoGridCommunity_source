classdef MASK_IPA < matlab.mixin.Copyable

    
    properties
        PARENT
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function mask = provide_PARA(mask)
            mask.PARA.IPA_map_path = [];
            mask.PARA.filename = [];
            mask.PARA.additive = [];
        end

        function mask = provide_STATVAR(mask)

        end
        
        function mask = provide_CONST(mask)
            
        end
        
        function mask = finalize_init(mask)
            
        end
        

        function mask = apply_mask(mask)
            load([mask.PARA.IPA_map_path mask.PARA.filename])
            a=struct2cell(IPA_mask);
            h_index=squeeze(cell2mat(a(1,1,:)));
            v_index=squeeze(cell2mat(a(2,1,:)));
            mask_temp = IPA_mask(find(mask.PARENT.PARA.horizontal==h_index(:,1) & mask.PARENT.PARA.vertical==v_index(:,1))).mask;

            if mask.PARA.additive
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask | mask_temp;
            else
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask & mask_temp;
            end
            
        end
            
            
    end
end

