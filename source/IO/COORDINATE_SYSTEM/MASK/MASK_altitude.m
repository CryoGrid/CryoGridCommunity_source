classdef MASK_altitude < matlab.mixin.Copyable

    
    properties
        PARENT
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function mask = provide_PARA(mask)
            mask.PARA.min_altitude = [];
            mask.PARA.max_altitude = [];
        end

        function mask = provide_STATVAR(mask)

        end
        
        function mask = provide_CONST(mask)
            
        end
        
        function mask = finalize_init(mask)
            
        end
        

        function mask = apply_mask(mask)

            mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask & mask.PARENT.STATVAR.altitude <= mask.PARA.max_altitude & mask.PARENT.STATVAR.altitude >= mask.PARA.min_altitude;

        end
     
            
    end
end

