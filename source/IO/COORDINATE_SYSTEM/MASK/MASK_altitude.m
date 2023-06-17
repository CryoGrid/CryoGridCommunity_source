%========================================================================
% CryoGrid DATA_MASK class MASK_altitude
% selects the region of interest as the target points in an altitudinal
% range
%
% S. Westermann, Dec 2022
%========================================================================

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
        
        
        
        %-------------param file generation-----
        function mask = param_file_info(mask)
            mask = provide_PARA(mask);
            
            mask.PARA.STATVAR = [];
            mask.PARA.class_category = 'DATA_MASK';
            mask.PARA.default_value = [];
            mask.PARA.options = [];
            
            mask.PARA.comment.min_altitude = {'minimum altitude of target points'};
            
            mask.PARA.comment.max_altitude = {'maximum altitude of target points'};
        end
     
            
    end
end

