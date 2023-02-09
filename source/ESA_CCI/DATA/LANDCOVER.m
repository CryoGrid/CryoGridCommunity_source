classdef LANDCOVER < matlab.mixin.Copyable

    properties
        PARENT
        CONST
        STATVAR
        PARA
    end

    
    methods
        
        function lc = provide_PARA(lc)
            lc.PARA.landcover_file = [];
            lc.PARA.landcover_path = [];
        end
        
        function lc = provide_STATVAR(lc)

        end
        
        function lc = provide_CONST(lc)
            
        end
        
        function lc = finalize_init(lc)
            
        end
        
        function lc = load_data(lc)

            lc.PARENT.STATVAR.landcover = [lc.PARENT.STATVAR.latitude .* 0 + 1 repmat(lc.PARENT.STATVAR.latitude .* 0, 1,9)]; 
            
       end


    end
end

