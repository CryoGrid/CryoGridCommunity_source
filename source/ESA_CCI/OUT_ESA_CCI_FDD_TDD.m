classdef OUT_ESA_CCI_FDD_TDD < matlab.mixin.Copyable


    properties
        CONST
        PARA
        STATVAR
        TEMP

    end
    
    methods
        
        
        function out = provide_PARA(out)

        end
        
        function out = provide_CONST(out)
        end
        
        function out = provide_STATVAR(out)

        end
        
        
        function out = finalize_init(out, tile)
                      
            out.TEMP.FDD=0;
            out.TEMP.TDD=0;
            out.TEMP.N=0;
            

        end
        
        function out = store_OUT(out, tile) 
            
            
            out.TEMP.FDD = out.TEMP.FDD + tile.SUBSURFACE_CLASS.STATVAR.T(5,:).*double(tile.SUBSURFACE_CLASS.STATVAR.T(5,:)<=0); 
            out.TEMP.TDD = out.TEMP.TDD + tile.SUBSURFACE_CLASS.STATVAR.T(5,:).*double(tile.SUBSURFACE_CLASS.STATVAR.T(5,:)>0);
            out.TEMP.N = out.TEMP.N + 1;

        end
        


    end
end

