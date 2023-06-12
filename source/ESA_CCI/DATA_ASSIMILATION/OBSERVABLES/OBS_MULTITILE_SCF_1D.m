classdef OBS_MULTITILE_SCF_1D < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        STATVAR
    end
    
    methods

        function obs = provide_PARA(obs)
          end
        
        function obs = provide_CONST(obs)
            
        end
        
        function obs = provide_STATVAR(obs)
            
        end
        
        function obs = finalize_init(obs, tile)
        end
        
        function result = observable_operator(obs, tile) 

            result = double(sum(tile.SUBSURFACE_CLASS.STATVAR.waterIce_snow,1)>1e-5);

        end
        
    end
end

