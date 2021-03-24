
classdef PREPROC_ESA_CCI_ERA_ONLY  < PREPROC_ESA_CCI
    

    
    
    methods
        
        %-----initialize-----------------
        
        function preproc = provide_PARA(preproc)
            preproc.PARA.timestep = []; %[days]
            preproc.PARA.threshold_T_snowfall = []; %[degree C]
            preproc.PARA.threshold_T_snowmelt = [];

        end
        
        

        
        function preproc = finalize_init(preproc, tile)

        end
        
        %-----mandatory functions------------------------
  
        function preproc = get_boundary_condition_u(preproc, tile) %get_MODIS
            
        
        end
      
        function preproc = get_derivatives_prognostic(preproc, tile)  %calculate snowfall and melt
                      
            preproc = get_derivatives_prognostic@PREPROC_ESA_CCI(preproc, tile);
            
            preproc.STATVAR.final_av_T = preproc.STATVAR.ERA_T_downscaled;
            preproc.STATVAR.final_MODIS_weight = preproc.STATVAR.ERA_T_downscaled.*0;
        end
        
        function preproc = compute_diagnostic(preproc, tile)  
            
        end
        
        
        
        
    end
    
end

