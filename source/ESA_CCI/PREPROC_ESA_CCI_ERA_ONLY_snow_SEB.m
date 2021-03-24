
classdef PREPROC_ESA_CCI_ERA_ONLY_snow_SEB  < PREPROC_ESA_CCI_snow_SEB
    

    
    
    methods
        
        %-----initialize-----------------
        
%         function preproc = provide_PARA(preproc)
%             preproc.PARA.timestep = []; %[days]
%             preproc.PARA.threshold_T_snowfall = []; %[degree C]
%             preproc.PARA.threshold_T_snowmelt = [];
% 
%         end
        
        

        
        function preproc = finalize_init(preproc, tile)

            preproc.PARA.Lupwelling = preproc.PARA.emissivity_snow.*preproc.CONST.sigma.*preproc.CONST.Tmfw.^4; % upwelling longwave radiation for T=273.15K
                        
            preproc.STATVAR.albedo_bare = repmat((preproc.PARA.albsmax_bare + preproc.PARA.albsmin_bare)./2, size(tile.RUN_INFO.STATVAR.key,1), 1);
            preproc.STATVAR.albedo_forest = repmat((preproc.PARA.albsmax_forest + preproc.PARA.albsmin_forest)./2, size(tile.RUN_INFO.STATVAR.key,1), 1);
            
        end
        
        %-----mandatory functions------------------------
  
        function preproc = get_boundary_condition_u(preproc, tile) %get_MODIS
            
        
        end
      
        function preproc = get_derivatives_prognostic(preproc, tile)  %calculate snowfall and melt
                      
            preproc = get_derivatives_prognostic@PREPROC_ESA_CCI_snow_SEB(preproc, tile);
            
            preproc.STATVAR.final_av_T = preproc.STATVAR.ERA_T_downscaled;
            preproc.STATVAR.final_MODIS_weight = preproc.STATVAR.ERA_T_downscaled.*0;
        end
        
        function preproc = compute_diagnostic(preproc, tile)  
            
        end
        
        
        
        
    end
    
end

