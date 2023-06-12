classdef ENSEMBLE_MULTITILE_general < matlab.mixin.Copyable


    properties
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function ensemble = provide_PARA(ensemble) 

            ensemble.PARA.gaussian_variables_name= [];
            ensemble.PARA.gaussian_variables_center= [];
            ensemble.PARA.gaussian_variables_width= [];
            
            ensemble.PARA.boxcar_variables_name= [];
            ensemble.PARA.boxcar_variables_center= [];
            ensemble.PARA.boxcar_variables_width= [];
            
            ensemble.PARA.modify_class_name = [];
            ensemble.PARA.modify_class_index = [];
            ensemble.PARA.variable = [];
        end
        
        function ensemble = provide_CONST(ensemble)

        end
        
        function ensemble = provide_STATVAR(ensemble)

        end 
        
        
        function ensemble = finalize_init(ensemble, tile)
            
           
           for i=1:size(ensemble.PARA.gaussian_variables_name,1)
               ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = ensemble.PARA.gaussian_variables_center(i,1) + randn(1,tile.PARA.ensemble_size) .* ensemble.PARA.gaussian_variables_width(i,1);
           end
           for i=1:size(ensemble.PARA.boxcar_variables_name,1)
               ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1}) = ensemble.PARA.boxcar_variables_center(i,1) + (rand(1,tile.PARA.ensemble_size)-0.5) .* ensemble.PARA.boxcar_variables_width(i,1);
           end

            %write variables in PROVIDER
            for i=1:size(ensemble.PARA.modify_class_name,1)
                tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable{i,1}) = ...
                    ensemble.STATVAR.(ensemble.PARA.variable{i,1});
                tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.ensemble_size = tile.PARA.ensemble_size; 
            end
            

            
%             ensemble.STATVAR.stratigraphy_ID_list = [];
%             for i=1:size(tile.STRATIGRAPHY,2)
%                 ensemble.STATVAR.stratigraphy_ID_list = [ensemble.STATVAR.stratigraphy_ID_list; tile.STRATIGRAPHY{1,i}.PARA.ID];
%             end
%             
%             subgrid_index = (ensemble.STATVAR.stratigraphy - floor(ensemble.STATVAR.stratigraphy)).*10;  %zero to one
%             id_list = floor(ensemble.STATVAR.stratigraphy);
            

            
%                         %assign snowfall_factor and landcover class
%             ensemble.STATVAR.snowfall_factor = zeros(tile.PARA.number_of_realizations, tile.PARA.ensemble_size);   %FORCING + TTOP
%             ensemble.STATVAR.stratigraphy = ensemble.STATVAR.snowfall_factor; %reshape later      %STRATIGRAPHY
%             ensemble.STATVAR.wind_compaction_timescale = ensemble.STATVAR.snowfall_factor;    %GROUND
%             ensemble.STATVAR.water_table_depth = ensemble.STATVAR.snowfall_factor;   %GROUND
%             ensemble.STATVAR.wind_speed_class = ensemble.STATVAR.snowfall_factor;   %GROUND
%             ensemble.STATVAR.bare_forest_fraction = ensemble.STATVAR.snowfall_factor;   %FORCING + TTOP
%             ensemble.STATVAR.rk_init = ensemble.STATVAR.snowfall_factor;    %TTOP
            
            
%             variables = fieldnames(tile.STRATIGRAPHY{1,1}.STATVAR);
%             for index = 1:size(variables,1)
%                 tile.GRID.STATVAR.(variables{index,1}) = [];
%                 for i=1:size(id_list,2)
%                     pos = find(ensemble.STATVAR.stratigraphy_ID_list(:,1) == id_list(1,i));
%                     tile.GRID.STATVAR.(variables{index,1}) = [tile.GRID.STATVAR.(variables{index,1}) ...
%                         tile.STRATIGRAPHY{1,pos}.STATVAR.(variables{index,1})(:,1) .* subgrid_index(1,i) + tile.STRATIGRAPHY{1,pos}.STATVAR.(variables{index,1})(:,2) .* (1-subgrid_index(1,i)) ];
%                 end
%             end
            ensemble.STATVAR = [];
            
            
        end
        
        %non-mandatory
        
        function snowfall_factor = get_snowfall_factor(ensemble, percentage, CV)

            av_snow = 1; %100%
            
            psi_square = log(1 + CV.^2);
            lambda = log(av_snow) - 0.5.*psi_square;
            
            snowfall_factor = exp(erfinv(2.*percentage - 1).*sqrt(2.*psi_square) + lambda);
            
        end

    end
end

