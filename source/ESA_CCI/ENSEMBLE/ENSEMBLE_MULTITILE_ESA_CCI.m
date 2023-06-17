classdef ENSEMBLE_MULTITILE_ESA_CCI < matlab.mixin.Copyable


    properties
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function ensemble = provide_PARA(ensemble) 
            ensemble.PARA.ensemble_folder = [];
            ensemble.PARA.write_ensemble_info = [];
            ensemble.PARA.ensemble_information = [];

            ensemble.PARA.modify_class_name = [];
            ensemble.PARA.modify_class_index = [];
            ensemble.PARA.variable = [];
        end
        
        function ensemble = provide_CONST(ensemble)

        end
        
        function ensemble = provide_STATVAR(ensemble)

        end 
        
        
        function ensemble = finalize_init(ensemble, tile)
            
           
            landcover_list = tile.PARA.landcover;
            %assume landcoverlist = ensemble list of base_classes, this
            %could be change din future updates
            landcover_list(sum(landcover_list,2)==0,1) = 1; %cells with 100% lake assigned the open class
            landcover_list = landcover_list ./ repmat(sum(landcover_list,2),1, size(landcover_list,2));  %normalize to 1
            landcover_list = landcover_list .* tile.PARA.ensemble_size;
            ensemble_size_per_class = floor(landcover_list);
            for i = 1:size(ensemble_size_per_class,1)
                while sum(ensemble_size_per_class(i,:),2) < tile.PARA.ensemble_size
                    [maxi, posi] = max(landcover_list(i,:)-ensemble_size_per_class(i,:));
                    ensemble_size_per_class(i,posi) = ensemble_size_per_class(i,posi) + 1;
                end
            end
        
            ensemble.STATVAR.ensemble_size_per_class = ensemble_size_per_class;
            
            %assign snowfall_factor and landcover class
            ensemble.STATVAR.snowfall_factor = zeros(tile.PARA.number_of_realizations, tile.PARA.ensemble_size);   %FORCING + TTOP
            ensemble.STATVAR.stratigraphy = ensemble.STATVAR.snowfall_factor; %reshape later      %STRATIGRAPHY
            ensemble.STATVAR.wind_compaction_timescale = ensemble.STATVAR.snowfall_factor;    %GROUND
            ensemble.STATVAR.water_table_depth = ensemble.STATVAR.snowfall_factor;   %GROUND
            ensemble.STATVAR.wind_speed_class = ensemble.STATVAR.snowfall_factor;   %GROUND
            ensemble.STATVAR.bare_forest_fraction = ensemble.STATVAR.snowfall_factor;   %FORCING + TTOP
            ensemble.STATVAR.rk_init = ensemble.STATVAR.snowfall_factor;    %TTOP
            
            CV_list = ensemble.PARA.ensemble_information.CV';
            stratigraphy_slope = ensemble.PARA.ensemble_information.stratigraphy_slope';
            stratigraphy_intercept = ensemble.PARA.ensemble_information.stratigraphy_intercept';
            water_table_slope = ensemble.PARA.ensemble_information.water_table_slope';
            water_table_intercept = ensemble.PARA.ensemble_information.water_table_intercept';
            wind_compaction_timescale = ensemble.PARA.ensemble_information.wind_compaction_timescale'; %[in days];
            wind_speed_per_class = ensemble.PARA.ensemble_information.wind_speed_per_class';
            bare_forest_fraction  = ensemble.PARA.ensemble_information.bare_forest_fraction';
            rk_init = ensemble.PARA.ensemble_information.rk_init';
            
            for i=1:size(ensemble_size_per_class,1) %loop over cells
                index=0;
                for j=1:size(ensemble_size_per_class,2) %loop over the different classes
                    
                    if j==4 || j==5 || j==9
                        
                        random_list=rand(ensemble_size_per_class(i,j),1);%CHANGED
                        [~, order] = sort(random_list);%CHANGED
                        
                    else
                        
                        order=[1:ensemble_size_per_class(i,j)]';  %REVISED 2019, no randomization of stratigraphies
                    end
                    
                    for k=1:ensemble_size_per_class(i,j) %insert
                        fraction = (0.5.*double(ensemble_size_per_class(i,j) == 1) +  (order(k,1)-1)./(max(1, ensemble_size_per_class(i,j)-1)));
                        stratigraphy_fraction =  max(0, min(1, stratigraphy_intercept (1,j) + fraction .* stratigraphy_slope(1,j)));
                        ensemble.STATVAR.stratigraphy(i, index+k) = j + 0.1.* stratigraphy_fraction;
                        
                        
                        %waterTable_fraction = max(0, min(1, water_table_intercept(1,j) + fraction .* water_table_slope(1,j)));
                        ensemble.STATVAR.water_table_depth(i,index+k) = max(0, min(1, water_table_intercept(1,j) + fraction .* water_table_slope(1,j)));
                        
                        ensemble.STATVAR.snowfall_factor(i,index+k) = get_snowfall_factor(ensemble, k./(ensemble_size_per_class(i,j)+1) , CV_list(1,j));
                        ensemble.STATVAR.wind_compaction_timescale(i, index+k) = wind_compaction_timescale(1,j);
                        ensemble.STATVAR.wind_speed_class(i, index+k) = wind_speed_per_class(1,j);
                        ensemble.STATVAR.bare_forest_fraction(i, index+k) = bare_forest_fraction(1,j);
                        ensemble.STATVAR.rk_init(i, index+k) = rk_init(1,j);
                    end
                    index = index + ensemble_size_per_class(i,j);
                end
            end
            ensemble.STATVAR.stratigraphy = ensemble.STATVAR.stratigraphy(:)';
            ensemble.STATVAR.snowfall_factor = ensemble.STATVAR.snowfall_factor(:)';
            ensemble.STATVAR.wind_compaction_timescale = ensemble.STATVAR.wind_compaction_timescale(:)';
            ensemble.STATVAR.water_table_depth = ensemble.STATVAR.water_table_depth(:)';
            ensemble.STATVAR.wind_speed_class = ensemble.STATVAR.wind_speed_class(:)';
            ensemble.STATVAR.bare_forest_fraction = ensemble.STATVAR.bare_forest_fraction(:)';
            ensemble.STATVAR.rk_init = ensemble.STATVAR.rk_init(:)';
            ensemble.STATVAR.stratigraphy_index = floor(ensemble.STATVAR.stratigraphy);
            ensemble.STATVAR.stratigraphy_fraction = (ensemble.STATVAR.stratigraphy - floor(ensemble.STATVAR.stratigraphy)).*10;
           
            
            %write variables in PROVIDER
            for i=1:size(ensemble.PARA.modify_class_name,1)
                tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable{i,1}) = ...
                    ensemble.STATVAR.(ensemble.PARA.variable{i,1});
            end
            tile.PARA.geothermal = repmat(tile.PARA.geothermal', 1, tile.PARA.ensemble_size);

            if ensemble.PARA.write_ensemble_info
                ensemble.PARA.ensemble_folder = [ensemble.PARA.ensemble_folder tile.RUN_INFO.PARA.run_name '/'];
                if ~(exist(ensemble.PARA.ensemble_folder)==7)
                    mkdir(ensemble.PARA.ensemble_folder);
                end
                ensemble_info = ensemble.STATVAR;
                ensemble_info.geothermal = tile.PARA.geothermal;
                save([ensemble.PARA.ensemble_folder 'ensemble_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.mat'], 'ensemble_info')
                
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

