classdef ENSEMBLE_ESA_CCI < matlab.mixin.Copyable


    properties
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function ensemble = provide_PARA(ensemble) 

            ensemble.PARA.ensemble_size = [];
            ensemble.PARA.ensemble_information = [];

        end
        
        function ensemble = provide_CONST(ensemble)

        end
        
        function ensemble = provide_STATVAR(ensemble)

        end 
        
        
        function ensemble = finalize_init(ensemble, tile)
            
            landcover_list = tile.RUN_INFO.STATVAR.landcover;
            %assume landcoverlist = ensemble list of base_classes, this
            %could be change din future updates
            landcover_list(sum(landcover_list,2)==0,1) = 1; %cells with 100% lake assigned the open class
            landcover_list = landcover_list ./ repmat(sum(landcover_list,2),1, size(landcover_list,2));  %normalize to 1
            landcover_list = landcover_list .* ensemble.PARA.ensemble_size;
            ensemble_size_per_class = floor(landcover_list);
            for i = 1:size(ensemble_size_per_class,1)
                while sum(ensemble_size_per_class(i,:),2) < ensemble.PARA.ensemble_size
                    [maxi, posi] = max(landcover_list(i,:)-ensemble_size_per_class(i,:));
                    ensemble_size_per_class(i,posi) = ensemble_size_per_class(i,posi) + 1;
                end
            end
        
            ensemble.STATVAR.ensemble_size_per_class = ensemble_size_per_class;
            
            %assign snowfall_factor and landcover class
            ensemble.STATVAR.snowfall_factor = zeros(tile.PARA.number_of_realizations, ensemble.PARA.ensemble_size);
            ensemble.STATVAR.stratigraphy = ensemble.STATVAR.snowfall_factor; %reshape later
            ensemble.STATVAR.wind_compaction_timescale = ensemble.STATVAR.snowfall_factor;
            ensemble.STATVAR.water_table_depth = ensemble.STATVAR.snowfall_factor;
            ensemble.STATVAR.wind_speed_class = ensemble.STATVAR.snowfall_factor;
            
            
            CV_list = ensemble.PARA.ensemble_information.CV';
            stratigraphy_slope = ensemble.PARA.ensemble_information.stratigraphy_slope';
            stratigraphy_intercept = ensemble.PARA.ensemble_information.stratigraphy_intercept';
            water_table_slope = ensemble.PARA.ensemble_information.water_table_slope';
            water_table_intercept = ensemble.PARA.ensemble_information.water_table_intercept';
            wind_compaction_timescale = ensemble.PARA.ensemble_information.wind_compaction_timescale' .*24.*3600; %[in seconds];
            wind_speed_per_class = ensemble.PARA.ensemble_information.wind_speed_per_class';
            
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
                        
                        %PROFILE.stratigraphy(i, index+k) = j + 0.1.* (0.5.*double(ensemble_size_per_class(i,j) == 1) +  (order(k,1)-1)./(max(0.1, ensemble_size_per_class(i,j)-1))); %CHANGED
                        %PROFILE.stratigraphy(i, index+k) = j;
                        
                        %waterTable_fraction = max(0, min(1, water_table_intercept(1,j) + fraction .* water_table_slope(1,j)));
                        ensemble.STATVAR.water_table_depth(i,index+k) = max(0, min(1, water_table_intercept(1,j) + fraction .* water_table_slope(1,j)));
                        
                        ensemble.STATVAR.snowfall_factor(i,index+k) = get_snowfall_factor(ensemble, k./(ensemble_size_per_class(i,j)+1) , CV_list(1,j));
                        ensemble.STATVAR.wind_compaction_timescale(i, index+k) = wind_compaction_timescale(1,j);
                        ensemble.STATVAR.wind_speed_class(i, index+k) = wind_speed_per_class(1,j);
                        
                    end
                    index = index + ensemble_size_per_class(i,j);
                end
            end
            ensemble.STATVAR.stratigraphy = ensemble.STATVAR.stratigraphy(:)';
            ensemble.STATVAR.snowfall_factor = ensemble.STATVAR.snowfall_factor(:)';
            ensemble.STATVAR.wind_compaction_timescale = ensemble.STATVAR.wind_compaction_timescale(:)';
            ensemble.STATVAR.water_table_depth = ensemble.STATVAR.water_table_depth(:)';
            ensemble.STATVAR.wind_speed_class = ensemble.STATVAR.wind_speed_class(:)';
            
            ensemble.STATVAR.stratigraphy_ID_list = [];
            for i=1:size(tile.STRATIGRAPHY,2)
                ensemble.STATVAR.stratigraphy_ID_list = [ensemble.STATVAR.stratigraphy_ID_list; tile.STRATIGRAPHY{1,i}.PARA.ID];
            end
            
            subgrid_index = (ensemble.STATVAR.stratigraphy - floor(ensemble.STATVAR.stratigraphy)).*10;  %zero to one
            id_list = floor(ensemble.STATVAR.stratigraphy);
            

            variables = fieldnames(tile.STRATIGRAPHY{1,1}.STATVAR);
            for index = 1:size(variables,1)
                tile.GRID.STATVAR.(variables{index,1}) = [];
                for i=1:size(id_list,2)
                    pos = find(ensemble.STATVAR.stratigraphy_ID_list(:,1) == id_list(1,i));
                    tile.GRID.STATVAR.(variables{index,1}) = [tile.GRID.STATVAR.(variables{index,1}) ...
                        tile.STRATIGRAPHY{1,pos}.STATVAR.(variables{index,1})(:,1) .* subgrid_index(1,i) + tile.STRATIGRAPHY{1,pos}.STATVAR.(variables{index,1})(:,2) .* (1-subgrid_index(1,i)) ];
                end
            end
            
        
            
        
        end
        
        %non-mandatory
        
        function snowfall_factor = get_snowfall_factor(ensemble, percentage, CV)

            av_snow = 1; %100%
            
            psi_square = log(1 + CV.^2);
            lambda = log(av_snow) - 0.5.*psi_square;
            
            snowfall_factor = exp(erfinv(2.*percentage - 1).*sqrt(2.*psi_square) + lambda);
            
        end
        
        

%         function ensemble = generate_ensemble(ensemble, tile) %delivers snowfall_factors and stratigraphy
% 
%             
%             
%             for i=1:size(ensemble_size_per_class,1) %loop over cells
%                 index=0;
%                 for j=1:size(ensemble_size_per_class,2) %loop over the different classes
%                     
%                     if j==4 || j==5 || j==9
%                         
%                         random_list=rand(ensemble_size_per_class(i,j),1);%CHANGED
%                         [dummy, order] = sort(random_list);%CHANGED
%                         
%                     else
%                         
%                         order=[1:ensemble_size_per_class(i,j)]';  %REVISED 2019, no randomization of stratigraphies
%                     end
%                     
%                     for k=1:ensemble_size_per_class(i,j) %insert
%                         fraction = (0.5.*double(ensemble_size_per_class(i,j) == 1) +  (order(k,1)-1)./(max(1, ensemble_size_per_class(i,j)-1)));
%                         stratigraphy_fraction =  max(0, min(1, stratigraphy_intercept (1,j) + fraction .* stratigraphy_slope(1,j)));
%                         PROFILE.stratigraphy(i, index+k) = j + 0.1.* stratigraphy_fraction;
%                         
%                         %PROFILE.stratigraphy(i, index+k) = j + 0.1.* (0.5.*double(ensemble_size_per_class(i,j) == 1) +  (order(k,1)-1)./(max(0.1, ensemble_size_per_class(i,j)-1))); %CHANGED
%                         %PROFILE.stratigraphy(i, index+k) = j;
%                         
%                         %waterTable_fraction = max(0, min(1, water_table_intercept(1,j) + fraction .* water_table_slope(1,j)));
%                         PROFILE.water_table_depth(i,index+k) = max(0, min(1, water_table_intercept(1,j) + fraction .* water_table_slope(1,j)));
%                         
%                         PROFILE.snowfall_factor(i,index+k) = get_snowfall_factor(k./(ensemble_size_per_class(i,j)+1) , CV_list(1,j));
%                         PROFILE.z0(i, index+k) =  z0_list(1,j);
%                         PROFILE.wind_compaction_timescale(i, index+k) = wind_compaction_timescale_list(1,j);
%                         PROFILE.wind_speed_class(i, index+k) = wind_speed_per_class(1,j);
%                         
%                     end
%                     index = index + ensemble_size_per_class(i,j);
%                 end
%             end
% 
%             %wind_compaction_timescale = PROFILE.wind_compaction_timescale;
%             %stratigraphy = PROFILE.stratigraphy;
%             %sf_factor = PROFILE.snowfall_factor;
%             %water_table_depth = PROFILE.water_table_depth;
%             %save('save_stratigraphy.mat', 'ensemble_size_per_class', 'stratigraphy', 'wind_compaction_timescale', 'sf_factor', 'water_table_depth')
%             PROFILE.stratigraphy = PROFILE.stratigraphy(:)';
%             PROFILE.snowfall_factor = PROFILE.snowfall_factor(:)';
%             PROFILE.z0 = PROFILE.z0(:)';
%             PROFILE.wind_compaction_timescale = PROFILE.wind_compaction_timescale(:)';
%             PROFILE.water_table_depth = PROFILE.water_table_depth(:)';
%             PROFILE.wind_speed_class = PROFILE.wind_speed_class(:)';
%         end
        
    end
end

