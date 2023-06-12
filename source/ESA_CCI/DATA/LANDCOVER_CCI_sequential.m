classdef LANDCOVER_CCI_sequential < matlab.mixin.Copyable

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
            
            lc.PARA.deg_tile_list_folder = [];
            lc.PARA.deg_tile_list_file = [];

            lc.PARA.yedoma_file = [];
            lc.PARA.yedoma_path = [];
            
            lc.PARA.accumulated{1,1} = [140 150 152 153]; %sparse vegetation
            lc.PARA.accumulated{2,1} = [10  11  12  20 130]; %grasslands and croplands
            lc.PARA.accumulated{3,1} = [30 40 100 110 120 121 122]; %shrubs
            lc.PARA.accumulated{4,1} = [50 60 61 62 80 81 82 90]; %deciduous forest
            lc.PARA.accumulated{5,1} = [70 71 72]; %evergreen forest
            lc.PARA.accumulated{6,1} = [160 170 180]; %wetlands
            lc.PARA.accumulated{7,1} = [190 200 201 202 220]; %bare areas and urban
            
            %Yedoma
            lc.PARA.accumulated{8,1} = [140 150 152 153]; %sparse vegetation
            lc.PARA.accumulated{9,1} = [50 60 61 62 80 81 82 90]; %deciduous forest
            lc.PARA.accumulated{10,1} = [160 170 180]; %wetlands
            lc.PARA.urban = 190;
        end
        
        function lc = provide_STATVAR(lc)

        end
        
        function lc = provide_CONST(lc)
            
        end
        
        function lc = finalize_init(lc)
            
        end
        
        function lc = load_data(lc)
            
            load([lc.PARA.deg_tile_list_folder lc.PARA.deg_tile_list_file], 'MODIS_deg_list');
            
            tic
            lc.PARENT.STATVAR.landcover = repmat(lc.PARENT.STATVAR.latitude .* 0, 1,10);
            
            lc.PARENT.STATVAR.urban = lc.PARENT.STATVAR.latitude .* 0;
            for i=1:size(lc.PARENT.STATVAR.latitude,1)
                
                max_lat = lc.PARENT.STATVAR.latitude(i)+0.1;
                min_lat = lc.PARENT.STATVAR.latitude(i)-0.1;
                max_lon = lc.PARENT.STATVAR.longitude(i)+0.1;
                min_lon = lc.PARENT.STATVAR.longitude(i)-0.1;
                delta_1km_lat  = 0.01;
                delta_1km_lon  = MODIS_deg_list(lc.PARENT.STATVAR.tile_number(i), 5);
                
                delta_300m_lat = 1/360;
                delta_300m_lon = 1/360;
                
                roi_start_index_lat = max(1, round((90 - (max_lat + delta_1km_lat/2))./delta_300m_lat + 1)-1);
                roi_start_index_lon = max(1, round((180 + (min_lon - delta_1km_lon/2))./delta_300m_lon + 1)-1);
                
                %no check for right boundary performed
                number_of_elements_lat = min(180/delta_300m_lat, round((max_lat - min_lat + delta_1km_lat) ./ delta_300m_lat)+2);
                number_of_elements_lon = min(360/delta_300m_lon, round((max_lon - min_lon + delta_1km_lon) ./ delta_300m_lon)+2);
                
                landcover = double(ncread([lc.PARA.landcover_path lc.PARA.landcover_file], 'lc', [roi_start_index_lat roi_start_index_lon], [number_of_elements_lat number_of_elements_lon], [1 1]));
                
                target_latitude = lc.PARENT.STATVAR.latitude(i,1);
                target_longitude = lc.PARENT.STATVAR.longitude(i,1); %make sure longitude is indeed -180->180
                
                index_1km_lat = round((90 - target_latitude + delta_1km_lat/2) ./ delta_1km_lat);
                index_1km_lon = round((180 + target_longitude + delta_1km_lon/2) ./ delta_1km_lon);
                
                start_index_300m_lat = round((index_1km_lat - 1) .*  delta_1km_lat ./ delta_300m_lat.*1e6) ./1e6; %0 ; 3.6
                end_index_300m_lat = round(index_1km_lat .*  delta_1km_lat ./ delta_300m_lat .* 1e6) ./1e6; % 3.6  ; 7.2
                start_index_300m_lon = round((index_1km_lon - 1) .*  delta_1km_lon ./ delta_300m_lon .*1e6) ./1e6;
                end_index_300m_lon = round(index_1km_lon .*  delta_1km_lon ./ delta_300m_lon .* 1e6) ./1e6;
                
                start_index_300m_lat2 = floor(start_index_300m_lat + 1); %1 ; 4
                end_index_300m_lat2 = floor(end_index_300m_lat + 1); %4 ; 8
                start_index_300m_lon2 = floor(start_index_300m_lon + 1);
                end_index_300m_lon2 = floor(end_index_300m_lon + 1);
                
                weight_300m = ones(end_index_300m_lat2 - start_index_300m_lat2 + 1, end_index_300m_lon2 - start_index_300m_lon2 + 1);
                weight_300m(1,:) = weight_300m(1,:) .*(start_index_300m_lat2 - start_index_300m_lat);
                weight_300m(end,:) = weight_300m(end,:) .*(-(end_index_300m_lat2 - 1) + end_index_300m_lat);
                weight_300m(:,1) = weight_300m(:,1) .*(start_index_300m_lon2 - start_index_300m_lon);
                weight_300m(:,end) = weight_300m(:,end) .*(-(end_index_300m_lon2 - 1) + end_index_300m_lon);
                
                weight_300m = weight_300m ./ (delta_1km_lat ./ delta_300m_lat) ./ (delta_1km_lon ./ delta_300m_lon); %normalize
                
                start_index_300m_lat3 = start_index_300m_lat2 - roi_start_index_lat + 1;
                end_index_300m_lat3 = end_index_300m_lat2 - roi_start_index_lat + 1;
                start_index_300m_lon3 = start_index_300m_lon2 - roi_start_index_lon + 1;
                end_index_300m_lon3 = end_index_300m_lon2 - roi_start_index_lon + 1;
                
                if start_index_300m_lat3 <1 || start_index_300m_lon3 <1
                    a=1;
                end
                
                for k = 1:7%size(lc.PARA.accumulated,1) %no yedoma yet
                    land_cover_fraction = weight_300m.* 0;
                    for j=1:size(lc.PARA.accumulated{k,1},2)
                        land_cover_fraction = land_cover_fraction + double(landcover(start_index_300m_lat3:end_index_300m_lat3,start_index_300m_lon3:end_index_300m_lon3) == lc.PARA.accumulated{k,1}(1,j));
                        
                    end
                    lc.PARENT.STATVAR.landcover(i,k) = sum(weight_300m(:) .* land_cover_fraction(:));
                end
                urban_fraction = double(landcover(start_index_300m_lat3:end_index_300m_lat3,start_index_300m_lon3:end_index_300m_lon3) == lc.PARA.urban);
                lc.PARENT.STATVAR.urban(i,1) = sum(weight_300m(:) .* urban_fraction(:));
                
            end
        
            %normalize to account for water bodies
            lc.PARENT.STATVAR.landcover = lc.PARENT.STATVAR.landcover ./ repmat(sum(lc.PARENT.STATVAR.landcover,2), 1, size(lc.PARENT.STATVAR.landcover,2));
            toc
            
            
            %add grassland class to tundra class for areas out of Central Asia
            grasslandCentralAsia = lc.PARENT.STATVAR.landcover(:,2) .* double(lc.PARENT.STATVAR.longitude > 50 & lc.PARENT.STATVAR.latitude < 55);
            grasslandOther = lc.PARENT.STATVAR.landcover(:,2) - grasslandCentralAsia;
            lc.PARENT.STATVAR.landcover(:,1) = lc.PARENT.STATVAR.landcover(:,1) + grasslandOther;
            lc.PARENT.STATVAR.landcover(:,2) = grasslandCentralAsia;
            
            
            %Yedoma map         
            delta_yedoma = double(ncread([lc.PARA.yedoma_path lc.PARA.yedoma_file], 'delta_lat_lon'));
            
            max_lat = max(lc.PARENT.STATVAR.latitude);
            min_lat = min(lc.PARENT.STATVAR.latitude);
            max_lon = max(lc.PARENT.STATVAR.longitude);
            min_lon = min(lc.PARENT.STATVAR.longitude);
            
            roi_start_index_lat = max(1, round((90 - (max_lat + delta_1km_lat/2))./delta_yedoma + 1)-1);
            roi_start_index_lon = max(1, round((180 + (min_lon - delta_1km_lon/2))./delta_yedoma + 1)-1);
            
            number_of_elements_lat = min(180/delta_yedoma, round((max_lat - min_lat + delta_1km_lat) ./ delta_yedoma)+2);
            number_of_elements_lon = min(360/delta_yedoma, round((max_lon - min_lon + delta_1km_lon) ./ delta_yedoma)+2);
            yedoma = double(ncread([lc.PARA.yedoma_path lc.PARA.yedoma_file], 'yedoma', [roi_start_index_lat roi_start_index_lon], [number_of_elements_lat number_of_elements_lon], [1 1]));
            
            lat_yedoma_start = 90-delta_yedoma/2 - (roi_start_index_lat-1) .* delta_yedoma; %latitude of 1st Yedoma pixel in ROI
            lon_yedoma_start = -180 + delta_yedoma/2 + (roi_start_index_lon-1) .* delta_yedoma; %longitude of 1st Yedoma pixel in ROI
            
            factor_lat = delta_1km_lat ./ delta_yedoma;
            factor_lon = delta_1km_lon ./ delta_yedoma;
            yedoma_lc = lc.PARENT.STATVAR.latitude.*0;
            delta_index_lat = (factor_lat-1)/2; 
            delta_index_lon = (factor_lon-1)/2;
            for i=1:size(lc.PARENT.STATVAR.latitude,1)
                center_pixel_lat = (lat_yedoma_start - lc.PARENT.STATVAR.latitude(i,1)) ./ delta_yedoma + 1;
                center_pixel_lon = (-lon_yedoma_start + lc.PARENT.STATVAR.longitude(i,1)) ./ delta_yedoma +1;
                offset_lat = delta_index_lat;
                offset_lon = delta_index_lon;
                yedoma_lc(i,1) = yedoma_lc(i,1) + sum(sum(yedoma(round(center_pixel_lat-offset_lat):round(center_pixel_lat+offset_lat), round(center_pixel_lon-offset_lon):round(center_pixel_lon+offset_lon))));
            end
            yedoma_lc = yedoma_lc ./ factor_lon ./ factor_lat; %check that normalized
            
            lc.PARENT.STATVAR.landcover(:,8) = lc.PARENT.STATVAR.landcover(:,1) .* yedoma_lc;
            lc.PARENT.STATVAR.landcover(:,1) = lc.PARENT.STATVAR.landcover(:,1) .* (1-yedoma_lc);
            lc.PARENT.STATVAR.landcover(:,9) = lc.PARENT.STATVAR.landcover(:,4) .* yedoma_lc;
            lc.PARENT.STATVAR.landcover(:,4) = lc.PARENT.STATVAR.landcover(:,4) .* (1-yedoma_lc);
            % wetland Yedoma not taken into account now
            %INFO.landcover(:,6) = INFO.landcover(:,6) .* (1-Yedoma_yes_no);
%             lc_list(:,10) = lc_list(:,10) .* 0; %no effect on wetland
            %INFO.landcover(:,6) = INFO.landcover(:,6) .* (1-Yedoma_yes_no);
            %INFO.landcover(:,10) = INFO.landcover(:,10) .* Yedoma_yes_no;

            lc.PARENT.STATVAR.landcover(lc.PARENT.STATVAR.landcover<0) = 0; %eliminate rounding errors
       end


    end
end

