classdef LANDCOVER_CCI < matlab.mixin.Copyable

    
    properties
        RUN_INFO
        PARA
        CONST
        STATVAR
        
    end
    
    methods
        
        function lc = provide_PARA(lc)
            lc.PARA.landcover_path = [];
            lc.PARA.filename = [];
            lc.PARA.yedoma_path = [];
            lc.PARA.yedoma_filename = [];
        end
        
        function lc = provide_STATVAR(lc)

        end
        
        function lc = provide_CONST(lc)
            
        end
        
        function lc = finalize_init(lc)

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
        end

        
        function [landcover_list, urban] = get_landcover(lc, run_info) %landcover_list = get_landcover(lc, run_info)
           
            landcover_list = [];     
            urban = []; %CHANGED
           
            for index = 1:size(run_info.STATVAR.list_of_MODIS_tiles,1)
                
                range = run_info.STATVAR.list_of_MODIS_tiles(index,3):run_info.STATVAR.list_of_MODIS_tiles(index,4);
                
                h = run_info.STATVAR.list_of_MODIS_tiles(index,1);
                h=h+100;
                h=num2str(h);
                h=h(2:3);
                v = run_info.STATVAR.list_of_MODIS_tiles(index,2);
                v=v+100;
                v=num2str(v);
                v=v(2:3);
                lc.PARA.filename(2:3) = h;
                lc.PARA.filename(5:6) = v;
                filename = lc.PARA.filename(1:16);
                extension = lc.PARA.filename(19:end);
                
                
                final=[];
                classes=[];
                for i=1:220
                    if exist([lc.PARA.landcover_path filename num2str(i) extension])==2
                        load([lc.PARA.landcover_path filename num2str(i) extension]);
                        
                        classes=[classes; i];
                        final = [final subcellstat(run_info.STATVAR.key(range))];
                        
                        %CHANGED
                        if i==190
                            urban = [urban; subcellstat(run_info.STATVAR.key(range))];
                        end
                    end
                end
                
                
                
                lc_list =[];
                
                for i=1:size(lc.PARA.accumulated,1)
                    acc_percentage = 0;
                    for j=1:size(lc.PARA.accumulated{i,1},2)
                        acc_percentage = acc_percentage + final(:,find(classes(:,1) == lc.PARA.accumulated{i,1}(1,j)));
                    end
                    lc_list = [lc_list acc_percentage];
                end
                
                
                %add grassland class to bare class for areas out of Central Asia
                grasslandCentralAsia = lc_list(:,2) .* double(run_info.STATVAR.longitude(range) > 50 & run_info.STATVAR.latitude(range) < 55);
                grasslandOther = lc_list(:,2) - grasslandCentralAsia;
                lc_list(:,1) = lc_list(:,1) + grasslandOther;
                lc_list(:,2) = grasslandCentralAsia;
                
                
                yedoma_filename = lc.PARA.yedoma_filename;
                yedoma_filename(2:3) = h;
                yedoma_filename(5:6) = v;
                %Yedoma
                if exist([lc.PARA.yedoma_path yedoma_filename])==2
                    load([lc.PARA.yedoma_path yedoma_filename]);
                else
                    subcellstat = zeros(1200,1200);
                end
                Yedoma_yes_no = subcellstat(run_info.STATVAR.key(range));
                lc_list(:,1) = lc_list(:,1) .* (1-Yedoma_yes_no);
                lc_list(:,8) = lc_list(:,8) .* Yedoma_yes_no;
                lc_list(:,4) = lc_list(:,4) .* (1-Yedoma_yes_no);
                lc_list(:,9) = lc_list(:,9) .* Yedoma_yes_no;
                %INFO.landcover(:,6) = INFO.landcover(:,6) .* (1-Yedoma_yes_no);
                lc_list(:,10) = lc_list(:,10) .* 0; %no effect on wetland
                %INFO.landcover(:,6) = INFO.landcover(:,6) .* (1-Yedoma_yes_no);
                %INFO.landcover(:,10) = INFO.landcover(:,10) .* Yedoma_yes_no;
                
                
                lc_list(lc_list<0) = 0; %eliminate rounding errors
                    
                
                landcover_list = [landcover_list; lc_list];
            end
            
        end
        
        
    end
end

