classdef DEM_gmted < matlab.mixin.Copyable

    properties
        PARA
        CONST
        STATVAR
        TEMP
    end

    
    methods
        
        function dem = provide_PARA(dem)
            dem.PARA.DEM_file = [];
            dem.PARA.DEM_path = [];
        end
        
        function dem = provide_STATVAR(dem)

        end
        
        function dem = provide_CONST(dem)
            
        end
        
        function dem = finalize_init(dem)
            
        end
        

        function altitudes = get_altitudes(dem, run_info)
             altitudes=[];
            for i=1:size(run_info.STATVAR.list_of_MODIS_tiles, 1)
                
                %get corner coordinates of ROI
                minLat=min(run_info.STATVAR.latitude(run_info.STATVAR.list_of_MODIS_tiles(i,3):run_info.STATVAR.list_of_MODIS_tiles(i,4), 1));
                maxLat=max(run_info.STATVAR.latitude(run_info.STATVAR.list_of_MODIS_tiles(i,3):run_info.STATVAR.list_of_MODIS_tiles(i,4), 1));
                minLong=min(run_info.STATVAR.longitude(run_info.STATVAR.list_of_MODIS_tiles(i,3):run_info.STATVAR.list_of_MODIS_tiles(i,4), 1));
                maxLong=max(run_info.STATVAR.longitude(run_info.STATVAR.list_of_MODIS_tiles(i,3):run_info.STATVAR.list_of_MODIS_tiles(i,4), 1));
                
                %read global DEM - this can be problematic, 17GB file size!!
                [DEM, R]= geotiffread([dem.PARA.DEM_path dem.PARA.DEM_file]);
                DEM_latitudes= R.LatitudeLimits;
                DEM_longitudes= R.LongitudeLimits;
                
                %make a separate list of longitudes and latitudes pixels from DEM
                DEM_latitude_list=linspace(max(DEM_latitudes), min(DEM_latitudes), size(DEM,1))';
                DEM_longitude_list=linspace(min(DEM_longitudes), max(DEM_longitudes), size(DEM,2));
                
                %Find the indeces of coordinates in the list that corresponds with MODIS extent coordinates
                [mini, maxPosLat]=min(abs(maxLat-DEM_latitude_list)); %finding an index of an value where where the difference between MODIS latitude and DEM latitude is the smallest
                maxPosLat=max(1,maxPosLat-5); %expanding the border by 5 pixels unless out of a list range
                
                [mini, minPosLat]=min(abs(minLat-DEM_latitude_list));
                minPosLat=min(size(DEM_latitude_list,1),minPosLat+5);
                
                [mini, minPosLon]=min(abs(minLong-DEM_longitude_list));
                minPosLon=max(1,minPosLon-5);
                
                [mini, maxPosLon]=min(abs(maxLong-DEM_longitude_list));
                maxPosLon=min(size(DEM_longitude_list,2),maxPosLon+5);
                
                %subset the DEM based on indexes from coordinate lists to MODIS extent
                DEM=DEM(maxPosLat:minPosLat, minPosLon:maxPosLon);
                DEM=double(DEM);
                
                %subset DEM coordinate vectors
                DEM_latitude_list = DEM_latitude_list(maxPosLat:minPosLat);
                DEM_longitude_list = DEM_longitude_list(minPosLon:maxPosLon);
                
                %create coordinate matrices that are used for interpolation reference
                [DEM_latitudes, DEM_longitudes]=meshgrid(DEM_latitude_list, DEM_longitude_list);
                DEM_latitudes = DEM_latitudes';
                DEM_longitudes = DEM_longitudes';
                
                %interpolation to a MODIS grid usind cubic interpolation
                altitudes = [altitudes; interp2(DEM_latitudes', DEM_longitudes', DEM', run_info.STATVAR.latitude(run_info.STATVAR.list_of_MODIS_tiles(i,3):run_info.STATVAR.list_of_MODIS_tiles(i,4), 1), run_info.STATVAR.longitude(run_info.STATVAR.list_of_MODIS_tiles(i,3):run_info.STATVAR.list_of_MODIS_tiles(i,4), 1), 'cubic')];
                
            end
            altitudes(altitudes < 0) = 0;
            
        end
        
        
    end
end

