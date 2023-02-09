%========================================================================
% CryoGrid DATA_PROVIDER class DEM
% DEM class deriving information from a digital 
% elevation model (DEM),including slope, aspect (terrain shading still 
% lacking). The class can ingest Copernicus 30m DEMs downloaded from 
% https://portal.opentopography.org/   
%
% S. Westermann, Dec 2022
%========================================================================

classdef DEM_CCI < DEM_BASE

    properties
        PARENT
    end
    
    methods
        function dem = provide_PARA(dem)            
%             dem.PARA.variables = []; % altitude or altitude, slope_angle, aspect or altitude, slope_angle, aspect, horizon_angles
%                         
%             dem.PARA.DEM_folder = [];
%             dem.PARA.DEM_filename = [];
%             
%             dem.PARA.reproject2utm = [];

        end
        
        function dem = provide_STATVAR(dem)

        end
        
        function dem = provide_CONST(dem)
            
        end
        
        function dem = finalize_init(dem)
            
%             dem = read_DEM_raster(dem);
            
%             dem.STATVAR.Z = readgeoraster([dem.PARA.DEM_folder dem.PARA.DEM_filename]);
%             info = georasterinfo([dem.PARA.DEM_folder dem.PARA.DEM_filename]);
%             dem.STATVAR.info = info;
%             dem.STATVAR.Z(dem.STATVAR.Z==single(info.MissingDataIndicator)) = NaN;
%             xlim = info.RasterReference.XWorldLimits; 
%             ylim = info.RasterReference.YWorldLimits;
%             dx = info.RasterReference.CellExtentInWorldX; 
%             dy = info.RasterReference.CellExtentInWorldY;
%             x=((xlim(1)+dx/2):dx:(xlim(end)-dx/2))';
%             y=flipud(((ylim(1)+dy/2):dy:(ylim(end)-dy/2))');
%             [dem.STATVAR.X,dem.STATVAR.Y] = meshgrid(x,y);  
%             
        end
        
        function dem = load_data(dem)
%             if ~strcmp(dem.TEMP.info.RasterReference.CoordinateSystemType, 'geographic') 
%                 %[dem.TEMP.X_target, dem.TEMP.Y_target] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARA.latitude, dem.PARA.longitude);
%                 [dem.TEMP.X_target, dem.TEMP.Y_target] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARENT.STATVAR.latitude, dem.PARENT.STATVAR.longitude);
%             else
%                 if dem.PARA.reproject2utm
%                     [dem.TEMP.X_target,dem.TEMP.Y_target] = ll2utm(dem, dem.PARENT.STATVAR.latitude, dem.PARENT.STATVAR.longitude, dem.TEMP.utm_zone);
%                 else
%                     dem.TEMP.X_target = dem.PARENT.STATVAR.longitude;
%                     dem.TEMP.Y_target = dem.PARENT.STATVAR.latitude;
%                     %[LAT,LON]=projinv(info.CoordinateReferenceSystem,X,Y);
%                 end
%             end
%             
%             
%             %[dem.TEMP.X_target, dem.TEMP.Y_target] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARENT.STATVAR.latitude, dem.PARENT.STATVAR.longitude);
%             
%             dem.PARENT.STATVAR.altitude = dem.PARENT.STATVAR.longitude.*0;
%             dem.PARENT.STATVAR.slope_angle = dem.PARENT.STATVAR.longitude.*0;
%             dem.PARENT.STATVAR.aspect = dem.PARENT.STATVAR.longitude.*0;
%             dem.PARENT.STATVAR.skyview_factor = dem.PARENT.STATVAR.longitude.*0;
%             dem.PARENT.STATVAR.horizon_bins = dem.PARENT.STATVAR.longitude.*0;
%             dem.PARENT.STATVAR.horizon_angles = dem.PARENT.STATVAR.longitude.*0;
%             
%             if ~strcmp(dem.TEMP.info.RasterReference.CoordinateSystemType, 'geographic')
%                 [X1, Y1] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARENT.STATVAR.latitude - 0.1, dem.PARENT.STATVAR.longitude);
%                 [X2, Y2] = projfwd(dem.TEMP.info.CoordinateReferenceSystem, dem.PARENT.STATVAR.latitude + 0.1, dem.PARENT.STATVAR.longitude);
%                 dem.TEMP.offset_angle_trueNorth = atand((X1-X2)./(Y1-Y2)); %offset grid North and true North
%             else
%                 if dem.PARA.reproject2utm
%                     [X1,Y1] = ll2utm(dem, dem.PARENT.STATVAR.latitude-0.1, dem.PARENT.STATVAR.longitude, dem.TEMP.utm_zone);
%                     [X2,Y2] = ll2utm(dem, dem.PARENT.STATVAR.latitude+0.1, dem.PARENT.STATVAR.longitude, dem.TEMP.utm_zone);
%                     dem.TEMP.offset_angle_trueNorth = atand((X1-X2)./(Y1-Y2)); %offset grid North and true North
%                 else
%                     dem.TEMP.offset_angle_trueNorth = 0;
%                 end
%             end
%             
%             for i=1:size(dem.PARA.variables,1)
%                 a = str2func(['get_' dem.PARA.variables{i,1}]);
%                 dem = a(dem);
%             end
%             dem.STATVAR = []; 
%             dem.TEMP=[];
            
            dem.PARENT.STATVAR.altitude = dem.PARENT.STATVAR.latitude .* 0;
        end
        
        %service functions using DEM_base
        
        function dem = get_altitude(dem)
            
            dem = get_altitude_base(dem);
            dem.PARENT.STATVAR.altitude=  double(dem.STATVAR.altitude);
           % dem.PARENT.STATVAR.altitude = double(interp2(dem.STATVAR.X, dem.STATVAR.Y, dem.STATVAR.Z, dem.STATVAR.X_target, dem.STATVAR.Y_target));
            
        end
       
        
        

        
%         %-------------param file generation-----
%         function point = param_file_info(point)
%             point = provide_PARA(point);
%             
%             point.PARA.STATVAR = [];
%             point.PARA.class_category = 'DATA_PROVIDER';
%             
%             point.PARA.comment.variables = {'properties calculated from DEM: altitude OR altitude, slope_angle, aspect'};
%             point.PARA.options.variables.name = 'H_LIST';
%             point.PARA.options.variables.entries_x = {'altitude' 'slope_angle' 'aspect'};
%                         
% %             point.PARA.comment.number_of_horizon_bins = {'number of angular points for which horizon is calculated; must be multiple of 4'};
% %             point.PARA.default_value.number_of_horizon_bins = {24};
%             
%             point.PARA.comment.DEM_folder = {'folder in which DEM file is located'}; 
%             
%             point.PARA.comment.DEM_filename = {'name of DEM file'}; 
%             
%             point.PARA.comment.reproject2utm = {'select 1 when using a DEM in geographic coordinates (or similar) and computing more than just altitude; select 0 to speed up altitde computation in big DEMs'};
%             point.PARA.default_value.reproject2utm = {1};
%         end
        
    end
end

