%========================================================================
% CryoGrid SPATIAL_REFERENCE class POINT_3D_SIMPLE
% 3D POINT class designed to provide information for several 
% coupled TILE classes that are run in parallel with the RUN_INFO class 
% RUN_3D_POINT. The topological relationships between the TILE classes are
% provided in this class
%
% S. Westermann, Dec 2022
%========================================================================

classdef POINT_3D_SIMPLE < matlab.mixin.Copyable

    properties
        RUN_INFO
        PARA
        CONST
        STATVAR
        TEMP
        ACTION
    end
    
    methods
        function point = provide_PARA(point)
            point.PARA.number_of_tiles = []; %3;
            point.PARA.latitude = [];
            point.PARA.longitude = [];
            point.PARA.altitude = []; 
            point.PARA.area = [];
            point.PARA.param_file_number = []; %[1;2;3];

            point.PARA.connected = [];
            point.PARA.contact_length = [];
            point.PARA.distance = [];
            
            %can be provided optionally
%             %provide default values
%             point.PARA.slope_angle = 0;     
%             point.PARA.aspect = 0;     
%             point.PARA.skyview_factor = 0;     
%             point.PARA.horizon_bins = 0;
%             point.PARA.horizon_angles = 0;
            
        end
        
        function point = provide_STATVAR(point)

        end
        
        function point = provide_CONST(point)
            
        end
        
        function point = finalize_init(point)
            %provide default values
            point.PARA.slope_angle = 0;     
            point.PARA.aspect = 0;     
            point.PARA.skyview_factor = 0;     
            point.PARA.horizon_bins = 0;
            point.PARA.horizon_angles = 0;
            
            if size(point.PARA.latitude,1) == 1
                point.STATVAR.latitude = repmat(point.PARA.latitude, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.latitude = point.PARA.latitude;
            end
            if size(point.PARA.longitude,1) == 1
                point.STATVAR.longitude = repmat(point.PARA.longitude, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.longitude = point.PARA.longitude;
            end
            if size(point.PARA.altitude,1) == 1
                point.STATVAR.altitude = repmat(point.PARA.altitude, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.altitude = point.PARA.altitude;
            end
            if size(point.PARA.area,1) == 1
                point.STATVAR.area = repmat(point.PARA.area, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.area = point.PARA.area;
            end
            if size(point.PARA.slope_angle,1) == 1
                point.STATVAR.slope_angle = repmat(point.PARA.slope_angle, point.PARA.number_of_tiles, 1);
            else
                point.PARA.slope_angle = point.PARA.slope_angle;
            end
            if size(point.PARA.slope_angle,1) == 1
                point.STATVAR.slope_angle = repmat(point.PARA.slope_angle, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.slope_angle = point.PARA.slope_angle;
            end
            if size(point.PARA.aspect,1) == 1
                point.STATVAR.aspect = repmat(point.PARA.aspect, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.aspect = point.PARA.aspect;
            end
            if size(point.PARA.skyview_factor,1) == 1
                point.STATVAR.skyview_factor = repmat(point.PARA.skyview_factor, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.skyview_factor = point.PARA.skyview_factor;
            end
            if size(point.PARA.horizon_angles,1) == 1
                point.STATVAR.horizon_angles = repmat(point.PARA.horizon_angles, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.horizon_angles = point.PARA.horizon_angles;
            end
            if size(point.PARA.horizon_bins,1) == 1
                point.STATVAR.horizon_bins = repmat(point.PARA.horizon_bins, point.PARA.number_of_tiles, 1);
            else
                point.STATVAR.horizon_bins = point.PARA.horizon_bins;
            end
        end
        
        
        
        %-------------param file generation-----
        function point = param_file_info(point)
            point = provide_PARA(point);
            
            point.PARA.STATVAR = [];
            point.PARA.class_category = 'SPATIAL_REFERENCE';
            point.PARA.default_value = [];
            
            point.PARA.default_value.number_of_tiles = {3};
            point.PARA.comment.number_of_tiles = {'number of tiles/cores'};
            
            point.PARA.comment.latitude = {'latitude in decimal degrees, either single value assumed for all tiles, or list for all TILEs'};
            
            point.PARA.comment.longitude = {'longitude in decimal degrees, either single value assumed for all tiles, or list for all TILEs'};
            
            point.PARA.comment.altitude = {'altitude in m, either single value assumed for all tiles, or list for all TILEs'};
            point.PARA.options.altitude.name = 'H_LIST';
            point.PARA.options.altitude.entries_x = {31 32.4 33};
            
            point.PARA.comment.area = {'area in m2, either single value assumed for all tiles, or list for all TILEs'};
            point.PARA.default_value.area = {1};
            
            point.PARA.comment.param_file_number = {'index of parameter file used for each TILE'};
            point.PARA.options.param_file_number.name = 'H_LIST';
            point.PARA.options.param_file_number.entries_x = {'1' '2' '3'};

            point.PARA.comment.connected = {'matrix of conectivity between TILES, 1: connected; 0: not connected'};
            point.PARA.options.connected.name = 'MATRIX';
            point.PARA.options.connected.entries_matrix = {'0' '1' '0'; '1' '0' '1'; '0' '1' '0'};

            point.PARA.comment.contact_length = {'contact length in m between TILES'};
            point.PARA.options.contact_length.name = 'MATRIX';
            point.PARA.options.contact_length.entries_matrix = {'0' '1' '0'; '1' '0' '1'; '0' '1' '0'};

            point.PARA.comment.distance = {'distance in m between TILES'};
            point.PARA.options.distance.name = 'MATRIX';
            point.PARA.options.distance.entries_matrix = {'0' '1' '0'; '1' '0' '1'; '0' '1' '0'};
        end
        
    end
end

