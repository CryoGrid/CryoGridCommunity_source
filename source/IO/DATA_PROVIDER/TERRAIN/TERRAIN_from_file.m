%========================================================================
% TERRAIN class for simple and predescribed terrain information.
%
% Thomas Ingeman-Nielsen, October 2022
%========================================================================

classdef TERRAIN_from_file < TERRAIN_base
    
    methods
        
        function terrain = provide_PARA(terrain)
           terrain = provide_PARA@TERRAIN_base(terrain);
           terrain.PARA.terrain_path = [];
           terrain.PARA.terrain_file = [];
           terrain.PARA.tp_number = [];
        end
        
        function terrain = provide_CONST(terrain)
            % empty; CONSTs are provided in the relevant terrain classes
        end
        
        function terrain = provide_STATVAR(terrain)
            % empty; STATVARs are provided in the relevant terrain classes
        end
        
        %%%%%% ----- finalize init functions ----- %%%%%%
        
        function terrain = initialize_TEMP(terrain)
            % The TEMP variables required for all classes
        end
        
        function terrain = finalize_init(terrain, tile)
            % All angles in degrees and clockwise from N (maybe change)  
            i = terrain.PARA.tp_number;
            path = terrain.PARA.terrain_path;
            tp_file = terrain.PARA.terrain_file;
            load([path tp_file])
            
%             terrain.TEMP.azimuth = 0;
            terrain.PARA.altitude = tp.z(i,:);
            terrain.PARA.h = rad2deg(tp.h(i,:));
            terrain.PARA.hbins = 180 - 180/pi*tp.hbins;
            terrain.PARA.slope_angle = rad2deg(tp.slp(i));
            terrain.PARA.aspect = 180 - 180/pi*tp.asp(i);
            terrain.PARA.skyview_factor = tp.svf(i);
        end
        
    end
end